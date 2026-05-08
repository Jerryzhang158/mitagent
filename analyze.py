"""
analyze.py — 统计分析模块
- 单方法 Top-K 排名
- 多方法交叉排名（Borda Count 融合）
- 年度趋势分析
- CSV 导出
"""

import csv
import logging
import os
from collections import Counter, defaultdict
from pathlib import Path

import config

log = logging.getLogger(__name__)


# ── Top-K 排名 ────────────────────────────────────────────────────────────────

def top_k(counter: Counter, k: int = None) -> list[tuple[str, int]]:
    k = k or config.TOP_K_DISPLAY
    return counter.most_common(k)


# ── Borda Count 融合多方法排名 ─────────────────────────────────────────────────

def borda_merge(method_counters: dict[str, Counter], top_n: int = 500) -> Counter:
    """
    对多个方法各取 top_n，用 Borda Count 融合排名。
    每个词在某方法中排第 i 名 → 得到 (top_n - i) 分。
    """
    borda_scores: dict[str, float] = defaultdict(float)

    for method, counter in method_counters.items():
        ranked = counter.most_common(top_n)
        n = len(ranked)
        for rank, (kw, _) in enumerate(ranked):
            borda_scores[kw] += (n - rank)

    return Counter(borda_scores)


# ── 年度趋势分析 ──────────────────────────────────────────────────────────────

def yearly_trend(
    records: list[dict],
    top_keywords: list[str],
    method: str = "author_kw",
) -> dict[str, dict[str, int]]:
    """
    返回 {keyword: {year: count}} 的嵌套字典
    method: 使用哪个字段来统计（author_kw / mesh / abstract_text）
    """
    # 按年份分组记录
    by_year: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        year = r.get("year", "")
        if year and year.isdigit():
            by_year[year].append(r)

    trend = {kw: {} for kw in top_keywords}

    for year, yr_records in sorted(by_year.items()):
        # 对该年份重新统计
        if method == "author_kw":
            from extract import extract_author_keywords
            counter = extract_author_keywords(yr_records)
        elif method == "mesh":
            from extract import extract_mesh_keywords
            counter = extract_mesh_keywords(yr_records)
        else:
            from extract import extract_tfidf
            counter = extract_tfidf(yr_records)

        for kw in top_keywords:
            trend[kw][year] = counter.get(kw, 0)

    return trend


# ── 导出 CSV ──────────────────────────────────────────────────────────────────

def export_csv(
    method_counters: dict[str, Counter],
    borda: Counter,
    output_dir: str,
):
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # 1. 各方法分别导出
    for method, counter in method_counters.items():
        path = os.path.join(output_dir, f"keywords_{method}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["rank", "keyword", "count"])
            for rank, (kw, cnt) in enumerate(counter.most_common(config.TOP_K_DISPLAY), 1):
                writer.writerow([rank, kw, cnt])
        log.info(f"导出: {path}")

    # 2. Borda 融合排名导出
    path = os.path.join(output_dir, "keywords_borda_merged.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "keyword", "borda_score"])
        for rank, (kw, score) in enumerate(borda.most_common(config.TOP_K_DISPLAY), 1):
            writer.writerow([rank, kw, f"{score:.1f}"])
    log.info(f"导出: {path}")

    # 3. 综合对照表（多方法同一个文件）
    all_kws = set()
    for counter in method_counters.values():
        all_kws.update(kw for kw, _ in counter.most_common(config.TOP_K_DISPLAY))

    path = os.path.join(output_dir, "keywords_comparison.csv")
    methods = list(method_counters.keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["keyword", "borda_score"] + methods)
        for kw, bscore in borda.most_common(config.TOP_K_DISPLAY):
            row = [kw, f"{bscore:.1f}"]
            for m in methods:
                row.append(method_counters[m].get(kw, 0))
            writer.writerow(row)
    log.info(f"导出对照表: {path}")


# ── 打印摘要到 console ─────────────────────────────────────────────────────────

def print_summary(borda: Counter, method_counters: dict[str, Counter], n: int = 30):
    methods = list(method_counters.keys())
    print("\n" + "═" * 70)
    print(f"  Top {n} Functional Keywords (Borda Count 融合排名)")
    print("═" * 70)
    header = f"{'Rank':>4}  {'Keyword':<35}" + "".join(f"{m:>12}" for m in methods)
    print(header)
    print("─" * 70)
    for rank, (kw, bscore) in enumerate(borda.most_common(n), 1):
        row = f"{rank:>4}  {kw:<35}"
        for m in methods:
            cnt = method_counters[m].get(kw, 0)
            row += f"{cnt:>12}"
        print(row)
    print("═" * 70 + "\n")

"""
main.py — 主程序入口
运行方式：
    python main.py                    # 完整流程（抓取 + 提取 + 分析 + 可视化）
    python main.py --skip-fetch       # 跳过抓取，直接用缓存
    python main.py --query "miRNA"    # 临时覆盖搜索词
    python main.py --max 2000         # 临时覆盖最大篇数
"""

import argparse
import logging
import os
import sys

# ── 依赖检查 ───────────────────────────────────────────────────────────────────
REQUIRED = {
    "Bio":          "biopython",
    "sklearn":      "scikit-learn",
    "matplotlib":   "matplotlib",
    "yake":         "yake",
    "wordcloud":    "wordcloud",
}

def check_deps():
    missing = []
    for mod, pkg in REQUIRED.items():
        try:
            __import__(mod)
        except ImportError:
            missing.append(pkg)
    if missing:
        print("以下依赖未安装，请运行：")
        print(f"  pip install {' '.join(missing)}")
        sys.exit(1)

check_deps()

import config
from fetch    import search_pmids, fetch_abstracts, load_all_cached
from extract  import run_extraction
from analyze  import borda_merge, export_csv, print_summary, yearly_trend
from visualize import plot_all_barplots, plot_wordcloud, plot_trend

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="PubMed Functional Keyword Analyzer")
    parser.add_argument("--skip-fetch", action="store_true",
                        help="跳过抓取步骤，直接从缓存加载")
    parser.add_argument("--query", type=str, default=None,
                        help="临时覆盖 config.SEARCH_QUERY")
    parser.add_argument("--max", type=int, default=None,
                        help="临时覆盖 config.MAX_RECORDS")
    parser.add_argument("--method", type=str, default=None,
                        help="临时覆盖 config.EXTRACTION_METHOD")
    parser.add_argument("--no-trend", action="store_true",
                        help="跳过年度趋势分析（节省时间）")
    return parser.parse_args()


def main():
    args = parse_args()

    # 临时覆盖配置
    if args.query:
        config.SEARCH_QUERY = args.query
    if args.max:
        config.MAX_RECORDS = args.max
    if args.method:
        config.EXTRACTION_METHOD = args.method

    os.makedirs("data",           exist_ok=True)
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    # ── Step 1: 抓取 ──────────────────────────────────────────────────────────
    if args.skip_fetch:
        log.info("跳过抓取，加载缓存...")
        records = load_all_cached(config.CACHE_FILE)
    else:
        pmids   = search_pmids(config.SEARCH_QUERY, config.MAX_RECORDS)
        records = fetch_abstracts(pmids, config.CACHE_FILE, config.BATCH_SIZE)

    if not records:
        log.error("没有抓取到任何记录，请检查网络或搜索词")
        sys.exit(1)

    log.info(f"总记录数: {len(records)}")

    # ── Step 2: 关键词提取 ────────────────────────────────────────────────────
    method_counters = run_extraction(records)

    # ── Step 3: 分析 ─────────────────────────────────────────────────────────
    # Borda Count 融合多方法排名
    borda = borda_merge(method_counters, top_n=500)

    # 打印到控制台
    print_summary(borda, method_counters, n=config.TOP_K_DISPLAY)

    # 导出 CSV
    if config.CSV_EXPORT:
        export_csv(method_counters, borda, config.OUTPUT_DIR)

    # ── Step 4: 可视化 ────────────────────────────────────────────────────────
    if config.BARPLOT_ENABLE:
        plot_all_barplots(method_counters, borda, config.OUTPUT_DIR)

    if config.WORDCLOUD_ENABLE:
        plot_wordcloud(borda, title="Functional Keywords — Borda Merged",
                       filename="wordcloud_borda.png", output_dir=config.OUTPUT_DIR)
        for method, counter in method_counters.items():
            plot_wordcloud(counter, title=f"Functional Keywords — {method.upper()}",
                           filename=f"wordcloud_{method}.png", output_dir=config.OUTPUT_DIR)

    # ── Step 5: 年度趋势（可选）──────────────────────────────────────────────
    if not args.no_trend:
        log.info("开始年度趋势分析...")
        # 取 Borda Top 20 做趋势图
        top_kws = [kw for kw, _ in borda.most_common(20)]
        # 趋势分析优先用 author_kw，其次 mesh
        trend_method = "author_kw" if "author_kw" in method_counters else list(method_counters.keys())[0]
        trend = yearly_trend(records, top_kws, method=trend_method)
        plot_trend(trend,
                   title=f"Top 20 Keywords Yearly Trend ({trend_method})",
                   filename="trend.png",
                   output_dir=config.OUTPUT_DIR)

    log.info(f"\n✅ 全部完成！输出文件在 ./{config.OUTPUT_DIR}/ 目录")


if __name__ == "__main__":
    main()
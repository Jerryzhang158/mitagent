"""
extract.py — 关键词提取模块
支持三种模式：author_kw / tfidf / yake / combined
"""

import logging
import re
from collections import Counter
from typing import Optional

import config

log = logging.getLogger(__name__)

# ── 通用清洗 ───────────────────────────────────────────────────────────────────

def _normalize(kw: str) -> str:
    """小写 + 去标点 + 合并空格"""
    kw = kw.lower().strip()
    kw = re.sub(r"[^a-z0-9\s\-]", "", kw)
    kw = re.sub(r"\s+", " ", kw)
    return kw.strip()

def _is_valid(kw: str) -> bool:
    """过滤掉太短、纯数字、stopword 的关键词"""
    if len(kw) < 3:
        return False
    if kw.replace("-", "").isnumeric():
        return False
    tokens = kw.split()
    # 如果所有词都是 stopword，丢弃
    if all(t in config.STOPWORDS_EXTRA for t in tokens):
        return False
    return True


# ── 方法 1：直接统计 Author Keywords ──────────────────────────────────────────

def extract_author_keywords(records: list[dict]) -> Counter:
    counter = Counter()
    for r in records:
        for kw in r.get("author_keywords", []):
            kw = _normalize(kw)
            if _is_valid(kw):
                counter[kw] += 1
    log.info(f"[author_kw] 提取到 {len(counter)} 个唯一关键词")
    return counter


# ── 方法 2：直接统计 MeSH Terms ───────────────────────────────────────────────

def extract_mesh_keywords(records: list[dict]) -> Counter:
    counter = Counter()
    for r in records:
        for kw in r.get("mesh_terms", []):
            # MeSH 可能含斜杠（subheading），只保留主词
            main = kw.split("/")[0]
            main = _normalize(main)
            if _is_valid(main):
                counter[main] += 1
    log.info(f"[mesh] 提取到 {len(counter)} 个唯一 MeSH 词")
    return counter


# ── 方法 3：TF-IDF（摘要全文）────────────────────────────────────────────────

def extract_tfidf(records: list[dict]) -> Counter:
    try:
        from sklearn.feature_extraction.text import TfidfVectorizer
        import numpy as np
    except ImportError:
        log.error("sklearn 未安装，请 pip install scikit-learn")
        return Counter()

    # 合并 title + abstract 作为文档
    docs = []
    for r in records:
        text = f"{r.get('title', '')} {r.get('abstract', '')}"
        docs.append(text.lower())

    docs = [d for d in docs if len(d.strip()) > 20]
    if not docs:
        return Counter()

    log.info(f"[tfidf] 对 {len(docs)} 篇文档建模...")
    vectorizer = TfidfVectorizer(
        max_features=config.TFIDF_MAX_FEATURES,
        ngram_range=config.TFIDF_NGRAM_RANGE,
        stop_words="english",
        min_df=config.MIN_FREQ,
        sublinear_tf=True,
    )
    tfidf_matrix = vectorizer.fit_transform(docs)
    feature_names = vectorizer.get_feature_names_out()

    # 对每个词，取所有文档的 TF-IDF 之和作为全局重要性
    global_scores = np.asarray(tfidf_matrix.sum(axis=0)).flatten()
    ranked_idx = global_scores.argsort()[::-1][: config.TFIDF_TOP_N * 3]

    counter = Counter()
    for idx in ranked_idx:
        kw = feature_names[idx]
        kw = _normalize(kw)
        if _is_valid(kw):
            # 用 document frequency 作为 counter 值（更直观）
            doc_freq = int((tfidf_matrix[:, idx] > 0).sum())
            counter[kw] = doc_freq

    log.info(f"[tfidf] 提取到 {len(counter)} 个关键词")
    return counter


# ── 方法 4：YAKE（无监督，逐篇提取）─────────────────────────────────────────

def extract_yake(records: list[dict]) -> Counter:
    try:
        import yake
    except ImportError:
        log.error("yake 未安装，请 pip install yake")
        return Counter()

    extractor = yake.KeywordExtractor(
        lan=config.YAKE_LANG,
        n=config.YAKE_MAX_NGRAM_SIZE,
        dedupLim=config.YAKE_DEDUP_THRESH,
        top=config.YAKE_TOP_N_PER_DOC,
        features=None,
    )

    counter = Counter()
    total = len(records)
    for i, r in enumerate(records):
        if i % 500 == 0:
            log.info(f"[yake] {i}/{total} 篇...")
        text = f"{r.get('title', '')} {r.get('abstract', '')}"
        if len(text.strip()) < 20:
            continue
        try:
            keywords = extractor.extract_keywords(text)
            for kw, score in keywords:
                kw = _normalize(kw)
                if _is_valid(kw):
                    counter[kw] += 1
        except Exception:
            pass

    log.info(f"[yake] 提取到 {len(counter)} 个唯一关键词")
    return counter


# ── Combined：多来源融合 ───────────────────────────────────────────────────────

def extract_combined(records: list[dict]) -> dict[str, Counter]:
    """
    返回各方法的 Counter，供后续分别分析或合并排名
    """
    results = {}
    log.info("=== 开始 combined 提取 ===")

    results["author_kw"] = extract_author_keywords(records)
    results["mesh"]      = extract_mesh_keywords(records)
    results["tfidf"]     = extract_tfidf(records)
    results["yake"]      = extract_yake(records)

    # 过滤低频词
    for method, counter in results.items():
        before = len(counter)
        results[method] = Counter({k: v for k, v in counter.items() if v >= config.MIN_FREQ})
        log.info(f"[{method}] 过滤后: {before} → {len(results[method])} 个关键词")

    return results


# ── 统一入口 ──────────────────────────────────────────────────────────────────

def run_extraction(records: list[dict]) -> dict[str, Counter]:
    method = config.EXTRACTION_METHOD
    log.info(f"提取方法: {method}")

    if method == "author_kw":
        return {"author_kw": extract_author_keywords(records)}
    elif method == "mesh":
        return {"mesh": extract_mesh_keywords(records)}
    elif method == "tfidf":
        return {"tfidf": extract_tfidf(records)}
    elif method == "yake":
        return {"yake": extract_yake(records)}
    elif method == "combined":
        return extract_combined(records)
    else:
        raise ValueError(f"未知方法: {method}")

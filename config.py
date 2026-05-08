"""
Configuration for PubMed Functional Keyword Analyzer
"""

# ── PubMed / NCBI ──────────────────────────────────────────────────────────────
ENTREZ_EMAIL = "your_email@example.com"   # NCBI 要求填邮箱，不验证但必须填
ENTREZ_API_KEY = ""                        # 可选：填了速率从 3 req/s → 10 req/s
                                           # 申请：https://www.ncbi.nlm.nih.gov/account/

# ── 搜索参数 ────────────────────────────────────────────────────────────────────
SEARCH_QUERY = (
    '("miRNA"[Title/Abstract] OR "microRNA"[Title/Abstract] OR '
    '"miR-"[Title/Abstract] OR "MicroRNAs"[MeSH Terms]) '
    'AND ("2020/01/01"[PDAT] : "2025/12/31"[PDAT]) '
    'AND ("journal article"[PT]) '
    'AND (hasabstract[text]) '
    'AND ("English"[Language])'
)

# ── 可选的更聚焦搜索式（按需替换上面的 SEARCH_QUERY）────────────────────────
# 专注 miRNA 靶点预测 / 验证
# SEARCH_QUERY = (
#     '("miRNA target"[Title/Abstract] OR "microRNA target"[Title/Abstract] OR '
#     '"miRNA-target interaction"[Title/Abstract]) '
#     'AND ("2020/01/01"[PDAT] : "2024/12/31"[PDAT]) '
#     'AND (hasabstract[text]) AND ("English"[Language])'
# )

# 专注 miRNA + 癌症
# SEARCH_QUERY = (
#     '("miRNA"[Title/Abstract] OR "microRNA"[Title/Abstract]) '
#     'AND "Neoplasms"[MeSH Terms] '
#     'AND ("2020/01/01"[PDAT] : "2024/12/31"[PDAT]) '
#     'AND (hasabstract[text]) AND ("English"[Language])'
# )

# 专注 miRNA + 实验验证（luciferase / qPCR）
# SEARCH_QUERY = (
#     '("miRNA"[Title/Abstract] OR "microRNA"[Title/Abstract]) '
#     'AND ("luciferase"[Title/Abstract] OR "qPCR"[Title/Abstract] '
#     '     OR "western blot"[Title/Abstract]) '
#     'AND ("2020/01/01"[PDAT] : "2024/12/31"[PDAT]) '
#     'AND (hasabstract[text]) AND ("English"[Language])'
# )

MAX_RECORDS = 10000      # 最多抓取论文数（建议先用 2000 测试，确认无误再跑大量）
BATCH_SIZE  = 200        # 每次请求的批大小（不超过 500）
CACHE_FILE  = "data/abstracts_cache.jsonl"   # 断点续抓缓存

# ── NLP 参数 ────────────────────────────────────────────────────────────────────
# 关键词提取方法：'author_kw' | 'tfidf' | 'yake' | 'combined'
EXTRACTION_METHOD = "combined"

# TF-IDF
TFIDF_MAX_FEATURES  = 5000
TFIDF_NGRAM_RANGE   = (1, 3)   # unigram + bigram + trigram
TFIDF_TOP_N         = 200      # 保留 top N 关键词

# YAKE（单篇提取，再汇总）
YAKE_LANG           = "en"
YAKE_MAX_NGRAM_SIZE = 3
YAKE_DEDUP_THRESH   = 0.9
YAKE_TOP_N_PER_DOC  = 10      # 每篇提 10 个

# 过滤设置
MIN_FREQ            = 5        # 出现少于 N 次的词丢弃
STOPWORDS_EXTRA     = {        # 生物领域无意义高频词
    # 通用方法学词
    "study", "result", "results", "analysis", "method", "methods", "data",
    "approach", "paper", "show", "suggest", "identify", "investigate",
    "significant", "significantly", "increase", "decrease", "high", "low",
    "using", "based", "via", "associated", "related", "found", "performed",
    "observed", "indicated", "demonstrated", "revealed", "evaluated",
    # miRNA 本身不是 functional keyword，过滤掉
    "mirna", "microrna", "mir", "mirnas", "micrornas",
    # 过于宽泛的生物词
    "expression", "gene", "genes", "protein", "proteins", "cell", "cells",
    "role", "function", "effect", "effects", "pathway", "pathways",
    "model", "level", "levels", "activity", "regulation", "target", "targets",
    "binding", "interaction", "interactions", "mechanism", "mechanisms",
    "patient", "patients", "group", "control", "sample", "samples",
    "tissue", "tissues", "type", "types", "line", "lines",
}

# ── 输出参数 ────────────────────────────────────────────────────────────────────
OUTPUT_DIR       = "output"
TOP_K_DISPLAY    = 50          # 最终展示 Top K 关键词
WORDCLOUD_ENABLE = True
BARPLOT_ENABLE   = True
CSV_EXPORT       = True
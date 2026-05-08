"""
fetch.py — PubMed 摘要抓取模块
支持断点续抓、自动重试、速率控制
"""

import json
import os
import time
import logging
from pathlib import Path
from typing import Generator

from Bio import Entrez, Medline

import config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

Entrez.email   = config.ENTREZ_EMAIL
Entrez.api_key = config.ENTREZ_API_KEY or None

# 有 API key → 10 req/s；否则 3 req/s（保守用 2.5）
REQUEST_INTERVAL = 0.11 if config.ENTREZ_API_KEY else 0.4


def _retry(fn, retries=3, wait=5):
    """简单重试包装器"""
    for attempt in range(retries):
        try:
            return fn()
        except Exception as e:
            if attempt == retries - 1:
                raise
            log.warning(f"请求失败（{e}），{wait}s 后重试 ({attempt+1}/{retries})...")
            time.sleep(wait)


def search_pmids(query: str, max_records: int) -> list[str]:
    """用 esearch 获取 PMID 列表"""
    log.info(f"搜索 PubMed: {query[:80]}...")
    handle = _retry(lambda: Entrez.esearch(
        db="pubmed", term=query,
        retmax=max_records, usehistory="y"
    ))
    record = Entrez.read(handle)
    handle.close()
    pmids = record["IdList"]
    total = int(record["Count"])
    log.info(f"命中 {total} 篇，将抓取前 {len(pmids)} 篇")
    return pmids


def _load_cache(cache_file: str) -> set[str]:
    """读取已缓存的 PMID"""
    cached = set()
    path = Path(cache_file)
    if path.exists():
        with open(path, encoding="utf-8") as f:
            for line in f:
                try:
                    obj = json.loads(line)
                    cached.add(obj["pmid"])
                except Exception:
                    pass
    log.info(f"缓存已有 {len(cached)} 篇")
    return cached


def _append_cache(records: list[dict], cache_file: str):
    Path(cache_file).parent.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def fetch_abstracts(
    pmids: list[str],
    cache_file: str,
    batch_size: int = 200,
) -> list[dict]:
    """
    批量抓取摘要，返回 list of dict：
    {pmid, title, abstract, author_keywords, mesh_terms, year}
    """
    cached_pmids = _load_cache(cache_file)
    todo = [p for p in pmids if p not in cached_pmids]
    log.info(f"需要抓取 {len(todo)} 篇（跳过已缓存 {len(cached_pmids)} 篇）")

    for start in range(0, len(todo), batch_size):
        batch = todo[start: start + batch_size]
        log.info(f"  抓取 {start+1}–{start+len(batch)} / {len(todo)}")
        ids_str = ",".join(batch)

        raw = _retry(lambda: Entrez.efetch(
            db="pubmed", id=ids_str,
            rettype="medline", retmode="text"
        ))
        records_raw = list(Medline.parse(raw))
        raw.close()

        parsed = []
        for r in records_raw:
            # 提取 author keywords（OT 字段）
            author_kws = list(r.get("OT", []))
            # 提取 MeSH terms（MH 字段）
            mesh = [m.lstrip("*") for m in r.get("MH", [])]
            # 发表年份
            dp = r.get("DP", "")
            year = dp[:4] if dp else ""

            parsed.append({
                "pmid":            r.get("PMID", ""),
                "title":           r.get("TI",   ""),
                "abstract":        r.get("AB",   ""),
                "author_keywords": author_kws,
                "mesh_terms":      mesh,
                "year":            year,
            })

        _append_cache(parsed, cache_file)
        time.sleep(REQUEST_INTERVAL)

    # 读取完整缓存并返回
    return load_all_cached(cache_file)


def load_all_cached(cache_file: str) -> list[dict]:
    records = []
    with open(cache_file, encoding="utf-8") as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except Exception:
                pass
    log.info(f"共加载 {len(records)} 篇缓存记录")
    return records
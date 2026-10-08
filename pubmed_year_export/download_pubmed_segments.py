#!/usr/bin/env python3
"""Download a large PubMed query in monthly publication-date segments.

PubMed ESearch/EFetch cannot page beyond the first 10,000 PubMed results from a
single query. Monthly date segments keep every History result set below that
limit. State is persisted per month, and output filenames are unique across
segments. Existing batches in the same export directory are preserved and are
deduplicated by PMID when build_local_index.py creates the final artifacts.
"""

from __future__ import annotations

import argparse
import calendar
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from pubmed_year_export.download_pubmed_year import (
    common_params,
    download_batch,
    parse_xml_gz,
    request_bytes,
    save_state,
    utc_now,
    write_gzip,
    write_jsonl_gz,
)


def month_segments(start_year: int, end_year: int):
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            last_day = calendar.monthrange(year, month)[1]
            key = f"{year:04d}-{month:02d}"
            date_clause = (
                f"{year:04d}/{month:02d}/01:"
                f"{year:04d}/{month:02d}/{last_day:02d}[dp]"
            )
            yield key, date_clause


def create_history(args: argparse.Namespace, query: str) -> tuple[int, str, str]:
    params: dict[str, Any] = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 0,
        "usehistory": "y",
        **common_params(args),
    }
    result = json.loads(request_bytes("esearch.fcgi", params))["esearchresult"]
    return int(result["count"]), result["querykey"], result["webenv"]


def load_segmented_state(path: Path, args: argparse.Namespace) -> dict[str, Any]:
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        expected = (args.start_year, args.end_year, args.base_query)
        actual = (
            state.get("start_year"),
            state.get("end_year"),
            state.get("base_query"),
        )
        if actual != expected:
            raise ValueError(f"Segment state belongs to different settings: {actual}")
        return state
    return {
        "start_year": args.start_year,
        "end_year": args.end_year,
        "base_query": args.base_query,
        "segments": {},
        "created_at": utc_now(),
    }


def save_segmented_state(path: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(state, ensure_ascii=False, sort_keys=True, indent=2),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> int:
    root = Path(args.output_dir).resolve()
    raw_dir = root / "raw"
    parsed_dir = root / "parsed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    parsed_dir.mkdir(parents=True, exist_ok=True)
    state_path = root / "segmented_state.json"
    state = load_segmented_state(state_path, args)
    delay = args.delay if args.delay is not None else (0.11 if args.api_key else 0.4)

    for key, date_clause in month_segments(args.start_year, args.end_year):
        query = f"({date_clause}) AND ({args.base_query})"
        segment = state["segments"].setdefault(
            key,
            {"query": query, "next_retstart": 0, "records_written": 0},
        )
        count, query_key, webenv = create_history(args, query)
        segment["pubmed_match_count"] = count
        if count >= 10000:
            raise RuntimeError(
                f"Segment {key} has {count:,} records; split it into smaller dates"
            )
        start = int(segment["next_retstart"])
        if start >= count:
            print(f"[{key}] already complete: {start:,}/{count:,}")
            continue
        print(f"[{key}] matches={count:,}, resume={start:,}")

        while start < count:
            size = min(args.batch_size, count - start)
            stem = f"segment_{key.replace('-', '_')}_batch_{start:07d}_{start + size - 1:07d}"
            raw_path = raw_dir / f"{stem}.xml.gz"
            parsed_path = parsed_dir / f"{stem}.jsonl.gz"
            if not raw_path.exists():
                data = download_batch(args, query_key, webenv, start, size)
                write_gzip(raw_path, data)
            records = parse_xml_gz(raw_path)
            if not records:
                raise RuntimeError(f"No records parsed for {key} retstart={start}")
            for record in records:
                record["query_segment"] = key
                record["fetched_at"] = utc_now()
            write_jsonl_gz(parsed_path, records)

            start += size
            segment["next_retstart"] = start
            segment["records_written"] = int(segment["records_written"]) + len(records)
            segment["batches_completed"] = int(segment.get("batches_completed", 0)) + 1
            segment["last_raw_file"] = str(raw_path.relative_to(root))
            segment["last_parsed_file"] = str(parsed_path.relative_to(root))
            save_segmented_state(state_path, state)
            print(f"[{key}] {start:,}/{count:,}, parsed={len(records):,}")
            if start < count:
                time.sleep(delay)

    total_matches = sum(
        int(segment.get("pubmed_match_count", 0))
        for segment in state["segments"].values()
    )
    total_written = sum(
        int(segment.get("records_written", 0))
        for segment in state["segments"].values()
    )
    state["complete"] = True
    state["total_segment_matches"] = total_matches
    state["total_records_written"] = total_written
    save_segmented_state(state_path, state)
    print(f"Complete: segment matches={total_matches:,}, parsed={total_written:,}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-year", type=int, default=2020)
    parser.add_argument("--end-year", type=int, default=2026)
    parser.add_argument("--base-query", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    parser.add_argument("--tool", default="mti_agent_local_pubmed_segments")
    parser.add_argument("--delay", type=float, default=None)
    args = parser.parse_args()
    if args.end_year < args.start_year:
        parser.error("--end-year must be >= --start-year")
    if not 1 <= args.batch_size <= 1000:
        parser.error("--batch-size must be between 1 and 1000")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())

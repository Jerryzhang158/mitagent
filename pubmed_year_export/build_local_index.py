#!/usr/bin/env python3
"""Build reproducible JSONL and a compact SQLite FTS5 index from parsed batches."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import sqlite3
from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO


SCHEMA_VERSION = 1
PARSER_FORMAT_VERSION = 1
YEAR_RE = re.compile(r"(?<!\d)(?:18|19|20|21)\d{2}(?!\d)")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def deterministic_gzip_writer(path: Path) -> tuple[TextIO, Any, Any]:
    """Return a UTF-8 writer plus resources that must be closed in order."""
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    return text, compressed, raw


def iter_records(parsed_dir: Path) -> Iterator[dict[str, Any]]:
    for path in sorted(parsed_dir.glob("*.jsonl.gz")):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"Invalid JSON in {path}:{line_number}: {exc}") from exc


def flatten_mesh(record: dict[str, Any]) -> str:
    values: list[str] = []
    for heading in record.get("mesh_headings", []):
        descriptor = heading.get("descriptor", "")
        if descriptor:
            values.append(descriptor)
        values.extend(
            qualifier.get("name", "")
            for qualifier in heading.get("qualifiers", [])
            if qualifier.get("name")
        )
    return " ".join(values)


def flatten_keywords(record: dict[str, Any]) -> str:
    return " ".join(
        keyword.get("text", "")
        for keyword in record.get("keywords", [])
        if keyword.get("text")
    )


def publication_year(record: dict[str, Any]) -> int | None:
    publication_date = record.get("journal", {}).get("publication_date", {})
    values = [
        publication_date.get("year", ""),
        publication_date.get("medline_date", ""),
    ]
    values.extend(item.get("year", "") for item in record.get("article_dates", []))
    for value in values:
        match = YEAR_RE.search(str(value))
        if match:
            return int(match.group(0))
    return None


def init_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE articles (
            pmid TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            abstract TEXT NOT NULL,
            publication_year INTEGER,
            publication_date_json TEXT NOT NULL,
            journal TEXT NOT NULL,
            doi TEXT,
            pmcid TEXT,
            languages_json TEXT NOT NULL,
            abstract_sections_json TEXT NOT NULL,
            other_abstracts_json TEXT NOT NULL,
            authors_json TEXT NOT NULL,
            publication_types_json TEXT NOT NULL,
            mesh_json TEXT NOT NULL,
            keywords_json TEXT NOT NULL,
            date_revised_json TEXT NOT NULL,
            source_file TEXT NOT NULL
        );

        CREATE INDEX idx_articles_year ON articles(publication_year);
        CREATE INDEX idx_articles_doi ON articles(doi);
        CREATE INDEX idx_articles_pmcid ON articles(pmcid);

        CREATE VIRTUAL TABLE articles_fts USING fts5(
            pmid UNINDEXED,
            title,
            abstract,
            mesh_terms,
            keywords,
            tokenize='unicode61 remove_diacritics 2'
        );
        """
    )
    return connection


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def article_row(record: dict[str, Any]) -> tuple[Any, ...]:
    identifiers = record.get("identifiers", {})
    journal = record.get("journal", {})
    return (
        record["pmid"],
        record.get("title", ""),
        record.get("abstract", ""),
        publication_year(record),
        compact_json(journal.get("publication_date", {})),
        journal.get("title", ""),
        identifiers.get("doi", ""),
        identifiers.get("pmc", ""),
        compact_json(record.get("languages", [])),
        compact_json(record.get("abstract_sections", [])),
        compact_json(record.get("other_abstracts", [])),
        compact_json(record.get("authors", [])),
        compact_json(record.get("publication_types", [])),
        compact_json(record.get("mesh_headings", [])),
        compact_json(record.get("keywords", [])),
        compact_json(record.get("date_revised", {})),
        record.get("source_file", ""),
    )


def fts_row(record: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        record["pmid"],
        record.get("title", ""),
        record.get("abstract", ""),
        flatten_mesh(record),
        flatten_keywords(record),
    )


def build(args: argparse.Namespace) -> int:
    root = args.export_dir.resolve()
    parsed_dir = root / "parsed"
    state_path = root / "state.json"
    segmented_state_path = root / "segmented_state.json"
    baseline_state_path = root / "baseline_state.json"
    if not parsed_dir.exists() or not (
        state_path.exists()
        or segmented_state_path.exists()
        or baseline_state_path.exists()
    ):
        raise SystemExit(f"Not a completed export directory: {root}")
    baseline: dict[str, Any] | None = None
    if baseline_state_path.exists():
        baseline = json.loads(baseline_state_path.read_text(encoding="utf-8"))
        if not baseline.get("complete"):
            raise SystemExit(
                "Baseline download/filter is incomplete: "
                f"{baseline.get('completed_file_count', 0)}/"
                f"{baseline.get('discovered_file_count', 0)}"
            )
        files = baseline.get("files", {})
        incomplete = [name for name, entry in files.items() if not entry.get("completed")]
        if incomplete:
            raise SystemExit(
                f"Baseline state contains incomplete source files: {incomplete[:3]}"
            )
        state = {
            "query": baseline.get("filter_description", ""),
            "pubmed_match_count": baseline.get("records_scanned", 0),
        }
    elif segmented_state_path.exists():
        segmented = json.loads(segmented_state_path.read_text(encoding="utf-8"))
        if not segmented.get("complete"):
            raise SystemExit("Segmented download is incomplete")
        state = {
            "query": (
                f"{segmented['start_year']}:{segmented['end_year']}[dp] AND "
                f"({segmented['base_query']})"
            ),
            "pubmed_match_count": segmented.get("total_segment_matches", 0),
        }
    else:
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("next_retstart", 0) < state.get("target_records", 0):
            raise SystemExit(
                f"Download incomplete: {state.get('next_retstart', 0)}/"
                f"{state.get('target_records', 0)}"
            )

    database_path = root / "pubmed_mirna.sqlite"
    database_tmp = root / "pubmed_mirna.sqlite.tmp"
    jsonl_path = root / "articles.jsonl.gz"
    jsonl_tmp = root / "articles.jsonl.gz.tmp"
    pmids_path = root / "pmids.txt.gz"
    pmids_tmp = root / "pmids.txt.gz.tmp"
    for temporary in (database_tmp, jsonl_tmp, pmids_tmp):
        if temporary.exists():
            temporary.unlink()

    connection = init_database(database_tmp)
    json_writer, json_compressed, json_raw = deterministic_gzip_writer(jsonl_tmp)
    seen: set[str] = set()
    count = 0
    batch_articles: list[tuple[Any, ...]] = []
    batch_fts: list[tuple[str, str, str, str, str]] = []

    try:
        for record in iter_records(parsed_dir):
            pmid = record.get("pmid", "")
            if not pmid or pmid in seen:
                continue
            seen.add(pmid)
            json_writer.write(compact_json(record))
            json_writer.write("\n")
            batch_articles.append(article_row(record))
            batch_fts.append(fts_row(record))
            count += 1

            if len(batch_articles) >= 1000:
                with connection:
                    connection.executemany(
                        "INSERT INTO articles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        batch_articles,
                    )
                    connection.executemany(
                        "INSERT INTO articles_fts VALUES (?,?,?,?,?)",
                        batch_fts,
                    )
                batch_articles.clear()
                batch_fts.clear()
                if count % 10000 == 0:
                    print(f"Indexed {count:,} records")

        if batch_articles:
            with connection:
                connection.executemany(
                    "INSERT INTO articles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch_articles,
                )
                connection.executemany(
                    "INSERT INTO articles_fts VALUES (?,?,?,?,?)",
                    batch_fts,
                )

        metadata = {
            "schema_version": str(SCHEMA_VERSION),
            "parser_format_version": str(PARSER_FORMAT_VERSION),
            "query": state.get("query", ""),
            "record_count": str(count),
        }
        with connection:
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES (?, ?)", metadata.items()
            )
            connection.execute("INSERT INTO articles_fts(articles_fts) VALUES ('optimize')")
        connection.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    finally:
        json_writer.flush()
        json_writer.close()
        json_compressed.close()
        json_raw.close()
        connection.close()

    pmid_writer, pmid_compressed, pmid_raw = deterministic_gzip_writer(pmids_tmp)
    try:
        for pmid in sorted(seen, key=int):
            pmid_writer.write(pmid)
            pmid_writer.write("\n")
    finally:
        pmid_writer.flush()
        pmid_writer.close()
        pmid_compressed.close()
        pmid_raw.close()

    os.replace(database_tmp, database_path)
    os.replace(jsonl_tmp, jsonl_path)
    os.replace(pmids_tmp, pmids_path)

    raw_files = sorted((root / "raw").glob("*.xml.gz"))
    parsed_files = sorted(parsed_dir.glob("*.jsonl.gz"))
    manifest = {
        "corpus_version": args.corpus_version,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "query": state.get("query", ""),
        "pubmed_match_count": state.get("pubmed_match_count"),
        "record_count": count,
        "schema_version": SCHEMA_VERSION,
        "parser_format_version": PARSER_FORMAT_VERSION,
        "sqlite_version": sqlite3.sqlite_version,
        "artifacts": {
            "articles.jsonl.gz": {
                "bytes": jsonl_path.stat().st_size,
                "sha256": sha256_file(jsonl_path),
            },
            "pmids.txt.gz": {
                "bytes": pmids_path.stat().st_size,
                "sha256": sha256_file(pmids_path),
            },
            "pubmed_mirna.sqlite": {
                "bytes": database_path.stat().st_size,
                "sha256": sha256_file(database_path),
            },
        },
        "source_batches": {
            "raw_xml_files": len(raw_files),
            "raw_xml_bytes": sum(path.stat().st_size for path in raw_files),
            "parsed_jsonl_files": len(parsed_files),
            "parsed_jsonl_bytes": sum(path.stat().st_size for path in parsed_files),
        },
    }
    if baseline is not None:
        manifest["baseline"] = {
            "year": baseline.get("baseline_year"),
            "base_url": baseline.get("base_url", ""),
            "filter_version": baseline.get("filter_version", ""),
            "filter_description": baseline.get("filter_description", ""),
            "start_year": baseline.get("start_year"),
            "end_year": baseline.get("end_year"),
            "source_files": [
                {
                    "source_file": name,
                    "official_md5": entry.get("official_md5", ""),
                    "records_scanned": entry.get("records_scanned", 0),
                    "records_selected": entry.get("records_selected", 0),
                    "filtered_sha256": entry.get("filtered_sha256", ""),
                }
                for name, entry in sorted(baseline.get("files", {}).items())
            ],
        }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2),
        encoding="utf-8",
    )
    print(f"Built {count:,} records")
    print(f"JSONL: {jsonl_path}")
    print(f"SQLite: {database_path}")
    print(f"Manifest: {manifest_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_dir", type=Path)
    parser.add_argument(
        "--corpus-version", default="pubmed-mirna-2020-2026-v1"
    )
    return build(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

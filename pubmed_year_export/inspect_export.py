#!/usr/bin/env python3
"""Inspect a PubMed year export and write a compact format report."""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
from collections import Counter
from pathlib import Path


def uncompressed_size(path: Path) -> int:
    total = 0
    with gzip.open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            total += len(chunk)
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_dir", type=Path)
    parser.add_argument("--sample-pmid", default="")
    args = parser.parse_args()

    root = args.export_dir.resolve()
    connection = sqlite3.connect(root / "pubmed_abstracts.sqlite")
    rows = connection.execute("SELECT record_json FROM articles").fetchall()
    records = [json.loads(row[0]) for row in rows]
    if not records:
        raise SystemExit("No records found")

    languages = Counter(
        language for record in records for language in record.get("languages", [])
    )
    publication_years = Counter(
        record.get("journal", {}).get("publication_date", {}).get("year", "unknown")
        for record in records
    )
    abstract_lengths = [len(record.get("abstract", "")) for record in records]
    raw_files = sorted((root / "raw").glob("*.xml.gz"))
    parsed_files = sorted((root / "parsed").glob("*.jsonl.gz"))

    selected = records[0]
    if args.sample_pmid:
        selected = next(
            (record for record in records if record["pmid"] == args.sample_pmid),
            selected,
        )

    report = {
        "record_count": len(records),
        "records_with_primary_abstract": sum(bool(x) for x in abstract_lengths),
        "records_with_multiple_abstract_sections": sum(
            len(record.get("abstract_sections", [])) > 1 for record in records
        ),
        "records_with_labeled_abstract_sections": sum(
            any(section.get("label") for section in record.get("abstract_sections", []))
            for record in records
        ),
        "records_with_other_abstracts": sum(
            bool(record.get("other_abstracts")) for record in records
        ),
        "records_with_doi": sum(
            bool(record.get("identifiers", {}).get("doi")) for record in records
        ),
        "records_with_pmcid": sum(
            bool(record.get("identifiers", {}).get("pmc")) for record in records
        ),
        "abstract_characters": {
            "minimum": min(abstract_lengths),
            "maximum": max(abstract_lengths),
            "mean": round(sum(abstract_lengths) / len(abstract_lengths), 1),
        },
        "languages": dict(languages.most_common()),
        "xml_publication_years": dict(publication_years.most_common()),
        "storage_bytes": {
            "raw_xml_gzip": sum(path.stat().st_size for path in raw_files),
            "raw_xml_uncompressed": sum(uncompressed_size(path) for path in raw_files),
            "jsonl_gzip": sum(path.stat().st_size for path in parsed_files),
            "jsonl_uncompressed": sum(uncompressed_size(path) for path in parsed_files),
            "sqlite": (root / "pubmed_abstracts.sqlite").stat().st_size,
        },
        "sample_record": selected,
    }
    report_path = root / "inspection_report.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in report.items() if key != "sample_record"}, ensure_ascii=False, indent=2))
    print(f"Sample PMID: {selected['pmid']}")
    print(f"Report: {report_path}")
    connection.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

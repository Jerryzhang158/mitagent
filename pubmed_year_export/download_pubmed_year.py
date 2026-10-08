#!/usr/bin/env python3
"""Download and normalize all PubMed records with abstracts for one year.

The downloader uses the official NCBI E-utilities History server. Every batch is
stored twice:

1. raw/*.xml.gz       exact XML returned by NCBI, compressed without changes
2. parsed/*.jsonl.gz  one normalized JSON object per PubMed record

The same records are also upserted into a queryable SQLite database. A state
file makes interrupted runs resumable; a fresh ESearch History session is
created on every invocation, so an expired WebEnv does not invalidate progress.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "mti-agent-pubmed-year-export/1.0"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean_text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def node_text(parent: ET.Element | None, path: str) -> str:
    return clean_text(parent.find(path)) if parent is not None else ""


def date_object(node: ET.Element | None) -> dict[str, str]:
    if node is None:
        return {}
    result = {
        "year": node_text(node, "Year"),
        "month": node_text(node, "Month"),
        "day": node_text(node, "Day"),
        "season": node_text(node, "Season"),
        "medline_date": node_text(node, "MedlineDate"),
    }
    return {key: value for key, value in result.items() if value}


def abstract_sections(container: ET.Element | None) -> list[dict[str, str]]:
    if container is None:
        return []
    result: list[dict[str, str]] = []
    for node in container.findall("AbstractText"):
        text = clean_text(node)
        if not text:
            continue
        section = {"text": text}
        if node.get("Label"):
            section["label"] = node.get("Label", "")
        if node.get("NlmCategory"):
            section["nlm_category"] = node.get("NlmCategory", "")
        result.append(section)
    return result


def flatten_sections(sections: Iterable[dict[str, str]]) -> str:
    chunks = []
    for section in sections:
        text = section["text"]
        label = section.get("label", "").strip()
        chunks.append(f"{label}: {text}" if label else text)
    return "\n".join(chunks)


def parse_authors(container: ET.Element | None) -> list[dict[str, Any]]:
    if container is None:
        return []
    authors = []
    for node in container.findall("Author"):
        author: dict[str, Any] = {}
        for xml_name, output_name in (
            ("LastName", "last_name"),
            ("ForeName", "fore_name"),
            ("Initials", "initials"),
            ("Suffix", "suffix"),
            ("CollectiveName", "collective_name"),
        ):
            value = node_text(node, xml_name)
            if value:
                author[output_name] = value
        affiliations = [
            clean_text(item)
            for item in node.findall("AffiliationInfo/Affiliation")
            if clean_text(item)
        ]
        if affiliations:
            author["affiliations"] = affiliations
        if node.get("ValidYN"):
            author["valid"] = node.get("ValidYN") == "Y"
        if author:
            authors.append(author)
    return authors


def parse_mesh(citation: ET.Element | None) -> list[dict[str, Any]]:
    if citation is None:
        return []
    headings = []
    for node in citation.findall("MeshHeadingList/MeshHeading"):
        descriptor = node.find("DescriptorName")
        if descriptor is None:
            continue
        qualifiers = []
        for qualifier in node.findall("QualifierName"):
            qualifiers.append(
                {
                    "name": clean_text(qualifier),
                    "ui": qualifier.get("UI", ""),
                    "major_topic": qualifier.get("MajorTopicYN") == "Y",
                }
            )
        headings.append(
            {
                "descriptor": clean_text(descriptor),
                "ui": descriptor.get("UI", ""),
                "major_topic": descriptor.get("MajorTopicYN") == "Y",
                "qualifiers": qualifiers,
            }
        )
    return headings


def parse_identifiers(record: ET.Element) -> dict[str, str]:
    identifiers: dict[str, str] = {}
    for node in record.findall(".//PubmedData/ArticleIdList/ArticleId"):
        value = clean_text(node)
        kind = node.get("IdType", "").lower()
        if value and kind:
            identifiers[kind] = value
    return identifiers


def parse_record(record: ET.Element, source_file: str) -> dict[str, Any] | None:
    citation = record.find("MedlineCitation")
    article = citation.find("Article") if citation is not None else None
    book_document = record.find("BookDocument")
    content = article if article is not None else book_document
    if content is None:
        return None

    pmid = node_text(citation, "PMID") or node_text(book_document, "PMID")
    if not pmid:
        return None

    primary_abstract = content.find("Abstract")
    sections = abstract_sections(primary_abstract)
    other_abstracts = []
    other_parent = citation if citation is not None else book_document
    if other_parent is not None:
        for other in other_parent.findall("OtherAbstract"):
            other_sections = abstract_sections(other)
            if not other_sections:
                continue
            item: dict[str, Any] = {
                "text": flatten_sections(other_sections),
                "sections": other_sections,
            }
            for attr, output_name in (
                ("Type", "type"),
                ("Language", "language"),
            ):
                if other.get(attr):
                    item[output_name] = other.get(attr)
            copyright_text = node_text(other, "CopyrightInformation")
            if copyright_text:
                item["copyright"] = copyright_text
            other_abstracts.append(item)

    journal = content.find("Journal")
    journal_issue = journal.find("JournalIssue") if journal is not None else None
    pub_date = journal_issue.find("PubDate") if journal_issue is not None else None
    identifiers = parse_identifiers(record)

    publication_types = []
    for node in content.findall("PublicationTypeList/PublicationType"):
        publication_types.append(
            {"name": clean_text(node), "ui": node.get("UI", "")}
        )

    keywords = []
    if citation is not None:
        for keyword_list in citation.findall("KeywordList"):
            owner = keyword_list.get("Owner", "")
            for node in keyword_list.findall("Keyword"):
                keywords.append(
                    {
                        "text": clean_text(node),
                        "owner": owner,
                        "major_topic": node.get("MajorTopicYN") == "Y",
                    }
                )

    languages = [
        clean_text(node) for node in content.findall("Language") if clean_text(node)
    ]
    article_dates = []
    for node in content.findall("ArticleDate"):
        value = date_object(node)
        if node.get("DateType"):
            value["date_type"] = node.get("DateType", "")
        if value:
            article_dates.append(value)

    result: dict[str, Any] = {
        "pmid": pmid,
        "title": node_text(content, "ArticleTitle")
        or node_text(content, "Book/BookTitle"),
        "vernacular_title": node_text(content, "VernacularTitle"),
        "abstract": flatten_sections(sections),
        "abstract_sections": sections,
        "other_abstracts": other_abstracts,
        "copyright": node_text(primary_abstract, "CopyrightInformation"),
        "authors": parse_authors(content.find("AuthorList")),
        "languages": languages,
        "journal": {
            "title": node_text(journal, "Title"),
            "iso_abbreviation": node_text(journal, "ISOAbbreviation"),
            "issn": node_text(journal, "ISSN"),
            "issn_type": (
                journal.find("ISSN").get("IssnType", "")
                if journal is not None and journal.find("ISSN") is not None
                else ""
            ),
            "volume": node_text(journal_issue, "Volume"),
            "issue": node_text(journal_issue, "Issue"),
            "publication_date": date_object(pub_date),
        },
        "pagination": node_text(content, "Pagination/MedlinePgn"),
        "electronic_location_ids": [
            {"type": node.get("EIdType", ""), "value": clean_text(node)}
            for node in content.findall("ELocationID")
            if clean_text(node)
        ],
        "article_dates": article_dates,
        "publication_types": publication_types,
        "keywords": keywords,
        "mesh_headings": parse_mesh(citation),
        "identifiers": identifiers,
        "citation_status": citation.get("Status", "") if citation is not None else "",
        "date_completed": date_object(citation.find("DateCompleted"))
        if citation is not None
        else {},
        "date_revised": date_object(citation.find("DateRevised"))
        if citation is not None
        else {},
        "source_file": source_file,
    }
    return result


def iter_xml_gz_records(path: Path) -> Iterator[dict[str, Any]]:
    """Yield normalized records without loading a complete XML file into memory."""
    with gzip.open(path, "rb") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            tag = elem.tag.rsplit("}", 1)[-1]
            if tag not in {"PubmedArticle", "PubmedBookArticle"}:
                continue
            record = parse_record(elem, path.name)
            if record is not None:
                yield record
            elem.clear()


def parse_xml_gz(path: Path) -> list[dict[str, Any]]:
    return list(iter_xml_gz_records(path))


def request_bytes(
    endpoint: str,
    params: dict[str, Any],
    retries: int = 6,
) -> bytes:
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            retryable = exc.code == 429 or 500 <= exc.code < 600
            if not retryable or attempt == retries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries - 1:
                raise
        time.sleep(min(60, 2**attempt))
    raise RuntimeError("unreachable")


def common_params(args: argparse.Namespace) -> dict[str, str]:
    params = {"tool": args.tool}
    if args.email:
        params["email"] = args.email
    if args.api_key:
        params["api_key"] = args.api_key
    return params


def create_history(args: argparse.Namespace) -> tuple[int, str, str, str]:
    query = args.query.strip() if args.query else f'{args.year}[dp] AND hasabstract[text]'
    params: dict[str, Any] = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 0,
        "usehistory": "y",
        **common_params(args),
    }
    data = json.loads(request_bytes("esearch.fcgi", params))
    result = data["esearchresult"]
    return int(result["count"]), result["querykey"], result["webenv"], query


def download_batch(
    args: argparse.Namespace,
    query_key: str,
    webenv: str,
    start: int,
    size: int,
) -> bytes:
    params: dict[str, Any] = {
        "db": "pubmed",
        "query_key": query_key,
        "WebEnv": webenv,
        "retstart": start,
        "retmax": size,
        "retmode": "xml",
        "rettype": "abstract",
        **common_params(args),
    }
    return request_bytes("efetch.fcgi", params)


def atomic_write_bytes(path: Path, data: bytes) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def write_gzip(path: Path, data: bytes) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temporary, "wb", compresslevel=6) as handle:
        handle.write(data)
    os.replace(temporary, path)


def write_jsonl_gz(path: Path, records: list[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=6) as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")))
            handle.write("\n")
    os.replace(temporary, path)


def init_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS articles (
            pmid TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            abstract TEXT NOT NULL,
            publication_year TEXT,
            journal TEXT,
            doi TEXT,
            pmcid TEXT,
            languages_json TEXT NOT NULL,
            abstract_sections_json TEXT NOT NULL,
            other_abstracts_json TEXT NOT NULL,
            authors_json TEXT NOT NULL,
            publication_types_json TEXT NOT NULL,
            keywords_json TEXT NOT NULL,
            mesh_headings_json TEXT NOT NULL,
            record_json TEXT NOT NULL,
            source_file TEXT NOT NULL
        )
        """
    )
    connection.execute(
        "CREATE INDEX IF NOT EXISTS idx_articles_year ON articles(publication_year)"
    )
    return connection


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def upsert_records(connection: sqlite3.Connection, records: list[dict[str, Any]]) -> None:
    rows = []
    for record in records:
        journal = record["journal"]
        identifiers = record["identifiers"]
        pub_year = journal.get("publication_date", {}).get("year", "")
        rows.append(
            (
                record["pmid"],
                record["title"],
                record["abstract"],
                pub_year,
                journal.get("title", ""),
                identifiers.get("doi", ""),
                identifiers.get("pmc", ""),
                compact_json(record["languages"]),
                compact_json(record["abstract_sections"]),
                compact_json(record["other_abstracts"]),
                compact_json(record["authors"]),
                compact_json(record["publication_types"]),
                compact_json(record["keywords"]),
                compact_json(record["mesh_headings"]),
                compact_json(record),
                record["source_file"],
            )
        )
    with connection:
        connection.executemany(
            """
            INSERT OR REPLACE INTO articles VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            )
            """,
            rows,
        )


def load_state(path: Path, year: int, query: str) -> dict[str, Any]:
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        if state.get("year") != year or state.get("query") != query:
            raise ValueError(
                f"State file belongs to a different query: {state.get('query')!r}"
            )
        return state
    return {
        "year": year,
        "query": query,
        "next_retstart": 0,
        "records_written": 0,
        "batches_completed": 0,
        "created_at": utc_now(),
    }


def save_state(path: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> int:
    output = Path(args.output_dir).resolve()
    raw_dir = output / "raw"
    parsed_dir = output / "parsed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    parsed_dir.mkdir(parents=True, exist_ok=True)

    count, query_key, webenv, query = create_history(args)
    target = count if args.max_records == 0 else min(count, args.max_records)
    state_path = output / "state.json"
    state = load_state(state_path, args.year, query)
    state["pubmed_match_count"] = count
    state["target_records"] = target
    state["batch_size"] = args.batch_size
    save_state(state_path, state)

    start = int(state["next_retstart"])
    if start >= target:
        print(f"Already complete: {start:,}/{target:,} records")
        return 0

    connection = (
        None
        if args.no_batch_sqlite
        else init_database(output / "pubmed_abstracts.sqlite")
    )
    delay = args.delay
    if delay is None:
        delay = 0.11 if args.api_key else 0.4

    print(f"Query: {query}")
    print(f"PubMed matches: {count:,}; this run target: {target:,}")
    print(f"Resuming at retstart={start:,}; output={output}")

    try:
        while start < target:
            size = min(args.batch_size, target - start)
            stem = f"batch_{start:09d}_{start + size - 1:09d}"
            raw_path = raw_dir / f"{stem}.xml.gz"
            parsed_path = parsed_dir / f"{stem}.jsonl.gz"

            if not raw_path.exists():
                xml_data = download_batch(args, query_key, webenv, start, size)
                write_gzip(raw_path, xml_data)
            records = parse_xml_gz(raw_path)
            if not records:
                raise RuntimeError(f"NCBI returned no records for retstart={start}")

            for record in records:
                record["query_year"] = args.year
                record["fetched_at"] = utc_now()
            write_jsonl_gz(parsed_path, records)
            if connection is not None:
                upsert_records(connection, records)

            state["next_retstart"] = start + size
            state["records_written"] = int(state["records_written"]) + len(records)
            state["batches_completed"] = int(state["batches_completed"]) + 1
            state["last_batch"] = {
                "retstart": start,
                "requested": size,
                "records_parsed": len(records),
                "raw_file": str(raw_path.relative_to(output)),
                "raw_gzip_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                "parsed_file": str(parsed_path.relative_to(output)),
            }
            save_state(state_path, state)
            start += size
            print(
                f"[{start:,}/{target:,}] parsed={len(records):,} "
                f"raw={raw_path.stat().st_size / 1024 / 1024:.1f} MiB"
            )
            if start < target:
                time.sleep(delay)
    finally:
        if connection is not None:
            connection.close()

    print(f"Complete: {state['records_written']:,} records written")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2025)
    parser.add_argument(
        "--query",
        default="",
        help="custom PubMed query; when set it replaces the default YEAR[dp] query",
    )
    parser.add_argument(
        "--output-dir",
        default="pubmed_year_export/data/pubmed_2025",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=1000,
        help="maximum records to download; 0 means the complete year (default: 1000)",
    )
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    parser.add_argument("--tool", default="mti_agent_pubmed_year_export")
    parser.add_argument("--delay", type=float, default=None)
    parser.add_argument(
        "--no-batch-sqlite",
        action="store_true",
        help="store raw/XML and JSONL batches only; build the compact FTS database later",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.batch_size < 1 or args.batch_size > 1000:
        print("--batch-size must be between 1 and 1000", file=sys.stderr)
        return 2
    if args.max_records < 0:
        print("--max-records must be >= 0", file=sys.stderr)
        return 2
    if not args.email:
        print(
            "Warning: set --email or NCBI_EMAIL before a long/full run so NCBI can "
            "contact you about the traffic.",
            file=sys.stderr,
        )
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())

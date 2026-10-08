#!/usr/bin/env python3
"""Build filtered PubMed shards from an official Annual Baseline.

The command discovers source files from the official directory, downloads only
a small bounded set concurrently, validates every file against NLM's MD5, and
streams matching records into deterministic JSONL gzip shards. By default raw
XML files are deleted after their filtered shard and state have been committed.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import http.client
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, TextIO

try:
    from .download_pubmed_year import iter_xml_gz_records
except ImportError:  # Allow direct execution: python pubmed_year_export/build_...py
    from download_pubmed_year import iter_xml_gz_records


DEFAULT_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
USER_AGENT = "mti-agent-pubmed-baseline/1.0"
STATE_SCHEMA_VERSION = 1
FILTER_VERSION = "pubmed-mirna-abstract-years-v1"
YEAR_RE = re.compile(r"(?<!\d)(?:18|19|20|21)\d{2}(?!\d)")
MIRNA_RE = re.compile(
    r"(?:\bmicro[\s-]*RNAs?\b|\bmiRNAs?\b|"
    r"\b(?:[a-z]{3}-)?miR-\d[\w.-]*\b|\blet-7[\w.-]*\b)",
    re.IGNORECASE,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def save_state(path: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    atomic_write_text(
        path,
        json.dumps(state, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
    )


def load_state(
    path: Path,
    baseline_year: int,
    base_url: str,
    start_year: int,
    end_year: int,
) -> dict[str, Any]:
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        expected = {
            "baseline_year": baseline_year,
            "base_url": base_url,
            "filter_version": FILTER_VERSION,
            "start_year": start_year,
            "end_year": end_year,
        }
        mismatches = {
            key: (state.get(key), value)
            for key, value in expected.items()
            if state.get(key) != value
        }
        if mismatches:
            raise ValueError(f"State configuration mismatch: {mismatches}")
        state.setdefault("files", {})
        return state
    return {
        "state_schema_version": STATE_SCHEMA_VERSION,
        "baseline_year": baseline_year,
        "base_url": base_url,
        "filter_version": FILTER_VERSION,
        "filter_description": (
            "has Abstract or OtherAbstract; publication year intersects "
            f"{start_year}-{end_year}; title/abstract/MeSH/keywords match miRNA"
        ),
        "start_year": start_year,
        "end_year": end_year,
        "created_at": utc_now(),
        "complete": False,
        "files": {},
    }


def request_bytes(url: str, timeout: float, retries: int) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            retryable = exc.code == 429 or 500 <= exc.code < 600
            if not retryable or attempt == retries - 1:
                raise
        except (
            urllib.error.URLError,
            TimeoutError,
            ConnectionError,
            http.client.IncompleteRead,
        ):
            if attempt == retries - 1:
                raise
        time.sleep(min(60, 2**attempt))
    raise RuntimeError("unreachable")


def discover_source_files(
    base_url: str,
    baseline_year: int,
    timeout: float = 120,
    retries: int = 6,
) -> list[str]:
    listing = request_bytes(base_url, timeout, retries).decode("utf-8", errors="replace")
    two_digit_year = baseline_year % 100
    pattern = re.compile(rf"pubmed{two_digit_year:02d}n\d{{4}}\.xml\.gz")
    files = sorted(set(pattern.findall(listing)))
    if not files:
        raise RuntimeError(
            f"No PubMed {baseline_year} Baseline XML files found at {base_url}"
        )
    return files


def select_source_files(
    discovered: list[str], start_file: str, limit_files: int
) -> list[str]:
    if start_file:
        try:
            start_index = discovered.index(start_file)
        except ValueError as exc:
            raise ValueError(f"--start-file was not discovered: {start_file}") from exc
    else:
        start_index = 0
    candidates = discovered[start_index:]
    return candidates[:limit_files] if limit_files else candidates


def parse_official_md5(text: str, source_file: str) -> str:
    match = re.search(r"(?i)(?<![0-9a-f])[0-9a-f]{32}(?![0-9a-f])", text)
    if match is None:
        raise ValueError(f"No MD5 found in checksum for {source_file}")
    return match.group(0).lower()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_official_md5(
    base_url: str,
    source_file: str,
    checksum_path: Path,
    timeout: float,
    retries: int,
) -> str:
    if checksum_path.exists():
        text = checksum_path.read_text(encoding="utf-8", errors="replace")
    else:
        url = urllib.parse.urljoin(base_url, source_file + ".md5")
        text = request_bytes(url, timeout, retries).decode("utf-8", errors="replace")
        atomic_write_text(checksum_path, text)
    return parse_official_md5(text, source_file)


def _part_is_complete(part_path: Path, expected_md5: str, final_path: Path) -> bool:
    if part_path.exists() and md5_file(part_path) == expected_md5:
        os.replace(part_path, final_path)
        return True
    return False


def download_with_resume(
    url: str,
    final_path: Path,
    expected_md5: str,
    timeout: float,
    retries: int,
) -> int:
    """Download one file with Range resume and atomically commit after MD5."""
    if final_path.exists():
        if md5_file(final_path) == expected_md5:
            return final_path.stat().st_size
        final_path.unlink()

    part_path = final_path.with_name(final_path.name + ".part")
    if _part_is_complete(part_path, expected_md5, final_path):
        return final_path.stat().st_size

    for attempt in range(retries):
        offset = part_path.stat().st_size if part_path.exists() else 0
        headers = {"User-Agent": USER_AGENT}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", None) or response.getcode()
                append = bool(offset and status == 206)
                mode = "ab" if append else "wb"
                with part_path.open(mode) as handle:
                    while chunk := response.read(1024 * 1024):
                        handle.write(chunk)
            actual_md5 = md5_file(part_path)
            if actual_md5 != expected_md5:
                if attempt == retries - 1:
                    raise ValueError(
                        f"MD5 mismatch for {final_path.name}: "
                        f"expected {expected_md5}, got {actual_md5}"
                    )
                part_path.unlink(missing_ok=True)
                time.sleep(min(60, 2**attempt))
                continue
            os.replace(part_path, final_path)
            return final_path.stat().st_size
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and _part_is_complete(part_path, expected_md5, final_path):
                return final_path.stat().st_size
            retryable = exc.code == 416 or exc.code == 429 or 500 <= exc.code < 600
            if not retryable or attempt == retries - 1:
                raise
            if exc.code == 416:
                part_path.unlink(missing_ok=True)
        except (
            urllib.error.URLError,
            TimeoutError,
            ConnectionError,
            http.client.IncompleteRead,
        ):
            if attempt == retries - 1:
                raise
        time.sleep(min(60, 2**attempt))
    raise RuntimeError("unreachable")


def publication_years(record: dict[str, Any]) -> set[int]:
    """Return publication-related years, excluding completion/revision dates."""
    values: list[Any] = []
    publication_date = record.get("journal", {}).get("publication_date", {})
    values.extend(
        [publication_date.get("year", ""), publication_date.get("medline_date", "")]
    )
    for article_date in record.get("article_dates", []):
        values.append(article_date.get("year", ""))

    years: set[int] = set()
    for value in values:
        years.update(int(match) for match in YEAR_RE.findall(str(value)))
    return years


def has_abstract(record: dict[str, Any]) -> bool:
    if str(record.get("abstract", "")).strip():
        return True
    return any(str(item.get("text", "")).strip() for item in record.get("other_abstracts", []))


def mirna_search_text(record: dict[str, Any]) -> str:
    values: list[str] = [record.get("title", ""), record.get("abstract", "")]
    values.extend(
        item.get("text", "") for item in record.get("other_abstracts", [])
    )
    for heading in record.get("mesh_headings", []):
        values.append(heading.get("descriptor", ""))
        values.extend(
            qualifier.get("name", "") for qualifier in heading.get("qualifiers", [])
        )
    values.extend(item.get("text", "") for item in record.get("keywords", []))
    return "\n".join(str(value) for value in values if value)


def record_matches(
    record: dict[str, Any], start_year: int = 2020, end_year: int = 2026
) -> tuple[bool, list[int]]:
    matched_years = sorted(
        year for year in publication_years(record) if start_year <= year <= end_year
    )
    matches = bool(
        has_abstract(record)
        and matched_years
        and MIRNA_RE.search(mirna_search_text(record))
    )
    return matches, matched_years


def deterministic_gzip_writer(path: Path) -> tuple[TextIO, Any, Any]:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    return text, compressed, raw


def process_source_file(
    raw_path: Path,
    shard_path: Path,
    start_year: int,
    end_year: int,
) -> tuple[int, int, str]:
    temporary = shard_path.with_name(shard_path.name + ".tmp")
    temporary.unlink(missing_ok=True)
    writer, compressed, raw = deterministic_gzip_writer(temporary)
    scanned = 0
    selected = 0
    try:
        for record in iter_xml_gz_records(raw_path):
            scanned += 1
            matches, matched_years = record_matches(record, start_year, end_year)
            if not matches:
                continue
            record["baseline_filter"] = {
                "version": FILTER_VERSION,
                "matched_publication_years": matched_years,
            }
            writer.write(
                json.dumps(
                    record,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
            )
            writer.write("\n")
            selected += 1
    except BaseException:
        writer.close()
        compressed.close()
        raw.close()
        temporary.unlink(missing_ok=True)
        raise
    else:
        writer.flush()
        writer.close()
        compressed.close()
        raw.close()
    os.replace(temporary, shard_path)
    return scanned, selected, sha256_file(shard_path)


def completed_shard_is_valid(root: Path, entry: dict[str, Any]) -> bool:
    if not entry.get("completed") or not entry.get("filtered_shard"):
        return False
    shard = root / entry["filtered_shard"]
    expected = entry.get("filtered_sha256", "")
    return bool(shard.exists() and expected and sha256_file(shard) == expected)


def download_source(
    source_file: str,
    base_url: str,
    raw_dir: Path,
    checksum_dir: Path,
    timeout: float,
    retries: int,
) -> tuple[str, Path, str, int]:
    checksum_path = checksum_dir / f"{source_file}.md5"
    official_md5 = fetch_official_md5(
        base_url, source_file, checksum_path, timeout, retries
    )
    raw_path = raw_dir / source_file
    source_url = urllib.parse.urljoin(base_url, source_file)
    size = download_with_resume(
        source_url, raw_path, official_md5, timeout, retries
    )
    return source_file, raw_path, official_md5, size


def process_download_result(
    result: tuple[str, Path, str, int],
    root: Path,
    parsed_dir: Path,
    state_path: Path,
    state: dict[str, Any],
    start_year: int,
    end_year: int,
    keep_raw: bool,
) -> None:
    source_file, raw_path, official_md5, download_bytes = result
    shard_path = parsed_dir / source_file.replace(".xml.gz", ".jsonl.gz")
    scanned, selected, filtered_sha256 = process_source_file(
        raw_path, shard_path, start_year, end_year
    )
    entry = {
        "source_file": source_file,
        "source_url": urllib.parse.urljoin(state["base_url"], source_file),
        "official_md5": official_md5,
        "download_bytes": download_bytes,
        "records_scanned": scanned,
        "records_selected": selected,
        "filtered_shard": shard_path.relative_to(root).as_posix(),
        "filtered_sha256": filtered_sha256,
        "completed": True,
        "completed_at": utc_now(),
    }
    state["files"][source_file] = entry
    save_state(state_path, state)
    if not keep_raw:
        raw_path.unlink(missing_ok=True)
    print(
        f"[{source_file}] scanned={scanned:,} selected={selected:,} "
        f"raw={download_bytes / 1024 / 1024:.1f} MiB"
    )


def run(args: argparse.Namespace) -> int:
    root = Path(args.output_dir).resolve()
    raw_dir = root / "raw"
    parsed_dir = root / "parsed"
    checksum_dir = root / "source_md5"
    for directory in (raw_dir, parsed_dir, checksum_dir):
        directory.mkdir(parents=True, exist_ok=True)

    base_url = args.base_url.rstrip("/") + "/"
    discovered = discover_source_files(
        base_url, args.baseline_year, args.timeout, args.retries
    )
    selected = select_source_files(discovered, args.start_file, args.limit_files)
    state_path = root / "baseline_state.json"
    state = load_state(
        state_path,
        args.baseline_year,
        base_url,
        args.start_year,
        args.end_year,
    )
    state["discovered_file_count"] = len(discovered)
    state["discovered_first_file"] = discovered[0]
    state["discovered_last_file"] = discovered[-1]

    pending: list[str] = []
    for source_file in selected:
        entry = state["files"].get(source_file, {})
        if completed_shard_is_valid(root, entry):
            continue
        if entry.get("completed"):
            entry["completed"] = False
            entry["error"] = "filtered shard missing or SHA-256 mismatch"
        pending.append(source_file)
    state["complete"] = False
    save_state(state_path, state)

    print(
        f"Discovered {len(discovered):,} files; selected {len(selected):,}; "
        f"pending {len(pending):,}; output={root}"
    )

    source_iter = iter(pending)
    futures: dict[Future[tuple[str, Path, str, int]], str] = {}
    executor = ThreadPoolExecutor(max_workers=args.workers)
    try:
        for source_file in source_iter:
            future = executor.submit(
                download_source,
                source_file,
                base_url,
                raw_dir,
                checksum_dir,
                args.timeout,
                args.retries,
            )
            futures[future] = source_file
            if len(futures) >= args.workers:
                break

        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                source_file = futures.pop(future)
                try:
                    process_download_result(
                        future.result(),
                        root,
                        parsed_dir,
                        state_path,
                        state,
                        args.start_year,
                        args.end_year,
                        args.keep_raw_baseline,
                    )
                except BaseException as exc:
                    entry = state["files"].setdefault(source_file, {})
                    entry.update(
                        {
                            "source_file": source_file,
                            "completed": False,
                            "error": f"{type(exc).__name__}: {exc}",
                        }
                    )
                    save_state(state_path, state)
                    raise
                try:
                    next_source = next(source_iter)
                except StopIteration:
                    continue
                next_future = executor.submit(
                    download_source,
                    next_source,
                    base_url,
                    raw_dir,
                    checksum_dir,
                    args.timeout,
                    args.retries,
                )
                futures[next_future] = next_source
    finally:
        executor.shutdown(wait=True, cancel_futures=True)

    completed_files = sum(
        completed_shard_is_valid(root, state["files"].get(name, {}))
        for name in discovered
    )
    state["completed_file_count"] = completed_files
    state["records_scanned"] = sum(
        int(entry.get("records_scanned", 0)) for entry in state["files"].values()
    )
    state["records_selected"] = sum(
        int(entry.get("records_selected", 0)) for entry in state["files"].values()
    )
    state["complete"] = completed_files == len(discovered)
    if state["complete"]:
        state["completed_at"] = utc_now()
    save_state(state_path, state)
    print(
        f"Run complete: {completed_files:,}/{len(discovered):,} source files; "
        f"selected records={state['records_selected']:,}; "
        f"baseline_complete={state['complete']}"
    )
    if args.build_local_index:
        if not state["complete"]:
            print("Skipping local index build because the complete Baseline is not ready")
        else:
            try:
                from .build_local_index import build as build_local_index
            except ImportError:
                from build_local_index import build as build_local_index
            print("Building deterministic merged corpus and SQLite FTS5 index")
            return build_local_index(
                argparse.Namespace(
                    export_dir=root,
                    corpus_version=args.corpus_version,
                )
            )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-year", type=int, default=2026)
    parser.add_argument("--start-year", type=int, default=2020)
    parser.add_argument("--end-year", type=int, default=2026)
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument(
        "--output-dir",
        default="pubmed_local_backend/data/pubmed_mirna_baseline_2026",
    )
    parser.add_argument(
        "--limit-files",
        type=int,
        default=0,
        help="process only the first N discovered files; 0 means all files",
    )
    parser.add_argument(
        "--start-file",
        default="",
        help="begin the selected range at this exact discovered source filename",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument(
        "--keep-raw-baseline",
        action="store_true",
        help="retain validated source XML.gz files after filtering",
    )
    parser.add_argument(
        "--build-local-index",
        action="store_true",
        help="after a complete Baseline, build merged JSONL, PMID list, and FTS5",
    )
    parser.add_argument(
        "--corpus-version",
        default="pubmed-mirna-baseline-2026-v1",
        help="corpus version recorded when --build-local-index runs",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.start_year > args.end_year:
        print("--start-year must be <= --end-year", file=sys.stderr)
        return 2
    if args.limit_files < 0:
        print("--limit-files must be >= 0", file=sys.stderr)
        return 2
    if args.workers < 1 or args.workers > 4:
        print("--workers must be between 1 and 4", file=sys.stderr)
        return 2
    if args.retries < 1:
        print("--retries must be >= 1", file=sys.stderr)
        return 2
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())

"""
Summarize GO/KEGG enrichment CSV files with a local Ollama model.

Expected results directory layout:

    results7/
      gokegg/
        *_GO_*.csv
        *_KEGG_*.csv
      summarized_keywords/
        *_summary.txt
        *_summary.csv
        all_gokegg_summaries.csv

If the gokegg folder is missing or empty, the script can copy matching
GO/KEGG CSV files from the selected results directory root into gokegg/.
Original files are not moved or deleted.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import shutil
import socket
import sys
import textwrap
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


QUESTION = (
    "Below are mutiple kegg terms, please summarize three most important "
    "biologocial function terms. Make sure to choose the functional terms "
    "from the list i provided"
)

TERM_COLUMN_CANDIDATES = (
    "Term",
    "term",
    "Description",
    "description",
    "Pathway",
    "pathway",
    "Pathway Name",
    "pathway_name",
    "Name",
    "name",
)

METRIC_COLUMN_CANDIDATES = (
    "Adjusted P-value",
    "P-value",
    "Combined Score",
    "Overlap",
    "NES",
    "FDR q-val",
    "NOM p-val",
    "fdr",
    "nes",
    "pval",
    "size",
)

PREFERRED_MODELS = (
    "qwen3.5:9b",
    "qwen3:8b",
    "qwen3:latest",
    "llama3.1:latest",
    "llama3.1",
    "llama3:latest",
    "llama3",
    "qwen2.5:7b",
    "qwen2.5:latest",
    "mistral:latest",
    "mistral",
)


@dataclass
class TermRecord:
    term: str
    metrics: dict[str, str]


@dataclass
class SummaryResult:
    source_file: str
    selected_terms: list[str]
    summary: str
    reason: str
    raw_response: str
    validation_status: str
    model: str
    term_count: int


def read_text(path: Path) -> str:
    for encoding in ("utf-8-sig", "utf-8", "utf-16", "gbk", "latin-1"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    return path.read_text(encoding="utf-8", errors="replace")


def dedupe_headers(headers: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    result: list[str] = []
    for header in headers:
        key = header.strip() or "unnamed"
        seen[key] = seen.get(key, 0) + 1
        if seen[key] == 1:
            result.append(key)
        else:
            result.append(f"{key}__{seen[key]}")
    return result


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    text = read_text(path)
    reader = csv.reader(io.StringIO(text))
    try:
        headers = dedupe_headers(next(reader))
    except StopIteration:
        return [], []

    rows: list[dict[str, str]] = []
    for raw_row in reader:
        if not any(cell.strip() for cell in raw_row):
            continue
        if len(raw_row) > len(headers):
            extra_headers = [f"extra_{i}" for i in range(len(raw_row) - len(headers))]
            local_headers = headers + extra_headers
        else:
            local_headers = headers
        rows.append({h: raw_row[i].strip() if i < len(raw_row) else "" for i, h in enumerate(local_headers)})
    return headers, rows


def normalized_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def choose_term_column(headers: list[str]) -> str | None:
    for candidate in TERM_COLUMN_CANDIDATES:
        if candidate in headers:
            return candidate

    normalized = {normalized_column_name(h): h for h in headers}
    for key in ("term", "description", "pathway", "pathwayname"):
        if key in normalized:
            return normalized[key]

    for header in headers:
        key = normalized_column_name(header)
        if "term" in key or "description" in key or "pathway" in key:
            return header
    return None


def read_terms_from_csv(path: Path, max_terms: int) -> list[TermRecord]:
    headers, rows = read_csv_rows(path)
    term_column = choose_term_column(headers)
    if not term_column:
        raise ValueError(f"No term-like column found in {path}")

    seen: set[str] = set()
    records: list[TermRecord] = []
    for row in rows:
        term = (row.get(term_column) or "").strip()
        if not term:
            continue
        term_key = normalize_for_match(term)
        if term_key in seen:
            continue
        seen.add(term_key)

        metrics: dict[str, str] = {}
        for col in METRIC_COLUMN_CANDIDATES:
            value = row.get(col, "")
            if value:
                metrics[col] = compact_metric(value)
        records.append(TermRecord(term=term, metrics=metrics))
        if len(records) >= max_terms:
            break
    return records


def compact_metric(value: str, max_chars: int = 90) -> str:
    value = re.sub(r"\s+", " ", value.strip())
    if len(value) <= max_chars:
        return value
    return value[: max_chars - 3] + "..."


def is_gokegg_csv(path: Path) -> bool:
    name = path.name.lower()
    if path.suffix.lower() != ".csv":
        return False
    excluded = ("differential_expression", "ranked_genes", "ranked_gene", "counts")
    if any(token in name for token in excluded):
        return False
    return (
        "kegg" in name
        or "gokegg" in name
        or bool(re.search(r"(^|[_\-.])go([_\-.]|$)", name))
        or "go_biological_process" in name
    )


def prepare_folders(results_dir: Path, input_subdir: str, output_subdir: str, organize_existing: bool) -> tuple[Path, Path]:
    input_dir = results_dir / input_subdir
    output_dir = results_dir / output_subdir
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    existing_input_files = sorted(input_dir.glob("*.csv"))
    if organize_existing and not existing_input_files:
        for source in sorted(results_dir.glob("*.csv")):
            if not is_gokegg_csv(source):
                continue
            dest = input_dir / source.name
            if not dest.exists():
                shutil.copy2(source, dest)

    return input_dir, output_dir


def read_extractor_keywords(path: Path, top_n: int) -> list[str]:
    if not path.exists():
        return []

    headers, rows = read_csv_rows(path)
    if not rows:
        return []

    keyword_col = None
    for candidate in ("keyword", "Keyword", "term", "Term"):
        if candidate in headers:
            keyword_col = candidate
            break
    if not keyword_col:
        keyword_col = headers[0]

    keywords: list[str] = []
    seen: set[str] = set()
    for row in rows:
        keyword = (row.get(keyword_col) or "").strip()
        key = normalize_for_match(keyword)
        if keyword and key not in seen:
            keywords.append(keyword)
            seen.add(key)
        if len(keywords) >= top_n:
            break
    return keywords


def format_term_list(records: list[TermRecord]) -> str:
    lines: list[str] = []
    for index, record in enumerate(records, 1):
        metric_bits = [f"{k}={v}" for k, v in record.metrics.items() if v]
        if metric_bits:
            lines.append(f"{index}. {record.term} | " + " | ".join(metric_bits))
        else:
            lines.append(f"{index}. {record.term}")
    return "\n".join(lines)


def build_prompt(source_file: str, records: list[TermRecord], extractor_keywords: list[str]) -> str:
    extractor_context = ", ".join(extractor_keywords) if extractor_keywords else "No extractor keywords provided."
    term_list = format_term_list(records)
    return textwrap.dedent(
        f"""
        You are analyzing GO/KEGG enrichment results for a biological study.

        User question:
        {QUESTION}

        Important constraints:
        - Select exactly three biological function terms.
        - Each selected term must be copied exactly from the provided GO/KEGG term list.
        - Do not invent, merge, rename, or paraphrase terms.
        - Use the extractor keywords only as background context, not as selectable terms unless they also appear in the GO/KEGG list.
        - Prefer terms that are biologically functional and important, not broad disease labels when a clearer pathway/process term is available.
        - Do not include chain-of-thought, thinking tags, Markdown, or any text outside the JSON object.
        - Keep the summary and reason concise.

        Extractor keyword context:
        {extractor_context}

        Source CSV:
        {source_file}

        Provided GO/KEGG terms:
        {term_list}

        Return only valid JSON with this schema:
        {{
          "selected_terms": ["exact term 1", "exact term 2", "exact term 3"],
          "summary": "one concise paragraph summarizing the shared biological functions",
          "reason": "brief reason for choosing these three exact terms"
        }}
        """
    ).strip()


def ollama_get_models(base_url: str, timeout: int) -> list[str]:
    url = base_url.rstrip("/") + "/api/tags"
    request = urllib.request.Request(url, method="GET")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return [model.get("name", "") for model in payload.get("models", []) if model.get("name")]


def choose_model(base_url: str, requested_model: str | None, timeout: int) -> str:
    if requested_model:
        return requested_model

    try:
        models = ollama_get_models(base_url, timeout)
    except Exception:
        return "qwen3:8b"

    model_set = set(models)
    for preferred in PREFERRED_MODELS:
        if preferred in model_set:
            return preferred
    return models[0] if models else "qwen3:8b"


def ollama_generate(base_url: str, model: str, prompt: str, timeout: int) -> str:
    url = base_url.rstrip("/") + "/api/generate"
    body = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "think": False,
        "options": {
            "temperature": 0,
            "top_p": 0.9,
            "repeat_penalty": 1.05,
            "num_ctx": 8192,
            "num_predict": 450,
        },
    }
    data = json.dumps(body).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
        raise RuntimeError(
            f"Could not connect to Ollama at {base_url}. Start Ollama or pass --ollama-url. Details: {exc}"
        ) from exc

    return str(payload.get("response", "")).strip()


def extract_json_object(text: str) -> dict | None:
    cleaned = text.strip()
    cleaned = re.sub(r"<think>.*?</think>", "", cleaned, flags=re.IGNORECASE | re.DOTALL).strip()
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned, flags=re.IGNORECASE)
        cleaned = re.sub(r"\s*```$", "", cleaned)

    candidates = [cleaned]
    start = cleaned.find("{")
    end = cleaned.rfind("}")
    if start != -1 and end != -1 and end > start:
        candidates.append(cleaned[start : end + 1])

    for candidate in candidates:
        try:
            obj = json.loads(candidate)
            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            continue
    return None


def normalize_for_match(value: str) -> str:
    value = re.sub(r"\(GO:\d+\)", "", value, flags=re.IGNORECASE)
    value = value.lower().strip()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def map_selected_terms(raw_terms: Iterable[str], allowed_terms: list[str]) -> tuple[list[str], list[str]]:
    exact_map = {normalize_for_match(term): term for term in allowed_terms}
    selected: list[str] = []
    invalid: list[str] = []

    for raw in raw_terms:
        raw = str(raw).strip()
        if not raw:
            continue
        key = normalize_for_match(raw)
        match = exact_map.get(key)
        if not match:
            match = find_loose_match(raw, allowed_terms)
        if match and match not in selected:
            selected.append(match)
        elif raw:
            invalid.append(raw)

    return selected[:3], invalid


def find_loose_match(raw: str, allowed_terms: list[str]) -> str | None:
    raw_key = normalize_for_match(raw)
    if not raw_key:
        return None

    candidates = []
    for term in allowed_terms:
        term_key = normalize_for_match(term)
        if raw_key == term_key or raw_key in term_key or term_key in raw_key:
            candidates.append(term)
    if len(candidates) == 1:
        return candidates[0]
    return None


def parse_llm_response(raw_response: str, allowed_terms: list[str]) -> tuple[list[str], str, str, str]:
    obj = extract_json_object(raw_response)
    summary = ""
    reason = ""
    raw_selected: list[str] = []

    if obj:
        value = obj.get("selected_terms", [])
        if isinstance(value, str):
            raw_selected = [value]
        elif isinstance(value, list):
            raw_selected = [str(item) for item in value]
        summary = str(obj.get("summary", "")).strip()
        reason = str(obj.get("reason", "")).strip()
    else:
        raw_lower = normalize_for_match(raw_response)
        for term in allowed_terms:
            if normalize_for_match(term) in raw_lower:
                raw_selected.append(term)
        summary = raw_response.strip()

    selected, invalid = map_selected_terms(raw_selected, allowed_terms)
    if len(selected) == 3 and not invalid:
        status = "ok"
    elif len(selected) == 3:
        status = "ok_with_invalid_extra"
    else:
        status = f"needs_review: selected {len(selected)} valid terms"
    return selected, summary, reason, status


def retry_prompt(original_prompt: str, raw_response: str, allowed_terms: list[str]) -> str:
    allowed = "\n".join(f"- {term}" for term in allowed_terms)
    return textwrap.dedent(
        f"""
        Your previous answer did not provide exactly three exact terms from the allowed list.

        Previous answer:
        {raw_response}

        Allowed terms:
        {allowed}

        Please answer again. Return only valid JSON:
        {{
          "selected_terms": ["exact term 1", "exact term 2", "exact term 3"],
          "summary": "one concise paragraph",
          "reason": "brief reason"
        }}

        Original task:
        {original_prompt}
        """
    ).strip()


def summarize_file(
    csv_path: Path,
    records: list[TermRecord],
    extractor_keywords: list[str],
    base_url: str,
    model: str,
    timeout: int,
    retries: int,
) -> SummaryResult:
    allowed_terms = [record.term for record in records]
    prompt = build_prompt(csv_path.name, records, extractor_keywords)

    raw_response = ollama_generate(base_url, model, prompt, timeout)
    selected, summary, reason, status = parse_llm_response(raw_response, allowed_terms)

    attempts = 0
    while len(selected) != 3 and attempts < retries:
        attempts += 1
        raw_response = ollama_generate(base_url, model, retry_prompt(prompt, raw_response, allowed_terms), timeout)
        selected, summary, reason, status = parse_llm_response(raw_response, allowed_terms)

    if len(selected) < 3:
        fill_terms = [term for term in allowed_terms if term not in selected]
        selected.extend(fill_terms[: 3 - len(selected)])
        status = status + "; fallback_filled_from_input_order"

    return SummaryResult(
        source_file=csv_path.name,
        selected_terms=selected[:3],
        summary=summary,
        reason=reason,
        raw_response=raw_response,
        validation_status=status,
        model=model,
        term_count=len(records),
    )


def write_result_files(result: SummaryResult, output_dir: Path) -> None:
    stem = Path(result.source_file).stem
    txt_path = output_dir / f"{stem}_summary.txt"
    csv_path = output_dir / f"{stem}_summary.csv"

    txt = textwrap.dedent(
        f"""
        Source file: {result.source_file}
        Model: {result.model}
        Validation: {result.validation_status}

        Selected terms:
        1. {result.selected_terms[0] if len(result.selected_terms) > 0 else ""}
        2. {result.selected_terms[1] if len(result.selected_terms) > 1 else ""}
        3. {result.selected_terms[2] if len(result.selected_terms) > 2 else ""}

        Summary:
        {result.summary}

        Reason:
        {result.reason}

        Raw response:
        {result.raw_response}
        """
    ).strip()
    txt_path.write_text(txt + "\n", encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fieldnames())
        writer.writeheader()
        writer.writerow(result_to_row(result))


def summary_fieldnames() -> list[str]:
    return [
        "source_file",
        "selected_term_1",
        "selected_term_2",
        "selected_term_3",
        "summary",
        "reason",
        "validation_status",
        "model",
        "term_count",
        "raw_response",
    ]


def result_to_row(result: SummaryResult) -> dict[str, str | int]:
    terms = result.selected_terms + ["", "", ""]
    return {
        "source_file": result.source_file,
        "selected_term_1": terms[0],
        "selected_term_2": terms[1],
        "selected_term_3": terms[2],
        "summary": result.summary,
        "reason": result.reason,
        "validation_status": result.validation_status,
        "model": result.model,
        "term_count": result.term_count,
        "raw_response": result.raw_response,
    }


def write_combined_summary(results: list[SummaryResult], output_dir: Path) -> Path:
    path = output_dir / "all_gokegg_summaries.csv"
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fieldnames())
        writer.writeheader()
        for result in results:
            writer.writerow(result_to_row(result))
    return path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Use Ollama to summarize the three most important GO/KEGG biological function terms.",
    )
    parser.add_argument("--results-dir", default="results7", help="Results folder containing gokegg/ and summarized_keywords/.")
    parser.add_argument(
        "--extractor-csv",
        default="extractor/output/keywords_comparison.csv",
        help="Extractor keyword CSV used as background context.",
    )
    parser.add_argument("--input-subdir", default="gokegg", help="Subfolder under results-dir containing GO/KEGG CSV files.")
    parser.add_argument(
        "--output-subdir",
        default="summarized_keywords",
        help="Subfolder under results-dir where summaries are saved.",
    )
    parser.add_argument("--model", default=None, help="Ollama model name. If omitted, the script picks an available model.")
    parser.add_argument("--ollama-url", default="http://localhost:11434", help="Ollama base URL.")
    parser.add_argument("--max-terms", type=int, default=50, help="Maximum terms from each CSV to send to the LLM.")
    parser.add_argument("--extractor-top-n", type=int, default=30, help="Number of extractor keywords to include as context.")
    parser.add_argument("--timeout", type=int, default=600, help="Ollama request timeout in seconds.")
    parser.add_argument("--retries", type=int, default=1, help="Retry count when the model returns invalid terms.")
    parser.add_argument(
        "--no-organize-existing",
        action="store_true",
        help="Do not copy matching root-level GO/KEGG CSV files into gokegg/ when gokegg/ is empty.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Parse files and show planned work without calling Ollama.")
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    extractor_csv = Path(args.extractor_csv).resolve()
    if not results_dir.exists():
        print(f"Results directory does not exist: {results_dir}", file=sys.stderr)
        return 2

    input_dir, output_dir = prepare_folders(
        results_dir=results_dir,
        input_subdir=args.input_subdir,
        output_subdir=args.output_subdir,
        organize_existing=not args.no_organize_existing,
    )

    csv_files = sorted(path for path in input_dir.glob("*.csv") if is_gokegg_csv(path))
    if not csv_files:
        print(f"No GO/KEGG CSV files found in {input_dir}", file=sys.stderr)
        return 2

    extractor_keywords = read_extractor_keywords(extractor_csv, args.extractor_top_n)

    planned: list[tuple[Path, list[TermRecord]]] = []
    for csv_path in csv_files:
        try:
            records = read_terms_from_csv(csv_path, args.max_terms)
        except Exception as exc:
            print(f"Skipping {csv_path.name}: {exc}", file=sys.stderr)
            continue
        if not records:
            print(f"Skipping {csv_path.name}: no terms found", file=sys.stderr)
            continue
        planned.append((csv_path, records))

    if not planned:
        print("No usable GO/KEGG CSV files found.", file=sys.stderr)
        return 2

    print(f"Results dir: {results_dir}")
    print(f"GO/KEGG input dir: {input_dir}")
    print(f"Summary output dir: {output_dir}")
    print(f"Extractor CSV: {extractor_csv if extractor_csv.exists() else 'not found'}")
    print(f"Extractor keywords loaded: {len(extractor_keywords)}")
    print(f"Files to summarize: {len(planned)}")
    for csv_path, records in planned:
        print(f"  - {csv_path.name}: {len(records)} terms")

    if args.dry_run:
        return 0

    model = choose_model(args.ollama_url, args.model, args.timeout)
    print(f"Using Ollama model: {model}")

    results: list[SummaryResult] = []
    for index, (csv_path, records) in enumerate(planned, 1):
        print(f"[{index}/{len(planned)}] Summarizing {csv_path.name}...")
        start = time.time()
        try:
            result = summarize_file(
                csv_path=csv_path,
                records=records,
                extractor_keywords=extractor_keywords,
                base_url=args.ollama_url,
                model=model,
                timeout=args.timeout,
                retries=args.retries,
            )
        except Exception as exc:
            result = SummaryResult(
                source_file=csv_path.name,
                selected_terms=[],
                summary="",
                reason="",
                raw_response=str(exc),
                validation_status=f"error: {exc}",
                model=model,
                term_count=len(records),
            )
        write_result_files(result, output_dir)
        results.append(result)
        elapsed = time.time() - start
        print(f"    selected: {', '.join(result.selected_terms) if result.selected_terms else '(none)'}")
        print(f"    validation: {result.validation_status}; elapsed: {elapsed:.1f}s")

    combined_path = write_combined_summary(results, output_dir)
    print(f"Combined summary saved: {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

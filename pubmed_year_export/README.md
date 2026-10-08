# PubMed abstract corpus

This package builds a filtered, provenance-bearing abstract index for MiTAgent. Source records come from the [NLM PubMed Annual Baseline](https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/).

## Corpus definition

The default configuration selects the 2026 baseline and publication years 2020–2026. A record must have Abstract or OtherAbstract text and match miRNA terminology in its title, abstracts, MeSH headings, or author keywords. Publication-date fields determine eligibility; revision and completion dates do not substitute for publication dates.

The generated corpus version is `pubmed-mirna-baseline-2026-v1`. This label identifies the filtering specification. Preserve the manifest and official source MD5 files to identify the exact source snapshot. The annual baseline directory can change over time; a different baseline year or source set must receive a different corpus version. The tool rejects an incompatible discovered baseline year.

## Build

Run from the repository root. First process two source files:

```bash
python -m pubmed_year_export.build_from_baseline --baseline-year 2026 --start-year 2020 --end-year 2026 --limit-files 2 --workers 2 --output-dir outputs/pubmed-baseline
```

Inspect the filtered records and state before continuing the full build:

```bash
python -m pubmed_year_export.build_from_baseline --baseline-year 2026 --start-year 2020 --end-year 2026 --workers 2 --output-dir outputs/pubmed-baseline --build-local-index --corpus-version pubmed-mirna-baseline-2026-v1
```

Completed source shards are recorded in `baseline_state.json`. Downloads use `.part` files and HTTP Range when supported. Official MD5 checksums are verified before parsing. After the filtered shard and state are saved, source XML is deleted by default; `--keep-raw-baseline` retains it. A full build transfers substantial data even when the final corpus is small. `--limit-files` and `--start-file` permit bounded inspection; a partial build is not a complete baseline snapshot.

If the required baseline is no longer in the live NLM directory, use an archived copy of the same verified source set and retain its checksums. Do not silently substitute a newer baseline for a paper reproduction. `--base-url` can identify a compatible source directory.

## Outputs and configuration

A complete indexed build produces:

- `articles.jsonl.gz`: normalized records, including structured Abstract and OtherAbstract fields.
- `pmids.txt.gz`: retained identifiers.
- `pubmed_mirna.sqlite`: articles and SQLite FTS5 index.
- `manifest.json`: corpus version, source metadata, filtering information, and record counts.
- `baseline_state.json`: source-shard progress and completion state.

Set `MITAGENT_PUBMED_LOCAL_DB` to the resulting SQLite file. In PowerShell:

```powershell
$env:MITAGENT_PUBMED_LOCAL_DB = (Resolve-Path "outputs/pubmed-baseline/pubmed_mirna.sqlite").Path
```

In Bash:

```bash
export MITAGENT_PUBMED_LOCAL_DB="$PWD/outputs/pubmed-baseline/pubmed_mirna.sqlite"
```

`local_pubmed_backend.py` opens the index through a read-only SQLite URI. Mature miRNA matches are retrieved first; the default `exact_then_family` policy adds family matches until the requested limit. The Python backend also supports `exact_only`. The main pipeline uses the default policy; record this choice in study methods.

Each matching pair produces compatible TXT abstracts and JSONL provenance, including PMID, DOI, publication year, BM25 score, and match type. Full abstract fields are not article full text. Daily Updates are not automatically added to an Annual Baseline corpus.

## Other utilities

`download_pubmed_year.py` and `download_pubmed_segments.py` provide E-utilities exports for bounded retrieval and format checks. They are distinct from the Annual Baseline workflow. Do not mix independently queried and baseline-filtered files into one corpus without an explicit provenance plan. `inspect_export.py` supports export inspection; use each module's `--help` for its interface.

## Tests

```bash
python -m unittest -v pubmed_year_export.test_parser pubmed_year_export.test_local_index pubmed_year_export.test_baseline_builder test_local_pubmed_backend
```

These tests use synthetic records and mocked HTTP responses. They do not download the PubMed baseline.

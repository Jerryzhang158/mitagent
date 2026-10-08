# Functional keyword extraction

This independent utility retrieves PubMed abstracts and ranks functional terms using author keywords, MeSH, TF-IDF, YAKE, and Borda aggregation. Its `config.py` is separate from the MTI pipeline configuration in the repository root.

## Setup

From the repository root:

```bash
python -m pip install -r requirements-core.txt
python -m pip install -r requirements-optional.txt
```

Configure `NCBI_EMAIL` in the environment before retrieval. `NCBI_API_KEY` is optional; keep its value outside source files and version control. Set `SEARCH_QUERY`, `MAX_RECORDS`, `EXTRACTION_METHOD`, and ranking options in `extractor/config.py` for the study.

## Execution

Run from the extractor directory:

```bash
cd extractor
python main.py --max 2000
```

Use `--query` to supply a query, or `--skip-fetch` to process an existing cache. The analyzer reads the full configured cache; use a separate clean cache for each study query. Otherwise results can include records from earlier queries.

Outputs in `output/` include term tables, figures, `keywords_borda_merged.csv`, and `keywords_comparison.csv`. The comparison table reports `keyword,borda_score,author_kw,mesh,tfidf,yake` and can be supplied to the independent GO/KEGG summarizer:

```bash
cd ..
python summarize_gokegg_with_ollama.py --results-dir outputs/study --extractor-csv extractor/output/keywords_comparison.csv --model qwen3.5:9b
```

Record the query, retrieval date, cache provenance, method, and filtering parameters alongside results. The outputs are descriptive corpus statistics, not evidence of a specific miRNA–gene interaction.

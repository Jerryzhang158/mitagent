# MiTAgent

MiTAgent is research software for prioritizing microRNA–target interactions (MTIs) in a defined biological context. It connects candidate interactions from prediction resources or user-specified pairs with PubMed-derived abstract evidence, embedding-based evidence scores, optional language-model interpretation, and network exports.

The intended workflow starts with a biological question and an expression-derived candidate set. Literature scores help identify interactions for follow-up experiments; they are evidence-ranking measures, not probabilities of binding or experimental confirmation.

## Workflow

1. Analyze expression matrices and relevant GO, KEGG, or GSEA results to define the study context.
2. Prepare gene and miRNA lists, or an explicit pair table.
3. Generate candidates using database-guided or direct input modes.
4. Retrieve abstracts and retain their identifiers and matching provenance.
5. Score the literature evidence against the supplied biological functions.
6. Optionally interpret the results with Ollama and export a Cytoscape-compatible network.

`rnaseq_analyzer.py` and `interactive_launcher.py` support upstream expression analysis. `main_pipeline.py` begins at the prepared-list or pair-table stage. Candidate selection from enrichment results is study-specific and must be documented by the investigator. The main CLI does not perform FASTQ processing or automatically select candidates from enrichment results.

## Candidate modes

| Mode | Input | Candidate construction |
|---|---|---|
| `gene` | Gene list | Search for candidate regulatory miRNAs across TargetScan, miRDB, and miRWalk; integrate miRTarBase evidence. |
| `mirna` | miRNA list; optional gene list | Search for target genes and optionally restrict them to the supplied genes. |
| `combined` | Both lists | Merge and deduplicate the gene-led and miRNA-led searches. This is their union, not a strict intersection. |
| `direct` | Pair table, or both lists | Evaluate explicit pairs or the complete Cartesian product, without prediction-database filtering. |

`direct` forces the local SQLite abstract backend and skips LLM interpretation. Network export is optional. The `--min-databases` setting governs the gene-led prediction intersections; the miRNA-led path can include single-resource pairs. Strong miRTarBase-supported pairs can enter independently of prediction-resource counts.

## Installation

Use Python 3.11. From a terminal:

```bash
git clone https://github.com/Jerryzhang158/mitagent.git
cd mitagent
python -m venv .venv
```

Activate with `.venv\Scripts\Activate.ps1` in PowerShell, or `source .venv/bin/activate` on Linux/macOS. Then:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
python main_pipeline.py --help
```

Alternatively, run `conda env create -f environment.yml` followed by `conda activate mitagent` from the repository root. The Python environment does not install external datasets, embedding weights, R, or Ollama. GPU users should install a compatible [PyTorch build](https://pytorch.org/get-started/locally/) for their hardware before the requirements file.

For upstream expression analysis, LLM interpretation, and keyword extraction:

```bash
python -m pip install -r requirements-optional.txt
```

For DESeq2, install R and use [Bioconductor](https://bioconductor.org/install/) from an R session:

```r
install.packages(c("BiocManager", "jsonlite", "dplyr", "ggplot2"))
BiocManager::install("DESeq2")
```

Make `Rscript` available on PATH, or set `RNA_ANALYZER_RSCRIPT` to its executable. Raw counts are required for DESeq2; TPM/FPKM are not interchangeable inputs. The expression module has a Python statistical fallback, which is not DESeq2 and must be identified separately in study methods. Inspect the generated configuration with `python rnaseq_analyzer.py --create-config` and the available options with `--help` before using study data.

## Data and model preparation

### Local abstract corpus

The retrieval backend uses SQLite FTS5. Build a miRNA-focused corpus using [pubmed_year_export](pubmed_year_export/README.md), or point to a compatible index through `MITAGENT_PUBMED_LOCAL_DB`.

The corpus specification `pubmed-mirna-baseline-2026-v1` filters the 2026 PubMed Annual Baseline for records with Abstract or OtherAbstract, publication dates intersecting 2020–2026, and miRNA terms in titles, abstracts, MeSH, or keywords. It is a topical abstract collection, not all PubMed or full-text articles. Annual Baseline and subsequent Daily Updates are distinct snapshots; the builder does not automatically incorporate Daily Updates.

Retrieval first matches the mature miRNA name and then supplements with family matches up to `--max-articles`. JSONL records retain the match type, PMID, DOI, year, and abstract fields. The default CLI uses `exact_then_family`; strict mature-arm comparisons require an explicitly controlled `exact_only` backend workflow. Do not pool these retrieval policies in a benchmark.

Downloaded corpora and indexes are excluded from Git. A full baseline build downloads substantial source data; use a bounded pilot before a full build and retain the source checksums and generated manifest.

### Prediction resources

The `gene`, `mirna`, and `combined` modes require prepared TargetScan, miRDB, and miRWalk tables. miRTarBase is optional. Obtain data from the resource providers and record their releases, species, filtering rules, and checksums. Prepare columns used by the loaders:

| Resource | Required identifying fields |
|---|---|
| TargetScan | `miRNA`, `Gene Symbol` |
| miRDB | `miRNA`, `RefSeq`, `Score` |
| miRWalk | `miRNA`, `Genesymbol`, `binding_probability` |
| miRTarBase | `miRNA`, `Target Gene`; evidence metadata includes `Support Type`, `Experiments`, `References (PMID)` |

Raw provider exports may need conversion to these schemas. RefSeq and network identifier mapping can contact MyGene. Database-guided execution therefore may require internet access even when abstract retrieval is local.

### Embeddings and LLMs

The default embedding model is [`NeuML/pubmedbert-base-embeddings`](https://huggingface.co/NeuML/pubmedbert-base-embeddings). It is selected with `--bert-model`, which also accepts a local SentenceTransformer model directory. The first use of a model ID may download weights. Pin and archive the exact model revision for a study.

The manuscript analyses report GeneralBERT-Mini (`sentence-transformers/all-MiniLM-L6-v2`). To select that model explicitly, pass `--bert-model sentence-transformers/all-MiniLM-L6-v2`. Selecting the name alone does not reconstruct the historical corpus, inputs, thresholds, or software environment.

LLM interpretation uses a local [Ollama installation](https://ollama.com/download), with `qwen3.5:9b` as the default (`--llm-model`). Run `ollama pull qwen3.5:9b` and ensure the service is reachable at `localhost:11434`. It is unnecessary for direct mode. LLM summaries are hypotheses and interpretations that require verification against the cited evidence.

## Synthetic example

The files in `examples/` contain a fictional pair and synthetic evidence text. They demonstrate software execution and do not represent a biological result or PubMed citation.

Create the fixture once:

```bash
python examples/create_demo_corpus.py --output outputs/demo/pubmed.sqlite
```

PowerShell:

```powershell
$env:MITAGENT_PUBMED_LOCAL_DB = (Resolve-Path "outputs/demo/pubmed.sqlite").Path
python main_pipeline.py --mode direct --pairs examples/pairs.csv --functions "oxidative stress" --output outputs/pairs
```

Bash:

```bash
export MITAGENT_PUBMED_LOCAL_DB="$PWD/outputs/demo/pubmed.sqlite"
python main_pipeline.py --mode direct --pairs examples/pairs.csv --functions "oxidative stress" --output outputs/pairs
```

To evaluate the two synthetic genes against one miRNA:

```bash
python main_pipeline.py --mode direct --genes examples/genes.txt --mirnas examples/mirnas.txt --functions "oxidative stress" --output outputs/cartesian
```

This creates two candidate pairs; only the pair with matching fixture evidence proceeds to scoring. For a retrieval-only check without embedding weights, append `--skip-bert`. For candidate construction alone, append `--skip-literature --skip-bert`. Skipped scoring assigns no numerical evidence scores.

## Study execution

Direct pair tables accept CSV, TSV, and XLSX with `Gene` and `miRNA` columns. Common aliases such as `gene_name`, `Target_Gene`, `miRNA_name`, and `miRNA_norm` are accepted. List files contain one identifier per line.

The compatibility output format uses family-level filenames and scoring keys. A batch containing distinct mature miRNAs that map to the same gene/family key is rejected before retrieval to prevent evidence overwriting. Run those pairs separately and retain their mature identities when combining results. This safeguard does not change the retrieval or scoring rules.

Database-guided example (replace paths with prepared study inputs):

```bash
python main_pipeline.py --mode combined --genes data/genes.txt --mirnas data/mirnas.txt --targetscan data/targetscan.tsv --mirdb data/mirdb.tsv --mirwalk data/mirwalk.tsv --mirtarbase data/mirtarbase.csv --functions "oxidative stress,ferroptosis" --skip-llm --output outputs/study
```

Omit `--skip-llm` to enable Ollama interpretation in a database-guided run. `--analysis-mode` accepts `standard`, `functional`, or `both`. `--existing-step1` reuses a candidate table; it does not load previous literature or BERT outputs and is not a general checkpoint-resume mechanism.

For expression-annotated network export, add `--generate-network --mirna-deseq data/mirna_deseq.csv --gene-deseq data/gene_deseq.csv --network-score-threshold 60`. Both expression files must be supplied together. Gene tables use `mRNA`; miRNA tables use `miRNA`; DESeq2-style numeric fields include `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. Inspect identifier mapping in the network log.

For topology-only export without measured expression, explicitly use `--generate-network --allow-placeholder-expression`. This creates labeled placeholder expression values (zero log2 fold change and p/padj of one), not inferred or measured differential expression. Do not use those values as experimental evidence. Network export cannot be combined with `--skip-bert`.

### Reproducibility record

For each analysis retain the repository commit, full command, input checksums, resource releases, corpus manifest and retrieval policy, model revision, biological functions, score thresholds, and Python/R package versions. Runs write `run_parameters.json` with CLI arguments and selected installed package versions; external resource provenance must also be retained.

The default BERT validity threshold is 30 and the default network inclusion threshold is 50. Study-specific thresholds, such as 60, must be supplied explicitly. Scores are not calibrated probabilities. Changing the embedding model, sentence segmentation resources, retrieval coverage, or function terms can change ranking. Archive NLTK resources as well as model weights; if NLTK sentence resources are unavailable, the scorer uses its regex sentence splitter.

The repository supplies source and a synthetic execution example. It does not bundle private expression data, manuscript drafts, full prediction databases, model weights, or a complete archive of the paper's experimental runs. The example and unit tests establish software behavior, not full manuscript-result reproduction.

## Outputs

Each run is written under `<output>/run_<timestamp>/`:

| Path | Contents |
|---|---|
| `mirna_selection/` | Candidate table, including pair origin and database annotations |
| `pubmed_articles/mining_summary.csv` | Retrieval coverage per pair |
| `pubmed_articles/*.txt`, `*.jsonl` | Abstract text and record provenance |
| `bert_validation/mti_validation_results.csv` | Overall, relation, function, and evidence component scores |
| `llm_summaries/` | Optional interpretation outputs |
| `cytoscape_network/` | Optional nodes, edges, GraphML/SIF, and style files |
| `run_parameters.json` | CLI parameters and package versions |
| `pipeline_log.txt`, `pipeline_report.txt` | Execution details, stage status, and existing output files |

Only literature-covered pairs are scored. A requested stage that fails returns a nonzero exit status and is recorded as failed; insufficient evidence for a requested scoring stage is not reported as successful scoring. Expression support, database support, literature evidence, and direct binding experiments are different evidence levels.

## Additional research utilities

- [Functional keyword extraction](extractor/README.md): author keywords, MeSH, TF-IDF, YAKE, and Borda aggregation from PubMed abstracts.
- `summarize_gokegg_with_ollama.py`: selects three functional terms from supplied enrichment tables and validates that they occur in the input list. It is independent of the MTI CLI.
- `rnaseq_analyzer.py`: expression-matrix analysis and enrichment; `interactive_launcher.py` provides an interactive entry point.

## Tests

The deterministic suite uses temporary synthetic databases and mocked external services:

```bash
python -m pip install -r requirements-core.txt
python -m unittest -v tests.test_pipeline test_local_pubmed_backend pubmed_year_export.test_parser pubmed_year_export.test_local_index pubmed_year_export.test_baseline_builder
```

CI runs this suite on Windows and Linux with Python 3.11. It does not download embedding weights, execute Ollama, rebuild PubMed, or reproduce biological experiments. Use the synthetic example separately to check embedding inference in the selected environment.

## License and citation

Source code is distributed under the [Apache License 2.0](LICENSE). External datasets, publications, model weights, and services retain their respective terms.

For software attribution, cite this repository and the exact commit used: https://github.com/Jerryzhang158/mitagent. A manuscript citation should be added when its public bibliographic record is available; no publication DOI or release DOI is assigned here.

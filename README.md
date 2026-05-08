# MiTAgent: Enhanced Modular miRNA Research Pipeline

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-success.svg)]()

> A comprehensive computational pipeline for systematic analysis of microRNA regulatory networks with enhanced modular architecture

## Overview

MiTAgent provides an end-to-end solution for miRNA-mRNA interaction analysis, featuring a completely redesigned modular architecture that integrates multiple prediction databases, machine learning validation, and advanced functional interpretation. The system supports flexible analytical workflows with improved error handling, centralized configuration, and seamless step continuity.

## Key Features

| Feature | Description |
|---------|-------------|
| **Modular Architecture** | Component-based design with centralized configuration management |
| **Multi-Database Integration** | TargetScan, miRDB, miRWalk, miRTarBase support |
| **Machine Learning Validation** | BERT-based interaction validation with enhanced processing |
| **Literature Mining** | Automated PubMed literature retrieval with validation |
| **Advanced LLM Analysis** | Dual-mode analysis (Standard + Functional) with Ollama integration |
| **PubMed Functional Keyword Extraction** | Standalone `extractor/` workflow for PubMed keyword ranking across author keywords, MeSH, TF-IDF, YAKE, and Borda fusion |
| **GO/KEGG LLM Summarization** | Standalone Ollama utility that reads GO/KEGG CSV files, chooses exactly three functional terms from the provided list, and saves validated summaries |
| **Network Visualization** | Cytoscape-compatible network generation with auto-fallback |
| **Flexible Execution** | Resume from any step, skip optional components |
| **Enhanced Error Handling** | Robust error management with detailed logging |

## Quick Start

### Prerequisites

- Python 3.8+
- 8GB+ RAM (16GB recommended)
- 2GB disk space
- Internet connection

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/Jerryzhang158/mitagent.git
   cd mitagent
   ```

2. **Create Python environment**
   ```bash
   conda env create -f environment.yml
   conda activate mitagent
   ```

   If Conda cannot solve the full environment on your machine, create a lighter environment and install the core packages manually:

   ```bash
   conda create -n mitagent python=3.11
   conda activate mitagent
   pip install pandas numpy requests biopython openpyxl scikit-learn matplotlib seaborn

   # Optional modules used by LLM / enrichment / network workflows
   pip install langchain langchain-community gseapy networkx
   ```

3. **Install Ollama** (for LLM analysis)
   ```bash
   curl -fsSL https://ollama.ai/install.sh | sh
   ollama serve

   # Main MTI LLM summarizer currently defaults to qwen3:8b
   ollama pull qwen3:8b

   # Standalone GO/KEGG summarizer examples use qwen3.5:9b
   ollama pull qwen3.5:9b
   ```

### Basic Usage

```bash
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv
```

---

## Standalone Functional Interpretation Utilities

In addition to the main MTI pipeline, this repository includes two standalone workflows that can be run independently.

### 1. PubMed Functional Keyword Extractor

The `extractor/` folder contains a PubMed keyword analyzer. It searches PubMed, caches abstracts, extracts functional keywords, and exports CSV files and figures.

```bash
cd extractor

# Full workflow: fetch PubMed records, extract keywords, export CSV/figures
python main.py

# Reuse the existing cache and only rerun extraction/analysis
python main.py --skip-fetch

# Test with fewer records
python main.py --max 2000
```

Important configuration lives in `extractor/config.py`:

| Setting | Purpose |
|---------|---------|
| `SEARCH_QUERY` | PubMed query and date/language filters |
| `MAX_RECORDS` | Maximum PubMed records to fetch |
| `EXTRACTION_METHOD` | `author_kw`, `mesh`, `tfidf`, `yake`, or `combined` |
| `MIN_FREQ` | Minimum frequency filter for extracted terms |
| `TOP_K_DISPLAY` | Number of top terms exported |

Key extractor outputs:

```text
extractor/output/
├── keywords_author_kw.csv
├── keywords_mesh.csv
├── keywords_tfidf.csv
├── keywords_yake.csv
├── keywords_borda_merged.csv
└── keywords_comparison.csv
```

`keywords_comparison.csv` is the recommended background keyword file for downstream GO/KEGG LLM summarization because it includes the fused Borda score plus counts from all extraction methods:

```text
keyword,borda_score,author_kw,mesh,tfidf,yake
```

Note: the current extractor loads the full cache when analysis runs, so if `data/abstracts_cache.jsonl` contains records from older searches, the exported keyword CSVs may reflect the whole cache rather than only the current `SEARCH_QUERY`. Use a clean cache when strict query reproducibility is required.

### 2. GO/KEGG Term Summarization with Ollama

Use `summarize_gokegg_with_ollama.py` to summarize GO/KEGG enrichment CSV files with a local Ollama model. The script is independent from `main_pipeline.py`.

Expected result folder structure:

```text
results7/
├── gokegg/
│   ├── Gene_GO_Biological_Process_2021_enrichment.csv
│   ├── Gene_KEGG_2019_Human_enrichment.csv
│   └── Gene_KEGG_2021_Human_GSEA.csv
└── summarized_keywords/
    ├── *_summary.txt
    ├── *_summary.csv
    └── all_gokegg_summaries.csv
```

If `gokegg/` is missing or empty, the script automatically copies matching root-level GO/KEGG CSV files from the selected results folder into `gokegg/`. Original files are not moved or deleted.

Recommended command:

```bash
python summarize_gokegg_with_ollama.py \
    --results-dir results7 \
    --extractor-csv extractor/output/keywords_comparison.csv \
    --model qwen3.5:9b
```

Dry-run without calling Ollama:

```bash
python summarize_gokegg_with_ollama.py \
    --results-dir results7 \
    --extractor-csv extractor/output/keywords_comparison.csv \
    --model qwen3.5:9b \
    --dry-run
```

The fixed LLM question is:

```text
Below are mutiple kegg terms, please summarize three most important biologocial function terms. Make sure to choose the functional terms from the list i provided
```

The script enforces this constraint by validating the LLM output against the terms in the input CSV. Each selected term must come from the provided GO/KEGG term list.

Key options:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--results-dir` | `results7` | Folder containing GO/KEGG results |
| `--input-subdir` | `gokegg` | Subfolder containing GO/KEGG CSV files |
| `--output-subdir` | `summarized_keywords` | Subfolder where summaries are saved |
| `--extractor-csv` | `extractor/output/keywords_comparison.csv` | Background keyword context from the PubMed extractor |
| `--model` | auto-selected, prefers `qwen3.5:9b` | Ollama model |
| `--max-terms` | `50` | Number of top GO/KEGG terms sent to the LLM per CSV |
| `--dry-run` | off | Parse files and show planned work without calling Ollama |

Summary output columns:

```text
source_file,selected_term_1,selected_term_2,selected_term_3,summary,reason,validation_status,model,term_count,raw_response
```

`validation_status = ok` means all three selected terms were matched back to the original CSV term list.

---

## Enhanced Usage Options

### Resume from Existing Results

One of the major improvements is the ability to resume analysis from existing Step 1 results:

```bash
# Continue from previously generated MTI selection
python main_pipeline.py \
    --mode combined \
    --existing-step1 previous_results/mti_selection_results_combined.xlsx \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --analysis-mode both
```

### Flexible Step Control

The enhanced pipeline exposes skip flags for optional steps. In the current implementation, LLM analysis still depends on literature files and BERT validation output from the same run.

```bash
# Skip literature mining, BERT, and LLM; keep MTI selection and optional network/report generation
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --skip-literature \
    --skip-bert \
    --skip-llm \
    --generate-network
```

Current code note: the main LLM step expects `bert_validation/mti_validation_results.csv` and files in `pubmed_articles/` inside the current run directory. For LLM analysis, run literature mining and BERT first. If you only need GO/KEGG term summarization, use the standalone `summarize_gokegg_with_ollama.py` utility.

### Network Generation with Auto-Fallback

The pipeline now automatically generates default expression files if DESeq2 results aren't provided:

```bash
# Network generation with auto-generated expression files
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --generate-network \
    --network-score-threshold 60
    # DESeq2 files will be auto-generated if not provided
```

---

## Architecture Overview

### Core Components

```
mitagent/
├── main_pipeline.py          # Main controller with enhanced modularity
├── config.py                 # Centralized configuration management
├── utils.py                  # Utilities (Logger, FileUtils, ResultsIntegrator)
├── mti_selection.py          # MTI selection logic
├── literature_mining.py      # PubMed literature retrieval
├── mirna_matcher.py          # miRNA name matching and normalization
├── refseq_cache.py          # RefSeq identifier management
├── mti_llm_summarize.py     # LLM-based analysis (optional)
├── mti_cytoscape_network.py # Network generation (optional)
├── rnaseq_analyzer.py       # RNA-seq differential expression and GO/KEGG enrichment
├── summarize_gokegg_with_ollama.py # Standalone GO/KEGG term summarization with Ollama
└── extractor/               # Standalone PubMed functional keyword extractor
```

### Data Flow

```
Input Files → MTI Selection → Literature Mining → BERT Validation → LLM Analysis → Network Generation
     ↓              ↓              ↓              ↓              ↓              ↓
Config.py    ResultsIntegrator  TextProcessor  Enhanced Logging  Dual Analysis  Auto-Fallback
```

---

## Command Line Parameters

### New Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--existing-step1` | Resume from existing Step 1 results | `previous_results.xlsx` |
| `--analysis-mode` | LLM analysis type | `standard`, `functional`, `both` |
| `--network-score-threshold` | Network inclusion threshold | `60` |

### Enhanced Parameters

| Parameter | Enhancement | Default | New Behavior |
|-----------|-------------|---------|--------------|
| `--skip-*` | Skip optional steps | - | LLM requires literature + BERT outputs in the current run |
| `--generate-network` | Auto-generates DESeq2 files | `False` | Creates defaults when files missing |
| `--output` | Enhanced directory structure | `pipeline_results` | Improved organization |

---

## Pipeline Steps (Enhanced)

### Step 1: Enhanced MTI Selection
**New Features**:
- Load existing results with `--existing-step1`
- Improved database integration through `MTISelector` class
- Better error handling and validation
- Enhanced logging with step progress

### Step 2: Robust Literature Mining  
**New Features**:
- File validation after download
- Better error recovery
- Enhanced article processing through `LiteratureMiner` class
- Detailed download statistics

### Step 3: Advanced BERT Validation
**New Features**:
- Improved miRNA name handling via `TextProcessor`
- Better batch processing
- Enhanced result integration
- Graceful fallback when articles unavailable

### Step 4: Dual-Mode LLM Analysis
**New Features**:
- **Standard Analysis**: Relationship summarization and evidence assessment
- **Functional Analysis**: Pathway enrichment and systems-level interpretation  
- **Both Modes**: Comprehensive analysis with integrated results
- Enhanced Ollama connection management
- Improved text processing and result integration

### Step 5: Smart Network Generation
**New Features**:
- Auto-generation of default DESeq2 files
- Enhanced network statistics
- Improved error handling
- Better integration with expression data

---

## Output Structure (Enhanced)

```
pipeline_results/run_YYYYMMDD_HHMMSS/
├── pipeline_report.txt                    # 🆕 Enhanced with architectural notes
├── pipeline_log.txt                       # 🆕 Comprehensive logging
├── mirna_selection/
│   └── mti_selection_results_combined.xlsx # Can be reused with --existing-step1
├── pubmed_articles/                       # 🆕 Enhanced validation
│   ├── [literature files]
│   └── mining_summary.csv                # Literature mining summary
├── bert_validation/
│   ├── mti_validation_results.csv
│   ├── integrated_validation_results.csv  # 🆕 Enhanced integration
│   └── batch_input.csv                    # 🆕 Batch processing data
├── llm_summaries/                         # 🆕 Dual-mode analysis
│   ├── *_summary_YYYYMMDD_HHMMSS.json
│   ├── all_mti_summaries_YYYYMMDD_HHMMSS.json
│   ├── comprehensive_mti_report_YYYYMMDD_HHMMSS.txt
│   └── functional_analysis/              # Pathway and systems analysis
├── cytoscape_network/                     # 🆕 Enhanced with auto-fallback
│   └── cytoscape_YYYYMMDD_HHMMSS/
│       ├── nodes.csv
│       ├── edges.csv
│       ├── network.sif
│       ├── network_enhanced.graphml
│       ├── cytoscape_style_enhanced.xml
│       ├── network_summary.json
│       └── COMPREHENSIVE_USAGE_GUIDE.txt
├── default_mirna_deseq.csv              # 🆕 Auto-generated when needed
└── default_gene_deseq.csv               # 🆕 Auto-generated when needed
```

Standalone GO/KEGG summary output:

```text
results7/
├── gokegg/
│   └── [GO/KEGG enrichment CSV files]
└── summarized_keywords/
    ├── [source_name]_summary.txt
    ├── [source_name]_summary.csv
    └── all_gokegg_summaries.csv
```

---

## Advanced Examples

### Example 1: Complete Enhanced Analysis

```bash
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --mirna-deseq miRNA.csv \
    --gene-deseq mRNA.csv \
    --functions "cell proliferation,apoptosis,migration" \
    --max-articles 150 \
    --min-databases 3 \
    --analysis-mode both \
    --generate-network \
    --network-score-threshold 65 \
    --output results/enhanced_analysis
```

### Example 2: Resume and Complete Previous Analysis

```bash
# Step 1: Generate initial MTI selection
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --skip-literature --skip-bert --skip-llm

# Step 2: Resume with full analysis
python main_pipeline.py \
    --mode combined \
    --existing-step1 pipeline_results/run_*/mirna_selection/mti_selection_results_combined.xlsx \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --analysis-mode both \
    --generate-network
```

### Example 3: Fast Network/Report from Existing MTI Selection

```bash
# Reuse existing Step 1 MTIs, skip literature/BERT/LLM, and generate a network with default expression values
python main_pipeline.py \
    --mode combined \
    --existing-step1 previous_mti_results.xlsx \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --skip-literature \
    --skip-bert \
    --skip-llm \
    --generate-network
```

---

## Configuration Management

### Centralized Configuration

The enhanced pipeline uses a centralized configuration system via `config.py`:

```python
# Configuration is automatically managed
from config import PipelineConfig
config = PipelineConfig()

# Access default functions
config.DEFAULT_FUNCTIONS  # ['cell proliferation', 'apoptosis', ...]

# Directory structure auto-created
config.create_output_structure(results_dir)
```

### Runtime Configuration Notes

Most runtime defaults are currently defined in `config.py` and module constructors rather than environment variables.

| Setting | Current Code Location |
|---------|-----------------------|
| Main pipeline Ollama URL check | `utils.check_ollama_connection(base_url="http://localhost:11434")` |
| Main MTI LLM model | `mti_llm_summarize.MTILLMSummarizer(model_name="qwen3:8b")` |
| Standalone GO/KEGG Ollama URL | `summarize_gokegg_with_ollama.py --ollama-url` |
| Standalone GO/KEGG model | `summarize_gokegg_with_ollama.py --model qwen3.5:9b` |
| Default biological functions | `PipelineConfig.DEFAULT_FUNCTIONS` |

---

## Error Handling and Recovery

### Enhanced Error Management

The new architecture provides robust error handling:

```bash
# Pipeline continues even if individual steps fail
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    # ... database files ...
    --analysis-mode both
    
# Check pipeline_log.txt for detailed error information
# Final report includes success/failure status for each step
```

### Recovery Strategies

| Scenario | Solution | Command |
|----------|----------|---------|
| Step 1 completed, others failed | Use `--existing-step1` | `--existing-step1 results.xlsx` |
| Literature mining timeout | Reduce articles or skip | `--max-articles 25 --skip-literature` |
| BERT validation memory error | Skip BERT and LLM, or reduce the BERT batch externally | `--skip-bert --skip-llm` |
| Network generation failure | Check score threshold | `--network-score-threshold 40` |
| Main pipeline Ollama connection issues | Check service and install the current default model | `ollama serve && ollama pull qwen3:8b` |
| GO/KEGG summarizer Ollama connection issues | Check service and install the selected model | `ollama serve && ollama pull qwen3.5:9b` |
| Empty Qwen/Ollama response | Use the GO/KEGG summarizer or pass `think: false` in custom Ollama API calls | `python summarize_gokegg_with_ollama.py --dry-run` |

---

## Performance and Optimization

### Resource Management

The enhanced pipeline provides better resource management:

| Component | Memory Usage | Optimization |
|-----------|--------------|-------------|
| **MTI Selection** | Low (< 1GB) | Efficient database loading |
| **Literature Mining** | Medium (2-4GB) | Batched API requests |
| **BERT Validation** | High (4-8GB) | Configurable batch sizes |
| **LLM Analysis** | Medium (2-6GB) | Streaming processing |
| **Network Generation** | Low (< 2GB) | Efficient graph algorithms |

### Recommended Settings

| Dataset Size | Settings | Performance |
|--------------|----------|-------------|
| **Small** (<20 MTIs) | `--max-articles 100 --analysis-mode both` | ~30-60 min |
| **Medium** (20-100 MTIs) | `--max-articles 75 --analysis-mode both` | ~1-3 hours |
| **Large** (100+ MTIs) | `--max-articles 50 --analysis-mode standard` | ~3-8 hours |

---

## Troubleshooting (Updated)

### Component-Specific Issues

#### Module Import Errors
**Problem**: `ImportError: cannot import name 'MTILLMSummarizer'`

**Solution**:
```bash
# Check if all modules are present
ls -la *.py | grep -E "(config|utils|mti_|mirna_|literature_|refseq_)"

# Missing modules will be reported as warnings but pipeline continues
```

#### Configuration Issues
**Problem**: Pipeline configuration errors

**Solution**:
```bash
# Check configuration
python -c "from config import PipelineConfig; print('Config OK')"

# Reset configuration if needed
rm -f config.pyc __pycache__/config.*
```

#### Step Continuity Problems
**Problem**: Cannot resume from existing results

**Solution**:
```bash
# Verify file format and location
python -c "import pandas as pd; df = pd.read_excel('your_results.xlsx'); print(f'Shape: {df.shape}')"

# Use absolute paths
python main_pipeline.py --existing-step1 /full/path/to/results.xlsx
```

### Enhanced Logging

Check detailed logs for troubleshooting:

```bash
# Real-time log monitoring
tail -f pipeline_results/run_*/pipeline_log.txt

# Search for specific errors
grep -i "error\|failed\|exception" pipeline_results/run_*/pipeline_log.txt

# Check step completion status
grep -i "step.*completed\|step.*failed" pipeline_results/run_*/pipeline_log.txt
```

---

## Development and Extension

### Adding New Components

The modular architecture makes it easy to add new components:

```python
# Create new module: my_analysis.py
class MyAnalyzer:
    def __init__(self, logger):
        self.logger = logger
    
    def analyze(self, data):
        # Your analysis logic
        pass

# Integrate in main_pipeline.py
from my_analysis import MyAnalyzer

# Add to pipeline steps
analyzer = MyAnalyzer(self.logger)
results = analyzer.analyze(self.validation_results)
```

### Configuration Extension

Extend the configuration system:

```python
# In config.py
class PipelineConfig:
    def __init__(self):
        self.MY_NEW_SETTING = "default_value"
        self.MY_PARAMETERS = {
            'param1': 10,
            'param2': 'setting'
        }
```

---

## Migration from Previous Versions

### Key Changes

| Previous Version | Enhanced Version | Migration |
|------------------|------------------|-----------|
| Single file pipeline | Modular architecture | No code changes needed |
| Basic error handling | Comprehensive error recovery | Better stability |
| Fixed step execution | Flexible step control | New command options |
| Manual continuation | Automatic step resumption | `--existing-step1` option |
| Single LLM mode | Dual analysis modes | `--analysis-mode` parameter |

### Backward Compatibility

The enhanced pipeline maintains full backward compatibility:

```bash
# Old command still works
python main_pipeline.py -m combined -g genes.txt -r mirnas.txt -t target.txt -d mirdb.txt -w mirwalk.txt

# Enhanced features are optional
python main_pipeline.py \
    --mode combined \
    --genes genes.txt \
    --mirnas mirnas.txt \
    --targetscan target.txt \
    --mirdb mirdb.txt \
    --mirwalk mirwalk.txt \
    --analysis-mode both
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

```bibtex
@software{mitagent_enhanced_2024,
  title={MiTAgent: Enhanced Modular miRNA Research Pipeline},
  author={Jerry Zhang},
  year={2024},
  url={https://github.com/Jerryzhang158/mitagent},
  version={2.0-enhanced},
  note={Enhanced modular architecture with improved error handling and flexible execution}
}
```

## Acknowledgments

This enhanced version builds upon the original MiTAgent pipeline with significant architectural improvements including modular design, centralized configuration management, enhanced error handling, and flexible execution workflows.

**Key Dependencies**:
- **Database Resources**: TargetScan, miRDB, miRWalk, miRTarBase  
- **Machine Learning**: Hugging Face Transformers, BERT models
- **LLM Integration**: Ollama framework with local model deployment
- **Visualization**: Cytoscape network export capabilities
- **Development**: Enhanced Python architecture with component separation

---

## Support

- **Documentation**: This enhanced README with architectural details
- **Bug Reports**: [Open an issue](https://github.com/Jerryzhang158/mitagent/issues) with component information
- **Feature Requests**: [Submit enhancement requests](https://github.com/Jerryzhang158/mitagent/issues)
- **Technical Support**: Include pipeline logs and configuration details

---

**Enhanced Modular Architecture - More Reliable, Flexible, and Extensible!**

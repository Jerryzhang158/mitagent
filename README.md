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

2. **Create conda environment**
   ```bash
   conda env create -f environment.yml
   conda activate mirna-pipeline
   ```

3. **Install Ollama** (for LLM analysis)
   ```bash
   curl -fsSL https://ollama.ai/install.sh | sh
   ollama serve
   ollama pull llama3.1
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

The enhanced pipeline allows you to skip any combination of steps:

```bash
# Skip literature mining and BERT, proceed directly to LLM analysis
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --skip-literature \
    --skip-bert \
    --analysis-mode both \
    --generate-network
```

### Network Generation with Auto-Fallback

The pipeline now automatically generates default expression files if DESeq2 results aren't provided:

```bash
# Network generation with auto-generated expression files
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
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
└── rnaseq_analyzer.py       # RNA-seq differential expression
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
| `--skip-*` | Can skip any combination | - | Pipeline continues gracefully |
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
│   └── download_validation.json          # 🆕 File validation results  
├── bert_validation/
│   ├── mti_validation_results.csv
│   ├── integrated_validation_results.csv  # 🆕 Enhanced integration
│   └── batch_input.csv                    # 🆕 Batch processing data
├── llm_summaries/                         # 🆕 Dual-mode analysis
│   ├── standard_summaries/               # Standard relationship analysis
│   └── functional_analysis/              # Pathway and systems analysis
├── cytoscape_network/                     # 🆕 Enhanced with auto-fallback
│   ├── network.sif
│   ├── node_attributes.txt
│   ├── edge_attributes.txt
│   └── network_statistics.json          # 🆕 Detailed network stats
├── default_mirna_deseq.csv              # 🆕 Auto-generated when needed
└── default_gene_deseq.csv               # 🆕 Auto-generated when needed
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

### Example 3: LLM-Only Analysis Pipeline

```bash
# For when you have literature but want to focus on LLM analysis
python main_pipeline.py \
    --mode combined \
    --existing-step1 previous_mti_results.xlsx \
    --skip-literature \
    --skip-bert \
    --analysis-mode both \
    --functions "proliferation,apoptosis,invasion,metastasis" \
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

### Environment Variables

Set these environment variables for enhanced functionality:

```bash
export OLLAMA_HOST=http://localhost:11434  # Custom Ollama endpoint
export MAX_BERT_BATCH_SIZE=50             # Batch processing control  
export LOG_LEVEL=INFO                     # Logging verbosity
```

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
| BERT validation memory error | Skip BERT, use LLM only | `--skip-bert --analysis-mode both` |
| Network generation failure | Check score threshold | `--network-score-threshold 40` |
| Ollama connection issues | Check service status | `ollama serve && ollama pull llama3.1` |

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
python main_pipeline.py --mode combined --genes genes.txt --mirnas mirnas.txt --analysis-mode both
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

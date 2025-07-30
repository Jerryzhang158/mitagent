# mitagent ——— An Enhanced Modular miRNA Research agent

## Overview

This pipeline provides an end-to-end solution for miRNA-mRNA interaction analysis, integrating multiple prediction databases, machine learning validation, and advanced functional interpretation. The system supports flexible analytical workflows for gene-to-miRNA, miRNA-to-gene, and bidirectional analysis.

## Key Features

| Feature | Description |
|---------|-------------|
| **Multi-Database Integration** | TargetScan, miRDB, miRWalk, miRTarBase support |
| **Machine Learning Validation** | BERT-based interaction validation |
| **Literature Mining** | Automated PubMed literature retrieval |
| **LLM Analysis** | Advanced functional interpretation |
| **Network Visualization** | Cytoscape-compatible network generation |
| **Modular Architecture** | Component-based design for extensibility |

---

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
   cd mirna-pipeline
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
    --mode gene \
    --genes input/genes.txt \
    --targetscan data/targetscan.txt \
    --mirdb data/mirdb.txt \
    --mirwalk data/mirwalk.txt
```

---

## Table of Contents

- [Installation](#installation)
- [Input Files](#input-files)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output](#output)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)

---

## Installation

### System Requirements

| Component | Requirement |
|-----------|-------------|
| **Operating System** | Linux, macOS, Windows (WSL recommended) |
| **Python Version** | 3.8 or higher |
| **Memory** | 8GB minimum, 16GB recommended |
| **Storage** | 2GB available space |
| **Network** | Internet connection for downloads |

### Step-by-Step Installation

#### 1. Environment Setup

```bash
# Clone repository
git clone https://github.com/username/mirna-pipeline.git
cd mirna-pipeline

# Create and activate conda environment
conda env create -f environment.yml
conda activate mirna-pipeline
```

#### 2. External Dependencies

**Ollama (Required for LLM features)**
```bash
# Install Ollama
curl -fsSL https://ollama.ai/install.sh | sh

# Start service
ollama serve

# Download language model
ollama pull llama3.1
```

**Cytoscape (Optional for visualization)**
- Download from: https://cytoscape.org/

#### 3. Verify Installation

```bash
# Test pipeline
python main_pipeline.py --help

# Check Ollama status
ollama list
```

---

## Input Files

### Required Database Files

The pipeline requires miRNA target prediction databases in tab-delimited format:

| Database | Description | Required |
|----------|-------------|----------|
| **TargetScan** | Conserved target predictions | Yes |
| **miRDB** | ML-based predictions | Yes |
| **miRWalk** | Comprehensive database | Yes |
| **miRTarBase** | Experimental validation | Optional |

### Input Data Formats

#### Gene List Format
Create a text file with one gene symbol per line:

```
TP53
BRCA1
EGFR
MYC
KRAS
```

#### miRNA List Format
Create a text file with one miRNA identifier per line:

```
hsa-miR-21-5p
hsa-miR-155-5p
hsa-miR-200c-3p
hsa-miR-34a-5p
```

#### Differential Expression Format (Optional)
CSV format with standard DESeq2 output columns:

```csv
GENE,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
TP53,1500.2,2.1,0.3,7.0,1.2e-12,3.4e-11
BRCA1,890.5,-1.8,0.4,-4.5,6.7e-06,2.1e-05
```

---

## Usage

### Command Structure

```bash
python main_pipeline.py [REQUIRED OPTIONS] [OPTIONAL PARAMETERS]
```

### Required Parameters

| Parameter | Description | Values |
|-----------|-------------|--------|
| `-m, --mode` | Analysis mode | `gene`, `mirna`, `combined` |
| `-t, --targetscan` | TargetScan database file | File path |
| `-d, --mirdb` | miRDB database file | File path |
| `-w, --mirwalk` | miRWalk database file | File path |

### Input Parameters

| Parameter | Description | Required For |
|-----------|-------------|--------------|
| `-g, --genes` | Gene list file | `gene`, `combined` modes |
| `-r, --mirnas` | miRNA list file | `mirna`, `combined` modes |
| `-b, --mirtarbase` | miRTarBase file | Optional (all modes) |

### Analysis Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min-databases` | Minimum database support | `2` |
| `--max-articles` | Max articles per MTI | `50` |
| `-f, --functions` | Biological functions | `"proliferation,apoptosis,migration,invasion"` |
| `--analysis-mode` | LLM analysis type | `both` |

### Pipeline Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip-literature` | Skip literature mining | `False` |
| `--skip-bert` | Skip BERT validation | `False` |
| `--skip-llm` | Skip LLM analysis | `False` |
| `--generate-network` | Create network files | `False` |

### Network Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mirna-deseq` | miRNA expression file | Auto-generated |
| `--gene-deseq` | Gene expression file | Auto-generated |
| `--network-score-threshold` | Network inclusion threshold | `50` |
| `-o, --output` | Output directory | `pipeline_results` |

---

## Pipeline Steps

### Step 1: MTI Selection
**Purpose**: Identify miRNA-mRNA interactions from multiple databases
- Multi-database consensus filtering
- miRTarBase validation integration
- Configurable support thresholds

### Step 2: Literature Mining
**Purpose**: Retrieve relevant scientific literature
- Automated PubMed queries
- Rate-limited API requests
- Literature validation and curation

### Step 3: BERT Validation
**Purpose**: Machine learning-based interaction validation
- Pre-trained language model analysis
- Multi-functional relevance scoring
- Quantitative evidence assessment

### Step 4: LLM Analysis
**Purpose**: Advanced functional interpretation

#### Standard Analysis
- Relationship summarization
- Evidence strength evaluation
- Clinical relevance assessment
- Research gap identification

#### Functional Analysis
- Pathway enrichment
- Disease association mapping
- Therapeutic potential assessment
- Systems-level interpretation

### Step 5: Network Generation
**Purpose**: Create visualization-ready network files
- Cytoscape-compatible formats
- Differential expression integration
- Network topology analysis

---

## Output

### Directory Structure

```
pipeline_results/
└── run_YYYYMMDD_HHMMSS/
    ├── 📄 pipeline_report.txt
    ├── 📄 pipeline_log.txt
    ├── 📁 mirna_selection/
    │   └── mti_selection_results_[mode].xlsx
    ├── 📁 pubmed_articles/
    │   └── [literature files]
    ├── 📁 bert_validation/
    │   ├── mti_validation_results.csv
    │   └── integrated_validation_results.csv
    ├── 📁 llm_summaries/
    │   ├── standard_summaries/
    │   └── functional_analysis/
    └── 📁 cytoscape_network/
        ├── network.sif
        ├── node_attributes.txt
        └── edge_attributes.txt
```

### Key Output Files

| File | Description | Usage |
|------|-------------|-------|
| `integrated_validation_results.csv` | **Main results file** with all analysis data | Primary analysis output |
| `mti_selection_results_[mode].xlsx` | Initial MTI selection with database info | Database integration results |
| `pipeline_report.txt` | Comprehensive execution summary | Analysis overview |
| `network.sif` | Cytoscape network file | Network visualization |

---

## Examples

### Example 1: Gene-to-miRNA Analysis

```bash
python main_pipeline.py \
    --mode gene \
    --genes examples/cancer_genes.txt \
    --targetscan data/targetscan_v8.txt \
    --mirdb data/mirdb_v6.txt \
    --mirwalk data/mirwalk_v3.txt \
    --mirtarbase data/mirtarbase_v9.txt \
    --functions "proliferation,apoptosis,invasion" \
    --max-articles 30 \
    --output results/cancer_analysis
```

### Example 2: Comprehensive Analysis with Networks

```bash
python main_pipeline.py \
    --mode combined \
    --genes examples/target_genes.txt \
    --mirnas examples/candidate_mirnas.txt \
    --targetscan data/targetscan_v8.txt \
    --mirdb data/mirdb_v6.txt \
    --mirwalk data/mirwalk_v3.txt \
    --generate-network \
    --mirna-deseq data/mirna_deseq2_results.csv \
    --gene-deseq data/gene_deseq2_results.csv \
    --analysis-mode both \
    --min-databases 3 \
    --network-score-threshold 60 \
    --output results/comprehensive_study
```

### Example 3: Quick Analysis (Skip Time-Intensive Steps)

```bash
python main_pipeline.py \
    --mode mirna \
    --mirnas examples/mirna_list.txt \
    --targetscan data/targetscan_v8.txt \
    --mirdb data/mirdb_v6.txt \
    --mirwalk data/mirwalk_v3.txt \
    --skip-literature \
    --skip-llm \
    --max-articles 10 \
    --output results/quick_analysis
```

---

## Configuration

### Biological Functions

Specify custom biological functions for analysis:

```bash
# Cancer-related functions
--functions "proliferation,apoptosis,invasion,metastasis,angiogenesis"

# Metabolic functions
--functions "glucose metabolism,lipid metabolism,energy production"

# Developmental functions
--functions "differentiation,development,morphogenesis"
```

### Analysis Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `standard` | MTI relationship analysis | Basic interaction validation |
| `functional` | Pathway and systems analysis | Comprehensive functional interpretation |
| `both` | Combined analysis | Complete analysis (recommended) |

### Performance Tuning

| Parameter | Small Dataset | Large Dataset |
|-----------|---------------|---------------|
| `--max-articles` | 50 | 20 |
| `--min-databases` | 2 | 3 |
| Memory allocation | 8GB | 16GB+ |

---

## Troubleshooting

### Common Issues

#### 🔧 Ollama Connection Problems

**Problem**: `Cannot connect to Ollama`

**Solution**:
```bash
# Check status
ollama list

# Restart service
ollama serve

# Verify model
ollama pull llama3.1
```

#### 🔧 Memory Issues

**Problem**: `Out of memory error`

**Solutions**:
- Reduce `--max-articles` parameter
- Process smaller gene/miRNA lists
- Increase system swap space
- Use `--skip-llm` for memory-constrained systems

#### 🔧 Literature Retrieval Failures

**Problem**: `PubMed connection timeout`

**Solutions**:
- Check internet connection
- Verify PubMed service status
- Retry with smaller article limits
- Use `--skip-literature` if persistent

#### 🔧 Database File Errors

**Problem**: `Database file format error`

**Solutions**:
- Verify file format matches database specifications
- Check file encoding (UTF-8 recommended)
- Ensure proper column headers
- Validate file completeness

### Getting Help

1. **Check the logs**: `tail -f pipeline_results/run_*/pipeline_log.txt`
2. **Verify inputs**: Ensure all required files exist and are properly formatted
3. **Test components**: Run individual steps to isolate issues
4. **Resource monitoring**: Check system memory and disk space

---

## Performance

### Computational Complexity

| Step | Complexity | Typical Runtime |
|------|------------|-----------------|
| MTI Selection | O(n) | Minutes |
| Literature Mining | O(n×m) | Hours |
| BERT Validation | O(n×m×k) | Hours |
| LLM Analysis | O(n×k) | Hours |
| Network Generation | O(n²) | Minutes |

*Where n = MTI count, m = articles per MTI, k = analysis complexity*

### Optimization Tips

- **Parallel Processing**: Pipeline automatically uses available CPU cores
- **Caching**: Literature and model results are cached for reuse
- **Batch Processing**: BERT validation uses efficient batch processing
- **Memory Management**: Automatic garbage collection for large datasets

---

## Contributing

### Development Setup

```bash
# Fork and clone
git clone https://github.com/yourusername/mirna-pipeline.git
cd mirna-pipeline

# Create development environment
conda env create -f environment-dev.yml
conda activate mirna-pipeline-dev

# Install in development mode
pip install -e .
```

### Code Standards

- **Style**: Follow PEP 8 guidelines
- **Documentation**: Comprehensive docstrings required
- **Testing**: Unit tests for new functionality
- **Commits**: Descriptive commit messages

### Contribution Process

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/new-analysis`)
3. **Develop** with appropriate testing
4. **Document** changes and new features
5. **Submit** a pull request with detailed description

### Areas for Contribution

- Additional database integrations
- New validation methods
- Performance optimizations
- Visualization enhancements
- Documentation improvements

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

```bibtex
@software{mirna_pipeline_2024,
  title={Enhanced Modular miRNA Research Pipeline},
  author={[Author Names]},
  year={2024},
  url={https://github.com/username/mirna-pipeline},
  version={1.0}
}
```

## Acknowledgments

This pipeline builds upon several established resources:

- **TargetScan** - Conserved miRNA target predictions
- **miRDB** - Machine learning-based target predictions  
- **miRWalk** - Comprehensive miRNA target database
- **miRTarBase** - Experimentally validated interactions
- **Hugging Face** - Pre-trained language models
- **Ollama** - Local LLM deployment framework

---

## Support

- **📖 Documentation**: Check this README and inline documentation
- **🐛 Bug Reports**: [Open an issue](https://github.com/username/mirna-pipeline/issues)
- **💡 Feature Requests**: [Submit a feature request](https://github.com/username/mirna-pipeline/issues)
- **📧 Contact**: your.email@institution.edu

---

**⭐ Star this repository if it helps your research!**

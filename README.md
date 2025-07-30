# MiTAgent
> A comprehensive computational pipeline for systematic analysis of microRNA regulatory networks

## Overview

MiTAgent provides an end-to-end solution for miRNA-mRNA interaction analysis, integrating multiple prediction databases, machine learning validation, and advanced functional interpretation. The system supports flexible analytical workflows for gene-to-miRNA, miRNA-to-gene, and bidirectional analysis.

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

### Step-by-Step Installation

#### 1. Environment Setup

```bash
# Clone repository
git clone https://github.com/Jerryzhang158/mitagent.git
cd mitagent

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

Download the following database files and place them in your working directory:

| Database | Filename | Description | Source |
|----------|----------|-------------|--------|
| **TargetScan** | `Predicted_Targets_Context_Scores.default_predictions.txt` | Conserved target predictions | TargetScan v8.0 |
| **miRDB** | `miRDB_v6.0_prediction_result_fixed.txt` | ML-based predictions | miRDB v6.0 |
| **miRWalk** | `hsa_miRWalk_3UTR.txt` | Comprehensive database | miRWalk 3.0 |
| **miRTarBase** | `miRTarBase_MTI_fixed.csv` | Experimental validation | miRTarBase v9.0 |

### Input Data Formats

#### Gene List Format (`gene_list.txt`)
Create a text file with one gene symbol per line:

```
TP53
BRCA1
EGFR
MYC
KRAS
AKT1
PIK3CA
PTEN
```

#### miRNA List Format (`mirna_list.txt`)
Create a text file with one miRNA identifier per line:

```
hsa-miR-21-5p
hsa-miR-155-5p
hsa-miR-200c-3p
hsa-miR-34a-5p
hsa-miR-125b-5p
hsa-miR-let-7a-5p
```

#### Differential Expression Formats

**miRNA Expression** (`miRNA.csv`):
```csv
miRNA,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
hsa-miR-21-5p,2500.3,1.8,0.2,9.0,2.5e-19,4.1e-18
hsa-miR-155-5p,1800.7,2.3,0.3,7.7,1.3e-14,1.8e-13
hsa-miR-200c-3p,950.2,-1.5,0.4,-3.8,1.4e-04,8.2e-04
```

**Gene Expression** (`mRNA.csv`):
```csv
mRNA,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
TP53,1500.2,2.1,0.3,7.0,1.2e-12,3.4e-11
BRCA1,890.5,-1.8,0.4,-4.5,6.7e-06,2.1e-05
EGFR,2100.8,1.2,0.2,6.0,2.0e-09,3.5e-08
```

---

## Usage

### Command Structure

```bash
python main_pipeline.py [REQUIRED OPTIONS] [OPTIONAL PARAMETERS]
```

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-m, --mode` | Analysis mode | `gene`, `mirna`, `combined` |
| `-t, --targetscan` | TargetScan database file | `Predicted_Targets_Context_Scores.default_predictions.txt` |
| `-d, --mirdb` | miRDB database file | `miRDB_v6.0_prediction_result_fixed.txt` |
| `-w, --mirwalk` | miRWalk database file | `hsa_miRWalk_3UTR.txt` |

### Input Parameters

| Parameter | Description | Required For |
|-----------|-------------|--------------|
| `-g, --genes` | Gene list file | `gene`, `combined` modes |
| `-r, --mirnas` | miRNA list file | `mirna`, `combined` modes |
| `-b, --mirtarbase` | miRTarBase file | Optional (all modes) |

### Analysis Options

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `--min-databases` | Minimum database support | `2` | `3` for high confidence |
| `--max-articles` | Max articles per MTI | `50` | `100-150` for comprehensive analysis |
| `-f, --functions` | Biological functions | `"proliferation,apoptosis,migration,invasion"` | Customize based on research |
| `--analysis-mode` | LLM analysis type | `both` | `both` for complete analysis |

### Pipeline Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip-literature` | Skip literature mining | `False` |
| `--skip-bert` | Skip BERT validation | `False` |
| `--skip-llm` | Skip LLM analysis | `False` |
| `--generate-network` | Create network files | `False` |

### Network Options

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--mirna-deseq` | miRNA expression file | Auto-generated | `miRNA.csv` |
| `--gene-deseq` | Gene expression file | Auto-generated | `mRNA.csv` |
| `--network-score-threshold` | Network inclusion threshold | `50` | `60` for high confidence |
| `-o, --output` | Output directory | `pipeline_results` | Custom path |

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
    ├── pipeline_report.txt
    ├── pipeline_log.txt
    ├── mirna_selection/
    │   └── mti_selection_results_combined.xlsx
    ├── pubmed_articles/
    │   └── [literature files]
    ├── bert_validation/
    │   ├── mti_validation_results.csv
    │   └── integrated_validation_results.csv
    ├── llm_summaries/
    │   ├── standard_summaries/
    │   └── functional_analysis/
    └── cytoscape_network/
        ├── network.sif
        ├── node_attributes.txt
        └── edge_attributes.txt
```

### Key Output Files

| File | Description | Usage |
|------|-------------|-------|
| `integrated_validation_results.csv` | **Main results file** with all analysis data | Primary analysis output |
| `mti_selection_results_combined.xlsx` | Initial MTI selection with database info | Database integration results |
| `pipeline_report.txt` | Comprehensive execution summary | Analysis overview |
| `network.sif` | Cytoscape network file | Network visualization |

---

## Examples

### Example 1: Complete Analysis (Recommended)

Based on the actual usage pattern, here's the comprehensive analysis command:

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
    --functions "cell proliferation,apoptosis" \
    --max-articles 150 \
    --min-databases 3 \
    --analysis-mode both \
    --generate-network \
    --network-score-threshold 60
```

### Example 2: Gene-to-miRNA Analysis

```bash
python main_pipeline.py \
    --mode gene \
    --genes gene_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --functions "proliferation,apoptosis,invasion,metastasis" \
    --max-articles 100 \
    --min-databases 2 \
    --output results/gene_to_mirna_analysis
```

### Example 3: High-Throughput Analysis

For large-scale studies with extensive gene/miRNA lists:

```bash
python main_pipeline.py \
    --mode combined \
    --genes large_gene_list.txt \
    --mirnas large_mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --mirtarbase miRTarBase_MTI_fixed.csv \
    --functions "cell proliferation,apoptosis,migration,invasion,angiogenesis" \
    --max-articles 200 \
    --min-databases 3 \
    --analysis-mode both \
    --generate-network \
    --network-score-threshold 70 \
    --output results/comprehensive_study
```

### Example 4: Quick Analysis (Time-Limited)

For rapid preliminary analysis:

```bash
python main_pipeline.py \
    --mode combined \
    --genes gene_list.txt \
    --mirnas mirna_list.txt \
    --targetscan Predicted_Targets_Context_Scores.default_predictions.txt \
    --mirdb miRDB_v6.0_prediction_result_fixed.txt \
    --mirwalk hsa_miRWalk_3UTR.txt \
    --skip-literature \
    --skip-llm \
    --min-databases 2 \
    --generate-network \
    --output results/quick_analysis
```

---

## Configuration

### Biological Functions

Customize biological functions based on your research focus:

```bash
# Cancer research
--functions "cell proliferation,apoptosis,invasion,metastasis,angiogenesis,drug resistance"

# Metabolic research
--functions "glucose metabolism,lipid metabolism,energy production,insulin signaling"

# Neurological research
--functions "neurogenesis,synaptic plasticity,neurodegeneration,memory formation"

# Cardiovascular research
--functions "cardiac development,angiogenesis,hypertrophy,fibrosis"
```

### Analysis Modes

| Mode | Description | Processing Time | Use Case |
|------|-------------|-----------------|----------|
| `standard` | MTI relationship analysis | Fast | Basic interaction validation |
| `functional` | Pathway and systems analysis | Medium | Comprehensive functional interpretation |
| `both` | Combined analysis | Slow | Complete analysis (recommended) |

### Performance Tuning

| Parameter | Small Dataset (<50 MTIs) | Medium Dataset (50-200 MTIs) | Large Dataset (>200 MTIs) |
|-----------|---------------------------|------------------------------|---------------------------|
| `--max-articles` | 150 | 100 | 50 |
| `--min-databases` | 2 | 3 | 3 |
| `--network-score-threshold` | 50 | 60 | 70 |

---

## Troubleshooting

### Common Issues

#### Ollama Connection Problems

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

#### Memory Issues

**Problem**: `Out of memory error during BERT validation`

**Solutions**:
- Reduce `--max-articles` to 50-100
- Process smaller gene/miRNA lists (10-20 items per run)
- Increase system swap space
- Use `--skip-llm` for memory-constrained systems

#### Database File Issues

**Problem**: `Database file format error`

**Solutions**:
- Verify file names match exactly:
  - `Predicted_Targets_Context_Scores.default_predictions.txt`
  - `miRDB_v6.0_prediction_result_fixed.txt`
  - `hsa_miRWalk_3UTR.txt`
  - `miRTarBase_MTI_fixed.csv`
- Check file encoding (UTF-8 recommended)
- Ensure files are not corrupted during download

#### Literature Retrieval Failures

**Problem**: `PubMed connection timeout or rate limiting`

**Solutions**:
- Reduce `--max-articles` to 50-100
- Check internet connection stability
- Retry analysis after waiting period
- Use `--skip-literature` if persistent issues

#### Network Generation Failures

**Problem**: `No network generated or empty network file`

**Solutions**:
- Lower `--network-score-threshold` (try 40-50)
- Verify differential expression files format
- Ensure sufficient validated MTIs for network construction
- Check that `--generate-network` flag is included

### Performance Monitoring

Monitor system resources during analysis:

```bash
# Check memory usage
free -h

# Monitor CPU usage
top

# Check disk space
df -h

# View pipeline logs in real-time
tail -f pipeline_results/run_*/pipeline_log.txt
```

---

## Data Sources and Versions

### Database Information

| Database | Version | Last Updated | Download Source |
|----------|---------|--------------|-----------------|
| **TargetScan** | v8.0 | 2021 | targetscan.org |
| **miRDB** | v6.0 | 2020 | mirdb.org |
| **miRWalk** | v3.0 | 2019 | mirwalk.umm.uni-heidelberg.de |
| **miRTarBase** | v9.0 | 2022 | mirtarbase.cuhk.edu.cn |

### File Specifications

The pipeline expects specific file formats. Contact the repository maintainer if you need assistance with file format conversion.

---

## Performance Benchmarks

### Typical Runtime (Combined Mode)

| Dataset Size | Literature Mining | BERT Validation | LLM Analysis | Total Time |
|--------------|-------------------|-----------------|--------------|------------|
| 10 genes + 10 miRNAs | 30 min | 45 min | 60 min | ~2.5 hours |
| 25 genes + 25 miRNAs | 90 min | 120 min | 180 min | ~6.5 hours |
| 50 genes + 50 miRNAs | 180 min | 240 min | 360 min | ~13 hours |

*Times based on max-articles=150, min-databases=3, analysis-mode=both*

---

## Contributing

### Development Setup

```bash
# Fork and clone
git clone https://github.com/yourusername/mitagent.git
cd mitagent

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

### Contribution Areas

- Additional database integrations
- Performance optimizations
- New analysis methods
- Visualization enhancements
- Documentation improvements

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

```bibtex
@software{mitagent_2024,
  title={MiTAgent: Enhanced Modular miRNA Research Pipeline},
  author={Jerry Zhang},
  year={2024},
  url={https://github.com/Jerryzhang158/mitagent},
  version={1.0}
}
```

## Acknowledgments

This pipeline builds upon several established resources:

- **TargetScan** - Conserved miRNA target predictions (Agarwal et al., 2015)
- **miRDB** - Machine learning-based target predictions (Chen & Wang, 2020)
- **miRWalk** - Comprehensive miRNA target database (Sticht et al., 2018)
- **miRTarBase** - Experimentally validated interactions (Huang et al., 2020)
- **Hugging Face** - Pre-trained language models
- **Ollama** - Local LLM deployment framework

---

## Support

- **Documentation**: Check this README and inline documentation
- **Bug Reports**: [Open an issue](https://github.com/Jerryzhang158/mitagent/issues)
- **Feature Requests**: [Submit a feature request](https://github.com/Jerryzhang158/mitagent/issues)
- **Contact**: Open an issue for questions and support

---

**Star this repository if it helps your research!**

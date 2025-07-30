# mitagent
A comprehensive computational pipeline for systematic analysis of microRNA regulatory networks. This modular system integrates multiple prediction databases, implements machine learning-based validation through pre-trained language models, and provides advanced functional analysis capabilities for miRNA-mRNA interactions.
Features

Multi-modal Analysis: Supports gene-to-miRNA, miRNA-to-gene, and bidirectional combined analysis workflows
Database Integration: Incorporates TargetScan, miRDB, miRWalk, and miRTarBase prediction databases
Literature Mining: Automated PubMed literature retrieval and analysis
BERT Validation: Machine learning-based interaction validation using pre-trained language models
LLM Analysis: Advanced functional interpretation using large language models
Network Visualization: Cytoscape-compatible network generation with differential expression integration
Modular Architecture: Component-based design for enhanced maintainability and extensibility

Table of Contents

Prerequisites
Installation
Input Requirements
Usage
Pipeline Workflow
Output Structure
Configuration
Troubleshooting
Contributing

Prerequisites
System Requirements

Python 3.8 or higher
Minimum 8GB RAM (16GB recommended for large datasets)
2GB available disk space
Internet connectivity for literature retrieval and model downloads

External Dependencies
Ollama Framework (Required for LLM analysis)
bashcurl -fsSL https://ollama.ai/install.sh | sh
ollama serve
ollama pull llama3.1
Cytoscape (Optional, for network visualization)
Download from: https://cytoscape.org/
Installation
Environment Setup
Create and activate the conda environment:
bash# Clone the repository
git clone https://github.com/username/mirna-pipeline.git
cd mirna-pipeline

# Create conda environment from provided file
conda env create -f environment.yml

# Activate the environment
conda activate mirna-pipeline
Verify Installation
bash# Test basic functionality
python main_pipeline.py --help

# Verify Ollama connection (if using LLM features)
ollama list
Input Requirements
Database Files
The pipeline requires access to miRNA target prediction databases in tab-delimited format:

TargetScan: Conserved miRNA target predictions
miRDB: Machine learning-based predictions
miRWalk: Comprehensive target database
miRTarBase: Experimentally validated interactions (optional)

Input File Formats
Gene List (gene_list.txt):
TP53
BRCA1
EGFR
MYC
KRAS
miRNA List (mirna_list.txt):
hsa-miR-21-5p
hsa-miR-155-5p
hsa-miR-200c-3p
hsa-miR-34a-5p
Differential Expression Files (optional, for network analysis):
GENE,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
TP53,1500.2,2.1,0.3,7.0,1.2e-12,3.4e-11
BRCA1,890.5,-1.8,0.4,-4.5,6.7e-06,2.1e-05
Usage
Basic Command Structure
bashpython main_pipeline.py [OPTIONS]
Required Parameters
ParameterDescriptionExample-m, --modeAnalysis mode: gene, mirna, or combined-m gene-t, --targetscanTargetScan database file path-t data/targetscan.txt-d, --mirdbmiRDB database file path-d data/mirdb.txt-w, --mirwalkmiRWalk database file path-w data/mirwalk.txt
Input File Parameters
ParameterDescriptionRequired For-g, --genesGene identifier list filegene, combined modes-r, --mirnasmiRNA identifier list filemirna, combined modes-b, --mirtarbasemiRTarBase validation fileOptional for all modes
Analysis Configuration
ParameterDescriptionDefault--min-databasesMinimum database support required2--max-articlesMaximum literature articles per MTI50-f, --functionsBiological functions for analysis"cell proliferation,apoptosis,migration,invasion"--analysis-modeLLM analysis typeboth
Pipeline Control
ParameterDescriptionDefault--skip-literatureBypass literature miningFalse--skip-bertBypass BERT validationFalse--skip-llmBypass LLM analysisFalse--generate-networkEnable network generationFalse
Network Analysis Options
ParameterDescriptionDefault--mirna-deseqmiRNA differential expression fileAuto-generated--gene-deseqGene differential expression fileAuto-generated--network-score-thresholdMinimum score for network inclusion50-o, --outputOutput directorypipeline_results
Example Commands
Gene-to-miRNA Analysis:
bashpython main_pipeline.py \
    --mode gene \
    --genes input/target_genes.txt \
    --targetscan data/targetscan_predictions.txt \
    --mirdb data/mirdb_predictions.txt \
    --mirwalk data/mirwalk_predictions.txt \
    --mirtarbase data/mirtarbase_validated.txt \
    --functions "proliferation,apoptosis,invasion" \
    --max-articles 30
Comprehensive Analysis with Network Generation:
bashpython main_pipeline.py \
    --mode combined \
    --genes input/genes.txt \
    --mirnas input/mirnas.txt \
    --targetscan data/targetscan.txt \
    --mirdb data/mirdb.txt \
    --mirwalk data/mirwalk.txt \
    --generate-network \
    --mirna-deseq input/mirna_deseq.csv \
    --gene-deseq input/gene_deseq.csv \
    --analysis-mode both \
    --min-databases 3 \
    --output results/comprehensive_analysis
Pipeline Workflow
Step 1: MTI Selection
Identifies miRNA-mRNA interactions from multiple databases using consensus-based filtering. Interactions are retained based on minimum database support requirements and integrated with experimental validation data from miRTarBase when available.
Step 2: Literature Mining
Performs automated literature retrieval from PubMed using structured queries combining miRNA and target gene identifiers. Implements rate limiting and error handling for robust data acquisition.
Step 3: BERT Validation
Utilizes pre-trained BERT models to analyze retrieved literature and assess biological relevance of predicted interactions. Evaluates multiple functional categories simultaneously with quantitative scoring.
Step 4: LLM Analysis
Employs large language models for comprehensive functional analysis:

Standard Analysis: Relationship summarization, evidence assessment, clinical relevance
Functional Analysis: Pathway interpretation, disease associations, therapeutic potential

Step 5: Network Generation
Constructs Cytoscape-compatible network files integrating differential expression data. Provides network topology analysis and identification of regulatory hubs.
Output Structure
pipeline_results/
└── run_YYYYMMDD_HHMMSS/
    ├── pipeline_report.txt                    # Comprehensive analysis report
    ├── pipeline_log.txt                       # Execution log
    ├── mirna_selection/
    │   └── mti_selection_results_[mode].xlsx  # Initial MTI selection
    ├── pubmed_articles/                       # Literature abstracts
    ├── bert_validation/
    │   ├── mti_validation_results.csv         # BERT validation scores
    │   └── integrated_validation_results.csv  # Complete analysis results
    ├── llm_summaries/
    │   ├── standard_summaries/                # Standard LLM analysis
    │   └── functional_analysis/               # Functional LLM analysis
    └── cytoscape_network/                     # Network visualization files
        ├── network.sif                        # Network structure
        ├── node_attributes.txt                # Node properties
        └── edge_attributes.txt                # Edge properties
Key Output Files
FileDescriptionintegrated_validation_results.csvComplete analysis results with all validation scores and LLM summariesmti_selection_results_[mode].xlsxInitial MTI selection with database support informationpipeline_report.txtComprehensive execution summary and statisticsnetwork.sifCytoscape-compatible network file for visualization
Configuration
Biological Functions
Default functions analyzed include cell proliferation, apoptosis, migration, and invasion. Custom functions can be specified:
bashpython main_pipeline.py -f "angiogenesis,metastasis,drug resistance,differentiation" ...
LLM Analysis Modes

standard: Focuses on MTI relationship summarization and evidence assessment
functional: Emphasizes pathway analysis and functional interpretation
both: Performs comprehensive analysis combining both approaches (recommended)

Database Support Thresholds
The --min-databases parameter controls the minimum number of databases required to support an interaction for inclusion in downstream analysis. Higher values increase specificity but may reduce sensitivity.
Troubleshooting
Common Issues and Solutions
Ollama Connection Failures
bash# Check service status
ollama list

# Restart Ollama service
ollama serve

# Verify model availability
ollama pull llama3.1
Memory Limitations

Reduce --max-articles parameter for literature mining
Process large gene/miRNA lists in smaller batches
Ensure sufficient system memory allocation

Literature Retrieval Issues

Verify internet connectivity and PubMed accessibility
Check for rate limiting or temporary service unavailability
Review query formatting in pipeline logs

BERT Model Loading Problems
bash# Set Hugging Face endpoint if needed
export HF_ENDPOINT=https://hf-mirror.com

# Clear model cache
rm -rf ~/.cache/huggingface/transformers/
Network Generation Failures

Ensure differential expression files follow required format
Verify adequate validation results for network construction
Check score threshold settings

Log Analysis
Detailed execution information is available in the pipeline log:
bashtail -f pipeline_results/run_*/pipeline_log.txt
Performance Considerations
Computational Complexity

Literature Mining: O(n) where n = number of MTIs
BERT Validation: O(n×m) where m = average articles per MTI
LLM Analysis: O(n×k) where k = analysis complexity factor
Network Generation: O(n²) for network construction

Optimization Strategies

Parallel processing implementation for literature retrieval
Batch processing for BERT validation to optimize GPU utilization
Caching mechanisms for repeated analyses
Memory-efficient data structures for large-scale datasets

Contributing
Development Guidelines

Follow PEP 8 coding standards
Implement comprehensive documentation for new functions
Include unit tests for new functionality
Maintain backward compatibility when possible

Contribution Process

Fork the repository
Create a feature branch (git checkout -b feature/new-analysis)
Implement changes with appropriate testing
Submit a pull request with detailed description

Code Structure
The pipeline follows a modular architecture with separate components for each analysis step:

config.py: Configuration management
utils.py: Utility functions and helper classes
mirna_matcher.py: MTI identification and matching
mti_selection.py: Database integration and filtering
literature_mining.py: PubMed retrieval and processing
mirna_bert.py: BERT-based validation
mti_llm_summarize.py: LLM analysis components
mti_cytoscape_network.py: Network generation

License
This project is licensed under the MIT License. See the LICENSE file for details.
Acknowledgments
This pipeline integrates data and methodologies from several established resources:

TargetScan (Lewis et al., Cell 2005): Conserved miRNA target predictions
miRDB (Chen and Wang, Nucleic Acids Research 2020): Machine learning-based predictions
miRWalk (Sticht et al., Nucleic Acids Research 2018): Comprehensive target database
miRTarBase (Huang et al., Nucleic Acids Research 2020): Experimentally validated interactions
Hugging Face Transformers: Pre-trained language model framework
Ollama: Local large language model deployment platform

Citation
If you use this pipeline in your research, please cite:
bibtex@software{mirna_pipeline_2024,
  title={Enhanced Modular miRNA Research Pipeline},
  author={[Author Names]},
  year={2024},
  url={https://github.com/username/mirna-pipeline},
  version={1.0}

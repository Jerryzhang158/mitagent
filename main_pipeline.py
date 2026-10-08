"""
miRNA Research Pipeline - 主控制器
整合所有模块，执行完整的分析流程
"""

import os
import json
import platform
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
import sys
import argparse
import subprocess
import pandas as pd
from datetime import datetime
from typing import Optional, List, Dict, Any

# 导入自定义模块
from config import PipelineConfig
from utils import Logger, FileUtils, ResultsIntegrator, TextProcessor, check_ollama_connection, create_timestamp
from mirna_matcher import MiRNAMatcher
from refseq_cache import RefSeqManager
from mti_selection import MTISelector
from literature_mining import LiteratureMiner

# 导入已存在的模块
try:
    from mti_llm_summarize import MTILLMSummarizer
except ImportError:
    print("Optional LLM dependencies are unavailable; install requirements-optional.txt to enable interpretation.")
    MTILLMSummarizer = None

try:
    from mti_cytoscape_network import MTICytoscapeGenerator
except ImportError:
    print("Warning: mti_cytoscape_network.py not found. Network generation will be unavailable.")
    MTICytoscapeGenerator = None

class MiRNAPipeline:
    """miRNA研究流水线主控制器"""
    
    def __init__(self, output_dir: str = "pipeline_results", min_databases: int = 2,
                 force_local_pubmed: bool = False,
                 bert_model: str = 'NeuML/pubmedbert-base-embeddings',
                 llm_model: str = 'qwen3.5:9b'):
        self.config = PipelineConfig()
        self.output_dir = output_dir
        self.min_databases = min_databases
        self.force_local_pubmed = force_local_pubmed
        self.bert_model = bert_model
        self.llm_model = llm_model
        self.run_status = 'running'
        self.step_status = {}
        self.timestamp = create_timestamp()
        self.results_dir = os.path.join(output_dir, f"run_{self.timestamp}")
        
        # 创建结果目录结构
        self.config.create_output_structure(self.results_dir)
        
        # 初始化日志
        log_file = os.path.join(self.results_dir, "pipeline_log.txt")
        self.logger = Logger(log_file)
        
        # 初始化组件
        self.mti_selector = MTISelector(min_databases, self.logger)
        self.literature_miner = None  # 稍后初始化
        
        # 记录初始化信息
        self.logger.info(f"Pipeline initialized at {self.timestamp}")
        self.logger.info(f"Results directory: {self.results_dir}")
        self.logger.info(f"Minimum database requirement: {self.min_databases}")
        self.logger.info(f"Force local PubMed: {self.force_local_pubmed}")
        
        # 存储流程数据
        self.databases = {}
        self.mti_results = []
        self.mining_results = []
        self.validation_results = None
        self.llm_summaries = None
        self.network_results = None
    
    def load_databases(self, targetscan_file: str, mirdb_file: str, 
                      mirwalk_file: str, mirtarbase_file: Optional[str] = None):
        """加载数据库文件"""
        self.logger.info("Loading database files...")
        
        try:
            self.databases = self.mti_selector.load_databases(
                targetscan_file, mirdb_file, mirwalk_file, mirtarbase_file
            )
            self.logger.info("All database files loaded successfully")
            return True
        except Exception as e:
            self.logger.error(f"Failed to load database files: {e}")
            return False
    
    def step1_mti_selection(self, mode: str, gene_list_file: Optional[str] = None, 
                        mirna_list_file: Optional[str] = None, 
                        existing_results_file: Optional[str] = None,
                        direct_pairs_file: Optional[str] = None) -> bool:
        """Step 1: MTI选择，支持加载已有结果"""
        try:
            # 如果提供了已有结果文件，直接加载
            if existing_results_file and os.path.exists(existing_results_file):
                self.logger.info(f"Loading existing Step 1 results from: {existing_results_file}")
                
                if existing_results_file.endswith('.xlsx'):
                    df = pd.read_excel(existing_results_file)
                else:
                    df = pd.read_csv(existing_results_file)
                
                # 转换为internal format
                self.mti_results = df.to_dict('records')
                self.logger.info(f"Loaded {len(self.mti_results)} MTI results")
                return True
            
            # 原有的MTI选择逻辑（仅在没有现有结果时执行）
            if mode == 'gene':
                if not gene_list_file:
                    raise ValueError("Gene list file required for gene mode")
                self.mti_results = self.mti_selector.gene_to_mirna_selection(
                    gene_list_file, self.databases
                )
            elif mode == 'mirna':
                if not mirna_list_file:
                    raise ValueError("miRNA list file required for mirna mode")
                self.mti_results = self.mti_selector.mirna_to_gene_selection(
                    mirna_list_file, self.databases, gene_list_file
                )
            elif mode == 'combined':
                if not gene_list_file or not mirna_list_file:
                    raise ValueError("Both gene and miRNA list files required for combined mode")
                self.mti_results = self.mti_selector.combined_selection(
                    gene_list_file, mirna_list_file, self.databases
                )
            elif mode == 'direct':
                if direct_pairs_file:
                    self.mti_results = self._load_direct_pairs(direct_pairs_file)
                else:
                    self.mti_results = self._build_direct_pairs(
                        gene_list_file, mirna_list_file
                    )
            else:
                raise ValueError(f"Invalid mode: {mode}")
            
            if not self.mti_results:
                self.logger.warning("No MTIs found in Step 1")
                return False
            
            # 保存结果
            output_file = os.path.join(self.results_dir, "mirna_selection", f"mti_selection_results_{mode}.xlsx")
            self.mti_selector.save_results(self.mti_results, output_file, mode)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Step 1 failed: {e}")
            return False

    def _load_direct_pairs(self, pairs_file: Optional[str]) -> List[Dict[str, Any]]:
        """Load explicit Gene-miRNA pairs without prediction-database filtering."""
        if not pairs_file:
            raise ValueError("Direct mode requires --pairs")
        if not os.path.exists(pairs_file):
            raise FileNotFoundError(f"Direct pair file not found: {pairs_file}")

        extension = os.path.splitext(pairs_file)[1].lower()
        if extension in {'.xlsx', '.xls'}:
            df = pd.read_excel(pairs_file)
        else:
            df, _ = FileUtils.smart_read_csv(pairs_file)

        normalized_columns = {
            ''.join(character for character in str(column).lower() if character.isalnum()): column
            for column in df.columns
        }

        def find_column(candidates: List[str]) -> Optional[str]:
            for candidate in candidates:
                if candidate in normalized_columns:
                    return normalized_columns[candidate]
            return None

        gene_column = find_column([
            'gene', 'targetgene', 'targetgenesymbol', 'genesymbol', 'genename'
        ])
        mirna_column = find_column([
            'mirna', 'microrna', 'mirnaname', 'mirnanorm', 'maturemirna'
        ])
        if not gene_column or not mirna_column:
            raise ValueError(
                "Direct pair file must contain Gene and miRNA columns "
                f"(available columns: {list(df.columns)})"
            )

        rename_map = {}
        if gene_column != 'Gene':
            rename_map[gene_column] = 'Gene'
        if mirna_column != 'miRNA':
            rename_map[mirna_column] = 'miRNA'
        direct_df = df.rename(columns=rename_map).copy()
        direct_df = direct_df.dropna(subset=['Gene', 'miRNA'])
        direct_df['Gene'] = direct_df['Gene'].astype(str).str.strip()
        direct_df['miRNA'] = direct_df['miRNA'].astype(str).str.strip()
        direct_df = direct_df[
            direct_df['Gene'].ne('') & direct_df['miRNA'].ne('')
        ]
        direct_df = direct_df.drop_duplicates(subset=['Gene', 'miRNA'], keep='first')

        defaults = {
            'Direction': 'DirectInput',
            'Databases': 0,
            'Priority': 'Unfiltered',
            'Database_Sources': 'DirectInput',
            'miRTarBase_Validation': False,
            'miRTarBase_Support_Type': 'Not assessed',
            'Validation_Strength': 'Not assessed',
        }
        for column, default in defaults.items():
            if column not in direct_df.columns:
                direct_df[column] = default
            else:
                direct_df[column] = direct_df[column].fillna(default)

        self.logger.info(
            f"Loaded {len(direct_df)} unique direct Gene-miRNA pairs from {pairs_file}"
        )
        return direct_df.to_dict('records')

    def _build_direct_pairs(self, gene_list_file: Optional[str],
                            mirna_list_file: Optional[str]) -> List[Dict[str, Any]]:
        """Build the Cartesian product of explicit gene and miRNA lists."""
        if not gene_list_file or not mirna_list_file:
            raise ValueError("Direct mode requires --pairs or both --genes and --mirnas")

        def load_ordered_values(path: str, label: str) -> List[str]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"{label} list file not found: {path}")
            with open(path, 'r', encoding='utf-8-sig') as handle:
                values = [line.strip() for line in handle if line.strip()]
            return list(dict.fromkeys(values))

        genes = load_ordered_values(gene_list_file, 'Gene')
        mirnas = load_ordered_values(mirna_list_file, 'miRNA')
        if not genes or not mirnas:
            raise ValueError("Direct mode gene and miRNA lists must not be empty")

        pairs = [
            {
                'Gene': gene,
                'miRNA': mirna,
                'Direction': 'DirectInput',
                'Databases': 0,
                'Priority': 'Unfiltered',
                'Database_Sources': 'DirectInput',
                'miRTarBase_Validation': False,
                'miRTarBase_Support_Type': 'Not assessed',
                'Validation_Strength': 'Not assessed',
            }
            for gene in genes
            for mirna in mirnas
        ]
        self.logger.info(
            f"Built {len(pairs)} direct pairs from {len(genes)} genes x "
            f"{len(mirnas)} miRNAs"
        )
        return pairs
        
    def step2_literature_mining(self, max_articles: int = None) -> bool:
        """Step 2: 文献挖掘"""
        try:
            if not self.mti_results:
                self.logger.error("No MTI results available for literature mining")
                return False
            
            # 初始化文献挖掘器
            abstracts_dir = os.path.join(self.results_dir, "pubmed_articles")
            backend = 'local' if self.force_local_pubmed else None
            self.literature_miner = LiteratureMiner(
                abstracts_dir, self.logger, backend=backend
            )
            # Family-based filenames and the downstream mapping must be unique.
            # Reject ambiguous batches instead of overwriting one arm's evidence.
            stems = {}
            for pair in self.mti_results:
                term = self.literature_miner.pubmed_miner.extract_mirna_for_search(pair['miRNA'])
                key = (pair['Gene'], term)
                previous = stems.setdefault(key, pair['miRNA'])
                if previous != pair['miRNA']:
                    raise ValueError(
                        f"Ambiguous family-level output for {pair['Gene']}: "
                        f"{previous} and {pair['miRNA']}. Run these mature miRNAs separately."
                    )
            
            # 执行文献挖掘
            self.mining_results = self.literature_miner.mine_literature(
                self.mti_results, max_articles
            )
            if any(result.get('Error') for result in self.mining_results):
                self.logger.error('Literature retrieval contains failed queries; inspect mining_summary.csv')
                return False
            
            # 验证下载的文件
            validation = self.literature_miner.validate_downloaded_files()
            self.logger.info(f"File validation results:")
            self.logger.info(f"  Valid files: {len(validation['valid_files'])}")
            self.logger.info(f"  Empty files: {len(validation['empty_files'])}")
            self.logger.info(f"  Missing files: {len(validation['missing_files'])}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Step 2 failed: {e}")
            return False
    
    def step3_bert_validation(self, functions: List[str] = None, valid_threshold: float = 30) -> bool:
        """Step 3: BERT验证"""
        try:
            if not self.mining_results:
                self.logger.error("No mining results available for BERT validation")
                return False
            
            if functions is None:
                functions = self.config.DEFAULT_FUNCTIONS
            
            # 创建批量输入文件
            batch_data = []
            mirna_mapping = {}
            
            for result in self.mining_results:
                if result['Articles_Found'] > 0:
                    mirna_simplified = result['miRNA_Search_Term']
                    mirna_full = result['miRNA']
                    
                    batch_data.append({
                        'mirna': mirna_simplified,
                        'gene': result['Gene'],
                        'file_path': os.path.join(
                            self.results_dir, "pubmed_articles", result['File']
                        )
                    })
                    
                    mirna_mapping[f"{mirna_simplified}_{result['Gene']}"] = mirna_full
            
            if not batch_data:
                self.logger.warning("No articles found for BERT validation")
                return False
            
            # 保存批量输入文件
            batch_file = os.path.join(self.results_dir, "bert_validation", "batch_input.csv")
            pd.DataFrame(batch_data).to_csv(batch_file, index=False)
            
            # 运行BERT验证
            output_file = os.path.join(self.results_dir, "bert_validation", "mti_validation_results.csv")
            
            bert_script = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), 'mirna_bert.py'
            )
            cmd = [
                sys.executable, bert_script,
                '--batch',
                '--batch-file', batch_file,
                '--model', self.bert_model,
                '-f', ','.join(functions),
                '--valid-threshold', str(valid_threshold),
                '-o', output_file
            ]
            
            self.logger.info(f"Running BERT validation for {len(batch_data)} MTIs...")
            self.logger.info(f"Functions to evaluate: {', '.join(functions)}")
            self.logger.info(f"BERT valid threshold: {valid_threshold}")
            
            subprocess.run(cmd, check=True)
            self.logger.info("BERT validation completed successfully")
            
            # 读取结果并处理miRNA名称映射
            validation_results = pd.read_csv(output_file)
            
            # 映射回全称miRNA名称
            validation_results['miRNA_Full'] = validation_results.apply(
                lambda row: mirna_mapping.get(f"{row['miRNA']}_{row['Target_Gene']}", row['miRNA']),
                axis=1
            )
            
            validation_results['miRNA_Simplified'] = validation_results['miRNA']
            validation_results['miRNA'] = validation_results['miRNA_Full']
            validation_results.drop('miRNA_Full', axis=1, inplace=True)
            
            # 添加其他信息
            info_map = {(r['Gene'], r['miRNA']): {
                'Direction': r['Direction'],
                'Databases': r['Databases'],
                'Database_Sources': r['Database_Sources'],
                'miRTarBase_Validation': r.get('miRTarBase_Validation', False),
                'Validation_Strength': r.get('Validation_Strength', 'Unknown')
            } for r in self.mining_results}
            
            for col, key in [
                ('Search_Direction', 'Direction'),
                ('Database_Count', 'Databases'),
                ('Database_Sources', 'Database_Sources'),
                ('miRTarBase_Validation', 'miRTarBase_Validation'),
                ('Validation_Strength', 'Validation_Strength')
            ]:
                validation_results[col] = validation_results.apply(
                    lambda row: info_map.get((row['Target_Gene'], row['miRNA']), {}).get(key, 'Unknown'),
                    axis=1
                )
            
            # 保存更新的结果
            validation_results.to_csv(output_file, index=False)
            self.validation_results = validation_results
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BERT validation failed: {e}")
            return False
        except Exception as e:
            self.logger.error(f"Step 3 failed: {e}")
            return False
    
    def step4_llm_analysis(self, functions: List[str] = None, analysis_mode: str = "both") -> bool:
        """Step 4: LLM分析"""
        try:
            if MTILLMSummarizer is None:
                self.logger.error("MTILLMSummarizer not available. Please ensure mti_llm_summarize.py exists.")
                return False
                
            if self.validation_results is None or len(self.validation_results) == 0:
                self.logger.warning("No validation results available for LLM analysis")
                return False
            
            if functions is None:
                functions = self.config.DEFAULT_FUNCTIONS
            
            # 检查Ollama连接
            is_connected, model_names = check_ollama_connection()
            if not is_connected:
                self.logger.error("Cannot connect to Ollama. Please start Ollama: ollama serve")
                self.logger.error(f"Install the selected model: ollama pull {self.llm_model}")
                return False
            
            self.logger.info(f"Ollama connected. Available models: {model_names}")
            
            # 创建LLM总结器
            summarizer = MTILLMSummarizer(model_name=self.llm_model)
            
            if summarizer.llm is None:
                self.logger.error("LLM not available")
                return False
            
            # 输出目录
            llm_output_dir = os.path.join(self.results_dir, "llm_summaries")
            
            # 准备文件路径
            abstracts_dir = os.path.join(self.results_dir, "pubmed_articles")
            validation_file = os.path.join(self.results_dir, "bert_validation", "mti_validation_results.csv")
            
            summaries = {}
            
            # Standard Analysis
            if analysis_mode in ["standard", "both"]:
                self.logger.info("Generating standard MTI relationship summaries...")
                try:
                    standard_summaries = summarizer.batch_summarize(
                        validation_results_file=validation_file,
                        abstracts_dir=abstracts_dir,
                        functions=functions,
                        output_dir=llm_output_dir
                    )
                    summaries['standard'] = standard_summaries
                    self.logger.info(f"Standard summarization completed. Generated {len(standard_summaries)} summaries.")
                except Exception as e:
                    self.logger.error(f"Error in standard summarization: {e}")
                    summaries['standard'] = {}
            
            # Functional Analysis
            if analysis_mode in ["functional", "both"]:
                self.logger.info("Generating functional pathway analysis...")
                try:
                    functional_output_dir = os.path.join(llm_output_dir, "functional_analysis")
                    functional_summaries = summarizer.batch_functional_analysis(
                        validation_results_file=validation_file,
                        abstracts_dir=abstracts_dir,
                        functions=functions,
                        output_dir=functional_output_dir
                    )
                    summaries['functional'] = functional_summaries
                    self.logger.info(f"Functional analysis completed. Generated {len(functional_summaries)} analyses.")
                except Exception as e:
                    self.logger.error(f"Error in functional analysis: {e}")
                    summaries['functional'] = {}
            
            self.llm_summaries = summaries
            def contains_error(value):
                if isinstance(value, dict):
                    return bool(value.get('error')) or any(contains_error(item) for item in value.values())
                if isinstance(value, list):
                    return any(contains_error(item) for item in value)
                return False
            if any(not results or contains_error(results) for results in summaries.values()):
                self.logger.error('A requested LLM analysis produced no summaries or returned an error')
                return False
            
            # 整合结果到CSV
            integrated_results = self._integrate_llm_results(analysis_mode)
            self.validation_results = integrated_results
            
            total_summaries = sum(len(s) for s in summaries.values())
            self.logger.info(f"LLM analysis completed. Total summaries generated: {total_summaries}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Step 4 failed: {e}")
            return False
    
    def step5_network_generation(self, mirna_deseq_file: Optional[str] = None, 
                                gene_deseq_file: Optional[str] = None, 
                                score_threshold: float = 50) -> bool:
        """Step 5: 网络生成"""
        try:
            if MTICytoscapeGenerator is None:
                self.logger.error("MTICytoscapeGenerator not available. Please ensure mti_cytoscape_network.py exists.")
                return False
                
            if self.validation_results is None or len(self.validation_results) == 0:
                self.logger.warning("No validation results available for network generation")
                return False
            
            # 检查DESeq2文件
            if not mirna_deseq_file or not gene_deseq_file:
                self.logger.warning("DESeq2 files not provided, using default expression values")
                mirna_deseq_file, gene_deseq_file = self._create_default_deseq_files()
            elif not os.path.exists(mirna_deseq_file) or not os.path.exists(gene_deseq_file):
                self.logger.warning("DESeq2 files not found, using default expression values")
                mirna_deseq_file, gene_deseq_file = self._create_default_deseq_files()
            
            # 创建网络生成器
            network_output_dir = os.path.join(self.results_dir, "cytoscape_network")
            network_generator = MTICytoscapeGenerator(
                output_dir=network_output_dir,
                parent_log_func=self.logger.info
            )
            
            network_generator.score_threshold = score_threshold
            
            # 准备MTI数据文件
            mti_file = os.path.join(self.results_dir, "bert_validation", "integrated_validation_results.csv")
            if not os.path.exists(mti_file):
                mti_file = os.path.join(self.results_dir, "bert_validation", "mti_validation_results.csv")
            
            self.logger.info(f"Using MTI file: {mti_file}")
            self.logger.info(f"Using miRNA DESeq2 file: {mirna_deseq_file}")
            self.logger.info(f"Using gene DESeq2 file: {gene_deseq_file}")
            self.logger.info(f"Score threshold: {score_threshold}")
            
            # 运行网络生成
            result = network_generator.run_full_pipeline(
                mti_file=mti_file,
                mirna_deseq_file=mirna_deseq_file,
                gene_deseq_file=gene_deseq_file
            )
            
            if result['success']:
                self.network_results = result
                net_stats = result['summary']['network_statistics']
                self.logger.info("Network generation completed successfully!")
                self.logger.info(f"Generated {net_stats['total_nodes']} nodes and {net_stats['total_edges']} edges")
                return True
            else:
                self.logger.error(f"Network generation failed: {result.get('error', 'Unknown error')}")
                return False
                
        except Exception as e:
            self.logger.error(f"Step 5 failed: {e}")
            return False
    
    def _integrate_llm_results(self, analysis_mode: str) -> pd.DataFrame:
        """整合LLM结果到验证CSV中"""
        if self.validation_results is None or self.llm_summaries is None:
            return self.validation_results
        
        self.logger.info("Integrating LLM analysis results into validation CSV...")
        integrated_df = self.validation_results.copy()
        
        # Standard Analysis结果
        if analysis_mode in ["standard", "both"] and 'standard' in self.llm_summaries:
            standard_results = self.llm_summaries['standard']
            
            # 添加标准分析列
            for col in ['LLM_Relationship_Summary', 'LLM_Evidence_Strength', 
                       'LLM_Clinical_Relevance', 'LLM_Key_Findings', 
                       'LLM_Contradictions', 'LLM_Research_Gaps', 'LLM_Overall_Assessment']:
                integrated_df[col] = ''
            
            # 填充数据
            for _, row in integrated_df.iterrows():
                mirna_simplified = row.get('miRNA_Simplified', '')
                if not mirna_simplified:
                    mirna_simplified = TextProcessor.extract_mirna_core(row['miRNA'])
                
                gene = row.get('Gene', row.get('Target_Gene', ''))
                mti_key = f"{mirna_simplified}_{gene}"
                
                if mti_key in standard_results:
                    llm_data = standard_results[mti_key]
                    if isinstance(llm_data, dict) and 'standard_summary' in llm_data:
                        summary_data = llm_data['standard_summary']
                    else:
                        summary_data = llm_data
                    
                    if 'error' not in summary_data:
                        idx = row.name
                        integrated_df.loc[idx, 'LLM_Relationship_Summary'] = TextProcessor.safe_get_text(summary_data, 'relationship_summary')
                        integrated_df.loc[idx, 'LLM_Evidence_Strength'] = TextProcessor.safe_get_text(summary_data, 'evidence_strength')
                        integrated_df.loc[idx, 'LLM_Clinical_Relevance'] = TextProcessor.safe_get_text(summary_data, 'clinical_relevance')
                        integrated_df.loc[idx, 'LLM_Key_Findings'] = TextProcessor.safe_get_list(summary_data, 'key_findings')
                        integrated_df.loc[idx, 'LLM_Contradictions'] = TextProcessor.safe_get_list(summary_data, 'contradictions')
                        integrated_df.loc[idx, 'LLM_Research_Gaps'] = TextProcessor.safe_get_list(summary_data, 'research_gaps')
                        integrated_df.loc[idx, 'LLM_Overall_Assessment'] = TextProcessor.safe_get_text(summary_data, 'overall_assessment')
        
        # Functional Analysis结果
        if analysis_mode in ["functional", "both"] and 'functional' in self.llm_summaries:
            functional_results = self.llm_summaries['functional']
            
            # 添加功能分析列
            for col in ['Gene_Primary_Functions', 'Gene_Signaling_Pathways', 'Gene_Disease_Associations',
                       'miRNA_Regulatory_Functions', 'miRNA_Target_Pathways', 'miRNA_Disease_Relevance',
                       'Shared_Pathways', 'Interaction_Impact', 'Therapeutic_Potential', 
                       'Research_Priority', 'Validation_Level']:
                integrated_df[col] = ''
            
            # 填充数据
            for _, row in integrated_df.iterrows():
                mirna_simplified = row.get('miRNA_Simplified', '')
                if not mirna_simplified:
                    mirna_simplified = TextProcessor.extract_mirna_core(row['miRNA'])
                
                gene = row.get('Gene', row.get('Target_Gene', ''))
                mti_key = f"{mirna_simplified}_{gene}"
                
                if mti_key in functional_results:
                    func_data = functional_results[mti_key]
                    if isinstance(func_data, dict) and 'functional_analysis' in func_data:
                        functional_data = func_data['functional_analysis']
                    else:
                        functional_data = func_data
                    
                    if 'error' not in functional_data:
                        idx = row.name
                        
                        # 基因分析数据
                        gene_analysis = functional_data.get('gene_analysis', {})
                        integrated_df.loc[idx, 'Gene_Primary_Functions'] = TextProcessor.safe_get_list(gene_analysis, 'primary_functions')
                        integrated_df.loc[idx, 'Gene_Signaling_Pathways'] = TextProcessor.safe_get_list(gene_analysis, 'signaling_pathways')
                        integrated_df.loc[idx, 'Gene_Disease_Associations'] = TextProcessor.safe_get_list(gene_analysis, 'disease_associations')
                        
                        # miRNA分析数据
                        mirna_analysis = functional_data.get('mirna_analysis', {})
                        integrated_df.loc[idx, 'miRNA_Regulatory_Functions'] = TextProcessor.safe_get_list(mirna_analysis, 'regulatory_functions')
                        integrated_df.loc[idx, 'miRNA_Target_Pathways'] = TextProcessor.safe_get_list(mirna_analysis, 'target_pathways')
                        integrated_df.loc[idx, 'miRNA_Disease_Relevance'] = TextProcessor.safe_get_list(mirna_analysis, 'disease_relevance')
                        
                        # 通路分析
                        pathway_analysis = functional_data.get('pathway_analysis', {})
                        integrated_df.loc[idx, 'Shared_Pathways'] = TextProcessor.safe_get_list(pathway_analysis, 'shared_pathways')
                        
                        # 功能总结
                        functional_summary = functional_data.get('functional_summary', {})
                        integrated_df.loc[idx, 'Interaction_Impact'] = TextProcessor.safe_get_text(functional_summary, 'interaction_impact')
                        integrated_df.loc[idx, 'Therapeutic_Potential'] = TextProcessor.safe_get_text(functional_summary, 'therapeutic_potential')
                        integrated_df.loc[idx, 'Research_Priority'] = TextProcessor.safe_get_text(functional_summary, 'research_priority')
                        
                        # 证据质量
                        evidence_quality = functional_data.get('evidence_quality', {})
                        integrated_df.loc[idx, 'Validation_Level'] = TextProcessor.safe_get_text(evidence_quality, 'validation_level')
        
        # 保存整合后的结果
        integrated_file = os.path.join(self.results_dir, "bert_validation", "integrated_validation_results.csv")
        integrated_df.to_csv(integrated_file, index=False)
        
        self.logger.info(f"Integrated results saved to: {integrated_file}")
        
        return integrated_df
    
    def _create_default_deseq_files(self) -> tuple:
        """创建默认的DESeq2文件"""
        self.logger.info("Creating default DESeq2 files for network generation...")
        
        if self.validation_results is None:
            return None, None
        
        # 获取唯一的基因和miRNA
        genes = set()
        mirnas = set()
        
        for _, row in self.validation_results.iterrows():
            if 'Gene' in row:
                genes.add(row['Gene'])
            elif 'Target_Gene' in row:
                genes.add(row['Target_Gene'])
            
            if 'miRNA' in row:
                mirnas.add(row['miRNA'])
        
        # 创建默认的miRNA DESeq2文件
        mirna_data = []
        for mirna in mirnas:
            mirna_data.append({
                'miRNA': mirna,
                'baseMean': 100.0,
                'log2FoldChange': 0.0,
                'lfcSE': 0.1,
                'stat': 0.0,
                'pvalue': 1.0,
                'padj': 1.0
            })
        
        mirna_df = pd.DataFrame(mirna_data)
        mirna_file = os.path.join(self.results_dir, "default_mirna_deseq.csv")
        mirna_df.to_csv(mirna_file, index=False)
        
        # 创建默认的基因DESeq2文件
        gene_data = []
        for gene in genes:
            gene_data.append({
                'mRNA': gene,
                'baseMean': 100.0,
                'log2FoldChange': 0.0,
                'lfcSE': 0.1,
                'stat': 0.0,
                'pvalue': 1.0,
                'padj': 1.0
            })
        
        gene_df = pd.DataFrame(gene_data)
        gene_file = os.path.join(self.results_dir, "default_gene_deseq.csv")
        gene_df.to_csv(gene_file, index=False)
        
        self.logger.info(f"Created default DESeq2 files:")
        self.logger.info(f"  miRNA: {mirna_file} ({len(mirnas)} entries)")
        self.logger.info(f"  Gene: {gene_file} ({len(genes)} entries)")
        self.logger.info("Note: Default values used (no differential expression)")
        
        return mirna_file, gene_file
    
    def generate_final_report(self, mode: str, analysis_mode: str = "both"):
        """Write the execution status and the results actually produced."""
        report_file = os.path.join(self.results_dir, "pipeline_report.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("MiTAgent analysis report\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"Status: {self.run_status}\n")
            f.write(f"Candidate mode: {mode}\n")
            f.write(f"Analysis mode: {analysis_mode}\n")
            f.write(f"BERT model: {self.bert_model}\n")
            f.write(f"LLM model (when enabled): {self.llm_model}\n")
            for step, status in self.step_status.items():
                f.write(f"{step}: {status}\n")
            if os.path.exists(os.path.join(self.results_dir, 'default_gene_deseq.csv')):
                f.write("Expression annotation: PLACEHOLDER values; not measured expression.\n")
            self._write_step_summary(f, "STEP 1: MTI Selection", self.mti_results)
            self._write_step_summary(f, "STEP 2: Literature Mining", self.mining_results)
            self._write_step_summary(f, "STEP 3: BERT Validation", self.validation_results)
            self._write_step_summary(f, "STEP 4: LLM Analysis", self.llm_summaries)
            self._write_step_summary(f, "STEP 5: Network Generation", self.network_results)
            self._write_output_files(f)
        self.logger.info(f"Report saved to: {report_file}")

    def _write_step_summary(self, f, step_name: str, step_results):
        """写入步骤摘要"""
        f.write(f"\n{step_name}\n")
        f.write("-"*50 + "\n")
        
        if step_name == "STEP 1: MTI Selection" and self.mti_results:
            f.write(f"Total MTIs found: {len(self.mti_results)}\n")
            stats = ResultsIntegrator.create_summary_stats(self.mti_results)
            for key, value in stats.items():
                f.write(f"  {key}: {value}\n")
        elif step_name == "STEP 2: Literature Mining" and self.mining_results:
            total_articles = sum(r.get('Articles_Found', 0) for r in self.mining_results)
            mtis_with_articles = len([r for r in self.mining_results if r.get('Articles_Found', 0) > 0])
            f.write(f"Total articles found: {total_articles}\n")
            f.write(f"MTIs with articles: {mtis_with_articles}\n")
        elif step_name == "STEP 3: BERT Validation" and self.validation_results is not None:
            f.write(f"MTIs scored: {len(self.validation_results)}\n")
            if 'Overall_Score' in self.validation_results.columns:
                high_score = len(self.validation_results[self.validation_results['Overall_Score'] > 80])
                f.write(f"Pairs with Overall_Score > 80: {high_score}\n")
        elif step_name == "STEP 4: LLM Analysis" and self.llm_summaries:
            total_summaries = sum(len(s) for s in self.llm_summaries.values())
            f.write(f"Total summaries generated: {total_summaries}\n")
        elif step_name == "STEP 5: Network Generation" and self.network_results:
            if self.network_results.get('success'):
                stats = self.network_results['summary']['network_statistics']
                f.write(f"Network nodes: {stats['total_nodes']}\n")
                f.write(f"Network edges: {stats['total_edges']}\n")
            else:
                f.write("Network generation: Failed\n")
        else:
            f.write("Not completed or no results\n")
    
    def _write_output_files(self, f):
        """List existing output paths relative to the run directory."""
        f.write("\nOUTPUT FILES\n")
        root = Path(self.results_dir)
        for path in sorted(root.rglob('*')):
            if path.is_file():
                f.write(path.relative_to(root).as_posix() + "\n")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='MiTAgent: context-informed miRNA-target evidence prioritization')
    
    # Mode selection (required)
    parser.add_argument('-m', '--mode', choices=['gene', 'mirna', 'combined', 'direct'], required=True,
                       help='Pipeline mode: gene, mirna, combined, or direct (explicit pairs to local PubMed + BERT)')
    
    # Input files
    parser.add_argument('-g', '--genes', help='Gene list file (gene/combined, or direct Cartesian input)')
    parser.add_argument('-r', '--mirnas', help='miRNA list file (mirna/combined, or direct Cartesian input)')
    parser.add_argument('--pairs', help='CSV/TSV/XLSX with Gene and miRNA columns (direct mode)')
    parser.add_argument('-t', '--targetscan', help='TargetScan database file')
    parser.add_argument('-d', '--mirdb', help='miRDB database file')
    parser.add_argument('-w', '--mirwalk', help='miRWalk database file')
    parser.add_argument('-b', '--mirtarbase', help='miRTarBase database file (optional)')
    parser.add_argument('--existing-step1', help='Existing Step 1 results file to continue from')
    # Pipeline steps
    parser.add_argument('--skip-literature', action='store_true', help='Skip literature mining step')
    parser.add_argument('--skip-bert', action='store_true', help='Skip BERT validation step')
    parser.add_argument('--skip-llm', action='store_true', help='Skip LLM summarization step')
    
    # LLM analysis configuration
    parser.add_argument('--analysis-mode', choices=['standard', 'functional', 'both'], default='both',
                       help='LLM analysis mode: standard (MTI relationships), functional (pathways), or both')
    
    # Network generation
    parser.add_argument('--generate-network', action='store_true', 
                       help='Generate Cytoscape network visualization')
    parser.add_argument('--mirna-deseq', help='miRNA DESeq2 results file for network generation')
    parser.add_argument('--gene-deseq', help='Gene DESeq2 results file for network generation')
    parser.add_argument('--network-score-threshold', type=float, default=50,
                       help='Minimum overall score for network inclusion (default: 50)')
    
    # Analysis parameters
    parser.add_argument('-f', '--functions', default="cell proliferation,apoptosis,migration,invasion",
                       help='Comma-separated list of biological functions to analyze')
    parser.add_argument('--min-databases', type=int, default=2,
                       help='Minimum number of databases required for MTI inclusion (default: 2)')
    parser.add_argument('--max-articles', type=int, default=50,
                       help='Maximum articles to retrieve per MTI (default: 50)')
    parser.add_argument('--bert-valid-threshold', type=float, default=30,
                       help='Minimum BERT Overall_Score used for the Valid flag (default: 30)')
    
    # Output
    parser.add_argument('-o', '--output', default="pipeline_results",
                       help='Output directory (default: pipeline_results)')
    
    parser.add_argument('--bert-model', default='NeuML/pubmedbert-base-embeddings',
                        help='SentenceTransformer model ID or local model directory')
    parser.add_argument('--llm-model', default='qwen3.5:9b', help='Ollama model name')
    parser.add_argument('--allow-placeholder-expression', action='store_true',
                        help='Allow topology-only network export with labeled placeholder expression')
    args = parser.parse_args()
    
    # Existing Step 1 and direct mode do not require prediction databases.
    if args.existing_step1 and not os.path.exists(args.existing_step1):
        parser.error(f"Existing Step 1 file not found: {args.existing_step1}")

    if not args.existing_step1:
        if args.mode in ['gene', 'combined'] and not args.genes:
            parser.error("Gene list file required for gene and combined modes")
        if args.mode in ['mirna', 'combined'] and not args.mirnas:
            parser.error("miRNA list file required for mirna and combined modes")
        if args.mode == 'direct' and not args.pairs and not (args.genes and args.mirnas):
            parser.error("Direct mode requires --pairs or both --genes and --mirnas")
        if args.mode != 'direct':
            missing_databases = [
                option for option, value in [
                    ('--targetscan', args.targetscan),
                    ('--mirdb', args.mirdb),
                    ('--mirwalk', args.mirwalk),
                ] if not value
            ]
            if missing_databases:
                parser.error(
                    "Prediction database files required: " + ', '.join(missing_databases)
                )

    if args.mode == 'direct':
        args.skip_llm = True
    if args.skip_literature and not args.skip_bert:
        parser.error('--skip-literature requires --skip-bert; existing literature is not loaded')
    if args.skip_bert and (not args.skip_llm or args.generate_network):
        parser.error('--skip-bert requires --skip-llm and cannot be used with --generate-network')
    if args.max_articles < 1 or args.min_databases < 1:
        parser.error('--max-articles and --min-databases must be positive')
    if not 0 <= args.bert_valid_threshold <= 100 or not 0 <= args.network_score_threshold <= 100:
        parser.error('Score thresholds must be between 0 and 100')
    if args.generate_network:
        if bool(args.mirna_deseq) != bool(args.gene_deseq):
            parser.error('Provide both --mirna-deseq and --gene-deseq')
        if not args.mirna_deseq and not args.allow_placeholder_expression:
            parser.error('Network export needs both expression files, or explicit --allow-placeholder-expression')
        for path in (args.mirna_deseq, args.gene_deseq):
            if path and not os.path.isfile(path):
                parser.error(f'Expression input not found: {path}')

    functions = [f.strip() for f in args.functions.split(',') if f.strip()]
    effective_analysis_mode = 'bert-only' if args.mode == 'direct' else args.analysis_mode
    if args.skip_bert:
        effective_analysis_mode = 'candidate-only' if args.skip_literature else 'literature-only'
    print('MiTAgent: context-informed miRNA-target evidence prioritization')
    pipeline = MiRNAPipeline(
        output_dir=args.output, min_databases=args.min_databases,
        force_local_pubmed=(args.mode == 'direct'),
        bert_model=args.bert_model, llm_model=args.llm_model,
    )
    parameters = {'arguments': vars(args), 'python': platform.python_version(), 'packages': {}}
    for package in ('pandas', 'numpy', 'sentence-transformers', 'transformers', 'torch', 'nltk'):
        try:
            parameters['packages'][package] = version(package)
        except PackageNotFoundError:
            parameters['packages'][package] = None
    with open(os.path.join(pipeline.results_dir, 'run_parameters.json'), 'w', encoding='utf-8') as handle:
        json.dump(parameters, handle, indent=2)

    def finish(status):
        pipeline.run_status = status
        pipeline.generate_final_report(args.mode, effective_analysis_mode)
        print(f'Analysis {status}. Results: {pipeline.results_dir}')
        return 0 if status == 'completed' else 1

    def run_step(name, operation):
        try:
            success = operation()
        except Exception as exc:
            pipeline.logger.error(f'{name} failed: {exc}')
            success = False
        pipeline.step_status[name] = 'completed' if success else 'failed'
        return success

    if not args.existing_step1 and args.mode != 'direct':
        if not run_step('Database loading', lambda: pipeline.load_databases(
            args.targetscan, args.mirdb, args.mirwalk, args.mirtarbase)):
            return finish('failed')
    if not run_step('Candidate selection', lambda: pipeline.step1_mti_selection(
        args.mode, args.genes, args.mirnas, args.existing_step1, args.pairs)):
        return finish('failed')
    if args.skip_literature:
        pipeline.step_status['Literature retrieval'] = 'skipped'
    elif not run_step('Literature retrieval', lambda: pipeline.step2_literature_mining(args.max_articles)):
        return finish('failed')
    if args.skip_bert:
        pipeline.step_status['BERT scoring'] = 'skipped; no scores assigned'
    elif not run_step('BERT scoring', lambda: pipeline.step3_bert_validation(functions, args.bert_valid_threshold)):
        return finish('failed')
    if args.skip_llm:
        pipeline.step_status['LLM interpretation'] = 'skipped'
    elif not run_step('LLM interpretation', lambda: pipeline.step4_llm_analysis(functions, args.analysis_mode)):
        return finish('failed')
    if args.generate_network:
        if not run_step('Network export', lambda: pipeline.step5_network_generation(
            args.mirna_deseq, args.gene_deseq, args.network_score_threshold)):
            return finish('failed')
    else:
        pipeline.step_status['Network export'] = 'skipped'
    return finish('completed')

if __name__ == '__main__':
    raise SystemExit(main())

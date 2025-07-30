"""
miRNA Research Pipeline - 主控制器
整合所有模块，执行完整的分析流程
"""

import os
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
    print("Warning: mti_llm_summarize.py not found. LLM analysis will be unavailable.")
    MTILLMSummarizer = None

try:
    from mti_cytoscape_network import MTICytoscapeGenerator
except ImportError:
    print("Warning: mti_cytoscape_network.py not found. Network generation will be unavailable.")
    MTICytoscapeGenerator = None

class MiRNAPipeline:
    """miRNA研究流水线主控制器"""
    
    def __init__(self, output_dir: str = "pipeline_results", min_databases: int = 2):
        self.config = PipelineConfig()
        self.output_dir = output_dir
        self.min_databases = min_databases
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
                           mirna_list_file: Optional[str] = None) -> bool:
        """Step 1: MTI选择"""
        try:
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
    
    def step2_literature_mining(self, max_articles: int = None) -> bool:
        """Step 2: 文献挖掘"""
        try:
            if not self.mti_results:
                self.logger.error("No MTI results available for literature mining")
                return False
            
            # 初始化文献挖掘器
            abstracts_dir = os.path.join(self.results_dir, "pubmed_articles")
            self.literature_miner = LiteratureMiner(abstracts_dir, self.logger)
            
            # 执行文献挖掘
            self.mining_results = self.literature_miner.mine_literature(
                self.mti_results, max_articles
            )
            
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
    
    def step3_bert_validation(self, functions: List[str] = None) -> bool:
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
            
            cmd = [
                'python', 'mirna_bert.py',
                '--batch',
                '--batch-file', batch_file,
                '-f', ','.join(functions),
                '-o', output_file
            ]
            
            self.logger.info(f"Running BERT validation for {len(batch_data)} MTIs...")
            self.logger.info(f"Functions to evaluate: {', '.join(functions)}")
            
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
                self.logger.error("Also install a model: ollama pull llama3.1")
                return False
            
            self.logger.info(f"Ollama connected. Available models: {model_names}")
            
            # 创建LLM总结器
            summarizer = MTILLMSummarizer()
            
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
        """生成最终报告"""
        self.logger.info("Generating comprehensive final report...")
        
        report_file = os.path.join(self.results_dir, "pipeline_report.txt")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("ENHANCED miRNA RESEARCH PIPELINE REPORT (Fixed Version)\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Pipeline Mode: {mode}\n")
            f.write(f"LLM Analysis Mode: {analysis_mode}\n")
            f.write(f"Results Directory: {self.results_dir}\n")
            f.write(f"Pipeline Version: Enhanced with Modular Architecture\n\n")
            
            # 修复说明
            f.write("🔧 ARCHITECTURAL IMPROVEMENTS:\n")
            f.write("-"*50 + "\n")
            f.write("✅ Modular design with separate components\n")
            f.write("✅ Centralized configuration management\n")
            f.write("✅ Enhanced error handling and logging\n")
            f.write("✅ Improved code maintainability\n")
            f.write("✅ Flexible pipeline execution\n\n")
            
            # Step summaries
            self._write_step_summary(f, "STEP 1: MTI Selection", self.mti_results)
            self._write_step_summary(f, "STEP 2: Literature Mining", self.mining_results)
            self._write_step_summary(f, "STEP 3: BERT Validation", self.validation_results)
            self._write_step_summary(f, "STEP 4: LLM Analysis", self.llm_summaries)
            self._write_step_summary(f, "STEP 5: Network Generation", self.network_results)
            
            # 输出文件列表
            self._write_output_files(f)
            
            f.write(f"\nPipeline completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("Thank you for using the Enhanced Modular miRNA Research Pipeline!\n")
        
        self.logger.info(f"Comprehensive report saved to: {report_file}")
    
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
            f.write(f"MTIs validated: {len(self.validation_results)}\n")
            if 'Overall_Score' in self.validation_results.columns:
                high_score = len(self.validation_results[self.validation_results['Overall_Score'] > 80])
                f.write(f"High confidence MTIs (>80): {high_score}\n")
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
        """写入输出文件信息"""
        f.write("\n" + "="*80 + "\n")
        f.write("KEY OUTPUT FILES:\n")
        f.write("="*80 + "\n")
        
        # 检查文件存在性并列出
        key_files = [
            ("integrated_validation_results.csv", "Complete analysis with experimental validation"),
            ("mti_validation_results.csv", "BERT validation results"),
            ("mti_selection_results_*.xlsx", "Initial MTI selection with miRTarBase info"),
            ("mining_summary.csv", "Literature mining summary"),
            ("pipeline_log.txt", "Complete execution log"),
            ("pipeline_report.txt", "This comprehensive report")
        ]
        
        for filename, description in key_files:
            f.write(f"📊 {filename}: {description}\n")
        
        if self.llm_summaries:
            f.write("📁 llm_summaries/: Detailed LLM analysis reports\n")
        
        if self.network_results and self.network_results.get('success'):
            f.write("🌐 cytoscape_network/: Network visualization files\n")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Enhanced Modular miRNA Research Pipeline')
    
    # Mode selection (required)
    parser.add_argument('-m', '--mode', choices=['gene', 'mirna', 'combined'], required=True,
                       help='Pipeline mode: gene (gene→miRNA), mirna (miRNA→gene), or combined')
    
    # Input files
    parser.add_argument('-g', '--genes', help='Gene list file (required for gene and combined modes)')
    parser.add_argument('-r', '--mirnas', help='miRNA list file (required for mirna and combined modes)')
    parser.add_argument('-t', '--targetscan', required=True, help='TargetScan database file')
    parser.add_argument('-d', '--mirdb', required=True, help='miRDB database file')
    parser.add_argument('-w', '--mirwalk', required=True, help='miRWalk database file')
    parser.add_argument('-b', '--mirtarbase', help='miRTarBase database file (optional)')
    
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
    
    # Output
    parser.add_argument('-o', '--output', default="pipeline_results",
                       help='Output directory (default: pipeline_results)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.mode in ['gene', 'combined'] and not args.genes:
        parser.error("Gene list file required for gene and combined modes")
    
    if args.mode in ['mirna', 'combined'] and not args.mirnas:
        parser.error("miRNA list file required for mirna and combined modes")
    
    # Parse functions
    functions = [f.strip() for f in args.functions.split(',')]
    
    print("="*80)
    print("ENHANCED MODULAR miRNA RESEARCH PIPELINE")
    print("="*80)
    print(f"Mode: {args.mode}")
    print(f"Analysis mode: {args.analysis_mode}")
    print(f"Functions to analyze: {', '.join(functions)}")
    print(f"Minimum databases: {args.min_databases}")
    print(f"Max articles per MTI: {args.max_articles}")
    if args.generate_network:
        print(f"Network generation: Enabled (threshold: {args.network_score_threshold})")
    print()
    
    # Initialize pipeline
    pipeline = MiRNAPipeline(output_dir=args.output, min_databases=args.min_databases)
    
    # Load databases
    success = pipeline.load_databases(
        args.targetscan, args.mirdb, args.mirwalk, args.mirtarbase
    )
    if not success:
        print("❌ Failed to load database files. Exiting.")
        return 1
    
    # Step 1: MTI Selection
    print("\n🔍 Step 1: MTI Selection")
    success = pipeline.step1_mti_selection(args.mode, args.genes, args.mirnas)
    if not success:
        print("❌ Step 1 failed. Exiting.")
        return 1
    
    # Step 2: Literature Mining
    if not args.skip_literature:
        print("\n📚 Step 2: Literature Mining")
        success = pipeline.step2_literature_mining(args.max_articles)
        if not success:
            print("⚠️ Step 2 failed, but continuing...")
    else:
        print("\n⏭️ Step 2: Skipped (Literature Mining)")
    
    # Step 3: BERT Validation
    if not args.skip_bert and pipeline.mining_results:
        print("\n🤖 Step 3: BERT Validation")
        success = pipeline.step3_bert_validation(functions)
        if not success:
            print("⚠️ Step 3 failed, but continuing...")
    elif args.skip_bert:
        print("\n⏭️ Step 3: Skipped (BERT Validation)")
        # Create basic validation results for subsequent steps
        pipeline.validation_results = pd.DataFrame(pipeline.mti_results)
        if 'Target_Gene' not in pipeline.validation_results.columns:
            pipeline.validation_results['Target_Gene'] = pipeline.validation_results['Gene']
        pipeline.validation_results['Overall_Score'] = 75
    
    # Step 4: LLM Analysis
    if not args.skip_llm and pipeline.validation_results is not None:
        print("\n🧠 Step 4: LLM Analysis")
        success = pipeline.step4_llm_analysis(functions, args.analysis_mode)
        if not success:
            print("⚠️ Step 4 failed, but continuing...")
    elif args.skip_llm:
        print("\n⏭️ Step 4: Skipped (LLM Analysis)")
    
    # Step 5: Network Generation
    if args.generate_network and pipeline.validation_results is not None:
        print("\n🌐 Step 5: Network Generation")
        success = pipeline.step5_network_generation(
            args.mirna_deseq, args.gene_deseq, args.network_score_threshold
        )
        if not success:
            print("⚠️ Step 5 failed, but continuing...")
    
    # Generate Final Report
    print("\n📄 Generating Final Report")
    pipeline.generate_final_report(args.mode, args.analysis_mode)
    
    print("\n" + "="*80)
    print("🎉 PIPELINE EXECUTION COMPLETED!")
    print("="*80)
    print(f"Results saved in: {pipeline.results_dir}")
    print(f"Check pipeline_report.txt for detailed results")
    print("Thank you for using the Enhanced Modular miRNA Research Pipeline!")
    
    return 0

if __name__ == "__main__":
    exit(main())
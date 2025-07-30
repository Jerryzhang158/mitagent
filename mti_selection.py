"""
MTI选择模块 - miRNA Research Pipeline
Step 1: 从多个数据库中选择miRNA-Target相互作用
"""

import os
import pandas as pd
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict

from config import PipelineConfig
from mirna_matcher import MiRNAMatcher
from refseq_cache import RefSeqManager
from utils import Logger, FileUtils, ResultsIntegrator, TextProcessor, DataValidator

class MTISelector:
    """MTI选择器 - 整合多个数据库进行MTI发现"""
    
    def __init__(self, min_databases: int = 2, logger: Optional[Logger] = None):
        self.config = PipelineConfig()
        self.min_databases = min_databases
        self.logger = logger
        self.mirna_matcher = MiRNAMatcher()
        self.refseq_manager = RefSeqManager()
        self.filename_mapping = {}  # 从全称miRNA到简化名称的映射
        self.mirtarbase_df = None
    
    def load_databases(self, targetscan_file: str, mirdb_file: str, 
                      mirwalk_file: str, mirtarbase_file: Optional[str] = None) -> Dict[str, pd.DataFrame]:
        """加载所有数据库文件"""
        databases = {}
        
        # 加载主要数据库
        try:
            databases['targetscan'], _ = FileUtils.smart_read_csv(targetscan_file)
            if self.logger:
                self.logger.info(f"✅ TargetScan loaded: {len(databases['targetscan'])} entries")
        except Exception as e:
            raise Exception(f"Error loading TargetScan: {e}")
        
        try:
            databases['mirdb'], _ = FileUtils.smart_read_csv(mirdb_file)
            if self.logger:
                self.logger.info(f"✅ miRDB loaded: {len(databases['mirdb'])} entries")
        except Exception as e:
            raise Exception(f"Error loading miRDB: {e}")
        
        try:
            databases['mirwalk'], _ = FileUtils.smart_read_csv(mirwalk_file)
            if self.logger:
                self.logger.info(f"✅ miRWalk loaded: {len(databases['mirwalk'])} entries")
        except Exception as e:
            raise Exception(f"Error loading miRWalk: {e}")
        
        # 加载miRTarBase（可选）
        if mirtarbase_file:
            databases['mirtarbase'] = self._load_mirtarbase(mirtarbase_file)
            if databases['mirtarbase'] is not None:
                self.mirtarbase_df = databases['mirtarbase']
        
        # 构建miRNA匹配索引
        self.mirna_matcher.build_database_index(
            databases['targetscan'], 
            databases['mirdb'], 
            databases['mirwalk'], 
            databases.get('mirtarbase')
        )
        
        return databases
    
    def _load_mirtarbase(self, mirtarbase_file: str) -> Optional[pd.DataFrame]:
        """加载miRTarBase数据 - 修复版：智能分隔符检测和列名映射"""
        if not mirtarbase_file or not os.path.exists(mirtarbase_file):
            if self.logger:
                self.logger.warning("miRTarBase file not provided or not found, experimental validation will be skipped")
            return None
        
        try:
            if self.logger:
                self.logger.info(f"Loading miRTarBase from: {mirtarbase_file}")
            
            # 使用工具函数智能读取
            df, metadata = FileUtils.smart_read_csv(mirtarbase_file)
            
            if self.logger:
                self.logger.info(f"miRTarBase loaded: {len(df)} total entries")
                self.logger.info(f"Detected delimiter: '{metadata['delimiter']}'")
                self.logger.info(f"Actual columns: {list(df.columns)}")
            
            # 检查必需列
            required_columns = ['miRNA', 'Target Gene']
            missing_cols = [col for col in required_columns if col not in df.columns]
            
            if missing_cols:
                if self.logger:
                    self.logger.info(f"Missing required columns: {missing_cols}")
                
                # 应用列名映射
                df = DataValidator.apply_column_mapping(df, self.config.MIRTARBASE_COLUMN_MAPPING)
                
                # 最终检查
                missing_cols = [col for col in required_columns if col not in df.columns]
                if missing_cols:
                    if self.logger:
                        self.logger.error(f"Still missing required columns after mapping: {missing_cols}")
                        self.logger.info(f"Available columns: {list(df.columns)}")
                    return None
            
            # 数据清理
            df = df.dropna(subset=['miRNA', 'Target Gene'])
            if self.logger:
                self.logger.info(f"After removing NaN values: {len(df)} entries")
            
            # 分析物种分布（仅用于统计）
            if 'miRNA' in df.columns:
                species_counts = df['miRNA'].str[:3].value_counts()
                if self.logger:
                    self.logger.info("Species distribution in miRTarBase:")
                    for species, count in species_counts.head(10).items():
                        self.logger.info(f"  {species}: {count} entries")
                
                hsa_count = len(df[df['miRNA'].str.startswith('hsa-', na=False)])
                if self.logger:
                    self.logger.info(f"Human miRNA entries (hsa-): {hsa_count}")
            
            # 统计实验验证类型
            if 'Support Type' in df.columns:
                support_counts = df['Support Type'].value_counts()
                if self.logger:
                    self.logger.info("miRTarBase Support Types:")
                    for support_type, count in support_counts.head(10).items():
                        self.logger.info(f"  {support_type}: {count}")
            
            if self.logger:
                self.logger.info(f"Successfully loaded miRTarBase with multi-species support")
            return df
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error loading miRTarBase: {e}")
                import traceback
                self.logger.debug(f"Full traceback: {traceback.format_exc()}")
            return None
    
    def get_experimental_validation(self, mirna: str, gene: str) -> Dict[str, any]:
        """查询miRTarBase中的实验验证信息 - 增强版支持多物种"""
        if self.mirtarbase_df is None:
            return {
                'has_validation': False,
                'support_type': None,
                'validation_strength': 'No miRTarBase data'
            }
        
        # 查找匹配的记录
        matches = self.mirtarbase_df[
            (self.mirtarbase_df['miRNA'] == mirna) & 
            (self.mirtarbase_df['Target Gene'] == gene)
        ]
        
        if matches.empty:
            # 尝试模糊匹配miRNA名称
            matched_mirna = self.mirna_matcher.find_best_match(mirna, 'mirtarbase')
            if matched_mirna:
                matches = self.mirtarbase_df[
                    (self.mirtarbase_df['miRNA'] == matched_mirna) & 
                    (self.mirtarbase_df['Target Gene'] == gene)
                ]
        
        # 如果仍然没有匹配，尝试跨物种匹配
        if matches.empty and mirna.startswith('hsa-'):
            mirna_core = mirna[4:]  # 移除hsa-前缀
            potential_matches = self.mirtarbase_df[
                (self.mirtarbase_df['miRNA'].str.contains(mirna_core, na=False, regex=False)) & 
                (self.mirtarbase_df['Target Gene'] == gene)
            ]
            if not potential_matches.empty:
                matches = potential_matches
                if self.logger:
                    self.logger.debug(f"Found cross-species validation for {mirna} → {gene}")
        
        if matches.empty:
            return {
                'has_validation': False,
                'support_type': None,
                'validation_strength': 'No experimental validation found'
            }
        
        # 分析支持类型
        support_types = matches['Support Type'].unique()
        
        # 确定验证强度
        has_functional = any('Functional MTI' in st and 'Weak' not in st for st in support_types)
        has_functional_weak = any('Functional MTI (Weak)' in st for st in support_types)
        has_nonfunctional = any('Non-Functional MTI' in st for st in support_types)
        
        if has_functional:
            validation_strength = 'Strong experimental validation'
            has_validation = True
        elif has_functional_weak:
            validation_strength = 'Weak experimental validation'
            has_validation = True
        elif has_nonfunctional:
            validation_strength = 'Non-functional (experimental negative)'
            has_validation = False
        else:
            validation_strength = 'Experimental evidence (unclassified)'
            has_validation = True
        
        return {
            'has_validation': has_validation,
            'support_type': '; '.join(support_types),
            'validation_strength': validation_strength,
            'experiments': '; '.join(matches['Experiments'].unique()) if 'Experiments' in matches.columns else 'N/A',
            'pmids': '; '.join(matches['References (PMID)'].astype(str).unique()) if 'References (PMID)' in matches.columns else 'N/A'
        }
    
    def find_mirna_for_gene(self, gene_name: str, databases: Dict[str, pd.DataFrame]) -> Dict[str, any]:
        """增强版基因到miRNA搜索，支持模糊匹配，返回原始全称miRNA"""
        refseq_ids = self.refseq_manager.convert_gene_to_refseq(gene_name)
        if not refseq_ids:
            return None
        
        # 从各数据库搜索
        mirdb_df = databases['mirdb']
        targetscan_df = databases['targetscan']
        mirwalk_df = databases['mirwalk']
        
        # 直接搜索
        mirdb_relevant = mirdb_df[mirdb_df['RefSeq'].isin(refseq_ids)]
        mirdb_mirnas = mirdb_relevant[['miRNA', 'Score']].drop_duplicates().values.tolist()

        targetscan_relevant = targetscan_df[targetscan_df['Gene Symbol'] == gene_name]
        targetscan_mirnas = targetscan_relevant[['miRNA', 'context++ score percentile']].drop_duplicates().values.tolist()

        mirwalk_relevant = mirwalk_df[mirwalk_df['Genesymbol'] == gene_name]
        mirwalk_mirnas = mirwalk_relevant[['miRNA', 'binding_probability']].drop_duplicates().values.tolist()

        # 标准化miRNA名称用于交集计算，但保留原始全称
        mirwalk_normalized = set()
        targetscan_normalized = set()
        mirdb_normalized = set()
        
        # 创建标准化到原始名称的映射
        mirwalk_mapping = {}
        targetscan_mapping = {}
        mirdb_mapping = {}
        
        for miRNA, _ in mirwalk_mirnas:
            normalized = self.mirna_matcher.normalize_mirna_name(miRNA)
            mirwalk_normalized.add(normalized)
            mirwalk_mapping[normalized] = miRNA
        
        for miRNA, _ in targetscan_mirnas:
            normalized = self.mirna_matcher.normalize_mirna_name(miRNA)
            targetscan_normalized.add(normalized)
            targetscan_mapping[normalized] = miRNA
        
        for miRNA, _ in mirdb_mirnas:
            normalized = self.mirna_matcher.normalize_mirna_name(miRNA)
            mirdb_normalized.add(normalized)
            mirdb_mapping[normalized] = miRNA
        
        # 计算交集
        mirwalk_targetscan_intersection = mirwalk_normalized.intersection(targetscan_normalized)
        mirwalk_mirdb_intersection = mirwalk_normalized.intersection(mirdb_normalized)
        mirdb_targetscan_intersection = mirdb_normalized.intersection(targetscan_normalized)
        all_three_intersection = mirwalk_normalized.intersection(targetscan_normalized, mirdb_normalized)

        # 转换回原始名称
        def get_best_original_name(normalized_name):
            candidates = []
            if normalized_name in mirwalk_mapping:
                candidates.append(mirwalk_mapping[normalized_name])
            if normalized_name in targetscan_mapping:
                candidates.append(targetscan_mapping[normalized_name])
            if normalized_name in mirdb_mapping:
                candidates.append(mirdb_mapping[normalized_name])
            
            return max(candidates, key=len) if candidates else normalized_name
        
        # 转换交集结果
        mirwalk_targetscan_final = {get_best_original_name(n) for n in mirwalk_targetscan_intersection}
        mirwalk_mirdb_final = {get_best_original_name(n) for n in mirwalk_mirdb_intersection}
        mirdb_targetscan_final = {get_best_original_name(n) for n in mirdb_targetscan_intersection}
        all_three_final = {get_best_original_name(n) for n in all_three_intersection}

        return {
            'targetscan_mirnas': targetscan_mirnas,
            'mirdb_mirnas': mirdb_mirnas,
            'mirwalk_mirnas': mirwalk_mirnas,
            'mirwalk_targetscan_intersection': mirwalk_targetscan_final,
            'mirwalk_mirdb_intersection': mirwalk_mirdb_final,
            'mirdb_targetscan_intersection': mirdb_targetscan_final,
            'all_three_intersection': all_three_final
        }
    
    def find_genes_for_mirna(self, mirna_name: str, databases: Dict[str, pd.DataFrame], 
                           gene_filter: Optional[Set[str]] = None) -> List[Dict[str, any]]:
        """增强版miRNA到基因搜索，支持模糊匹配，使用原始全称miRNA"""
        all_genes_found = []
        
        # 获取在各数据库中的最佳匹配名称
        matches = self.mirna_matcher.find_matches_all_databases(mirna_name)
        
        mirdb_df = databases['mirdb']
        targetscan_df = databases['targetscan']
        mirwalk_df = databases['mirwalk']
        
        # Step 1: miRDB搜索
        if 'mirdb' in matches:
            mirdb_mirna = matches['mirdb']
            mirdb_relevant = mirdb_df[mirdb_df['miRNA'] == mirdb_mirna]
            
            if not mirdb_relevant.empty:
                refseq_ids = mirdb_relevant['RefSeq'].unique().tolist()
                refseq_to_gene = self.refseq_manager.convert_refseq_to_genes(refseq_ids)
                
                for _, row in mirdb_relevant.iterrows():
                    gene_symbol = refseq_to_gene.get(row['RefSeq'])
                    if gene_symbol:
                        if gene_filter is None or gene_symbol in gene_filter:
                            all_genes_found.append({
                                'gene': gene_symbol,
                                'score': row['Score'],
                                'database': 'miRDB',
                                'database_count': 1,
                                'matched_mirna': mirdb_mirna
                            })
        
        # Step 2: TargetScan搜索
        if 'targetscan' in matches:
            targetscan_mirna = matches['targetscan']
            targetscan_relevant = targetscan_df[targetscan_df['miRNA'] == targetscan_mirna]
            for _, row in targetscan_relevant.iterrows():
                gene_symbol = row['Gene Symbol']
                if gene_filter is None or gene_symbol in gene_filter:
                    all_genes_found.append({
                        'gene': gene_symbol,
                        'score': row['context++ score percentile'],
                        'database': 'TargetScan',
                        'database_count': 1,
                        'matched_mirna': targetscan_mirna
                    })
        
        # Step 3: miRWalk搜索
        if 'mirwalk' in matches:
            mirwalk_mirna = matches['mirwalk']
            mirwalk_relevant = mirwalk_df[mirwalk_df['miRNA'] == mirwalk_mirna]
            for _, row in mirwalk_relevant.iterrows():
                gene_symbol = row['Genesymbol']
                if gene_filter is None or gene_symbol in gene_filter:
                    all_genes_found.append({
                        'gene': gene_symbol,
                        'score': row['binding_probability'] * 100,
                        'database': 'miRWalk',
                        'database_count': 1,
                        'matched_mirna': mirwalk_mirna
                    })
        
        # 记录匹配信息
        if matches and self.logger:
            match_info = ", ".join([f"{db}: {mirna}" for db, mirna in matches.items()])
            self.logger.info(f"  miRNA matches found - {match_info}")
        elif self.logger:
            self.logger.info(f"  No matches found for {mirna_name}")
        
        # 合并相同基因的结果
        gene_db_mapping = defaultdict(set)
        gene_scores = defaultdict(list)
        gene_matched_mirnas = defaultdict(set)
        
        for result in all_genes_found:
            gene = result['gene']
            gene_db_mapping[gene].add(result['database'])
            gene_scores[gene].append((result['database'], result['score']))
            gene_matched_mirnas[gene].add(result['matched_mirna'])
        
        final_results = []
        for gene, databases in gene_db_mapping.items():
            db_count = len(databases)
            db_sources = ','.join(sorted(databases))
            
            # 获取最佳匹配的miRNA全称
            best_mirna = max(gene_matched_mirnas[gene], key=len) if gene_matched_mirnas[gene] else mirna_name
            
            if db_count >= 3:
                priority = 'High'
            elif db_count == 2:
                priority = 'Medium'
            else:
                priority = 'Low'
            
            final_results.append({
                'Gene': gene,
                'miRNA': best_mirna,
                'Direction': 'miRNA→Gene',
                'Databases': db_count,
                'Priority': priority,
                'Database_Sources': db_sources
            })
        
        return final_results
    
    def gene_to_mirna_selection(self, gene_list_file: str, databases: Dict[str, pd.DataFrame]) -> List[Dict[str, any]]:
        """Step 1A: 从基因找miRNA"""
        if self.logger:
            self.logger.info("="*50)
            self.logger.info("Step 1A: Gene to miRNA Selection")
            self.logger.info("="*50)
        
        # 读取基因列表
        genes = FileUtils.load_gene_list(gene_list_file)
        if not genes:
            raise ValueError(f"Could not load genes from {gene_list_file}")
        
        genes = list(genes)
        
        if self.logger:
            self.logger.info(f"Processing {len(genes)} genes...")
            self.logger.info(f"Minimum database requirement: {self.min_databases}")
        
        mti_results = []
        
        # 处理每个基因
        for gene_name in genes:
            if self.logger:
                self.logger.info(f"Processing gene: {gene_name}")
            
            results = self.find_mirna_for_gene(gene_name, databases)
            
            if results is None:
                continue
            
            # 提取交集结果
            all_inter = results['all_three_intersection']
            mw_ts = results['mirwalk_targetscan_intersection']
            mw_md = results['mirwalk_mirdb_intersection']
            md_ts = results['mirdb_targetscan_intersection']
            
            # 收集3个数据库都有的MTI
            for miRNA in all_inter:
                validation_info = self.get_experimental_validation(miRNA, gene_name)
                
                # 建立文件名映射
                simplified_mirna = TextProcessor.extract_mirna_core(miRNA)
                if simplified_mirna:
                    file_key = f"{gene_name}_{simplified_mirna}"
                    self.filename_mapping[f"{gene_name}_{miRNA}"] = file_key
                
                mti_results.append({
                    'Gene': gene_name,
                    'miRNA': miRNA,
                    'Direction': 'Gene→miRNA',
                    'Databases': 3,
                    'Priority': 'High',
                    'Database_Sources': 'TargetScan,miRDB,miRWalk',
                    'miRTarBase_Validation': validation_info['has_validation'],
                    'miRTarBase_Support_Type': validation_info['support_type'],
                    'Validation_Strength': validation_info['validation_strength']
                })
            
            # 如果min_databases <= 2，也收集2个数据库的交集
            if self.min_databases <= 2:
                # 处理2数据库交集
                for intersection, db_sources in [
                    (mw_ts - all_inter, 'TargetScan,miRWalk'),
                    (mw_md - all_inter, 'miRDB,miRWalk'),
                    (md_ts - all_inter, 'TargetScan,miRDB')
                ]:
                    for miRNA in intersection:
                        validation_info = self.get_experimental_validation(miRNA, gene_name)
                        
                        simplified_mirna = TextProcessor.extract_mirna_core(miRNA)
                        if simplified_mirna:
                            file_key = f"{gene_name}_{simplified_mirna}"
                            self.filename_mapping[f"{gene_name}_{miRNA}"] = file_key
                        
                        mti_results.append({
                            'Gene': gene_name,
                            'miRNA': miRNA,
                            'Direction': 'Gene→miRNA',
                            'Databases': 2,
                            'Priority': 'Medium',
                            'Database_Sources': db_sources,
                            'miRTarBase_Validation': validation_info['has_validation'],
                            'miRTarBase_Support_Type': validation_info['support_type'],
                            'Validation_Strength': validation_info['validation_strength']
                        })
        
        if self.logger:
            self.logger.info(f"Found {len(mti_results)} MTIs from gene-to-miRNA search")
            self.logger.info(f"  - In 3 databases: {sum(1 for m in mti_results if m['Databases'] == 3)}")
            self.logger.info(f"  - In 2 databases: {sum(1 for m in mti_results if m['Databases'] == 2)}")
            
            validated_count = sum(1 for m in mti_results if m['miRTarBase_Validation'])
            self.logger.info(f"  - With experimental validation: {validated_count}")
        
        return mti_results
    
    def mirna_to_gene_selection(self, mirna_list_file: str, databases: Dict[str, pd.DataFrame], 
                               gene_list_file: Optional[str] = None) -> List[Dict[str, any]]:
        """Step 1B: 从miRNA找基因"""
        if self.logger:
            self.logger.info("="*50)
            self.logger.info("Step 1B: miRNA to Gene Selection (Modified)")
            self.logger.info("="*50)
        
        # 读取miRNA列表
        mirnas = FileUtils.load_gene_list(mirna_list_file)  # 重用函数
        if not mirnas:
            raise ValueError(f"Could not load miRNAs from {mirna_list_file}")
        
        mirnas = list(mirnas)
        
        # 加载基因过滤列表
        gene_filter = None
        if gene_list_file:
            gene_filter = FileUtils.load_gene_list(gene_list_file)
            if gene_filter and self.logger:
                self.logger.info(f"Gene filter loaded: {len(gene_filter)} genes")
                self.logger.info("Only genes in the gene list will be included in results")
            elif self.logger:
                self.logger.warning("Gene list file specified but could not be loaded")
        elif self.logger:
            self.logger.info("No gene filter specified - all genes will be included")
        
        if self.logger:
            self.logger.info(f"Processing {len(mirnas)} miRNAs...")
            self.logger.info("Modified mode: Collecting MTIs from ANY database (minimum 1 database)")
        
        all_mti_results = []
        
        # 处理每个miRNA
        for mirna_name in mirnas:
            if self.logger:
                self.logger.info(f"Processing miRNA: {mirna_name}")
            
            mti_results = self.find_genes_for_mirna(mirna_name, databases, gene_filter)
            
            # 为每个MTI添加实验验证信息和文件名映射
            for mti in mti_results:
                validation_info = self.get_experimental_validation(mti['miRNA'], mti['Gene'])
                mti['miRTarBase_Validation'] = validation_info['has_validation']
                mti['miRTarBase_Support_Type'] = validation_info['support_type']
                mti['Validation_Strength'] = validation_info['validation_strength']
                
                # 建立文件名映射
                simplified_mirna = TextProcessor.extract_mirna_core(mti['miRNA'])
                if simplified_mirna:
                    file_key = f"{mti['Gene']}_{simplified_mirna}"
                    self.filename_mapping[f"{mti['Gene']}_{mti['miRNA']}"] = file_key
            
            all_mti_results.extend(mti_results)
        
        if self.logger:
            self.logger.info(f"Found {len(all_mti_results)} MTIs from miRNA-to-gene search")
            self.logger.info(f"  - In 3 databases: {sum(1 for m in all_mti_results if m['Databases'] == 3)}")
            self.logger.info(f"  - In 2 databases: {sum(1 for m in all_mti_results if m['Databases'] == 2)}")
            self.logger.info(f"  - In 1 database: {sum(1 for m in all_mti_results if m['Databases'] == 1)}")
            
            validated_count = sum(1 for m in all_mti_results if m['miRTarBase_Validation'])
            self.logger.info(f"  - With experimental validation: {validated_count}")
            
            if gene_filter:
                unique_genes = set(m['Gene'] for m in all_mti_results)
                self.logger.info(f"  - Unique genes found (all in gene list): {len(unique_genes)}")
        
        # 保存RefSeq缓存
        self.refseq_manager.cleanup()
        
        return all_mti_results
    
    def combined_selection(self, gene_list_file: str, mirna_list_file: str, 
                          databases: Dict[str, pd.DataFrame]) -> List[Dict[str, any]]:
        """Step 1C: 组合基因和miRNA的搜索结果"""
        if self.logger:
            self.logger.info("="*50)
            self.logger.info("Step 1C: Combined Selection (Gene + miRNA)")
            self.logger.info("="*50)
        
        # 执行两种搜索
        gene_to_mirna_results = []
        mirna_to_gene_results = []
        
        if gene_list_file and os.path.exists(gene_list_file):
            gene_to_mirna_results = self.gene_to_mirna_selection(gene_list_file, databases)
        
        if mirna_list_file and os.path.exists(mirna_list_file):
            # 传入gene_list_file作为过滤器
            mirna_to_gene_results = self.mirna_to_gene_selection(mirna_list_file, databases, gene_list_file)
        
        # 合并结果并去重
        all_mti_results = gene_to_mirna_results + mirna_to_gene_results
        mti_results = ResultsIntegrator.merge_mti_results(all_mti_results)
        
        if self.logger:
            self.logger.info(f"Total unique MTIs after combining: {len(mti_results)}")
            self.logger.info(f"  - From gene search only: {sum(1 for m in mti_results if m['Direction'] == 'Gene→miRNA')}")
            self.logger.info(f"  - From miRNA search only: {sum(1 for m in mti_results if m['Direction'] == 'miRNA→Gene')}")
            self.logger.info(f"  - Bidirectional (found by both): {sum(1 for m in mti_results if m['Direction'] == 'Bidirectional')}")
            self.logger.info(f"  - In 3 databases: {sum(1 for m in mti_results if m['Databases'] == 3)}")
            self.logger.info(f"  - In 2 databases: {sum(1 for m in mti_results if m['Databases'] == 2)}")
            self.logger.info(f"  - In 1 database: {sum(1 for m in mti_results if m['Databases'] == 1)}")
            
            validated_count = sum(1 for m in mti_results if m.get('miRTarBase_Validation', False))
            self.logger.info(f"  - With experimental validation: {validated_count}")
        
        return mti_results
    
    def save_results(self, mti_results: List[Dict[str, any]], output_file: str, mode: str):
        """保存MTI选择结果"""
        # 转换为DataFrame
        df = pd.DataFrame(mti_results)
        
        # 创建统计摘要
        summary_data = ResultsIntegrator.create_summary_stats(mti_results)
        
        # 创建验证统计
        validation_stats = {}
        if mti_results:
            for strength in set(m.get('Validation_Strength', 'Unknown') for m in mti_results):
                validation_stats[strength] = sum(1 for m in mti_results if m.get('Validation_Strength') == strength)
        
        # 保存到Excel
        FileUtils.save_excel_with_summary(df, output_file, summary_data, validation_stats)
        
        if self.logger:
            self.logger.info(f"MTI selection results saved to {output_file}")
        
        return df
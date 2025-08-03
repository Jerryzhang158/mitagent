"""
miRNA名称匹配器 - miRNA Research Pipeline
处理不同数据库中miRNA名称的标准化和模糊匹配
"""

import re
import pandas as pd
from difflib import SequenceMatcher
from config import PipelineConfig
from typing import Optional, Dict, Set, List

class MiRNAMatcher:
    """miRNA名称标准化和模糊匹配工具 - 增强版支持miRTarBase"""
    
    def __init__(self):
        self.mirna_cache = {}  # 缓存标准化结果
        self.database_mirnas = {}  # 存储各数据库的miRNA名称
        self.original_names_map = {}  # 保存原始全称映射
        self.config = PipelineConfig()
        
    def normalize_mirna_name(self, mirna_name: str) -> Optional[str]:
        """标准化miRNA名称，移除版本号和多余字符，但保留原始全称信息"""
        if not mirna_name:
            return None
            
        # 缓存检查
        if mirna_name in self.mirna_cache:
            return self.mirna_cache[mirna_name]
        
        # 保存原始名称
        original = mirna_name
        
        # 转换为小写并去除空格
        normalized = mirna_name.lower().strip()
        
        # 移除版本号 (.1, .2等)
        normalized = re.sub(self.config.MIRNA_PATTERNS['version_pattern'], '', normalized)
        
        # 统一前缀格式
        if not normalized.startswith('hsa-'):
            if normalized.startswith('mir-') or normalized.startswith('let-'):
                normalized = 'hsa-' + normalized
            elif re.match(r'^\d', normalized):  # 以数字开头
                normalized = 'hsa-mir-' + normalized
        
        # 标准化命名格式
        normalized = re.sub(r'hsa-mir-', 'hsa-mir-', normalized)
        normalized = re.sub(r'hsa-let-', 'hsa-let-', normalized)
        
        # 缓存结果
        self.mirna_cache[mirna_name] = normalized
        
        # 保存原始名称映射
        self.original_names_map[normalized] = original
        
        return normalized
    
    def get_original_name(self, normalized_name: str) -> str:
        """获取miRNA的原始全称"""
        return self.original_names_map.get(normalized_name, normalized_name)
    
    def extract_core_pattern(self, mirna_name: str) -> Optional[str]:
        """提取miRNA的核心模式，用于模糊匹配"""
        if not mirna_name:
            return None
            
        # 标准化
        normalized = self.normalize_mirna_name(mirna_name)
        
        # 提取核心部分 (如: hsa-mir-183-5p -> mir-183-5p)
        match = re.search(self.config.MIRNA_PATTERNS['core_pattern'], normalized)
        if match:
            return match.group(1)
        
        # 处理let家族
        match = re.search(self.config.MIRNA_PATTERNS['let_pattern'], normalized)
        if match:
            return match.group(1)
            
        return normalized
    
    def build_database_index(self, targetscan_df: pd.DataFrame, mirdb_df: pd.DataFrame, 
                           mirwalk_df: pd.DataFrame, mirtarbase_df: Optional[pd.DataFrame] = None):
        """构建数据库miRNA索引，包括miRTarBase"""
        print("Building miRNA database index...")
        
        self.database_mirnas = {
            'targetscan': {},
            'mirdb': {},
            'mirwalk': {},
            'mirtarbase': {}
        }
        
        # TargetScan
        self._index_database_mirnas('targetscan', targetscan_df, 
                                   self.config.DATABASE_COLUMNS['targetscan']['mirna'])
        
        # miRDB
        self._index_database_mirnas('mirdb', mirdb_df, 
                                   self.config.DATABASE_COLUMNS['mirdb']['mirna'])
        
        # miRWalk
        self._index_database_mirnas('mirwalk', mirwalk_df, 
                                   self.config.DATABASE_COLUMNS['mirwalk']['mirna'])
        
        # miRTarBase - 修复版：支持多物种
        if mirtarbase_df is not None:
            self._index_database_mirnas('mirtarbase', mirtarbase_df, 
                                       self.config.DATABASE_COLUMNS['mirtarbase']['mirna'])
        
        total_indexed = sum(len(db) for db in self.database_mirnas.values())
        print(f"Indexed {total_indexed} miRNA variants across all databases")
        for db_name, db_mirnas in self.database_mirnas.items():
            print(f"  - {db_name}: {len(db_mirnas)} variants")
    
    def _index_database_mirnas(self, db_name: str, df: pd.DataFrame, mirna_column: str):
        """为单个数据库建立miRNA索引"""
        if mirna_column not in df.columns:
            print(f"Warning: Column '{mirna_column}' not found in {db_name} database")
            return
        
        unique_mirnas = df[mirna_column].unique()
        for mirna in unique_mirnas:
            if pd.notna(mirna):
                normalized = self.normalize_mirna_name(str(mirna))
                core = self.extract_core_pattern(str(mirna))
                self.database_mirnas[db_name][normalized] = str(mirna)
                if core and core != normalized:
                    self.database_mirnas[db_name][core] = str(mirna)
    
    def find_best_match(self, query_mirna: str, database_name: str, 
                       similarity_threshold: float = None) -> Optional[str]:
        """在指定数据库中找到最佳匹配的miRNA，返回原始全称"""
        if similarity_threshold is None:
            similarity_threshold = self.config.DEFAULT_SIMILARITY_THRESHOLD
            
        if not query_mirna or database_name not in self.database_mirnas:
            return None
        
        # 1. 直接精确匹配（标准化后）
        normalized_query = self.normalize_mirna_name(query_mirna)
        if normalized_query in self.database_mirnas[database_name]:
            matched_mirna = self.database_mirnas[database_name][normalized_query]
            return matched_mirna
        
        # 2. 核心模式匹配
        core_query = self.extract_core_pattern(query_mirna)
        if core_query and core_query in self.database_mirnas[database_name]:
            matched_mirna = self.database_mirnas[database_name][core_query]
            return matched_mirna
        
        # 3. 模糊匹配
        best_match = None
        best_score = 0
        
        for db_mirna_key, original_mirna in self.database_mirnas[database_name].items():
            # 计算相似度
            similarity = SequenceMatcher(None, normalized_query, db_mirna_key).ratio()
            
            if similarity > best_score and similarity >= similarity_threshold:
                best_score = similarity
                best_match = original_mirna
        
        return best_match
    
    def find_matches_all_databases(self, query_mirna: str, 
                                 similarity_threshold: float = None) -> Dict[str, str]:
        """在所有数据库中查找匹配，返回原始全称"""
        if similarity_threshold is None:
            similarity_threshold = self.config.DEFAULT_SIMILARITY_THRESHOLD
            
        matches = {}
        for db_name in self.database_mirnas.keys():
            match = self.find_best_match(query_mirna, db_name, similarity_threshold)
            if match:
                matches[db_name] = match
        return matches
    
    def get_database_statistics(self) -> Dict[str, Dict[str, int]]:
        """获取数据库统计信息"""
        stats = {}
        for db_name, db_mirnas in self.database_mirnas.items():
            stats[db_name] = {
                'total_variants': len(db_mirnas),
                'unique_normalized': len(set(self.normalize_mirna_name(mirna) 
                                           for mirna in db_mirnas.values()))
            }
        return stats
    
    def validate_mirna_name(self, mirna_name: str) -> Dict[str, any]:
        """验证miRNA名称格式"""
        validation_result = {
            'is_valid': False,
            'normalized_name': None,
            'core_pattern': None,
            'has_species_prefix': False,
            'issues': []
        }
        
        if not mirna_name:
            validation_result['issues'].append('Empty miRNA name')
            return validation_result
        
        # 标准化
        normalized = self.normalize_mirna_name(mirna_name)
        validation_result['normalized_name'] = normalized
        
        # 提取核心模式
        core = self.extract_core_pattern(mirna_name)
        validation_result['core_pattern'] = core
        
        # 检查物种前缀
        validation_result['has_species_prefix'] = mirna_name.lower().startswith('hsa-')
        
        # 检查格式
        if not re.match(r'hsa-(mir|let)-\d+', normalized or ''):
            validation_result['issues'].append('Invalid miRNA format')
        else:
            validation_result['is_valid'] = True
        
        return validation_result
    
    def suggest_corrections(self, query_mirna: str, max_suggestions: int = 3) -> List[str]:
        """为不匹配的miRNA名称建议可能的修正"""
        suggestions = []
        
        # 从所有数据库中收集候选项
        all_mirnas = set()
        for db_mirnas in self.database_mirnas.values():
            all_mirnas.update(db_mirnas.values())
        
        # 计算相似度并排序
        similarities = []
        normalized_query = self.normalize_mirna_name(query_mirna)
        
        for mirna in all_mirnas:
            normalized_mirna = self.normalize_mirna_name(mirna)
            similarity = SequenceMatcher(None, normalized_query, normalized_mirna).ratio()
            if similarity > 0.5:  # 只考虑相似度超过50%的
                similarities.append((similarity, mirna))
        
        # 按相似度降序排序，返回前N个建议
        similarities.sort(key=lambda x: x[0], reverse=True)
        suggestions = [mirna for _, mirna in similarities[:max_suggestions]]
        
        return suggestions
    
    def clear_cache(self):
        """清空缓存"""
        self.mirna_cache.clear()
        self.original_names_map.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """获取缓存统计信息"""
        return {
            'cached_normalizations': len(self.mirna_cache),
            'original_name_mappings': len(self.original_names_map)
        }
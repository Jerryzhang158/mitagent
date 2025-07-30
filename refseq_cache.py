"""
RefSeq缓存管理器 - miRNA Research Pipeline
管理RefSeq ID到基因名称的转换和缓存
"""

import json
import os
import requests
import time
from typing import Dict, List, Optional, Any
from config import PipelineConfig

class RefSeqCache:
    """缓存RefSeq到基因名的转换结果"""
    
    def __init__(self, cache_file: Optional[str] = None):
        self.config = PipelineConfig()
        self.cache_file = cache_file or self.config.get_cache_path()
        self.cache = {}
        self.load_cache()
    
    def load_cache(self):
        """从文件加载缓存"""
        try:
            if os.path.exists(self.cache_file):
                with open(self.cache_file, 'r', encoding=self.config.DEFAULT_ENCODING) as f:
                    self.cache = json.load(f)
                print(f"Loaded {len(self.cache)} cached RefSeq mappings")
            else:
                print("No cache file found, starting fresh")
        except Exception as e:
            print(f"Error loading cache: {e}")
            self.cache = {}
    
    def save_cache(self):
        """保存缓存到文件"""
        try:
            os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
            with open(self.cache_file, 'w', encoding=self.config.DEFAULT_ENCODING) as f:
                json.dump(self.cache, f, ensure_ascii=False, indent=2)
            print(f"Saved {len(self.cache)} RefSeq mappings to cache")
        except Exception as e:
            print(f"Error saving cache: {e}")
    
    def get_gene_name(self, refseq_id: str) -> Optional[str]:
        """获取基因名，优先使用缓存"""
        if refseq_id in self.cache:
            return self.cache[refseq_id]
        return None
    
    def add_mapping(self, refseq_id: str, gene_name: Optional[str]):
        """添加映射到缓存"""
        self.cache[refseq_id] = gene_name
    
    def batch_add_mappings(self, mappings: Dict[str, Optional[str]]):
        """批量添加映射到缓存"""
        self.cache.update(mappings)
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """获取缓存统计信息"""
        total_mappings = len(self.cache)
        valid_mappings = sum(1 for v in self.cache.values() if v is not None)
        null_mappings = total_mappings - valid_mappings
        
        return {
            'total_mappings': total_mappings,
            'valid_mappings': valid_mappings,
            'null_mappings': null_mappings,
            'cache_file': self.cache_file,
            'file_exists': os.path.exists(self.cache_file)
        }
    
    def clear_cache(self):
        """清空缓存"""
        self.cache.clear()
        if os.path.exists(self.cache_file):
            os.remove(self.cache_file)
        print("Cache cleared")

class MyGeneAPI:
    """MyGene.info API接口"""
    
    def __init__(self, cache: RefSeqCache):
        self.config = PipelineConfig()
        self.cache = cache
        self.base_url = self.config.MYGENE_BASE_URL
    
    def gene_name_to_refseq(self, gene_name: str) -> Optional[List[str]]:
        """Convert a gene name to RefSeq ID using MyGene.info API."""
        query = f"symbol:{gene_name}"
        
        params = {
            'q': query,
            'fields': 'refseq.rna',
            'species': 'human',
            'size': 1
        }
        
        try:
            response = requests.get(
                self.base_url, 
                params=params, 
                timeout=self.config.MYGENE_TIMEOUT
            )
            
            if response.status_code == 200:
                data = response.json()
                
                if 'hits' in data and len(data['hits']) > 0:
                    hit = data['hits'][0]
                    rna_ids = hit.get('refseq', {}).get('rna', [])
                    if isinstance(rna_ids, str):
                        rna_ids = [rna_ids]
                    rna_ids = [rna_id.split('.')[0] for rna_id in rna_ids]
                    return rna_ids
                else:
                    print(f"No RefSeq ID found for gene name: {gene_name}")
                    return None
            else:
                print(f"Error getting RefSeq for {gene_name}: {response.status_code}")
                return None
        except Exception as e:
            print(f"Exception getting RefSeq for {gene_name}: {e}")
            return None
    
    def batch_refseq_to_gene_name(self, refseq_list: List[str], 
                                 batch_size: int = None) -> Dict[str, Optional[str]]:
        """批量转换RefSeq ID到基因名，使用缓存优化"""
        if batch_size is None:
            batch_size = self.config.DEFAULT_BATCH_SIZE
            
        results = {}
        uncached_refseqs = []
        
        # 首先检查缓存
        for refseq_id in refseq_list:
            cached_result = self.cache.get_gene_name(refseq_id)
            if cached_result is not None:
                results[refseq_id] = cached_result
            else:
                uncached_refseqs.append(refseq_id)
        
        print(f"Found {len(results)} cached mappings, need to query {len(uncached_refseqs)} RefSeq IDs")
        
        # 批量查询未缓存的RefSeq ID
        if uncached_refseqs:
            for i in range(0, len(uncached_refseqs), batch_size):
                batch = uncached_refseqs[i:i + batch_size]
                batch_results = self._query_mygene_batch(batch)
                
                for refseq_id, gene_name in batch_results.items():
                    results[refseq_id] = gene_name
                    self.cache.add_mapping(refseq_id, gene_name)
                
                # 添加延迟避免API限制
                if i + batch_size < len(uncached_refseqs):
                    time.sleep(self.config.MYGENE_DELAY)
                
                print(f"Processed batch {i//batch_size + 1}/{(len(uncached_refseqs)-1)//batch_size + 1}")
        
        # 保存缓存
        self.cache.save_cache()
        
        return results
    
    def _query_mygene_batch(self, refseq_list: List[str]) -> Dict[str, Optional[str]]:
        """批量查询MyGene.info API"""
        # 构造批量查询
        query_terms = [f"refseq.rna:{refseq}" for refseq in refseq_list]
        query = " OR ".join(query_terms)
        
        params = {
            'q': query,
            'fields': 'symbol,refseq.rna',
            'species': 'human',
            'size': len(refseq_list)
        }
        
        try:
            response = requests.get(
                self.base_url, 
                params=params, 
                timeout=self.config.DEFAULT_TIMEOUT
            )
            
            if response.status_code == 200:
                data = response.json()
                results = {}
                
                if 'hits' in data:
                    for hit in data['hits']:
                        symbol = hit.get('symbol')
                        refseq_data = hit.get('refseq', {})
                        rna_ids = refseq_data.get('rna', [])
                        
                        if symbol and rna_ids:
                            if isinstance(rna_ids, str):
                                rna_ids = [rna_ids]
                            
                            for rna_id in rna_ids:
                                clean_id = rna_id.split('.')[0]
                                if clean_id in refseq_list:
                                    results[clean_id] = symbol
                
                # 对于没有找到的RefSeq ID，标记为None
                for refseq_id in refseq_list:
                    if refseq_id not in results:
                        results[refseq_id] = None
                
                return results
            else:
                print(f"API error: {response.status_code}")
                return {refseq_id: None for refseq_id in refseq_list}
        
        except Exception as e:
            print(f"Error querying API: {e}")
            return {refseq_id: None for refseq_id in refseq_list}
    
    def validate_gene_symbols(self, gene_symbols: List[str]) -> Dict[str, bool]:
        """验证基因符号是否有效"""
        results = {}
        
        for gene_symbol in gene_symbols:
            query = f"symbol:{gene_symbol}"
            params = {
                'q': query,
                'species': 'human',
                'size': 1
            }
            
            try:
                response = requests.get(
                    self.base_url, 
                    params=params, 
                    timeout=self.config.MYGENE_TIMEOUT
                )
                
                if response.status_code == 200:
                    data = response.json()
                    results[gene_symbol] = len(data.get('hits', [])) > 0
                else:
                    results[gene_symbol] = False
                    
                # 添加短暂延迟
                time.sleep(0.1)
                
            except Exception as e:
                print(f"Error validating gene symbol {gene_symbol}: {e}")
                results[gene_symbol] = False
        
        return results

class RefSeqManager:
    """RefSeq管理器，整合缓存和API功能"""
    
    def __init__(self, cache_file: Optional[str] = None):
        self.cache = RefSeqCache(cache_file)
        self.api = MyGeneAPI(self.cache)
    
    def convert_refseq_to_genes(self, refseq_list: List[str], 
                               batch_size: int = None) -> Dict[str, Optional[str]]:
        """转换RefSeq列表到基因名称"""
        return self.api.batch_refseq_to_gene_name(refseq_list, batch_size)
    
    def convert_gene_to_refseq(self, gene_name: str) -> Optional[List[str]]:
        """转换基因名称到RefSeq列表"""
        return self.api.gene_name_to_refseq(gene_name)
    
    def get_cache_info(self) -> Dict[str, Any]:
        """获取缓存信息"""
        return self.cache.get_cache_stats()
    
    def validate_genes(self, gene_list: List[str]) -> Dict[str, bool]:
        """验证基因列表"""
        return self.api.validate_gene_symbols(gene_list)
    
    def cleanup(self):
        """清理资源"""
        self.cache.save_cache()
"""
文献挖掘模块 - miRNA Research Pipeline
Step 2: 从PubMed获取相关文献摘要
"""

import requests
import xml.etree.ElementTree as ET
import time
import os
import pandas as pd
from typing import List, Dict, Optional
from urllib.parse import quote

from config import PipelineConfig
from utils import Logger, TextProcessor

class PubMedMiner:
    """PubMed文献挖掘器"""
    
    def __init__(self, output_dir: str, logger: Optional[Logger] = None):
        self.config = PipelineConfig()
        self.output_dir = output_dir
        self.logger = logger
        self.base_url = self.config.PUBMED_BASE_URL
        
        # 确保输出目录存在
        os.makedirs(self.output_dir, exist_ok=True)
    
    def extract_mirna_for_search(self, mirna: str) -> Optional[str]:
        """提取miRNA名称的核心部分用于PubMed搜索"""
        return TextProcessor.extract_mirna_core(mirna)
    
    def query_pubmed(self, gene_name: str, mirna: str, max_articles: int = 50) -> int:
        """查询PubMed文献"""
        # 提取miRNA搜索词
        mirna_search = self.extract_mirna_for_search(mirna)
        if not mirna_search:
            if self.logger:
                self.logger.warning(f"Could not extract search term from miRNA: {mirna}")
            return 0
        
        # 构建查询
        query_text = f"{gene_name} AND {mirna_search}"
        db = "pubmed"
        query = f"{query_text} AND 1900:2025[pdat]"
        
        if self.logger:
            self.logger.info(f"  Querying: {query_text} (full miRNA: {mirna})")
        
        # 搜索URL
        search_url = f"{self.base_url}esearch.fcgi?db={db}&term={quote(query)}&retmode=json&retmax={max_articles}&usehistory=y"
        
        try:
            response = requests.get(search_url, timeout=self.config.PUBMED_TIMEOUT)
            result = response.json()

            total_count = int(result["esearchresult"].get("count", 0))
            if self.logger:
                self.logger.info(f"  Found {total_count} articles for {gene_name} and {mirna_search}")

            idlist = result["esearchresult"]["idlist"]
            
            if not idlist:
                return 0

            # 下载并保存文章
            article_count = self._download_articles(idlist, gene_name, mirna_search, mirna)
            
            return article_count
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"  Error querying PubMed: {e}")
            return 0
    
    def _download_articles(self, idlist: List[str], gene_name: str, 
                          mirna_search: str, mirna_full: str) -> int:
        """下载文章摘要"""
        batch_size = 100
        output_file = os.path.join(self.output_dir, f"{gene_name}_{mirna_search}.txt")
        
        try:
            with open(output_file, 'w', encoding=self.config.DEFAULT_ENCODING) as f:
                for i in range(0, len(idlist), batch_size):
                    batch_ids = idlist[i:i+batch_size]
                    id_string = ",".join(batch_ids)

                    fetch_url = f"{self.base_url}efetch.fcgi?db=pubmed&retmode=xml&id={id_string}&rettype=medline"
                    response = requests.get(fetch_url, timeout=self.config.PUBMED_TIMEOUT)

                    try:
                        root = ET.fromstring(response.text.encode(self.config.DEFAULT_ENCODING))

                        for article in root.findall(".//PubmedArticle"):
                            # Title
                            title = article.find(".//ArticleTitle")
                            if title is not None and title.text is not None:
                                f.write(title.text.strip() + "\n")
                            else:
                                f.write("[No title available]\n")

                            # Abstract
                            f.write("Abstract\n")
                            abstract_elements = article.findall(".//AbstractText")
                            if abstract_elements:
                                for abstract in abstract_elements:
                                    abstract_text = ''.join(abstract.itertext()).strip()
                                    f.write(abstract_text + "\n")
                            else:
                                f.write("[No abstract available]\n")

                            f.write("\n")

                    except ET.ParseError as e:
                        if self.logger:
                            self.logger.warning(f"  Parse error: {e}")
                        continue

                    time.sleep(self.config.API_DELAY)  # 避免过快请求
            
            return len(idlist)
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"  Error downloading articles: {e}")
            return 0
    
    def validate_query_parameters(self, gene_name: str, mirna: str) -> Dict[str, any]:
        """验证查询参数"""
        validation = {
            'is_valid': True,
            'issues': []
        }
        
        if not gene_name or not gene_name.strip():
            validation['is_valid'] = False
            validation['issues'].append('Empty gene name')
        
        if not mirna or not mirna.strip():
            validation['is_valid'] = False
            validation['issues'].append('Empty miRNA name')
        
        mirna_search = self.extract_mirna_for_search(mirna)
        if not mirna_search:
            validation['is_valid'] = False
            validation['issues'].append(f'Cannot extract search term from miRNA: {mirna}')
        
        return validation
    
    def get_article_statistics(self, gene_name: str, mirna_search: str) -> Dict[str, any]:
        """获取下载的文章统计信息"""
        output_file = os.path.join(self.output_dir, f"{gene_name}_{mirna_search}.txt")
        
        stats = {
            'file_exists': os.path.exists(output_file),
            'file_size': 0,
            'estimated_articles': 0
        }
        
        if stats['file_exists']:
            stats['file_size'] = os.path.getsize(output_file)
            
            # 估算文章数量（通过计算"Abstract"出现次数）
            try:
                with open(output_file, 'r', encoding=self.config.DEFAULT_ENCODING) as f:
                    content = f.read()
                    stats['estimated_articles'] = content.count('\nAbstract\n')
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"Error reading file for statistics: {e}")
        
        return stats

class LiteratureMiner:
    """文献挖掘主控制器"""
    
    def __init__(self, output_dir: str, logger: Optional[Logger] = None):
        self.config = PipelineConfig()
        self.output_dir = output_dir
        self.logger = logger
        self.pubmed_miner = PubMedMiner(output_dir, logger)
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
    
    def mine_literature(self, mti_results: List[Dict[str, any]], 
                       max_articles: int = None) -> List[Dict[str, any]]:
        """Step 2: 文献爬取，使用原始全称miRNA"""
        if max_articles is None:
            max_articles = self.config.DEFAULT_MAX_ARTICLES
        
        if self.logger:
            self.logger.info("="*50)
            self.logger.info("Step 2: Literature Mining")
            self.logger.info("="*50)
            self.logger.info(f"Mining literature for {len(mti_results)} MTIs...")
            self.logger.info(f"Max articles per MTI: {max_articles}")
        
        mining_results = []
        successful_queries = 0
        total_articles = 0
        
        for i, mti in enumerate(mti_results, 1):
            gene = mti['Gene']
            mirna_full = mti['miRNA']  # 使用原始全称
            
            # 验证查询参数
            validation = self.pubmed_miner.validate_query_parameters(gene, mirna_full)
            if not validation['is_valid']:
                if self.logger:
                    self.logger.warning(f"  Skipping invalid query for {gene}-{mirna_full}: {validation['issues']}")
                
                mining_results.append({
                    'Gene': gene,
                    'miRNA': mirna_full,
                    'miRNA_Search_Term': None,
                    'Direction': mti.get('Direction', 'Unknown'),
                    'Databases': mti.get('Databases', 0),
                    'Database_Sources': mti.get('Database_Sources', ''),
                    'miRTarBase_Validation': mti.get('miRTarBase_Validation', False),
                    'Validation_Strength': mti.get('Validation_Strength', 'Unknown'),
                    'Articles_Found': 0,
                    'File': None,
                    'Error': '; '.join(validation['issues'])
                })
                continue
            
            mirna_for_search = self.pubmed_miner.extract_mirna_for_search(mirna_full)
            
            if self.logger:
                self.logger.info(f"Processing MTI {i}/{len(mti_results)}: {gene} - {mirna_full}")
            
            # 查询PubMed
            article_count = self.pubmed_miner.query_pubmed(gene, mirna_full, max_articles)
            
            if article_count > 0:
                successful_queries += 1
                total_articles += article_count
            
            mining_results.append({
                'Gene': gene,
                'miRNA': mirna_full,  # 保存原始全称
                'miRNA_Search_Term': mirna_for_search,  # 保存搜索用的简化名称
                'Direction': mti.get('Direction', 'Unknown'),
                'Databases': mti.get('Databases', 0),
                'Database_Sources': mti.get('Database_Sources', ''),
                'miRTarBase_Validation': mti.get('miRTarBase_Validation', False),
                'Validation_Strength': mti.get('Validation_Strength', 'Unknown'),
                'Articles_Found': article_count,
                'File': f"{gene}_{mirna_for_search}.txt" if article_count > 0 else None
            })
            
            # 添加延迟避免过快请求
            time.sleep(self.config.PUBMED_DELAY)
        
        # 保存挖掘结果摘要
        self._save_mining_summary(mining_results)
        
        if self.logger:
            self.logger.info(f"Literature mining completed:")
            self.logger.info(f"  - Successful queries: {successful_queries}/{len(mti_results)}")
            self.logger.info(f"  - Total articles downloaded: {total_articles}")
            self.logger.info(f"  - Average articles per successful query: {total_articles/max(successful_queries, 1):.1f}")
        
        return mining_results
    
    def _save_mining_summary(self, mining_results: List[Dict[str, any]]):
        """保存挖掘结果摘要"""
        mining_summary = pd.DataFrame(mining_results)
        summary_file = os.path.join(self.output_dir, "mining_summary.csv")
        mining_summary.to_csv(summary_file, index=False, encoding=self.config.DEFAULT_ENCODING)
        
        if self.logger:
            self.logger.info(f"Mining summary saved to: {summary_file}")
    
    def get_mining_statistics(self) -> Dict[str, any]:
        """获取文献挖掘统计信息"""
        txt_files = [f for f in os.listdir(self.output_dir) if f.endswith('.txt')]
        
        stats = {
            'total_files': len(txt_files),
            'total_file_size': 0,
            'files_info': []
        }
        
        for filename in txt_files:
            filepath = os.path.join(self.output_dir, filename)
            try:
                file_size = os.path.getsize(filepath)
                stats['total_file_size'] += file_size
                
                # 估算文章数量
                with open(filepath, 'r', encoding=self.config.DEFAULT_ENCODING) as f:
                    content = f.read()
                    article_count = content.count('\nAbstract\n')
                
                stats['files_info'].append({
                    'filename': filename,
                    'size_bytes': file_size,
                    'estimated_articles': article_count
                })
                
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"Error processing file {filename}: {e}")
        
        return stats
    
    def validate_downloaded_files(self) -> Dict[str, any]:
        """验证下载的文件"""
        validation_results = {
            'valid_files': [],
            'empty_files': [],
            'missing_files': [],
            'corrupted_files': []
        }
        
        summary_file = os.path.join(self.output_dir, "mining_summary.csv")
        if not os.path.exists(summary_file):
            return validation_results
        
        try:
            summary_df = pd.read_csv(summary_file)
            
            for _, row in summary_df.iterrows():
                filename = row.get('File')
                if not filename:
                    continue
                
                filepath = os.path.join(self.output_dir, filename)
                
                if not os.path.exists(filepath):
                    validation_results['missing_files'].append(filename)
                    continue
                
                try:
                    file_size = os.path.getsize(filepath)
                    if file_size == 0:
                        validation_results['empty_files'].append(filename)
                        continue
                    
                    # 尝试读取文件内容
                    with open(filepath, 'r', encoding=self.config.DEFAULT_ENCODING) as f:
                        content = f.read(100)  # 只读取前100个字符进行验证
                        if content.strip():
                            validation_results['valid_files'].append(filename)
                        else:
                            validation_results['empty_files'].append(filename)
                
                except Exception as e:
                    validation_results['corrupted_files'].append((filename, str(e)))
                    if self.logger:
                        self.logger.warning(f"File corruption detected in {filename}: {e}")
        
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error validating files: {e}")
        
        return validation_results
    
    def cleanup_empty_files(self):
        """清理空文件"""
        validation = self.validate_downloaded_files()
        empty_files = validation['empty_files']
        
        if not empty_files:
            if self.logger:
                self.logger.info("No empty files found to clean up")
            return
        
        cleaned_count = 0
        for filename in empty_files:
            filepath = os.path.join(self.output_dir, filename)
            try:
                os.remove(filepath)
                cleaned_count += 1
                if self.logger:
                    self.logger.info(f"Removed empty file: {filename}")
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"Could not remove {filename}: {e}")
        
        if self.logger:
            self.logger.info(f"Cleaned up {cleaned_count} empty files")
    
    def retry_failed_queries(self, mining_results: List[Dict[str, any]], 
                           max_articles: int = None) -> List[Dict[str, any]]:
        """重试失败的查询"""
        if max_articles is None:
            max_articles = self.config.DEFAULT_MAX_ARTICLES
        
        failed_results = [r for r in mining_results if r['Articles_Found'] == 0]
        
        if not failed_results:
            if self.logger:
                self.logger.info("No failed queries to retry")
            return mining_results
        
        if self.logger:
            self.logger.info(f"Retrying {len(failed_results)} failed queries...")
        
        updated_results = []
        
        for result in mining_results:
            if result['Articles_Found'] > 0:
                # 保留成功的结果
                updated_results.append(result)
            else:
                # 重试失败的查询
                gene = result['Gene']
                mirna_full = result['miRNA']
                
                if self.logger:
                    self.logger.info(f"Retrying: {gene} - {mirna_full}")
                
                article_count = self.pubmed_miner.query_pubmed(gene, mirna_full, max_articles)
                
                # 更新结果
                result_copy = result.copy()
                result_copy['Articles_Found'] = article_count
                if article_count > 0:
                    mirna_search = self.pubmed_miner.extract_mirna_for_search(mirna_full)
                    result_copy['File'] = f"{gene}_{mirna_search}.txt"
                    result_copy.pop('Error', None)  # 移除错误信息
                
                updated_results.append(result_copy)
                
                time.sleep(self.config.PUBMED_DELAY)
        
        # 保存更新的摘要
        self._save_mining_summary(updated_results)
        
        successful_retries = sum(1 for r in updated_results if r['Articles_Found'] > 0 and r in failed_results)
        if self.logger:
            self.logger.info(f"Retry completed. Successful retries: {successful_retries}")
        
        return updated_results
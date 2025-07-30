#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
手动提取富集分析结果
专门用于从运行成功的富集分析中提取和保存CSV文件
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import pickle
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_enrichment_from_analyzer(output_dir='results2'):
    """从分析器对象中提取富集分析结果"""
    
    # 首先检查是否有分析器的pickle文件或类似缓存
    cache_files = [
        'analyzer_cache.pkl',
        'enrichment_cache.pkl', 
        '.enrichment_results.pkl'
    ]
    
    for cache_file in cache_files:
        cache_path = Path(output_dir) / cache_file
        if cache_path.exists():
            try:
                with open(cache_path, 'rb') as f:
                    data = pickle.load(f)
                logger.info(f"找到缓存文件: {cache_file}")
                return data
            except:
                logger.info(f"无法读取缓存文件: {cache_file}")
    
    return None

def manual_enrichment_analysis(gene_lists, output_dir='results2'):
    """手动重新运行富集分析并保存结果"""
    
    try:
        import gseapy as gp
    except ImportError:
        logger.error("gseapy未安装")
        return False
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # 从差异表达结果文件中提取基因列表
    if not gene_lists:
        logger.info("从差异表达结果文件中提取基因列表...")
        gene_lists = {}
        
        de_files = list(output_path.glob('*_gene_differential_expression_results.csv'))
        logger.info(f"找到 {len(de_files)} 个差异表达结果文件")
        
        for de_file in de_files:
            try:
                df = pd.read_csv(de_file, index_col=0)
                
                # 提取显著差异基因
                if 'padj' in df.columns and 'log2FoldChange' in df.columns:
                    significant = df[
                        (df['padj'] < 0.05) & 
                        (abs(df['log2FoldChange']) > 1.0)
                    ]
                    
                    if len(significant) > 0:
                        comparison_name = de_file.stem.replace('_differential_expression_results', '')
                        gene_lists[comparison_name] = significant.index.tolist()
                        logger.info(f"提取 {comparison_name}: {len(significant)} 个显著基因")
                    
            except Exception as e:
                logger.error(f"处理文件 {de_file.name} 失败: {e}")
    
    if not gene_lists:
        logger.error("没有找到可用的基因列表")
        return False
    
    # 转换基因ID
    def convert_genes(ensembl_ids):
        try:
            import mygene
            mg = mygene.MyGeneInfo()
            
            clean_ids = [gene_id.split('.')[0] for gene_id in ensembl_ids[:1000]]  # 限制数量避免超时
            
            results = mg.querymany(clean_ids, 
                                 scopes='ensembl.gene', 
                                 fields='symbol', 
                                 species='human',
                                 returnall=True)
            
            converted = {}
            for result in results['out']:
                if 'symbol' in result and 'query' in result:
                    converted[result['query']] = result['symbol']
            
            symbols = [converted.get(gene_id, gene_id) for gene_id in clean_ids]
            return [s for s in symbols if s and isinstance(s, str)]
            
        except Exception as e:
            logger.error(f"基因转换失败: {e}")
            return [gene_id.split('.')[0] for gene_id in ensembl_ids[:1000]]
    
    # 对每个基因列表进行富集分析
    success_count = 0
    
    for comparison, gene_list in gene_lists.items():
        logger.info(f"处理 {comparison}: {len(gene_list)} 个基因")
        
        if len(gene_list) < 5:
            logger.warning(f"基因数量太少，跳过: {comparison}")
            continue
        
        # 转换基因ID
        gene_symbols = convert_genes(gene_list)
        if len(gene_symbols) < 3:
            logger.warning(f"转换后基因数量太少，跳过: {comparison}")
            continue
        
        logger.info(f"转换后基因数量: {len(gene_symbols)}")
        
        # 进行富集分析
        try:
            # 获取可用的基因集
            available_libs = gp.get_library_name()
            go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
            kegg_sets = [lib for lib in available_libs if 'KEGG' in lib and 'Human' in lib]
            
            go_sets.sort(reverse=True)
            kegg_sets.sort(reverse=True)
            
            gene_sets = []
            if go_sets:
                gene_sets.append(go_sets[0])
            if kegg_sets:
                gene_sets.append(kegg_sets[0])
            
            if not gene_sets:
                logger.error("找不到可用的基因集")
                continue
            
            logger.info(f"使用基因集: {gene_sets}")
            
            # 运行enrichr
            enr = gp.enrichr(gene_list=gene_symbols,
                           gene_sets=gene_sets,
                           organism='Human',
                           outdir=None,
                           cutoff=0.05)
            
            # 保存结果
            if enr and hasattr(enr, 'results'):
                logger.info(f"富集分析成功，结果类型: {type(enr.results)}")
                
                if isinstance(enr.results, dict):
                    # 字典格式：每个基因集一个DataFrame
                    for gene_set_name, result_df in enr.results.items():
                        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
                            safe_gene_set = gene_set_name.replace(' ', '_').replace('/', '_')
                            filename = f"{comparison}_{safe_gene_set}_enrichment.csv"
                            filepath = output_path / filename
                            
                            result_df.to_csv(filepath, index=False)
                            logger.info(f"✅ 保存: {filename} ({len(result_df)} 条目)")
                            success_count += 1
                            
                elif isinstance(enr.results, pd.DataFrame) and not enr.results.empty:
                    # 单一DataFrame格式
                    if 'Gene_set' in enr.results.columns:
                        # 按基因集分组保存
                        for gene_set_name in enr.results['Gene_set'].unique():
                            subset = enr.results[enr.results['Gene_set'] == gene_set_name]
                            if not subset.empty:
                                safe_gene_set = str(gene_set_name).replace(' ', '_').replace('/', '_')
                                filename = f"{comparison}_{safe_gene_set}_enrichment.csv"
                                filepath = output_path / filename
                                
                                subset.to_csv(filepath, index=False)
                                logger.info(f"✅ 保存: {filename} ({len(subset)} 条目)")
                                success_count += 1
                    else:
                        # 保存为单一文件
                        filename = f"{comparison}_enrichment_results.csv"
                        filepath = output_path / filename
                        
                        enr.results.to_csv(filepath, index=False)
                        logger.info(f"✅ 保存: {filename} ({len(enr.results)} 条目)")
                        success_count += 1
                        
        except Exception as e:
            logger.error(f"富集分析失败 {comparison}: {e}")
            continue
    
    logger.info(f"手动富集分析完成，成功保存 {success_count} 个文件")
    return success_count > 0

def verify_saved_files(output_dir='results2'):
    """验证保存的文件"""
    output_path = Path(output_dir)
    
    enrichment_files = list(output_path.glob('*enrichment*.csv'))
    
    logger.info(f"在 {output_dir} 中找到 {len(enrichment_files)} 个富集分析文件:")
    
    for file in enrichment_files:
        try:
            df = pd.read_csv(file)
            logger.info(f"  ✅ {file.name}: {len(df)} 行 x {len(df.columns)} 列")
            
            # 显示文件内容预览
            if len(df) > 0:
                if 'Term' in df.columns:
                    logger.info(f"      示例通路: {df['Term'].iloc[0]}")
                elif len(df.columns) > 0:
                    logger.info(f"      第一列: {df.columns[0]}")
                    logger.info(f"      示例值: {df.iloc[0, 0]}")
                    
        except Exception as e:
            logger.error(f"  ❌ 读取失败 {file.name}: {e}")
    
    return len(enrichment_files)

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='手动提取富集分析结果')
    parser.add_argument('--output-dir', default='results2', help='输出目录')
    parser.add_argument('--force-rerun', action='store_true', help='强制重新运行富集分析')
    
    args = parser.parse_args()
    
    logger.info("🔍 手动提取富集分析结果")
    logger.info("="*50)
    
    # 1. 检查现有文件
    existing_files = verify_saved_files(args.output_dir)
    
    if existing_files > 0 and not args.force_rerun:
        logger.info("✅ 已存在富集分析文件")
        choice = input("是否重新运行富集分析？(y/N): ").strip().lower()
        if choice not in ['y', 'yes']:
            logger.info("跳过重新运行")
            return
    
    # 2. 手动运行富集分析
    logger.info("开始手动富集分析...")
    success = manual_enrichment_analysis(None, args.output_dir)
    
    if success:
        logger.info("✅ 手动富集分析完成")
        verify_saved_files(args.output_dir)
    else:
        logger.error("❌ 手动富集分析失败")
    
    logger.info("完成!")

if __name__ == "__main__":
    main()
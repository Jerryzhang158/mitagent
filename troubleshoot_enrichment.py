#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
富集分析故障排除脚本
帮助诊断GO/KEGG富集分析问题
"""

import pandas as pd
import numpy as np
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_gseapy_installation():
    """检查gseapy安装和基本功能"""
    try:
        import gseapy as gp
        logger.info("✅ gseapy已安装")
        
        # 检查版本
        logger.info(f"gseapy版本: {gp.__version__}")
        
        # 测试网络连接
        try:
            available_libs = gp.get_library_name()
            logger.info(f"✅ 网络连接正常，可用基因集数量: {len(available_libs)}")
            return True, available_libs
        except Exception as e:
            logger.error(f"❌ 网络连接问题: {e}")
            return False, []
            
    except ImportError:
        logger.error("❌ gseapy未安装，请运行: pip install gseapy")
        return False, []

def find_best_gene_sets():
    """查找最合适的GO和KEGG基因集"""
    try:
        import gseapy as gp
        available_libs = gp.get_library_name()
        
        # 查找GO基因集
        go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
        go_sets.sort(reverse=True)  # 按年份排序，最新的在前
        
        # 查找KEGG基因集
        kegg_sets = [lib for lib in available_libs if 'KEGG' in lib and 'Human' in lib]
        kegg_sets.sort(reverse=True)
        
        logger.info("推荐的基因集:")
        if go_sets:
            logger.info(f"  GO生物过程: {go_sets[0]}")
        if kegg_sets:
            logger.info(f"  KEGG通路: {kegg_sets[0]}")
            
        return go_sets[0] if go_sets else None, kegg_sets[0] if kegg_sets else None
        
    except Exception as e:
        logger.error(f"查找基因集失败: {e}")
        return None, None

def test_sample_enrichment():
    """使用样本基因测试富集分析"""
    try:
        import gseapy as gp
        
        # 使用一些已知的基因符号进行测试
        test_genes = ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51', 'BARD1']
        
        go_set, kegg_set = find_best_gene_sets()
        if not go_set or not kegg_set:
            logger.error("找不到合适的基因集")
            return False
            
        gene_sets = [go_set, kegg_set]
        
        logger.info(f"测试富集分析，使用基因: {test_genes}")
        logger.info(f"使用基因集: {gene_sets}")
        
        enr = gp.enrichr(gene_list=test_genes,
                        gene_sets=gene_sets,
                        organism='Human',
                        outdir=None,
                        cutoff=0.05)
        
        if enr and hasattr(enr, 'results'):
            logger.info("✅ 富集分析测试成功！")
            for gene_set, result_df in enr.results.items():
                logger.info(f"  {gene_set}: {len(result_df)} 个富集条目")
            return True
        else:
            logger.error("❌ 富集分析测试失败：没有返回结果")
            return False
            
    except Exception as e:
        logger.error(f"❌ 富集分析测试失败: {e}")
        return False

def convert_ensembl_to_symbols(ensembl_ids):
    """将ENSEMBL ID转换为基因符号"""
    try:
        # 尝试使用mygene
        import mygene
        mg = mygene.MyGeneInfo()
        
        # 清理ENSEMBL ID
        clean_ids = [gene_id.split('.')[0] for gene_id in ensembl_ids[:10]]  # 只测试前10个
        
        logger.info(f"测试ENSEMBL ID转换，样本: {clean_ids}")
        
        results = mg.querymany(clean_ids, 
                             scopes='ensembl.gene', 
                             fields='symbol', 
                             species='human')
        
        converted = {}
        for result in results:
            if 'symbol' in result:
                converted[result['query']] = result['symbol']
        
        logger.info(f"转换成功: {len(converted)}/{len(clean_ids)} 个基因")
        logger.info(f"转换示例: {dict(list(converted.items())[:3])}")
        
        return converted
        
    except ImportError:
        logger.error("mygene未安装，无法进行基因ID转换")
        logger.info("安装命令: pip install mygene")
        return {}
    except Exception as e:
        logger.error(f"基因ID转换失败: {e}")
        return {}

def check_enrichment_results_files(output_dir='result_final'):
    """检查富集分析结果文件"""
    from pathlib import Path
    
    output_path = Path(output_dir)
    if not output_path.exists():
        logger.error(f"结果目录不存在: {output_dir}")
        return
    
    # 查找富集分析文件
    enrichment_files = list(output_path.glob('*enrichment.csv'))
    
    if not enrichment_files:
        logger.error("❌ 没有找到富集分析结果文件")
        logger.info("可能的原因:")
        logger.info("1. 富集分析没有运行成功")
        logger.info("2. 基因ID格式不正确")
        logger.info("3. 网络连接问题")
        return
    
    logger.info(f"✅ 找到 {len(enrichment_files)} 个富集分析文件:")
    for file in enrichment_files:
        try:
            df = pd.read_csv(file)
            logger.info(f"  {file.name}: {len(df)} 行")
            if len(df) > 0:
                logger.info(f"    前3个条目: {df['Term'].head(3).tolist()}")
        except Exception as e:
            logger.error(f"  读取 {file.name} 失败: {e}")

def main():
    """主函数 - 完整的故障排除流程"""
    logger.info("🔍 开始富集分析故障排除...")
    logger.info("="*50)
    
    # 1. 检查gseapy安装
    logger.info("步骤1: 检查gseapy安装和网络连接")
    gseapy_ok, available_libs = check_gseapy_installation()
    if not gseapy_ok:
        logger.error("请先解决gseapy安装问题")
        return
    
    # 2. 查找合适的基因集
    logger.info("\n步骤2: 查找合适的基因集")
    go_set, kegg_set = find_best_gene_sets()
    
    # 3. 测试富集分析
    logger.info("\n步骤3: 测试富集分析功能")
    test_ok = test_sample_enrichment()
    
    # 4. 检查现有结果文件
    logger.info("\n步骤4: 检查现有结果文件")
    check_enrichment_results_files()
    
    # 5. 测试基因ID转换
    logger.info("\n步骤5: 测试基因ID转换")
    # 使用一些示例ENSEMBL ID
    sample_ensembl_ids = ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000139618']
    converted = convert_ensembl_to_symbols(sample_ensembl_ids)
    
    # 6. 总结和建议
    logger.info("\n" + "="*50)
    logger.info("📋 诊断总结:")
    
    if gseapy_ok:
        logger.info("✅ gseapy安装正常")
    else:
        logger.error("❌ gseapy有问题")
    
    if test_ok:
        logger.info("✅ 富集分析功能正常")
    else:
        logger.error("❌ 富集分析功能异常")
    
    if converted:
        logger.info("✅ 基因ID转换可用")
    else:
        logger.error("❌ 基因ID转换有问题")
    
    logger.info("\n💡 建议的解决方案:")
    
    if go_set and kegg_set:
        logger.info(f"1. 使用以下基因集进行分析:")
        logger.info(f"   GO: {go_set}")
        logger.info(f"   KEGG: {kegg_set}")
    
    if not converted:
        logger.info("2. 安装mygene进行基因ID转换:")
        logger.info("   pip install mygene")
    
    if not test_ok:
        logger.info("3. 检查网络连接和防火墙设置")
        logger.info("4. 更新gseapy到最新版本:")
        logger.info("   pip install --upgrade gseapy")
    
    logger.info("\n🚀 修复后重新运行富集分析:")
    logger.info("python rnaseq_analyzer.py --enrichment-only --output-dir result_final")

if __name__ == "__main__":
    main()
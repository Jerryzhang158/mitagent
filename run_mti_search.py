#!/usr/bin/env python3
"""
运行MTI搜索 - 为9个重合基因查找上游miRNA
使用用户提供的三个数据库文件
"""

import os
import sys
from mti_selection import MTISelector
from utils import Logger

def main():
    """主函数 - 执行MTI搜索"""
    
    print("🧬 MTI搜索：查找9个重合基因的上游miRNA")
    print("="*60)
    
    # 1. 创建基因列表文件
    overlap_genes = [
        'SLC7A11',  # 溶质载体家族7成员11
        'TXNRD1',   # 硫氧还蛋白还原酶1  
        'FTH1',     # 铁蛋白重链1
        'ACSL1',    # 长链脂肪酸辅酶A合成酶1
        'SAT1',     # 精胺/精胺乙酰转移酶1
        'AKR1C3',   # 醛酮还原酶家族1成员C3
        'HMOX1',    # 血红素加氧酶1
        'AKR1C2',   # 醛酮还原酶家族1成员C2
        'AKR1C1'    # 醛酮还原酶家族1成员C1
    ]
    
    # 创建基因列表文件
    gene_list_file = "overlap_genes.txt"
    with open(gene_list_file, 'w') as f:
        for gene in overlap_genes:
            f.write(f"{gene}\n")
    
    print(f"✅ 创建基因列表文件: {gene_list_file}")
    print(f"📊 包含 {len(overlap_genes)} 个重合基因")
    
    # 2. 定义数据库文件路径
    database_files = {
        'targetscan_file': 'Predicted_Targets_Context_Scores.default_predictions.txt',
        'mirdb_file': 'miRDB_v6.0_prediction_result_fixed.txt',
        'mirwalk_file': 'hsa_miRWalk_3UTR.txt',
        'mirtarbase_file': None  # 暂时不使用miRTarBase
    }
    
    # 3. 检查文件是否存在
    print("\n🔍 检查数据库文件...")
    missing_files = []
    for db_name, file_path in database_files.items():
        if file_path and not os.path.exists(file_path):
            missing_files.append(f"{db_name}: {file_path}")
            print(f"❌ 文件不存在: {file_path}")
        elif file_path:
            file_size = os.path.getsize(file_path) / (1024*1024)  # MB
            print(f"✅ {db_name}: {file_path} ({file_size:.1f} MB)")
    
    if missing_files:
        print("\n❌ 以下数据库文件缺失:")
        for missing in missing_files:
            print(f"   {missing}")
        print("\n请确保所有数据库文件都在当前目录中")
        return
    
    # 4. 初始化MTI选择器
    print("\n🚀 初始化MTI选择器...")
    try:
        logger = Logger()
        mti_selector = MTISelector(
            min_databases=2,  # 至少在2个数据库中找到
            logger=logger
        )
        print("✅ MTI选择器初始化成功")
    except Exception as e:
        print(f"❌ MTI选择器初始化失败: {e}")
        return
    
    # 5. 加载数据库
    print("\n📚 加载miRNA数据库...")
    try:
        databases = mti_selector.load_databases(
            database_files['targetscan_file'],
            database_files['mirdb_file'],
            database_files['mirwalk_file'],
            database_files['mirtarbase_file']
        )
        print("✅ 所有数据库加载成功")
    except Exception as e:
        print(f"❌ 数据库加载失败: {e}")
        return
    
    # 6. 执行基因到miRNA搜索
    print("\n🎯 开始基因→miRNA搜索...")
    print("搜索模式: 基因到miRNA")
    print("最小数据库要求: 2个数据库")
    print("这可能需要几分钟时间...")
    
    try:
        mti_results = mti_selector.gene_to_mirna_selection(
            gene_list_file, 
            databases
        )
        print(f"✅ 搜索完成，找到 {len(mti_results)} 个miRNA-基因关系")
    except Exception as e:
        print(f"❌ 搜索失败: {e}")
        return
    
    # 7. 保存结果
    output_file = "overlap_genes_mirna_results.xlsx"
    print(f"\n💾 保存结果到: {output_file}")
    
    try:
        mti_selector.save_results(mti_results, output_file, 'gene_to_mirna')
        print("✅ 结果保存成功")
    except Exception as e:
        print(f"❌ 结果保存失败: {e}")
        return
    
    # 8. 分析并显示结果统计
    print("\n" + "="*60)
    print("📊 MTI搜索结果分析")
    print("="*60)
    
    if mti_results:
        analyze_results(mti_results, overlap_genes)
    else:
        print("⚠️  未找到任何miRNA-基因关系")
        print("\n可能的原因:")
        print("- 基因名称在数据库中不匹配")
        print("- min_databases设置过高（当前为2）") 
        print("- 数据库文件格式不正确")
        print("\n建议:")
        print("- 尝试降低min_databases到1")
        print("- 检查数据库文件的列名和格式")
    
    print(f"\n🎉 分析完成！详细结果请查看: {output_file}")

def analyze_results(mti_results, overlap_genes):
    """分析MTI搜索结果"""
    
    # 按基因分组统计
    gene_mirna_count = {}
    gene_validated_count = {}
    all_mirnas = set()
    
    for mti in mti_results:
        gene = mti['Gene']
        mirna = mti['miRNA']
        
        if gene not in gene_mirna_count:
            gene_mirna_count[gene] = set()
            gene_validated_count[gene] = 0
        
        gene_mirna_count[gene].add(mirna)
        all_mirnas.add(mirna)
        
        if mti.get('miRTarBase_Validation', False):
            gene_validated_count[gene] += 1
    
    print(f"🎯 总共找到 {len(mti_results)} 个miRNA-基因关系")
    print(f"🧬 涉及 {len(all_mirnas)} 个不同的miRNA")
    print(f"📈 有miRNA调控的基因: {len(gene_mirna_count)}/{len(overlap_genes)}")
    
    # 各基因的miRNA数量
    print(f"\n📊 各基因的miRNA数量:")
    for gene in overlap_genes:
        mirna_set = gene_mirna_count.get(gene, set())
        mirna_count = len(mirna_set)
        validated_count = gene_validated_count.get(gene, 0)
        status = "✅" if mirna_count > 0 else "❌"
        print(f"  {status} {gene}: {mirna_count} miRNAs ({validated_count} 实验验证)")
        
        # 显示前5个miRNA
        if mirna_count > 0:
            top_mirnas = list(mirna_set)[:5]
            mirna_preview = ", ".join(top_mirnas)
            if mirna_count > 5:
                mirna_preview += f" ... (还有{mirna_count-5}个)"
            print(f"      └─ {mirna_preview}")
    
    # 数据库分布统计
    db_counts = {}
    for mti in mti_results:
        db_count = mti['Databases']
        db_counts[db_count] = db_counts.get(db_count, 0) + 1
    
    print(f"\n🔍 数据库支持分布:")
    for db_count in sorted(db_counts.keys(), reverse=True):
        print(f"  {db_count} 个数据库: {db_counts[db_count]} 个MTI关系")
    
    # 优先级分布
    priority_counts = {}
    for mti in mti_results:
        priority = mti.get('Priority', 'Unknown')
        priority_counts[priority] = priority_counts.get(priority, 0) + 1
    
    print(f"\n⭐ 优先级分布:")
    for priority in ['High', 'Medium', 'Low']:
        count = priority_counts.get(priority, 0)
        if count > 0:
            print(f"  {priority}: {count} 个MTI关系")
    
    # Top miRNA（调控基因数最多的）
    mirna_gene_count = {}
    for mti in mti_results:
        mirna = mti['miRNA']
        mirna_gene_count[mirna] = mirna_gene_count.get(mirna, 0) + 1
    
    top_mirnas = sorted(mirna_gene_count.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\n🏆 Top 10 miRNA（调控基因数最多）:")
    for i, (mirna, gene_count) in enumerate(top_mirnas, 1):
        print(f"  {i:2d}. {mirna}: 调控 {gene_count} 个基因")
    
    # 实验验证统计
    validated_count = sum(1 for mti in mti_results if mti.get('miRTarBase_Validation', False))
    if validated_count > 0:
        print(f"\n🧪 实验验证:")
        print(f"  已验证: {validated_count} 个MTI关系")
        print(f"  未验证: {len(mti_results) - validated_count} 个MTI关系")

if __name__ == "__main__":
    main()
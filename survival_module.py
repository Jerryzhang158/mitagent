#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
癌症生存分析模块
用于分析基因表达与患者生存期的关系

主要功能：
1. TCGA数据获取和处理
2. Kaplan-Meier生存曲线
3. Cox回归分析
4. 风险评分模型
5. ROC曲线分析

依赖安装：
pip install lifelines pandas numpy matplotlib seaborn scikit-learn requests
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import requests
import json
from typing import List, Dict, Optional, Tuple, Any
import warnings
warnings.filterwarnings('ignore')

# 生存分析库
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("警告: lifelines未安装，生存分析功能将受限")
    print("安装命令: pip install lifelines")

# ROC分析
try:
    from sklearn.metrics import roc_curve, auc, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("警告: scikit-learn未安装，ROC分析将不可用")

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class TCGADataFetcher:
    """TCGA数据获取类"""
    
    def __init__(self):
        """初始化TCGA数据获取器"""
        self.gdc_api = "https://api.gdc.cancer.gov"
        self.gepia_api = "http://gepia2.cancer-pku.cn/api"
        
        # TCGA癌症类型代码
        self.cancer_types = {
            'BRCA': 'Breast invasive carcinoma',
            'LUAD': 'Lung adenocarcinoma',
            'LUSC': 'Lung squamous cell carcinoma',
            'COAD': 'Colon adenocarcinoma',
            'READ': 'Rectum adenocarcinoma',
            'PRAD': 'Prostate adenocarcinoma',
            'THCA': 'Thyroid carcinoma',
            'LIHC': 'Liver hepatocellular carcinoma',
            'STAD': 'Stomach adenocarcinoma',
            'KIRC': 'Kidney renal clear cell carcinoma',
            'BLCA': 'Bladder urothelial carcinoma',
            'HNSC': 'Head and neck squamous cell carcinoma',
            'ESCA': 'Esophageal carcinoma',
            'UCEC': 'Uterine corpus endometrial carcinoma',
            'OV': 'Ovarian serous cystadenocarcinoma',
            'PAAD': 'Pancreatic adenocarcinoma',
            'SKCM': 'Skin cutaneous melanoma',
            'GBM': 'Glioblastoma multiforme',
            'LGG': 'Brain lower grade glioma'
        }
        
        logger.info("TCGA数据获取器初始化完成")
    
    def get_sample_clinical_data(self, cancer_type: str = 'BRCA') -> pd.DataFrame:
        """
        获取示例临床数据（使用预定义数据）
        实际应用中应该从TCGA下载真实数据
        
        Args:
            cancer_type: 癌症类型
            
        Returns:
            临床数据DataFrame
        """
        logger.info(f"生成 {cancer_type} 示例临床数据...")
        
        # 生成示例数据（实际应用中应从TCGA下载）
        np.random.seed(42)
        n_samples = 500
        
        # 生成生存数据
        clinical_data = pd.DataFrame({
            'sample_id': [f'TCGA-{i:04d}' for i in range(n_samples)],
            'days_to_death': np.random.exponential(1000, n_samples),
            'days_to_last_followup': np.random.exponential(1200, n_samples),
            'vital_status': np.random.choice(['Alive', 'Dead'], n_samples, p=[0.6, 0.4]),
            'age_at_diagnosis': np.random.normal(60, 15, n_samples).clip(20, 90),
            'gender': np.random.choice(['Male', 'Female'], n_samples),
            'tumor_stage': np.random.choice(['I', 'II', 'III', 'IV'], n_samples, p=[0.2, 0.3, 0.3, 0.2]),
            'cancer_type': cancer_type
        })
        
        # 处理生存时间
        clinical_data['survival_time'] = clinical_data.apply(
            lambda x: x['days_to_death'] if x['vital_status'] == 'Dead' else x['days_to_last_followup'],
            axis=1
        )
        clinical_data['survival_event'] = (clinical_data['vital_status'] == 'Dead').astype(int)
        
        # 转换为月
        clinical_data['survival_months'] = clinical_data['survival_time'] / 30
        
        logger.info(f"生成 {n_samples} 个样本的临床数据")
        return clinical_data
    
    def get_gene_expression_data(self, gene_list: List[str], 
                                cancer_type: str = 'BRCA') -> pd.DataFrame:
        """
        获取基因表达数据（模拟数据）
        
        Args:
            gene_list: 基因列表
            cancer_type: 癌症类型
            
        Returns:
            基因表达DataFrame
        """
        logger.info(f"生成 {len(gene_list)} 个基因的表达数据...")
        
        # 获取临床数据的样本ID
        clinical = self.get_sample_clinical_data(cancer_type)
        sample_ids = clinical['sample_id'].tolist()
        
        # 生成模拟表达数据
        np.random.seed(42)
        expression_data = {}
        
        for gene in gene_list[:50]:  # 限制基因数量
            # 生成与生存相关的表达值
            base_expression = np.random.normal(10, 3, len(sample_ids))
            
            # 添加与生存的相关性（模拟）
            survival_effect = clinical['survival_event'].values * np.random.normal(-1, 0.5, len(sample_ids))
            expression = base_expression + survival_effect
            
            expression_data[gene] = expression
        
        expr_df = pd.DataFrame(expression_data, index=sample_ids)
        
        logger.info(f"生成表达矩阵: {expr_df.shape}")
        return expr_df


class SurvivalAnalyzer:
    """生存分析类"""
    
    def __init__(self):
        """初始化生存分析器"""
        self.tcga_fetcher = TCGADataFetcher()
        self.km_fitter = KaplanMeierFitter() if HAS_LIFELINES else None
        self.cox_fitter = CoxPHFitter() if HAS_LIFELINES else None
        self.results = {}
        
        logger.info("生存分析器初始化完成")
    
    def prepare_survival_data(self, expression_df: pd.DataFrame, 
                             clinical_df: pd.DataFrame,
                             gene: str) -> pd.DataFrame:
        """
        准备生存分析数据
        
        Args:
            expression_df: 表达数据
            clinical_df: 临床数据
            gene: 目标基因
            
        Returns:
            合并后的数据
        """
        if gene not in expression_df.columns:
            logger.error(f"基因 {gene} 不在表达数据中")
            return pd.DataFrame()
        
        # 合并表达和临床数据
        expr_series = expression_df[gene]
        
        # 确保索引匹配
        common_samples = list(set(expr_series.index) & set(clinical_df['sample_id'].values))
        
        if len(common_samples) == 0:
            logger.error("没有匹配的样本")
            return pd.DataFrame()
        
        # 创建分析数据框
        survival_data = clinical_df[clinical_df['sample_id'].isin(common_samples)].copy()
        survival_data['gene_expression'] = survival_data['sample_id'].map(expr_series.to_dict())
        
        # 基于中位数分组
        median_expr = survival_data['gene_expression'].median()
        survival_data['expression_group'] = (survival_data['gene_expression'] > median_expr).astype(int)
        survival_data['expression_label'] = survival_data['expression_group'].map({0: 'Low', 1: 'High'})
        
        return survival_data
    
    def kaplan_meier_analysis(self, survival_data: pd.DataFrame, 
                            gene: str,
                            time_col: str = 'survival_months',
                            event_col: str = 'survival_event',
                            group_col: str = 'expression_label') -> Dict:
        """
        执行Kaplan-Meier生存分析
        
        Args:
            survival_data: 生存数据
            gene: 基因名
            time_col: 时间列名
            event_col: 事件列名
            group_col: 分组列名
            
        Returns:
            分析结果
        """
        if not HAS_LIFELINES:
            logger.error("lifelines未安装")
            return {}
        
        logger.info(f"执行 {gene} 的Kaplan-Meier分析...")
        
        results = {'gene': gene}
        
        try:
            # 分组数据
            groups = survival_data[group_col].unique()
            
            if len(groups) != 2:
                logger.error(f"需要恰好2个组，但找到 {len(groups)} 个")
                return results
            
            # 准备数据
            group1_data = survival_data[survival_data[group_col] == groups[0]]
            group2_data = survival_data[survival_data[group_col] == groups[1]]
            
            # Log-rank检验
            lr_result = logrank_test(
                group1_data[time_col],
                group2_data[time_col],
                group1_data[event_col],
                group2_data[event_col]
            )
            
            results['logrank_pvalue'] = lr_result.p_value
            results['logrank_statistic'] = lr_result.test_statistic
            
            # 计算中位生存期
            km1 = KaplanMeierFitter()
            km1.fit(group1_data[time_col], group1_data[event_col], label=groups[0])
            
            km2 = KaplanMeierFitter()
            km2.fit(group2_data[time_col], group2_data[event_col], label=groups[1])
            
            results['median_survival'] = {
                groups[0]: km1.median_survival_time_,
                groups[1]: km2.median_survival_time_
            }
            
            results['km_fitters'] = {groups[0]: km1, groups[1]: km2}
            
            logger.info(f"  Log-rank p-value: {lr_result.p_value:.4f}")
            logger.info(f"  中位生存期: {groups[0]}={km1.median_survival_time_:.1f}, "
                       f"{groups[1]}={km2.median_survival_time_:.1f}")
            
        except Exception as e:
            logger.error(f"KM分析失败: {e}")
        
        return results
    
    def cox_regression_analysis(self, survival_data: pd.DataFrame,
                              covariates: List[str],
                              time_col: str = 'survival_months',
                              event_col: str = 'survival_event') -> Dict:
        """
        执行Cox比例风险回归
        
        Args:
            survival_data: 生存数据
            covariates: 协变量列表
            time_col: 时间列名
            event_col: 事件列名
            
        Returns:
            Cox回归结果
        """
        if not HAS_LIFELINES:
            logger.error("lifelines未安装")
            return {}
        
        logger.info("执行Cox回归分析...")
        
        # 准备数据
        cox_data = survival_data[[time_col, event_col] + covariates].copy()
        cox_data = cox_data.dropna()
        
        # 数值化分类变量
        for col in covariates:
            if cox_data[col].dtype == 'object':
                cox_data[col] = pd.Categorical(cox_data[col]).codes
        
        try:
            # 拟合Cox模型
            cph = CoxPHFitter()
            cph.fit(cox_data, duration_col=time_col, event_col=event_col)
            
            # 提取结果
            results = {
                'summary': cph.summary,
                'hazard_ratios': np.exp(cph.params_),
                'confidence_intervals': np.exp(cph.confidence_intervals_),
                'p_values': cph.summary['p'].to_dict(),
                'concordance_index': cph.concordance_index_
            }
            
            logger.info(f"  C-index: {cph.concordance_index_:.3f}")
            
            # 显示显著的协变量
            sig_vars = cph.summary[cph.summary['p'] < 0.05]
            if not sig_vars.empty:
                logger.info("  显著协变量:")
                for var in sig_vars.index:
                    hr = np.exp(cph.params_[var])
                    p_val = cph.summary.loc[var, 'p']
                    logger.info(f"    {var}: HR={hr:.2f}, p={p_val:.4f}")
            
            return results
            
        except Exception as e:
            logger.error(f"Cox回归失败: {e}")
            return {}
    
    def plot_survival_curve(self, km_results: Dict, 
                           title: str = None,
                           output_file: str = None,
                           figsize: Tuple[int, int] = (10, 6)) -> None:
        """
        绘制生存曲线
        
        Args:
            km_results: KM分析结果
            title: 图标题
            output_file: 输出文件
            figsize: 图像大小
        """
        if not HAS_LIFELINES:
            logger.error("lifelines未安装，无法绘图")
            return
        
        if 'km_fitters' not in km_results:
            logger.error("没有KM拟合结果")
            return
        
        plt.figure(figsize=figsize)
        
        # 绘制每组的生存曲线
        colors = ['blue', 'red']
        for i, (group, kmf) in enumerate(km_results['km_fitters'].items()):
            kmf.plot_survival_function(color=colors[i], label=group, ci_show=True)
        
        # 添加p值
        p_value = km_results.get('logrank_pvalue', 1)
        plt.text(0.1, 0.1, f'Log-rank p = {p_value:.4f}',
                transform=plt.gca().transAxes,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.xlabel('Time (months)')
        plt.ylabel('Survival Probability')
        
        if title:
            plt.title(title)
        else:
            gene = km_results.get('gene', 'Gene')
            plt.title(f'Kaplan-Meier Survival Curve - {gene}')
        
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"生存曲线已保存: {output_file}")
        
        plt.show()
    
    def calculate_risk_score(self, expression_df: pd.DataFrame,
                           gene_weights: Dict[str, float]) -> pd.Series:
        """
        计算风险评分
        
        Args:
            expression_df: 表达数据
            gene_weights: 基因权重（通常是Cox回归系数）
            
        Returns:
            风险评分Series
        """
        logger.info("计算风险评分...")
        
        risk_scores = pd.Series(index=expression_df.index, dtype=float)
        
        for sample in expression_df.index:
            score = 0
            for gene, weight in gene_weights.items():
                if gene in expression_df.columns:
                    score += expression_df.loc[sample, gene] * weight
            risk_scores[sample] = score
        
        return risk_scores
    
    def plot_risk_score_analysis(self, risk_scores: pd.Series,
                                clinical_data: pd.DataFrame,
                                output_file: str = None) -> None:
        """
        绘制风险评分分析图
        
        Args:
            risk_scores: 风险评分
            clinical_data: 临床数据
            output_file: 输出文件
        """
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        
        # 准备数据
        merged_data = clinical_data.copy()
        merged_data['risk_score'] = merged_data['sample_id'].map(risk_scores.to_dict())
        merged_data = merged_data.dropna(subset=['risk_score'])
        merged_data = merged_data.sort_values('risk_score')
        merged_data['rank'] = range(len(merged_data))
        
        # 1. 风险评分分布
        ax1 = axes[0]
        colors = ['green' if x == 0 else 'red' for x in merged_data['survival_event']]
        ax1.scatter(merged_data['rank'], merged_data['risk_score'], c=colors, alpha=0.6, s=20)
        ax1.axhline(y=merged_data['risk_score'].median(), color='black', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Patients (sorted by risk score)')
        ax1.set_ylabel('Risk Score')
        ax1.set_title('Risk Score Distribution')
        
        # 2. 生存时间分布
        ax2 = axes[1]
        ax2.scatter(merged_data['rank'], merged_data['survival_months'], c=colors, alpha=0.6, s=20)
        ax2.set_xlabel('Patients (sorted by risk score)')
        ax2.set_ylabel('Survival Time (months)')
        ax2.set_title('Survival Time Distribution')
        
        # 3. 生存状态热图
        ax3 = axes[2]
        survival_status = merged_data['survival_event'].values.reshape(1, -1)
        im = ax3.imshow(survival_status, aspect='auto', cmap='RdYlGn_r', vmin=0, vmax=1)
        ax3.set_xlabel('Patients (sorted by risk score)')
        ax3.set_yticks([0])
        ax3.set_yticklabels(['Status'])
        ax3.set_title('Survival Status (Green=Alive, Red=Dead)')
        
        plt.colorbar(im, ax=ax3, orientation='horizontal', pad=0.1)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"风险评分图已保存: {output_file}")
        
        plt.show()
    
    def roc_analysis(self, risk_scores: pd.Series,
                    clinical_data: pd.DataFrame,
                    time_point: float = 36) -> Dict:
        """
        ROC曲线分析
        
        Args:
            risk_scores: 风险评分
            clinical_data: 临床数据
            time_point: 时间点（月）
            
        Returns:
            ROC分析结果
        """
        if not HAS_SKLEARN:
            logger.error("scikit-learn未安装")
            return {}
        
        logger.info(f"执行{time_point}个月的ROC分析...")
        
        # 准备数据
        merged_data = clinical_data.copy()
        merged_data['risk_score'] = merged_data['sample_id'].map(risk_scores.to_dict())
        merged_data = merged_data.dropna(subset=['risk_score'])
        
        # 定义时间点的生存状态
        merged_data['event_at_time'] = (
            (merged_data['survival_months'] <= time_point) & 
            (merged_data['survival_event'] == 1)
        ).astype(int)
        
        # 计算ROC
        fpr, tpr, thresholds = roc_curve(merged_data['event_at_time'], 
                                        merged_data['risk_score'])
        roc_auc = auc(fpr, tpr)
        
        # 找最佳阈值
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
        
        results = {
            'fpr': fpr,
            'tpr': tpr,
            'thresholds': thresholds,
            'auc': roc_auc,
            'optimal_threshold': optimal_threshold,
            'time_point': time_point
        }
        
        logger.info(f"  AUC: {roc_auc:.3f}")
        logger.info(f"  最佳阈值: {optimal_threshold:.3f}")
        
        return results
    
    def plot_roc_curve(self, roc_results: Dict,
                      output_file: str = None) -> None:
        """
        绘制ROC曲线
        
        Args:
            roc_results: ROC分析结果
            output_file: 输出文件
        """
        plt.figure(figsize=(8, 6))
        
        plt.plot(roc_results['fpr'], roc_results['tpr'], 
                color='darkorange', lw=2,
                label=f'ROC curve (AUC = {roc_results["auc"]:.3f})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        
        # 标记最佳点
        optimal_idx = np.argmax(roc_results['tpr'] - roc_results['fpr'])
        plt.scatter(roc_results['fpr'][optimal_idx], 
                   roc_results['tpr'][optimal_idx],
                   color='red', s=100, label='Optimal threshold')
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'ROC Curve - {roc_results["time_point"]} months')
        plt.legend(loc="lower right")
        plt.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"ROC曲线已保存: {output_file}")
        
        plt.show()
    
    def batch_survival_analysis(self, gene_list: List[str],
                               cancer_type: str = 'BRCA',
                               output_dir: str = 'survival_results') -> pd.DataFrame:
        """
        批量生存分析
        
        Args:
            gene_list: 基因列表
            cancer_type: 癌症类型
            output_dir: 输出目录
            
        Returns:
            结果汇总DataFrame
        """
        logger.info(f"批量分析 {len(gene_list)} 个基因...")
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 获取数据
        clinical_data = self.tcga_fetcher.get_sample_clinical_data(cancer_type)
        expression_data = self.tcga_fetcher.get_gene_expression_data(gene_list, cancer_type)
        
        results_list = []
        
        for i, gene in enumerate(gene_list[:20], 1):  # 限制分析数量
            logger.info(f"分析 {i}/{min(len(gene_list), 20)}: {gene}")
            
            if gene not in expression_data.columns:
                logger.warning(f"  跳过 {gene} (无表达数据)")
                continue
            
            # 准备数据
            survival_data = self.prepare_survival_data(expression_data, clinical_data, gene)
            
            if survival_data.empty:
                continue
            
            # KM分析
            km_results = self.kaplan_meier_analysis(survival_data, gene)
            
            if km_results:
                # 保存结果
                results_list.append({
                    'gene': gene,
                    'logrank_pvalue': km_results.get('logrank_pvalue', np.nan),
                    'median_survival_low': km_results.get('median_survival', {}).get('Low', np.nan),
                    'median_survival_high': km_results.get('median_survival', {}).get('High', np.nan),
                    'hazard_ratio': np.nan  # 需要Cox回归计算
                })
                
                # 绘制生存曲线
                if km_results.get('logrank_pvalue', 1) < 0.05:
                    self.plot_survival_curve(
                        km_results,
                        output_file=str(output_path / f'{gene}_survival_curve.png')
                    )
        
        # 汇总结果
        if results_list:
            results_df = pd.DataFrame(results_list)
            results_df = results_df.sort_values('logrank_pvalue')
            results_df.to_csv(output_path / 'survival_analysis_summary.csv', index=False)
            
            logger.info(f"分析完成，结果保存在: {output_path}")
            return results_df
        else:
            logger.warning("没有成功分析的基因")
            return pd.DataFrame()


def analyze_gene_survival(gene_list: List[str],
                         cancer_type: str = 'BRCA',
                         output_dir: str = 'survival_analysis') -> Dict:
    """
    便捷函数：分析基因列表的生存影响
    
    Args:
        gene_list: 基因列表
        cancer_type: 癌症类型
        output_dir: 输出目录
        
    Returns:
        分析结果
    """
    logger.info("="*60)
    logger.info("开始生存分析")
    logger.info("="*60)
    
    # 创建分析器
    analyzer = SurvivalAnalyzer()
    
    # 批量分析
    results_df = analyzer.batch_survival_analysis(
        gene_list,
        cancer_type,
        output_dir
    )
    
    if not results_df.empty:
        # 显示显著基因
        sig_genes = results_df[results_df['logrank_pvalue'] < 0.05]
        
        if not sig_genes.empty:
            logger.info(f"发现 {len(sig_genes)} 个预后相关基因:")
            for _, row in sig_genes.head(10).iterrows():
                logger.info(f"  {row['gene']}: p={row['logrank_pvalue']:.4f}")
        else:
            logger.info("没有发现显著的预后相关基因")
    
    logger.info("="*60)
    logger.info("生存分析完成")
    logger.info("="*60)
    
    return {
        'results_df': results_df,
        'analyzer': analyzer
    }


if __name__ == "__main__":
    # 测试代码
    import argparse
    
    parser = argparse.ArgumentParser(description='生存分析模块')
    parser.add_argument('--genes', nargs='+', required=True, help='基因列表')
    parser.add_argument('--cancer', default='BRCA', help='癌症类型')
    parser.add_argument('--output', default='survival_results', help='输出目录')
    
    args = parser.parse_args()
    
    # 运行分析
    results = analyze_gene_survival(
        args.genes,
        args.cancer,
        args.output
    )
    
    print(f"分析完成！结果保存在: {args.output}")
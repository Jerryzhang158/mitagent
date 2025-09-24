#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from typing import List, Dict, Optional, Tuple, Any
import warnings
import sys
import os

warnings.filterwarnings('ignore')

# 尝试导入自定义模块，如果失败则提供默认实现
try:
    from ppi_network_module import PPINetworkAnalyzer, analyze_de_genes_ppi
    PPI_AVAILABLE = True
except ImportError:
    logging.warning("PPI网络模块未找到，将跳过PPI分析")
    PPI_AVAILABLE = False

try:
    from survival_module import SurvivalAnalyzer, analyze_gene_survival
    SURVIVAL_AVAILABLE = True
except ImportError:
    logging.warning("生存分析模块未找到，将跳过生存分析")
    SURVIVAL_AVAILABLE = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置中文字体，增加更多备选项
def setup_chinese_fonts():
    """设置中文字体支持"""
    chinese_fonts = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS', 'Microsoft YaHei', 'WenQuanYi Micro Hei']
    
    for font in chinese_fonts:
        try:
            plt.rcParams['font.sans-serif'] = [font]
            plt.rcParams['axes.unicode_minus'] = False
            # 测试字体是否可用
            fig, ax = plt.subplots(figsize=(1, 1))
            ax.text(0.5, 0.5, '测试', fontfamily=font)
            plt.close(fig)
            logger.info(f"成功设置中文字体: {font}")
            break
        except:
            continue
    else:
        logger.warning("未找到合适的中文字体，可能影响图表显示")
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False

setup_chinese_fonts()


class BiomarkerIntegrator:
    """生物标志物整合分析类"""
    
    def __init__(self):
        """初始化整合分析器"""
        # 初始化logger
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # 初始化分析结果存储
        self.de_results = None
        self.ppi_results = None
        self.survival_results = None
        self.integrated_scores = None
        self.biomarkers = None
        
        self.logger.info("生物标志物整合分析器初始化完成")
    
    def _safe_standardize_scores(self, df: pd.DataFrame, score_columns: List[str]) -> pd.DataFrame:
        """
        安全地标准化得分，处理无穷大值和异常值
        
        Args:
            df: 输入DataFrame
            score_columns: 需要标准化的列名列表
            
        Returns:
            标准化后的DataFrame
        """
        try:
            from sklearn.preprocessing import StandardScaler
        except ImportError:
            self.logger.warning("sklearn未安装，使用手动标准化方法")
            return self._manual_standardize_scores(df, score_columns)
        
        self.logger.info("检查数据中的异常值...")
        
        # 检查和清理每列
        cleaned_df = df.copy()
        
        for col in score_columns:
            if col not in cleaned_df.columns:
                self.logger.warning(f"列 {col} 不存在，跳过")
                continue
                
            original_values = cleaned_df[col].copy()
            
            # 统计问题值
            inf_count = np.isinf(original_values).sum()
            nan_count = original_values.isna().sum()
            very_large_count = (np.abs(original_values) > 1e10).sum()
            
            if inf_count > 0 or nan_count > 0 or very_large_count > 0:
                self.logger.warning(f"列 {col}: inf={inf_count}, nan={nan_count}, very_large={very_large_count}")
                
                # 获取有限值的统计量
                finite_values = original_values[np.isfinite(original_values)]
                
                if len(finite_values) > 0:
                    # 使用99%和1%分位数作为替换边界
                    upper_bound = np.percentile(finite_values, 99)
                    lower_bound = np.percentile(finite_values, 1)
                    
                    # 替换异常值
                    cleaned_values = original_values.copy()
                    
                    # 正无穷 -> 上界的1.1倍
                    cleaned_values[np.isposinf(cleaned_values)] = upper_bound * 1.1
                    
                    # 负无穷 -> 下界的1.1倍  
                    cleaned_values[np.isneginf(cleaned_values)] = lower_bound * 1.1
                    
                    # 极大正值 -> 上界
                    cleaned_values[cleaned_values > 1e10] = upper_bound
                    
                    # 极大负值 -> 下界
                    cleaned_values[cleaned_values < -1e10] = lower_bound
                    
                    # NaN -> 中位数
                    median_val = finite_values.median() if len(finite_values) > 0 else 0
                    cleaned_values = cleaned_values.fillna(median_val)
                    
                    cleaned_df[col] = cleaned_values
                    
                    self.logger.info(f"  {col}: 替换边界 [{lower_bound:.3f}, {upper_bound:.3f}]")
                else:
                    # 如果没有有限值，全部设为0
                    cleaned_df[col] = 0
                    self.logger.warning(f"  {col}: 没有有限值，全部设为0")
        
        # 最终验证
        for col in score_columns:
            if col in cleaned_df.columns:
                final_check = np.isfinite(cleaned_df[col]).all()
                if not final_check:
                    self.logger.error(f"清理后 {col} 仍有异常值，强制设为0")
                    cleaned_df[col] = 0
        
        # 标准化
        try:
            scaler = StandardScaler()
            existing_cols = [col for col in score_columns if col in cleaned_df.columns]
            
            if not existing_cols:
                self.logger.warning("没有有效的评分列可供标准化")
                return cleaned_df
            
            score_data = cleaned_df[existing_cols]
            
            # 再次检查
            if not np.all(np.isfinite(score_data.values)):
                raise ValueError("数据仍包含非有限值")
            
            # 检查是否有足够的方差
            for col in existing_cols:
                if score_data[col].std() < 1e-10:
                    self.logger.warning(f"列 {col} 方差接近零，跳过标准化")
                    cleaned_df[f'{col}_normalized'] = 0
                    continue
            
            # 只对有方差的列进行标准化
            valid_cols = [col for col in existing_cols if score_data[col].std() >= 1e-10]
            
            if valid_cols:
                normalized_scores = scaler.fit_transform(score_data[valid_cols])
                
                # 更新DataFrame
                for i, col in enumerate(valid_cols):
                    cleaned_df[f'{col}_normalized'] = normalized_scores[:, i]
            
            self.logger.info("✅ 数据标准化成功")
            return cleaned_df
            
        except Exception as e:
            self.logger.error(f"StandardScaler失败: {e}")
            return self._manual_standardize_scores(cleaned_df, score_columns)
    
    def _manual_standardize_scores(self, df: pd.DataFrame, score_columns: List[str]) -> pd.DataFrame:
        """手动标准化方法（备用）"""
        self.logger.info("使用备用标准化方法...")
        cleaned_df = df.copy()
        
        for col in score_columns:
            if col in cleaned_df.columns:
                values = cleaned_df[col]
                mean_val = values.mean()
                std_val = values.std()
                
                if std_val > 1e-10:  # 避免除以接近0的数
                    normalized = (values - mean_val) / std_val
                else:
                    normalized = pd.Series(0, index=values.index)
                
                # 再次检查结果
                if not np.all(np.isfinite(normalized)):
                    normalized = pd.Series(0, index=values.index)
                    self.logger.warning(f"{col} 备用标准化也失败，设为0")
                
                cleaned_df[f'{col}_normalized'] = normalized
        
        self.logger.info("✅ 备用标准化完成")
        return cleaned_df
    
    def load_de_results(self, de_file: str, 
                       padj_threshold: float = 0.05,
                       log2fc_threshold: float = 1.0) -> pd.DataFrame:
        """
        加载差异表达结果
        
        Args:
            de_file: 差异表达结果文件
            padj_threshold: 校正p值阈值
            log2fc_threshold: log2 fold change阈值
            
        Returns:
            显著差异基因DataFrame
        """
        try:
            self.logger.info(f"加载差异表达结果: {de_file}")
            
            # 检查文件是否存在
            if not os.path.exists(de_file):
                raise FileNotFoundError(f"文件不存在: {de_file}")
            
            # 读取数据
            de_df = pd.read_csv(de_file, index_col=0)
            
            # 检查必需列
            required_cols = ['padj', 'log2FoldChange']
            missing_cols = [col for col in required_cols if col not in de_df.columns]
            if missing_cols:
                raise ValueError(f"缺少必需列: {missing_cols}")
            
            # 筛选显著基因
            significant = de_df[
                (de_df['padj'] < padj_threshold) & 
                (abs(de_df['log2FoldChange']) >= log2fc_threshold)
            ].copy()
            
            if significant.empty:
                self.logger.warning("未找到符合条件的显著差异基因，放宽条件重试")
                # 放宽条件重试
                significant = de_df[
                    (de_df['padj'] < padj_threshold * 2) & 
                    (abs(de_df['log2FoldChange']) >= log2fc_threshold * 0.5)
                ].copy()
            
            # 计算差异表达得分
            significant['de_score'] = abs(significant['log2FoldChange']) * (-np.log10(significant['padj'].clip(lower=1e-300)))
            
            self.de_results = significant
            self.logger.info(f"加载 {len(significant)} 个显著差异基因")
            
            return significant
            
        except Exception as e:
            self.logger.error(f"加载差异表达结果失败: {e}")
            return pd.DataFrame()
    
    def run_ppi_analysis(self, string_score: int = 400) -> Dict:
        """
        运行PPI网络分析
        
        Args:
            string_score: STRING置信度阈值
            
        Returns:
            PPI分析结果
        """
        if not PPI_AVAILABLE:
            self.logger.warning("PPI分析模块不可用，跳过PPI分析")
            self.ppi_results = {}
            return {}
        
        if self.de_results is None or self.de_results.empty:
            self.logger.error("需要先加载差异表达结果")
            return {}
        
        try:
            self.logger.info("运行PPI网络分析...")
            
            # 获取基因列表
            gene_list = self.de_results.index.tolist()
            
            # PPI分析
            ppi_analyzer = PPINetworkAnalyzer()
            
            # 获取PPI网络
            ppi_data = ppi_analyzer.fetch_ppi_network(gene_list, string_score)
            
            if ppi_data.empty:
                self.logger.warning("未获得PPI数据")
                self.ppi_results = {}
                return {}
            
            # 构建网络
            G = ppi_analyzer.build_network(ppi_data)
            
            # 计算中心性
            centrality_df = ppi_analyzer.calculate_centrality_measures(G)
            
            # 识别hub基因
            hub_genes = ppi_analyzer.identify_hub_genes(centrality_df, method='top_n', threshold=30)
            
            # 检测模块
            modules = ppi_analyzer.detect_modules(G)
            
            self.ppi_results = {
                'network': G,
                'centrality': centrality_df,
                'hub_genes': hub_genes,
                'modules': modules,
                'analyzer': ppi_analyzer
            }
            
            self.logger.info(f"PPI分析完成: {G.number_of_nodes()} 个节点, {len(hub_genes)} 个hub基因")
            
            return self.ppi_results
            
        except Exception as e:
            self.logger.error(f"PPI分析失败: {e}")
            self.ppi_results = {}
            return {}
    
    def run_survival_analysis(self, cancer_type: str = 'BRCA', 
                            top_n: int = 50) -> pd.DataFrame:
        """
        运行生存分析
        
        Args:
            cancer_type: 癌症类型
            top_n: 分析前n个基因
            
        Returns:
            生存分析结果
        """
        if not SURVIVAL_AVAILABLE:
            self.logger.warning("生存分析模块不可用，跳过生存分析")
            self.survival_results = pd.DataFrame()
            return pd.DataFrame()
        
        if self.de_results is None:
            self.logger.error("需要先加载差异表达结果")
            return pd.DataFrame()
        
        try:
            self.logger.info(f"运行生存分析 ({cancer_type})...")
            
            # 选择要分析的基因
            if self.ppi_results and 'hub_genes' in self.ppi_results:
                # 优先分析hub基因
                gene_list = self.ppi_results['hub_genes'][:top_n]
            else:
                # 否则选择差异最显著的基因
                gene_list = self.de_results.head(top_n).index.tolist()
            
            # 生存分析
            survival_analyzer = SurvivalAnalyzer()
            results_df = survival_analyzer.batch_survival_analysis(
                gene_list,
                cancer_type,
                output_dir='temp_survival'
            )
            
            self.survival_results = results_df
            
            if not results_df.empty:
                sig_survival = results_df[results_df['logrank_pvalue'] < 0.05]
                self.logger.info(f"生存分析完成: {len(sig_survival)}/{len(results_df)} 个基因与预后相关")
            
            return results_df
            
        except Exception as e:
            self.logger.error(f"生存分析失败: {e}")
            self.survival_results = pd.DataFrame()
            return pd.DataFrame()
    
    def calculate_integrated_scores(self) -> pd.DataFrame:
        """
        计算整合得分
        
        Returns:
            整合得分DataFrame
        """
        self.logger.info("计算整合得分...")
        
        # 收集所有基因
        all_genes = set()
        
        if self.de_results is not None:
            all_genes.update(self.de_results.index)
        
        if self.ppi_results and 'centrality' in self.ppi_results:
            all_genes.update(self.ppi_results['centrality'].index)
        
        if self.survival_results is not None and not self.survival_results.empty:
            all_genes.update(self.survival_results['gene'])
        
        if not all_genes:
            self.logger.error("没有基因数据可用于整合分析")
            return pd.DataFrame()
        
        # 创建整合数据框
        integrated_df = pd.DataFrame(index=list(all_genes))
        
        # 1. 差异表达得分
        if self.de_results is not None:
            integrated_df['de_score'] = self.de_results['de_score']
            integrated_df['log2FC'] = self.de_results['log2FoldChange']
            integrated_df['padj'] = self.de_results['padj']
        
        # 2. PPI网络得分
        if self.ppi_results and 'centrality' in self.ppi_results:
            centrality = self.ppi_results['centrality']
            integrated_df['hub_score'] = centrality['hub_score']
            integrated_df['degree'] = centrality['degree_count']
        
        # 3. 生存分析得分
        if self.survival_results is not None and not self.survival_results.empty:
            survival_dict = self.survival_results.set_index('gene')['logrank_pvalue'].to_dict()
            integrated_df['survival_pvalue'] = integrated_df.index.map(survival_dict)
            integrated_df['survival_score'] = -np.log10(integrated_df['survival_pvalue'].fillna(1).clip(lower=1e-10, upper=0.999))
        
        # 填充缺失值
        integrated_df = integrated_df.fillna(0)
        
        # 标准化得分
        score_columns = ['de_score', 'hub_score', 'survival_score']
        existing_score_cols = [col for col in score_columns if col in integrated_df.columns and not integrated_df[col].isna().all()]
        
        if existing_score_cols:
            # 使用安全的标准化方法
            integrated_df = self._safe_standardize_scores(integrated_df, existing_score_cols)
            
            # 获取标准化后的列名
            normalized_cols = [f'{col}_normalized' for col in existing_score_cols if f'{col}_normalized' in integrated_df.columns]
            
            # 计算综合得分
            weights = {
                'de_score_normalized': 0.4,      # 差异表达权重
                'hub_score_normalized': 0.3,     # 网络中心性权重
                'survival_score_normalized': 0.3  # 生存影响权重
            }
            
            integrated_df['integrated_score'] = 0
            total_weight = 0
            
            for col, weight in weights.items():
                if col in integrated_df.columns:
                    integrated_df['integrated_score'] += integrated_df[col] * weight
                    total_weight += weight
            
            # 归一化权重
            if total_weight > 0:
                integrated_df['integrated_score'] /= total_weight
            
            # 排序
            integrated_df = integrated_df.sort_values('integrated_score', ascending=False)
        else:
            self.logger.warning("没有有效的评分列可用于计算综合得分")
            integrated_df['integrated_score'] = 0
        
        self.integrated_scores = integrated_df
        self.logger.info(f"计算完成: {len(integrated_df)} 个基因的整合得分")
        
        return integrated_df
    
    def identify_key_biomarkers(self, top_n: int = 20) -> pd.DataFrame:
        """
        识别关键生物标志物
        
        Args:
            top_n: 选择前n个基因
            
        Returns:
            关键生物标志物DataFrame
        """
        if self.integrated_scores is None:
            self.logger.error("需要先计算整合得分")
            return pd.DataFrame()
        
        self.logger.info(f"识别前 {top_n} 个关键生物标志物...")
        
        # 选择top基因
        biomarkers = self.integrated_scores.head(top_n).copy()
        
        # 添加标签
        biomarkers['biomarker_rank'] = range(1, len(biomarkers) + 1)
        
        # 分类
        biomarkers['category'] = 'Candidate'
        
        # 基于多个标准分类
        for gene in biomarkers.index:
            criteria_met = 0
            
            # 标准1: 显著差异表达
            if 'padj' in biomarkers.columns and not pd.isna(biomarkers.loc[gene, 'padj']) and biomarkers.loc[gene, 'padj'] < 0.01:
                criteria_met += 1
            
            # 标准2: Hub基因
            if self.ppi_results and gene in self.ppi_results.get('hub_genes', []):
                criteria_met += 1
            
            # 标准3: 预后相关
            if 'survival_pvalue' in biomarkers.columns and not pd.isna(biomarkers.loc[gene, 'survival_pvalue']) and biomarkers.loc[gene, 'survival_pvalue'] < 0.05:
                criteria_met += 1
            
            # 分类
            if criteria_met >= 3:
                biomarkers.loc[gene, 'category'] = 'Key Biomarker'
            elif criteria_met >= 2:
                biomarkers.loc[gene, 'category'] = 'Strong Candidate'
        
        self.biomarkers = biomarkers
        
        # 统计
        category_counts = biomarkers['category'].value_counts()
        for cat, count in category_counts.items():
            self.logger.info(f"  {cat}: {count} 个基因")
        
        return biomarkers
    
    def plot_integration_heatmap(self, output_file: str = None) -> None:
        """
        绘制整合分析热图
        
        Args:
            output_file: 输出文件路径
        """
        if self.biomarkers is None:
            self.logger.error("需要先识别生物标志物")
            return
        
        try:
            self.logger.info("绘制整合分析热图...")
            
            # 准备数据
            score_cols = [col for col in ['de_score', 'hub_score', 'survival_score'] 
                         if col in self.biomarkers.columns and not self.biomarkers[col].isna().all()]
            
            if not score_cols:
                self.logger.error("没有可绘制的数据")
                return
            
            plot_data = self.biomarkers[score_cols].fillna(0)
            
            if plot_data.empty:
                self.logger.error("没有可绘制的数据")
                return
            
            # 创建热图
            plt.figure(figsize=(10, 12))
            
            # 添加分类颜色条
            category_colors = {
                'Key Biomarker': 'red',
                'Strong Candidate': 'orange',
                'Candidate': 'yellow'
            }
            
            row_colors = self.biomarkers['category'].map(category_colors)
            
            # 绘制聚类热图
            try:
                g = sns.clustermap(plot_data.T,
                                  col_colors=row_colors,
                                  cmap='RdBu_r',
                                  center=0,
                                  figsize=(12, 6),
                                  cbar_kws={'label': 'Normalized Score'},
                                  yticklabels=True,
                                  xticklabels=True)
                
                g.fig.suptitle('Integrated Biomarker Analysis', y=1.02)
                
                if output_file:
                    plt.savefig(output_file, dpi=300, bbox_inches='tight')
                    self.logger.info(f"热图已保存: {output_file}")
                
                plt.show()
                
            except Exception as e:
                self.logger.error(f"聚类热图绘制失败: {e}")
                # 绘制简单热图作为备选
                plt.figure(figsize=(10, 6))
                sns.heatmap(plot_data.T, cmap='RdBu_r', center=0, annot=False)
                plt.title('Integrated Biomarker Analysis')
                
                if output_file:
                    plt.savefig(output_file, dpi=300, bbox_inches='tight')
                    self.logger.info(f"简化热图已保存: {output_file}")
                
                plt.show()
                
        except Exception as e:
            self.logger.error(f"绘制热图失败: {e}")
    
    def plot_biomarker_radar(self, gene: str, output_file: str = None) -> None:
        """
        绘制生物标志物雷达图
        
        Args:
            gene: 基因名
            output_file: 输出文件
        """
        if self.integrated_scores is None or gene not in self.integrated_scores.index:
            self.logger.error(f"基因 {gene} 数据不可用")
            return
        
        try:
            # 准备数据
            categories = []
            values = []
            
            gene_data = self.integrated_scores.loc[gene]
            
            # 收集有效的评分维度
            if 'de_score' in gene_data and not pd.isna(gene_data['de_score']):
                categories.append('Differential\nExpression')
                values.append(max(0, min(1, gene_data['de_score'])))  # 限制在0-1范围
            
            if 'hub_score' in gene_data and not pd.isna(gene_data['hub_score']):
                categories.append('Network\nCentrality')
                values.append(max(0, min(1, gene_data['hub_score'])))
            
            if 'survival_score' in gene_data and not pd.isna(gene_data['survival_score']):
                categories.append('Survival\nImpact')
                values.append(max(0, min(1, gene_data['survival_score'])))
            
            if len(categories) < 2:
                self.logger.warning(f"基因 {gene} 数据不足，无法绘制雷达图")
                return
            
            # 补足至少3个维度用于雷达图
            while len(categories) < 3:
                categories.append('Other')
                values.append(0)
            
            # 雷达图
            angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
            values = values + [values[0]]  # 闭合
            angles = angles + [angles[0]]
            
            fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
            
            ax.plot(angles, values, 'o-', linewidth=2, color='red')
            ax.fill(angles, values, alpha=0.25, color='red')
            
            ax.set_xticks(angles[:-1])
            ax.set_xticklabels(categories)
            ax.set_ylim(0, 1)
            ax.set_title(f'{gene} - Multi-dimensional Score', size=20, y=1.08)
            
            # 添加网格
            ax.grid(True)
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                self.logger.info(f"雷达图已保存: {output_file}")
            
            plt.show()
            
        except Exception as e:
            self.logger.error(f"绘制雷达图失败: {e}")
    
    def generate_report(self, output_dir: str = 'biomarker_report') -> None:
        """
        生成综合报告
        
        Args:
            output_dir: 输出目录
        """
        try:
            output_path = Path(output_dir)
            output_path.mkdir(exist_ok=True)
            
            self.logger.info(f"生成综合报告: {output_path}")
            
            # 1. 保存生物标志物列表
            if self.biomarkers is not None:
                self.biomarkers.to_csv(output_path / 'key_biomarkers.csv')
                self.logger.info("  保存关键生物标志物列表")
            
            # 2. 保存整合得分
            if self.integrated_scores is not None:
                self.integrated_scores.to_csv(output_path / 'integrated_scores.csv')
                self.logger.info("  保存整合得分")
            
            # 3. 生成报告文本
            report_lines = [
                "="*60,
                "生物标志物整合分析报告",
                "="*60,
                ""
            ]
            
            # 差异表达统计
            if self.de_results is not None:
                up_genes = sum(self.de_results['log2FoldChange'] > 0)
                down_genes = sum(self.de_results['log2FoldChange'] < 0)
                report_lines.extend([
                    "1. 差异表达分析",
                    f"   - 显著差异基因数: {len(self.de_results)}",
                    f"   - 上调基因: {up_genes}",
                    f"   - 下调基因: {down_genes}",
                    ""
                ])
            
            # PPI网络统计
            if self.ppi_results and 'network' in self.ppi_results:
                G = self.ppi_results['network']
                report_lines.extend([
                    "2. PPI网络分析",
                    f"   - 网络节点数: {G.number_of_nodes()}",
                    f"   - 网络边数: {G.number_of_edges()}",
                    f"   - Hub基因数: {len(self.ppi_results.get('hub_genes', []))}",
                    f"   - 功能模块数: {len(self.ppi_results.get('modules', {}))}",
                    ""
                ])
            elif not PPI_AVAILABLE:
                report_lines.extend([
                    "2. PPI网络分析",
                    "   - 模块不可用，已跳过",
                    ""
                ])
            
            # 生存分析统计
            if self.survival_results is not None and not self.survival_results.empty:
                sig_survival = self.survival_results[self.survival_results['logrank_pvalue'] < 0.05]
                report_lines.extend([
                    "3. 生存分析",
                    f"   - 分析基因数: {len(self.survival_results)}",
                    f"   - 预后相关基因: {len(sig_survival)}",
                    ""
                ])
            elif not SURVIVAL_AVAILABLE:
                report_lines.extend([
                    "3. 生存分析",
                    "   - 模块不可用，已跳过",
                    ""
                ])
            
            # 生物标志物统计
            if self.biomarkers is not None:
                category_counts = self.biomarkers['category'].value_counts()
                report_lines.extend([
                    "4. 关键生物标志物",
                    f"   - 候选生物标志物总数: {len(self.biomarkers)}",
                ])
                for cat, count in category_counts.items():
                    report_lines.append(f"   - {cat}: {count}")
                
                # Top 5生物标志物
                report_lines.extend(["", "Top 5 生物标志物:"])
                for i, (gene, row) in enumerate(self.biomarkers.head(5).iterrows(), 1):
                    score = row.get('integrated_score', 0)
                    report_lines.append(f"   {i}. {gene} (Score: {score:.3f})")
            
            report_lines.extend(["", "="*60])
            
            # 保存报告
            with open(output_path / 'analysis_report.txt', 'w', encoding='utf-8') as f:
                f.write('\n'.join(report_lines))
            
            self.logger.info("  生成文本报告")
            
            # 4. 生成可视化
            if self.biomarkers is not None and len(self.biomarkers) > 0:
                # 整合热图
                try:
                    self.plot_integration_heatmap(
                        output_file=str(output_path / 'integration_heatmap.png')
                    )
                except Exception as e:
                    self.logger.error(f"生成热图失败: {e}")
                
                # Top基因的雷达图
                try:
                    for gene in self.biomarkers.head(3).index:
                        self.plot_biomarker_radar(
                            gene,
                            output_file=str(output_path / f'{gene}_radar.png')
                        )
                except Exception as e:
                    self.logger.error(f"生成雷达图失败: {e}")
            
            self.logger.info(f"报告生成完成: {output_path}")
            
        except Exception as e:
            self.logger.error(f"生成报告失败: {e}")


def run_integrated_analysis(de_file: str,
                           cancer_type: str = 'BRCA',
                           output_dir: str = 'integrated_analysis',
                           padj_threshold: float = 0.05,
                           log2fc_threshold: float = 1.0,
                           string_score: int = 400,
                           top_biomarkers: int = 20) -> Dict:
    """
    运行完整的整合分析流程
    
    Args:
        de_file: 差异表达结果文件
        cancer_type: 癌症类型
        output_dir: 输出目录
        padj_threshold: 校正p值阈值
        log2fc_threshold: log2FC阈值
        string_score: STRING置信度
        top_biomarkers: 选择的生物标志物数量
        
    Returns:
        分析结果字典
    """
    logger.info("="*60)
    logger.info("开始生物标志物整合分析")
    logger.info("="*60)
    
    try:
        # 创建输出目录
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 创建整合分析器
        integrator = BiomarkerIntegrator()
        
        # 1. 加载差异表达结果
        logger.info("步骤 1/5: 加载差异表达结果")
        de_results = integrator.load_de_results(de_file, padj_threshold, log2fc_threshold)
        
        if de_results.empty:
            logger.error("差异表达数据加载失败，终止分析")
            return {'error': '差异表达数据加载失败'}
        
        # 2. PPI网络分析
        logger.info("步骤 2/5: PPI网络分析")
        integrator.run_ppi_analysis(string_score)
        
        # 3. 生存分析
        logger.info("步骤 3/5: 生存分析")
        integrator.run_survival_analysis(cancer_type)
        
        # 4. 计算整合得分
        logger.info("步骤 4/5: 计算整合得分")
        integrated_scores = integrator.calculate_integrated_scores()
        
        if integrated_scores.empty:
            logger.error("整合得分计算失败")
            return {'error': '整合得分计算失败'}
        
        # 5. 识别关键生物标志物
        logger.info("步骤 5/5: 识别关键生物标志物")
        biomarkers = integrator.identify_key_biomarkers(top_biomarkers)
        
        # 生成报告
        integrator.generate_report(output_dir)
        
        logger.info("="*60)
        logger.info("整合分析完成")
        logger.info("="*60)
        
        return {
            'integrator': integrator,
            'biomarkers': biomarkers,
            'de_results': integrator.de_results,
            'ppi_results': integrator.ppi_results,
            'survival_results': integrator.survival_results,
            'success': True
        }
        
    except Exception as e:
        logger.error(f"整合分析失败: {e}")
        return {'error': str(e), 'success': False}


# 模拟的PPI和生存分析模块（当真实模块不可用时）
class MockPPINetworkAnalyzer:
    """模拟PPI网络分析器"""
    
    def fetch_ppi_network(self, gene_list, string_score=400):
        """模拟获取PPI网络"""
        logger.info("使用模拟PPI网络分析")
        # 返回空DataFrame，表示没有PPI数据
        return pd.DataFrame()
    
    def build_network(self, ppi_data):
        """模拟构建网络"""
        import networkx as nx
        return nx.Graph()
    
    def calculate_centrality_measures(self, G):
        """模拟计算中心性"""
        return pd.DataFrame()
    
    def identify_hub_genes(self, centrality_df, method='top_n', threshold=30):
        """模拟识别hub基因"""
        return []
    
    def detect_modules(self, G):
        """模拟检测模块"""
        return {}


class MockSurvivalAnalyzer:
    """模拟生存分析器"""
    
    def batch_survival_analysis(self, gene_list, cancer_type, output_dir='temp_survival'):
        """模拟批量生存分析"""
        logger.info("使用模拟生存分析")
        # 返回空DataFrame，表示没有生存数据
        return pd.DataFrame()


# 如果原始模块不可用，使用模拟模块
if not PPI_AVAILABLE:
    PPINetworkAnalyzer = MockPPINetworkAnalyzer
    PPI_AVAILABLE = True  # 设置为可用，但使用模拟版本

if not SURVIVAL_AVAILABLE:
    SurvivalAnalyzer = MockSurvivalAnalyzer
    SURVIVAL_AVAILABLE = True  # 设置为可用，但使用模拟版本


if __name__ == "__main__":
    # 测试代码
    import argparse
    
    parser = argparse.ArgumentParser(description='生物标志物整合分析')
    parser.add_argument('--de-file', required=True, help='差异表达结果文件')
    parser.add_argument('--cancer', default='BRCA', help='癌症类型')
    parser.add_argument('--output', default='integrated_results', help='输出目录')
    parser.add_argument('--padj', type=float, default=0.05, help='padj阈值')
    parser.add_argument('--log2fc', type=float, default=1.0, help='log2FC阈值')
    parser.add_argument('--string-score', type=int, default=400, help='STRING置信度')
    parser.add_argument('--top-n', type=int, default=20, help='Top生物标志物数量')
    
    args = parser.parse_args()
    
    # 运行整合分析
    results = run_integrated_analysis(
        args.de_file,
        args.cancer,
        args.output,
        args.padj,
        args.log2fc,
        args.string_score,
        args.top_n
    )
    
    if results.get('success', False):
        print(f"分析完成！结果保存在: {args.output}")
        if 'biomarkers' in results and not results['biomarkers'].empty:
            print("\nTop 5 生物标志物:")
            for i, (gene, row) in enumerate(results['biomarkers'].head(5).iterrows(), 1):
                score = row.get('integrated_score', 0)
                category = row.get('category', 'Unknown')
                print(f"  {i}. {gene} (Score: {score:.3f}, Category: {category})")
    else:
        print(f"分析失败: {results.get('error', '未知错误')}")
        sys.exit(1)
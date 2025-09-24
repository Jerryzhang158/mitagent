#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import argparse
import sys
import warnings
import yaml
import json
import subprocess
import tempfile
import os
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try:
    from ppi_network_module import analyze_de_genes_ppi
    from survival_module import analyze_gene_survival
    from biomarker_integration import run_integrated_analysis
    HAS_BIOMARKER_MODULES = True
    logger.info("✅ 生物标志物分析模块加载成功")
except ImportError:
    logger.warning("⚠️ 生物标志物分析模块未找到，相关功能将不可用")
    logger.warning("   请确保 ppi_network_module.py, survival_module.py 和 biomarker_integration.py 在同一目录")
    HAS_BIOMARKER_MODULES = False

# 检查可选依赖
try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    logger.warning("gseapy未安装，GO/KEGG分析和GSEA将不可用")
    HAS_GSEAPY = False

try:
    from sklearn.preprocessing import quantile_transform
    HAS_SKLEARN = True
except ImportError:
    logger.warning("sklearn未安装，分位数标准化将不可用")
    HAS_SKLEARN = False

try:
    import mygene
    HAS_MYGENE = True
except ImportError:
    logger.warning("mygene未安装，基因ID转换将不可用")
    HAS_MYGENE = False

# 检查R是否可用
HAS_R = False
R_PATH = None
R_ERROR = None

def check_r_installation():
    """检查R是否安装并可用"""
    global HAS_R, R_PATH, R_ERROR
    
    # 尝试多个可能的R路径
    possible_r_paths = [
        'R',  # 系统PATH中的R
        'Rscript',  # 直接使用Rscript
        r'C:\Program Files\R\R-4.3.0\bin\R.exe',  # Windows常见路径
        r'C:\Program Files\R\R-4.2.0\bin\R.exe',
        r'C:\Program Files\R\R-4.1.0\bin\R.exe',
        '/usr/bin/R',  # Linux常见路径
        '/usr/local/bin/R',  # macOS常见路径
    ]
    
    for r_path in possible_r_paths:
        try:
            # 测试R是否可用
            result = subprocess.run([r_path, '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                HAS_R = True
                R_PATH = r_path
                logger.info(f"找到R: {r_path}")
                
                # 安全地获取R版本信息
                try:
                    stdout_lines = result.stdout.strip().splitlines()
                    stderr_lines = result.stderr.strip().splitlines()
                    
                    version_info = ""
                    if stdout_lines:
                        version_info = stdout_lines[0]
                    elif stderr_lines:
                        version_info = stderr_lines[0]
                    
                    if version_info:
                        logger.info(f"R版本: {version_info}")
                    else:
                        logger.info("R版本信息获取失败，但R可用")
                        
                except Exception as e:
                    logger.info(f"R可用，但版本信息解析失败: {e}")
                
                return True
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
            continue
    
    # 如果都没找到
    R_ERROR = "未找到R安装"
    logger.warning("⚠️ 未找到R安装，将使用Python备用方法")
    logger.warning("请安装R: https://cran.r-project.org/")
    return False

# 检查R
check_r_installation()


def create_deseq2_script():
    """创建DESeq2分析的R脚本"""
    r_script = '''
# DESeq2差异表达分析脚本
suppressPackageStartupMessages({
    library(DESeq2)
    library(utils)
})

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    cat("Usage: Rscript deseq2_analysis.R <count_file> <sample_info_file> <control_group> <treatment_group> <output_file> <config_file>\\n")
    quit(status = 1)
}

count_file <- args[1]
sample_info_file <- args[2]
control_group <- args[3]
treatment_group <- args[4]
output_file <- args[5]
config_file <- args[6]

# 读取配置
if (file.exists(config_file)) {
    config <- jsonlite::fromJSON(config_file)
    fit_type <- config$deseq2$fit_type
    test_type <- config$deseq2$test
    alpha <- config$deseq2$alpha
    independent_filtering <- config$deseq2$independent_filtering
    shrink_lfc <- config$deseq2$shrink_lfc
} else {
    # 默认配置
    fit_type <- "parametric"
    test_type <- "Wald"
    alpha <- 0.05
    independent_filtering <- TRUE
    shrink_lfc <- TRUE
}

tryCatch({
    cat("开始DESeq2分析...\\n")
    
    # 读取计数数据
    cat("读取计数数据...\\n")
    count_data <- read.csv(count_file, row.names = 1, stringsAsFactors = FALSE)
    
    # 读取样本信息
    cat("读取样本信息...\\n")
    sample_info <- read.csv(sample_info_file, row.names = 1, stringsAsFactors = FALSE)
    
    # 确保样本顺序一致
    common_samples <- intersect(colnames(count_data), rownames(sample_info))
    count_data <- count_data[, common_samples]
    sample_info <- sample_info[common_samples, ]
    
    # 转换为整数（DESeq2要求）
    count_data <- round(count_data)
    count_data[count_data < 0] <- 0
    
    # 设置条件因子
    sample_info$condition <- factor(sample_info$treatment, levels = c(control_group, treatment_group))
    
    cat(sprintf("对照组: %s (%d个样本)\\n", control_group, sum(sample_info$condition == control_group)))
    cat(sprintf("处理组: %s (%d个样本)\\n", treatment_group, sum(sample_info$condition == treatment_group)))
    cat(sprintf("分析基因数: %d\\n", nrow(count_data)))
    
    # 创建DESeqDataSet
    cat("创建DESeqDataSet...\\n")
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info,
        design = ~ condition
    )
    
    # 过滤低表达基因
    cat("过滤低表达基因...\\n")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat(sprintf("过滤后基因数: %d\\n", nrow(dds)))
    
    # 运行DESeq2
    cat("运行DESeq2分析...\\n")
    dds <- DESeq(dds, fitType = fit_type, test = test_type)
    
    # 获取结果
    cat("获取差异表达结果...\\n")
    res <- results(dds, 
                   contrast = c("condition", treatment_group, control_group),
                   alpha = alpha,
                   independentFiltering = independent_filtering)
    
    # LFC shrinkage（如果可用）
    if (shrink_lfc) {
        cat("应用LFC收缩...\\n")
        res <- lfcShrink(dds, 
                        contrast = c("condition", treatment_group, control_group), 
                        res = res, 
                        type = "normal")
    }
    
    # 转换为数据框
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    
    # 添加标准化计数均值
    normalized_counts <- counts(dds, normalized = TRUE)
    control_samples <- rownames(sample_info)[sample_info$condition == control_group]
    treatment_samples <- rownames(sample_info)[sample_info$condition == treatment_group]
    
    if (length(control_samples) > 0) {
        res_df$control_mean <- rowMeans(normalized_counts[, control_samples, drop = FALSE])
    } else {
        res_df$control_mean <- 0
    }
    
    if (length(treatment_samples) > 0) {
        res_df$treatment_mean <- rowMeans(normalized_counts[, treatment_samples, drop = FALSE])
    } else {
        res_df$treatment_mean <- 0
    }
    
    # 添加fold change（非log）
    res_df$foldChange <- 2^res_df$log2FoldChange
    
    # 重新排列列顺序
    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "foldChange", 
                        "lfcSE", "stat", "pvalue", "padj", 
                        "control_mean", "treatment_mean")]
    
    # 按padj排序
    res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
    
    # 保存结果
    cat("保存结果...\\n")
    write.csv(res_df, output_file, row.names = FALSE)
    
    # 统计结果
    total_genes <- nrow(res_df)
    significant_genes <- sum(res_df$padj < alpha, na.rm = TRUE)
    upregulated <- sum(res_df$padj < alpha & res_df$log2FoldChange >= 1, na.rm = TRUE)
    downregulated <- sum(res_df$padj < alpha & res_df$log2FoldChange <= -1, na.rm = TRUE)
    
    cat("DESeq2分析完成!\\n")
    cat(sprintf("总基因数: %d\\n", total_genes))
    cat(sprintf("显著差异基因: %d\\n", significant_genes))
    cat(sprintf("上调基因: %d\\n", upregulated))
    cat(sprintf("下调基因: %d\\n", downregulated))
    
    # 保存分析摘要  
    summary_file <- paste0(substr(output_file, 1, nchar(output_file)-4), "_summary.txt")
    writeLines(
        c(
            paste("DESeq2 Analysis Summary"),
            paste("======================"),
            paste("Control group:", control_group),
            paste("Treatment group:", treatment_group),
            paste("Total genes:", total_genes),
            paste("Significant genes (padj < ", alpha, "):", significant_genes),
            paste("Upregulated genes (log2FC >= 1):", upregulated),
            paste("Downregulated genes (log2FC <= -1):", downregulated),
            paste("Analysis completed at:", Sys.time())
        ),
        summary_file
    )
    
}, error = function(e) {
    cat("DESeq2分析失败:\\n")
    cat(paste("错误:", e$message, "\\n"))
    quit(status = 1)
})
'''
    return r_script


class HybridRNASeqAnalyzer:
    """混合Python+R的RNA-seq和miRNA-seq数据分析类"""
    
    def __init__(self, config_file=None):
        self.gene_data = None
        self.mirna_data = None
        self.sample_info = None
        self.de_results = {}
        self.enrichment_results = {}
        self.gsea_results = {}
        
        # 加载配置
        self.config = self.load_config(config_file) if config_file else self.get_default_config()
        
        # 创建临时目录用于R脚本
        self.temp_dir = Path(tempfile.mkdtemp(prefix='rnaseq_'))
        logger.info(f"临时文件目录: {self.temp_dir}")
        
        # 检查R和DESeq2状态
        if HAS_R:
            self.check_deseq2_installation()
        else:
            logger.warning("⚠️ R不可用，将使用Python备用方法进行差异表达分析")
    
    def get_default_config(self):
        """获取默认配置"""
        return {
            'sample_groups': {
                # 新格式：CT vs BAI
                'Control': ['CT_', 'CT.', 'C_', 'C.'],
                'Treatment': ['BAI_', 'BAI.', 'B_', 'B.'],
                # 旧格式：保持向后兼容
                'G_Treatment': ['G.L.', 'G.M.', 'G.H.'],
                'D_Treatment': ['D.L.', 'D.M.', 'D.H.']
            },
            'analysis': {
                'min_count': 10,      # DESeq2推荐的最小计数阈值
                'min_samples': 2,     # DESeq2推荐至少2个样本
                'padj_thresh': 0.05,
                'log2fc_thresh': 1.0,
                'normalization_method': 'deseq2'  # 使用DESeq2标准化
            },
            'deseq2': {
                'fit_type': 'parametric',    # parametric, local, mean
                'test': 'Wald',              # Wald, LRT
                'shrink_lfc': True,          # 是否shrink log2 fold changes
                'alpha': 0.05,               # 显著性水平
                'independent_filtering': True # 是否进行独立过滤
            },
            'gsea': {
                'permutation_num': 5000,     # 减少默认置换次数，避免卡死
                'min_size': 15,
                'max_size': 500,
                'ranking_method': 'wald_stat',  # 默认使用Wald统计量
                'random_seed': 42               # 固定随机种子
            }
        }
    
    def load_config(self, config_file):
        """加载配置文件"""
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    return yaml.safe_load(f)
                else:
                    return json.load(f)
        except Exception as e:
            logger.warning(f"加载配置文件失败: {e}，使用默认配置")
            return self.get_default_config()
    
    def check_deseq2_installation(self):
        """检查DESeq2是否在R中安装"""
        try:
            r_check_script = '''
            if (!requireNamespace("DESeq2", quietly = TRUE)) {
                cat("DESeq2_NOT_INSTALLED\\n")
                quit(status = 1)
            } else {
                cat("DESeq2_AVAILABLE\\n")
                tryCatch({
                    cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\\n")
                }, error = function(e) {
                    cat("DESeq2 version: unknown\\n")
                })
            }
            '''
            
            # 写临时R脚本
            check_script_path = self.temp_dir / 'check_deseq2.R'
            with open(check_script_path, 'w', encoding='utf-8') as f:
                f.write(r_check_script)
            
            # 运行检查
            result = subprocess.run([R_PATH, '--slave', '--no-restore', '--file=' + str(check_script_path)], 
                                  capture_output=True, text=True, timeout=30)
            
            # 检查输出（可能在stdout或stderr中）
            output_text = result.stdout + result.stderr
            
            if result.returncode == 0 and 'DESeq2_AVAILABLE' in output_text:
                logger.info("DESeq2在R中可用")
                
                # 尝试提取版本信息
                for line in output_text.split('\n'):
                    if 'DESeq2 version:' in line:
                        logger.info(f"版本信息: {line.strip()}")
                        break
                
                return True
            else:
                logger.warning("DESeq2未在R中安装")
                logger.warning("请在R中运行: install.packages('BiocManager'); BiocManager::install('DESeq2')")
                
                # 显示详细错误信息
                if output_text.strip():
                    logger.debug(f"R输出: {output_text.strip()}")
                
                return False
                
        except subprocess.TimeoutExpired:
            logger.warning("检查DESeq2超时")
            return False
        except Exception as e:
            logger.error(f"检查DESeq2失败: {e}")
            return False
    
    def validate_expression_data(self, data, data_type="expression"):
        """验证表达数据格式和质量"""
        try:
            # 检查数据类型
            if not isinstance(data, pd.DataFrame):
                raise ValueError("数据必须是pandas DataFrame")
            
            # 检查是否为空
            if data.empty:
                raise ValueError("数据为空")
            
            # 检查缺失值
            missing_count = data.isnull().sum().sum()
            if missing_count > 0:
                logger.warning(f"{data_type}数据中存在 {missing_count} 个缺失值")
            
            # 检查负值
            if (data < 0).any().any():
                negative_count = (data < 0).sum().sum()
                logger.warning(f"{data_type}数据中存在 {negative_count} 个负值")
                logger.info("负值将在DESeq2分析前被设为0")
            
            # 对于基因数据，检查是否适合计数分析
            if data_type == 'gene':
                # 检查数据范围，推断是否为计数数据
                max_value = data.max().max()
                min_value = data.min().min()
                
                # 检查是否包含小数
                has_decimals = not data.apply(lambda x: x.apply(lambda y: float(y).is_integer() if pd.notna(y) and y != 0 else True)).all().all()
                
                if max_value < 50 and has_decimals:
                    logger.warning("⚠️ 数据似乎已经标准化过(FPKM/TPM)，DESeq2需要原始计数数据")
                    logger.warning("   建议使用原始计数数据以获得最佳结果")
                elif max_value > 1000000:
                    logger.warning("⚠️ 数据值很大，请确认这是原始计数数据而非标准化数据")
                
                if has_decimals:
                    logger.info("🔧 检测到小数值，将四舍五入为整数用于DESeq2分析")
            
            # 基本统计信息
            logger.info(f"{data_type}数据验证通过: {data.shape[0]} 个基因/特征, {data.shape[1]} 个样本")
            logger.info(f"表达量范围: {data.min().min():.2f} - {data.max().max():.2f}")
            
            return True
            
        except Exception as e:
            logger.error(f"{data_type}数据验证失败: {e}")
            return False
    
    def load_data(self, gene_file=None, mirna_file=None):
        """加载基因表达和miRNA表达数据，支持不同格式和编码"""
        if gene_file:
            try:
                # 尝试多种编码格式
                encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252', 'gbk', 'utf-16']
                
                for encoding in encodings:
                    try:
                        logger.info(f"尝试使用编码 {encoding} 加载基因数据...")
                        
                        # 检测文件格式和分隔符
                        if gene_file.endswith('.csv'):
                            self.gene_data = pd.read_csv(gene_file, index_col=0, encoding=encoding)
                        else:
                            # 检测分隔符
                            with open(gene_file, 'r', encoding=encoding) as f:
                                first_line = f.readline().strip()
                            
                            if '\t' in first_line:
                                sep = '\t'
                            elif ',' in first_line:
                                sep = ','
                            else:
                                sep = '\t'  # 默认
                            
                            self.gene_data = pd.read_csv(gene_file, sep=sep, index_col=0, encoding=encoding)
                        
                        logger.info(f"成功加载基因数据: {self.gene_data.shape}")
                        break
                        
                    except UnicodeDecodeError:
                        if encoding == encodings[-1]:
                            raise
                        continue
                    except Exception as e:
                        if encoding == encodings[-1]:
                            raise
                        continue
                
                # 数据清理和验证
                if self.validate_expression_data(self.gene_data, "基因表达"):
                    # 清理样本名称
                    self.gene_data.columns = [str(col).replace('\x00', '').replace('\ufeff', '').replace('Ã¿Ã¾', '').strip() 
                                             for col in self.gene_data.columns]
                    
                    # 转换为数值类型
                    for col in self.gene_data.columns:
                        try:
                            self.gene_data[col] = pd.to_numeric(self.gene_data[col], errors='coerce')
                        except:
                            logger.warning(f"列 {col} 转换为数字类型失败")
                    
                    # 处理缺失值
                    self.gene_data = self.gene_data.fillna(0)
                    
                    # 确保非负值（DESeq2要求）
                    self.gene_data = self.gene_data.clip(lower=0)
                    
                    logger.info(f"基因表达数据加载完成: {self.gene_data.shape}")
                    logger.info(f"样本列表: {list(self.gene_data.columns)}")
                    logger.info(f"基因ID示例: {list(self.gene_data.index[:5])}")
                    logger.info(f"表达量统计: min={self.gene_data.min().min():.2f}, max={self.gene_data.max().max():.2f}, mean={self.gene_data.mean().mean():.2f}")
                else:
                    self.gene_data = None
                    
            except Exception as e:
                logger.error(f"加载基因表达数据失败: {e}")
                self.gene_data = None
        
        # miRNA数据加载（保持简单的统计方法）
        if mirna_file:
            try:
                encodings = ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252', 'gbk', 'utf-16']
                
                for encoding in encodings:
                    try:
                        if mirna_file.endswith('.txt'):
                            self.mirna_data = pd.read_csv(mirna_file, sep='\t', index_col=0, encoding=encoding)
                        else:
                            self.mirna_data = pd.read_csv(mirna_file, index_col=0, encoding=encoding)
                        break
                    except UnicodeDecodeError:
                        if encoding == encodings[-1]:
                            raise
                        continue
                
                if self.validate_expression_data(self.mirna_data, "miRNA表达"):
                    # 数据清理
                    for col in self.mirna_data.columns:
                        try:
                            self.mirna_data[col] = pd.to_numeric(self.mirna_data[col], errors='coerce')
                        except:
                            pass
                    
                    self.mirna_data = self.mirna_data.fillna(0)
                    logger.info(f"miRNA表达数据加载完成: {self.mirna_data.shape}")
                    logger.info(f"样本列表: {list(self.mirna_data.columns)}")
                    logger.info(f"miRNA ID示例: {list(self.mirna_data.index[:5])}")
                    logger.info(f"表达量统计: min={self.mirna_data.min().min():.2f}, max={self.mirna_data.max().max():.2f}, mean={self.mirna_data.mean().mean():.2f}")
                else:
                    self.mirna_data = None
                    
            except Exception as e:
                logger.error(f"加载miRNA表达数据失败: {e}")
                self.mirna_data = None
    
    def create_sample_info(self, data, data_type='gene'):
        """增强版样本分组方法 - 根据数据类型选择不同的分组策略，支持剂量分组"""
        samples = data.columns.tolist()
        sample_info = []
        
        logger.info(f"开始分组{data_type}数据样本: {samples}")
        
        for sample in samples:
            # 清理样本名称
            clean_sample = str(sample).replace('\x00', '').replace('\ufeff', '').replace('Ã¿Ã¾', '').strip()
            
            group = 'Unknown'
            treatment = 'Unknown'
            
            # 优先尝试剂量模式匹配
            group, treatment = self._match_dose_sample(clean_sample)
            
            # 如果剂量模式失败，使用原有的匹配方法
            if group == 'Unknown':
                # 根据数据类型选择不同的分组策略
                if data_type == 'gene':
                    # 基因数据：优先匹配CT_Rep和BAI_Rep模式
                    group, treatment = self._match_gene_sample(clean_sample)
                elif data_type == 'mirna':
                    # miRNA数据：优先匹配C1-3和B1-3模式
                    group, treatment = self._match_mirna_sample(clean_sample)
                
                # 如果特定匹配失败，使用通用匹配
                if group == 'Unknown':
                    group, treatment = self._match_general_sample(clean_sample)
            
            sample_info.append({
                'sample': sample,
                'clean_sample': clean_sample,
                'group': group,
                'treatment': treatment
            })
        
        sample_df = pd.DataFrame(sample_info).set_index('sample')
        
        # 记录分组信息
        group_counts = sample_df['treatment'].value_counts()
        logger.info(f"{data_type}数据样本分组信息:")
        for group, count in group_counts.items():
            group_samples = sample_df[sample_df['treatment'] == group]['clean_sample'].tolist()
            logger.info(f"  {group}: {count} 个样本 - {group_samples}")
        
        # 检查未分组样本
        unknown_samples = sample_df[sample_df['treatment'] == 'Unknown']['clean_sample'].tolist()
        if unknown_samples:
            logger.warning(f"以下{data_type}样本未能自动分组: {unknown_samples}")
            logger.warning("请检查配置文件中的样本分组设置")
        
        return sample_df
    
    def _match_gene_sample(self, sample_name):
        """匹配基因数据样本（CT_Rep, BAI_Rep等）"""
        sample_lower = sample_name.lower()
        
        # 基因数据的特定模式
        if any(pattern in sample_lower for pattern in ['ct_rep', 'ct_', 'control']):
            return 'Control', 'Control'
        elif any(pattern in sample_lower for pattern in ['bai_rep', 'bai_', 'treatment']):
            return 'Treatment', 'Treatment'
        
        return 'Unknown', 'Unknown'
    
    def _match_mirna_sample(self, sample_name):
        """匹配miRNA数据样本（C1, C2, C3, B1, B2, B3等）"""
        # miRNA数据的特定模式
        if sample_name in ['C1', 'C2', 'C3', 'c1', 'c2', 'c3']:
            return 'Control', 'Control'
        elif sample_name in ['B1', 'B2', 'B3', 'b1', 'b2', 'b3']:
            return 'Treatment', 'Treatment'
        
        # 更通用的模式
        if sample_name.startswith('C') and len(sample_name) <= 2:
            return 'Control', 'Control'
        elif sample_name.startswith('B') and len(sample_name) <= 2:
            return 'Treatment', 'Treatment'
        
        return 'Unknown', 'Unknown'
    
    def _match_dose_sample(self, sample_name):
        """匹配剂量相关样本命名（Low/Medium/High） - 修复版：优先匹配完整词汇"""
        sample_lower = sample_name.lower()
        
        # 🔧 修复：按优先级顺序匹配，完整词汇优先
        # 1. 优先匹配完整的关键词（避免部分匹配导致的错误）
        if any(pattern in sample_lower for pattern in ['control', 'ctrl']):
            return 'Control', 'Control'
        elif any(pattern in sample_lower for pattern in ['medium', 'mediuem']):  # 包含拼写错误的情况
            return 'Medium', 'Medium'
        elif sample_lower.startswith('low') or any(pattern in sample_lower for pattern in ['low_', 'low1', 'low2', 'low3']):
            return 'Low', 'Low'
        elif sample_lower.startswith('high') or any(pattern in sample_lower for pattern in ['high_', 'high1', 'high2', 'high3']):
            return 'High', 'High'
        
        # 2. 然后匹配简短模式（但要避免与control冲突）
        elif sample_lower.startswith('l_') or sample_lower.startswith('l.') or sample_lower.endswith('_l') or sample_lower in ['l1', 'l2', 'l3']:
            return 'Low', 'Low'
        elif sample_lower.startswith('m_') or sample_lower.startswith('m.') or sample_lower.endswith('_m') or sample_lower in ['m1', 'm2', 'm3', 'mid']:
            return 'Medium', 'Medium'
        elif sample_lower.startswith('h_') or sample_lower.startswith('h.') or sample_lower.endswith('_h') or sample_lower in ['h1', 'h2', 'h3']:
            return 'High', 'High'
        
        # 3. 最后匹配其他control模式（但不包括可能冲突的简短模式）
        elif any(pattern in sample_lower for pattern in ['con', 'nc', 'neg']) and not any(x in sample_lower for x in ['control', 'ctrl']):
            return 'Control', 'Control'
        
        return 'Unknown', 'Unknown'
    
    def _order_samples_custom(self, sample_info, custom_order):
        """根据自定义顺序排列样本"""
        ordered_samples = []
        
        for group in custom_order:
            group_samples = sample_info[sample_info['treatment'] == group].index.tolist()
            # 在组内按样本名排序以保证一致性
            group_samples.sort()
            ordered_samples.extend(group_samples)
        
        # 添加未在custom_order中指定的样本
        remaining_samples = [s for s in sample_info.index if s not in ordered_samples]
        remaining_samples.sort()
        ordered_samples.extend(remaining_samples)
        
        return ordered_samples

    def _order_samples_default(self, sample_info):
        """默认样本排序：Control -> Low -> Medium -> High"""
        # 定义期望的组顺序 - 修改为Control优先
        preferred_order = ['Control', 'Low', 'Medium', 'High', 'Treatment']
        
        ordered_samples = []
        
        # 按预定义顺序添加样本
        for group in preferred_order:
            group_samples = sample_info[sample_info['treatment'] == group].index.tolist()
            if group_samples:
                group_samples.sort()  # 组内排序
                ordered_samples.extend(group_samples)
        
        # 添加其他未分类的组
        remaining_groups = [g for g in sample_info['treatment'].unique() 
                           if g not in preferred_order]
        
        for group in sorted(remaining_groups):
            group_samples = sample_info[sample_info['treatment'] == group].index.tolist()
            group_samples.sort()
            ordered_samples.extend(group_samples)
        
        return ordered_samples

    def _add_group_separators(self, clustermap_obj, sample_info_ordered):
        """在热图中添加组间分割线"""
        try:
            ax = clustermap_obj.ax_heatmap
            
            # 计算组边界
            group_boundaries = []
            current_group = None
            
            for i, treatment in enumerate(sample_info_ordered['treatment']):
                if current_group != treatment:
                    if current_group is not None:
                        group_boundaries.append(i)
                    current_group = treatment
            
            # 绘制分割线
            for boundary in group_boundaries:
                ax.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)
            
            logger.info(f"添加了 {len(group_boundaries)} 条组间分割线")
            
        except Exception as e:
            logger.warning(f"添加组间分割线失败: {e}")
    
    def _match_general_sample(self, sample_name):
        """通用样本匹配（使用配置文件）"""
        for group_name, group_config in self.config['sample_groups'].items():
            matched = False
            
            if isinstance(group_config, dict):
                # 新格式：字典配置
                if 'patterns' in group_config:
                    patterns = group_config['patterns']
                    match_mode = group_config.get('match_mode', 'prefix')
                    
                    for pattern in patterns:
                        pattern = str(pattern).strip()
                        
                        if match_mode == 'exact':
                            if sample_name == pattern:
                                matched = True
                                break
                        elif match_mode == 'prefix':
                            if sample_name.startswith(pattern):
                                matched = True
                                break
                        elif match_mode == 'suffix':
                            if sample_name.endswith(pattern):
                                matched = True
                                break
                        elif match_mode == 'contains':
                            if pattern.lower() in sample_name.lower():
                                matched = True
                                break
                    
                    if matched:
                        return group_name, group_config.get('alias', group_name)
                        
            elif isinstance(group_config, list):
                # 旧格式：列表配置
                for pattern in group_config:
                    if sample_name.startswith(str(pattern)):
                        return group_name, group_name
        
        # 最后的通用规则
        clean_lower = sample_name.lower()
        
        if any(keyword in clean_lower for keyword in ['control', 'ctrl', 'baseline', 'normal']):
            return 'Control', 'Control'
        elif any(keyword in clean_lower for keyword in ['treatment', 'treat', 'drug', 'compound']):
            return 'Treatment', 'Treatment'
        # 数字模式 (如: 1-1., 3-2. 等)
        elif any(char.isdigit() for char in sample_name) and '-' in sample_name:
            return 'Treatment', 'Treatment'
        
        return 'Unknown', 'Unknown'
    
    def run_r_deseq2_analysis(self, count_data, sample_info, contrast, data_type='gene'):
        """使用R脚本运行DESeq2分析 - 修复编码问题版本"""
        if not HAS_R:
            logger.error("❌ R不可用，无法运行DESeq2分析")
            return None
        
        try:
            logger.info("开始DESeq2分析...")
            
            # 获取对比组信息
            control_group, treatment_group = contrast
            
            # 准备数据文件
            count_file = self.temp_dir / f'count_data_{control_group}_vs_{treatment_group}.csv'
            sample_file = self.temp_dir / f'sample_info_{control_group}_vs_{treatment_group}.csv'
            config_file = self.temp_dir / f'config_{control_group}_vs_{treatment_group}.json'
            output_file = self.temp_dir / f'deseq2_results_{control_group}_vs_{treatment_group}.csv'
            r_script_file = self.temp_dir / 'deseq2_analysis.R'
            
            # 筛选相关样本
            relevant_samples = sample_info[
                sample_info['treatment'].isin([control_group, treatment_group])
            ]
            count_subset = count_data[relevant_samples.index]
            
            logger.info(f"对照组: {control_group} ({sum(relevant_samples['treatment'] == control_group)} 个样本)")
            logger.info(f"处理组: {treatment_group} ({sum(relevant_samples['treatment'] == treatment_group)} 个样本)")
            logger.info(f"分析基因数: {count_subset.shape[0]}")
            
            # 保存计数数据（整数化）
            count_subset_int = count_subset.round().astype(int)
            count_subset_int.to_csv(count_file)
            
            # 保存样本信息
            relevant_samples.to_csv(sample_file)
            
            # 保存配置
            with open(config_file, 'w') as f:
                json.dump(self.config, f)
            
            # 创建R脚本
            with open(r_script_file, 'w', encoding='utf-8') as f:
                f.write(create_deseq2_script())
            
            # 运行R脚本 - 修复编码问题
            logger.info("执行R脚本...")
            cmd = [
                R_PATH, '--slave', '--no-restore',
                '--file=' + str(r_script_file),
                '--args',
                str(count_file),
                str(sample_file),
                control_group,
                treatment_group,  
                str(output_file),
                str(config_file)
            ]
            
            # 关键修复：设置正确的编码和环境变量
            env = os.environ.copy()
            env['LC_ALL'] = 'C'  # 使用C locale避免编码问题
            env['LANG'] = 'C'
            
            try:
                # 首先尝试UTF-8编码
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True, 
                    timeout=600,
                    encoding='utf-8',  # 强制使用UTF-8
                    errors='replace',  # 替换无法解码的字符
                    env=env
                )
            except UnicodeDecodeError:
                # 如果UTF-8失败，尝试系统默认编码
                logger.warning("UTF-8编码失败，尝试系统默认编码...")
                try:
                    result = subprocess.run(
                        cmd, 
                        capture_output=True, 
                        text=True, 
                        timeout=600,
                        encoding='gbk',  # Windows系统编码
                        errors='replace',
                        env=env
                    )
                except:
                    # 最后尝试不指定编码
                    result = subprocess.run(
                        cmd, 
                        capture_output=True, 
                        timeout=600,
                        env=env
                    )
                    # 手动解码
                    try:
                        stdout_text = result.stdout.decode('utf-8', errors='replace') if result.stdout else ''
                        stderr_text = result.stderr.decode('utf-8', errors='replace') if result.stderr else ''
                    except:
                        stdout_text = result.stdout.decode('gbk', errors='replace') if result.stdout else ''
                        stderr_text = result.stderr.decode('gbk', errors='replace') if result.stderr else ''
                    
                    # 创建一个简单的结果对象
                    class SimpleResult:
                        def __init__(self, returncode, stdout, stderr):
                            self.returncode = returncode
                            self.stdout = stdout
                            self.stderr = stderr
                    
                    result = SimpleResult(result.returncode, stdout_text, stderr_text)
            
            # 显示R的输出 - 安全处理编码
            if hasattr(result, 'stdout') and result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        # 清理可能的编码问题字符
                        try:
                            clean_line = line.encode('utf-8', errors='replace').decode('utf-8')
                            logger.info(f"R> {clean_line}")
                        except:
                            logger.info(f"R> {line}")
            
            if hasattr(result, 'stderr') and result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip() and 'Loading required package' not in line:
                        try:
                            clean_line = line.encode('utf-8', errors='replace').decode('utf-8')
                            logger.warning(f"R> {clean_line}")
                        except:
                            logger.warning(f"R> {line}")
            
            # 检查是否成功
            if result.returncode != 0:
                logger.error(f"❌ R脚本执行失败 (返回码: {result.returncode})")
                return None
            
            if not output_file.exists():
                logger.error("DESeq2结果文件未生成")
                return None
            
            # 读取结果
            logger.info("读取DESeq2分析结果...")
            results_df = pd.read_csv(output_file)
            
            if 'gene_id' in results_df.columns:
                results_df = results_df.set_index('gene_id')
            
            logger.info(f"DESeq2分析完成，获得 {len(results_df)} 个基因的结果")
            
            # 读取分析摘要 - 安全处理编码
            summary_file = str(output_file).replace('.csv', '_summary.txt')
            if os.path.exists(summary_file):
                try:
                    # 尝试多种编码读取摘要文件
                    encodings = ['utf-8', 'gbk', 'latin-1']
                    summary_content = None
                    
                    for encoding in encodings:
                        try:
                            with open(summary_file, 'r', encoding=encoding) as f:
                                summary_content = f.read()
                                break
                        except UnicodeDecodeError:
                            continue
                    
                    if summary_content:
                        logger.info("分析摘要:")
                        for line in summary_content.split('\n')[2:]:  # 跳过标题行
                            if line.strip():
                                logger.info(f"  {line}")
                    else:
                        logger.warning("无法读取摘要文件（编码问题）")
                        
                except Exception as e:
                    logger.warning(f"读取摘要文件时出错: {e}")
            
            return results_df
            
        except subprocess.TimeoutExpired:
            logger.error("R脚本执行超时")
            return None
        except Exception as e:
            logger.error(f"R DESeq2分析失败: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def run_python_fallback_analysis(self, count_data, sample_info, contrast, data_type='gene', method='ttest'):
        """Python备用差异表达分析"""
        logger.warning("使用Python备用方法进行差异表达分析")
        
        try:
            # 标准化数据
            normalized_data = self.normalize_counts(count_data, method='log2_cpm')
            
            # 获取两组样本
            control_group, treatment_group = contrast
            control_samples = sample_info[sample_info['treatment'] == control_group].index
            treatment_samples = sample_info[sample_info['treatment'] == treatment_group].index
            
            control_data = normalized_data[control_samples]
            treatment_data = normalized_data[treatment_samples]
            
            logger.info(f"对照组样本数: {len(control_samples)}")
            logger.info(f"处理组样本数: {len(treatment_samples)}")
            
            results = []
            epsilon = 1e-8
            
            for gene in normalized_data.index:
                control_expr = control_data.loc[gene].values
                treatment_expr = treatment_data.loc[gene].values
                
                # 基本统计量
                control_mean = np.mean(control_expr)
                treatment_mean = np.mean(treatment_expr)
                
                # Fold change计算
                fold_change = (treatment_mean + epsilon) / (control_mean + epsilon)
                log2_fold_change = np.log2(fold_change)
                
                # 统计检验
                if method == 'ttest' and len(control_expr) >= 2 and len(treatment_expr) >= 2:
                    stat, pvalue = stats.ttest_ind(control_expr, treatment_expr, equal_var=False)
                elif method == 'mannwhitney' and len(control_expr) >= 1 and len(treatment_expr) >= 1:
                    try:
                        stat, pvalue = stats.mannwhitneyu(control_expr, treatment_expr, alternative='two-sided')
                    except:
                        pvalue = 1.0
                        stat = 0.0
                else:
                    pvalue = 1.0
                    stat = 0.0
                
                results.append({
                    'gene_id': gene,
                    'baseMean': (control_mean + treatment_mean) / 2,
                    'log2FoldChange': log2_fold_change,
                    'foldChange': fold_change,
                    'control_mean': control_mean,
                    'treatment_mean': treatment_mean,
                    'stat': stat,
                    'pvalue': pvalue
                })
            
            # 转换为DataFrame
            result_df = pd.DataFrame(results).set_index('gene_id')
            
            # 多重检验校正
            valid_pvalues = result_df['pvalue'].replace([np.inf, -np.inf], np.nan).dropna()
            if len(valid_pvalues) > 0:
                rejected, padj, alpha_sidak, alpha_bonf = multipletests(
                    valid_pvalues, alpha=0.05, method='fdr_bh'
                )
                result_df['padj'] = np.nan
                result_df.loc[valid_pvalues.index, 'padj'] = padj
            else:
                result_df['padj'] = 1.0
            
            result_df = result_df.sort_values('pvalue')
            
            logger.warning(f"备用方法完成，分析了 {len(result_df)} 个基因")
            logger.warning("注意：建议安装R和DESeq2以获得更准确的结果")
            
            return result_df
            
        except Exception as e:
            logger.error(f"备用分析方法失败: {e}")
            return None
    
    def normalize_counts(self, count_data, method='log2_cpm'):
        """标准化计数数据（用于可视化和备用分析）"""
        count_data = count_data.astype(float)
        
        if method == 'log2_cpm':
            # CPM + log2 transformation
            cpm = count_data.div(count_data.sum(axis=0), axis=1) * 1e6
            normalized = np.log2(cpm + 1)
        elif method == 'log2':
            normalized = np.log2(count_data + 1)
        elif method == 'quantile' and HAS_SKLEARN:
            normalized = pd.DataFrame(
                quantile_transform(count_data.T).T,
                index=count_data.index,
                columns=count_data.columns
            )
        else:
            normalized = count_data
            
        return normalized
    
    def filter_low_expression(self, count_data, min_count=None, min_samples=None):
        """过滤低表达基因"""
        if min_count is None:
            min_count = self.config['analysis']['min_count']
        if min_samples is None:
            min_samples = self.config['analysis']['min_samples']
        
        keep_genes = (count_data >= min_count).sum(axis=1) >= min_samples
        filtered_data = count_data[keep_genes]
        
        logger.info(f"过滤前: {count_data.shape[0]} 个基因")
        logger.info(f"过滤后: {filtered_data.shape[0]} 个基因") 
        logger.info(f"过滤标准: 至少{min_samples}个样本中表达量>={min_count}")
        
        return filtered_data
    
    def run_differential_expression(self, count_data, sample_info, contrast, data_type='gene'):
        """运行差异表达分析 - 优先使用R+DESeq2，否则使用Python备用方法"""
        try:
            # 先过滤低表达基因
            filtered_data = self.filter_low_expression(count_data)
            
            # 对于基因和miRNA数据，都优先尝试R+DESeq2（修改这里）
            if HAS_R:  # 移除data_type == 'gene'的限制
                logger.info(f"尝试使用R+DESeq2进行{data_type}差异表达分析...")
                result = self.run_r_deseq2_analysis(filtered_data, sample_info, contrast, data_type)
                
                if result is not None:
                    return result
                else:
                    logger.warning("R+DESeq2分析失败，切换到Python备用方法...")
            
            # 使用Python备用方法
            method = 'ttest'  # 对基因和miRNA都使用t-test
            logger.info(f"使用Python {method}方法分析{data_type}数据...")
            
            return self.run_python_fallback_analysis(filtered_data, sample_info, contrast, data_type, method)
            
        except Exception as e:
            logger.error(f"差异表达分析失败: {e}")
            return None
        
    # 以下方法保持不变，用于可视化和富集分析
    def filter_significant_genes(self, de_result, padj_thresh=None, log2fc_thresh=None):
        """筛选显著差异表达的基因"""
        if de_result is None or de_result.empty:
            return None
        
        if padj_thresh is None:
            padj_thresh = self.config['analysis']['padj_thresh']
        if log2fc_thresh is None:
            log2fc_thresh = self.config['analysis']['log2fc_thresh']
        
        # 处理无穷大值
        valid_results = de_result.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2FoldChange', 'padj'])
        
        significant = valid_results[
            (valid_results['padj'] < padj_thresh) & 
            (abs(valid_results['log2FoldChange']) >= log2fc_thresh)  # 改为 >= 包含等于
        ]
        
        return significant
    
    def plot_volcano(self, de_result, title="", 
                    padj_thresh=None, log2fc_thresh=None, figsize=(10, 8)):
        """绘制火山图 - 移除标题"""
        if de_result is None or de_result.empty:
            logger.error("没有差异表达结果用于绘制火山图")
            return None
        
        if padj_thresh is None:
            padj_thresh = self.config['analysis']['padj_thresh']
        if log2fc_thresh is None:
            log2fc_thresh = self.config['analysis']['log2fc_thresh']
            
        plt.figure(figsize=figsize)
        
        # 准备数据
        data = de_result.copy()
        data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2FoldChange', 'padj'])
        
        if len(data) == 0:
            logger.error("没有有效数据用于绘制火山图")
            return None
            
        data['-log10(padj)'] = -np.log10(data['padj'].clip(lower=1e-300))
        
        # 分类基因
        data['significant'] = 'Not Significant'
        data.loc[(data['padj'] < padj_thresh) & (data['log2FoldChange'] >= log2fc_thresh), 'significant'] = 'Up-regulated'
        data.loc[(data['padj'] < padj_thresh) & (data['log2FoldChange'] <= -log2fc_thresh), 'significant'] = 'Down-regulated'
        
        # 绘制散点图
        colors = {'Not Significant': 'gray', 'Up-regulated': 'red', 'Down-regulated': 'blue'}
        for sig_type, color in colors.items():
            subset = data[data['significant'] == sig_type]
            plt.scatter(subset['log2FoldChange'], subset['-log10(padj)'], 
                       c=color, label=f"{sig_type} ({len(subset)})", alpha=0.6, s=20)
        
        # 添加阈值线
        plt.axhline(y=-np.log10(padj_thresh), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=log2fc_thresh, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-log2fc_thresh, color='black', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10(adjusted P-value)')
        # 移除标题设置
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        return plt.gcf()
    
    def plot_ma(self, de_result, title="", 
                padj_thresh=None, log2fc_thresh=None, figsize=(10, 8)):
        """绘制MA图 - 移除标题"""
        if de_result is None or de_result.empty:
            logger.error("没有差异表达结果用于绘制MA图")
            return None
        
        if padj_thresh is None:
            padj_thresh = self.config['analysis']['padj_thresh']
        if log2fc_thresh is None:
            log2fc_thresh = self.config['analysis']['log2fc_thresh']
            
        plt.figure(figsize=figsize)
        
        # 准备数据
        data = de_result.copy()
        data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2FoldChange', 'baseMean'])
        
        if len(data) == 0:
            logger.error("没有有效数据用于绘制MA图")
            return None
            
        data['log10_baseMean'] = np.log10(data['baseMean'].clip(lower=1))
        
        # 分类基因
        data['significant'] = 'Not Significant'
        data.loc[(data['padj'] < padj_thresh) & (data['log2FoldChange'] >= log2fc_thresh), 'significant'] = 'Up-regulated'
        data.loc[(data['padj'] < padj_thresh) & (data['log2FoldChange'] <= -log2fc_thresh), 'significant'] = 'Down-regulated'
        
        # 绘制散点图
        colors = {'Not Significant': 'gray', 'Up-regulated': 'red', 'Down-regulated': 'blue'}
        for sig_type, color in colors.items():
            subset = data[data['significant'] == sig_type]
            plt.scatter(subset['log10_baseMean'], subset['log2FoldChange'], 
                       c=color, label=f"{sig_type} ({len(subset)})", alpha=0.6, s=20)
        
        # 添加阈值线
        plt.axhline(y=log2fc_thresh, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=-log2fc_thresh, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        plt.xlabel('Log10(Base Mean)')
        plt.ylabel('Log2 Fold Change')
        # 移除标题设置
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        return plt.gcf()
    
    def plot_heatmap(self, de_result, data_type='gene', title="", 
                    top_n=50, figsize=(12, 10), custom_order=None):
        """绘制表达热图 - 支持自定义样本顺序，移除标题"""
        if de_result is None or de_result.empty:
            logger.error("没有差异表达结果用于绘制热图")
            return None
            
        # 选择数据源
        if data_type == 'gene' and self.gene_data is not None:
            expr_data = self.gene_data
        elif data_type == 'mirna' and self.mirna_data is not None:
            expr_data = self.mirna_data
        else:
            logger.error(f"没有可用的{data_type}数据")
            return None
        
        # 获取显著差异基因
        significant_genes = self.filter_significant_genes(de_result)
        if significant_genes is None or len(significant_genes) == 0:
            logger.warning("没有显著差异表达的基因用于绘制热图")
            return None
        
        # 选择top_n个基因
        genes_to_plot = significant_genes.head(top_n).index
        genes_to_plot = genes_to_plot.intersection(expr_data.index)
        
        if len(genes_to_plot) == 0:
            logger.warning("没有匹配的基因用于绘制热图")
            return None
        
        logger.info(f"热图将展示 {len(genes_to_plot)} 个显著差异表达的{data_type}")
        
        plot_data = expr_data.loc[genes_to_plot]
        
        # 标准化数据
        plot_data_norm = self.normalize_counts(plot_data, method='log2_cpm')
        
        # Z-score标准化
        try:
            plot_data_scaled = plot_data_norm.T
            means = plot_data_scaled.mean()
            stds = plot_data_scaled.std()
            stds = stds.replace(0, 1).fillna(1)
            stds[stds < 1e-8] = 1
            
            plot_data_scaled = (plot_data_scaled - means) / stds
            plot_data_scaled = plot_data_scaled.T
            
            plot_data_scaled = plot_data_scaled.replace([np.inf, -np.inf], np.nan).fillna(0)
            
            if not np.all(np.isfinite(plot_data_scaled.values)):
                logger.error("数据标准化后仍包含非有限值")
                return None
                
        except Exception as e:
            logger.error(f"数据标准化失败: {e}")
            return None
        
        # 创建样本信息用于颜色标注
        sample_info = self.create_sample_info(expr_data, data_type)
        
        # 自定义样本排序逻辑
        if custom_order:
            ordered_samples = self._order_samples_custom(sample_info, custom_order)
        else:
            # 默认排序：Control -> Low -> Medium -> High
            ordered_samples = self._order_samples_default(sample_info)
        
        # 重新排列数据和样本信息
        plot_data_ordered = plot_data_scaled[ordered_samples]
        sample_info_ordered = sample_info.loc[ordered_samples]
        
        try:
            plt.figure(figsize=figsize)
            
            # 创建颜色映射
            unique_treatments = sample_info_ordered['treatment'].unique()
            colors = plt.cm.Set1(np.linspace(0, 1, len(unique_treatments)))
            treatment_colors = dict(zip(unique_treatments, colors))
            col_colors = sample_info_ordered['treatment'].map(treatment_colors)
            
            # 绘制聚类热图 - 关闭列聚类保持自定义顺序
            g = sns.clustermap(plot_data_ordered, 
                              col_colors=col_colors,
                              cmap='RdBu_r', 
                              center=0,
                              figsize=figsize,
                              cbar_kws={'label': 'Z-score'},
                              yticklabels=True,
                              xticklabels=True,
                              row_cluster=True,    # 保持基因聚类
                              col_cluster=False,   # 关闭样本聚类，保持自定义顺序
                              method='average',
                              metric='euclidean')
            
            # 移除标题设置
            # g.fig.suptitle(title, y=1.02)
            
            # 添加分组分割线
            self._add_group_separators(g, sample_info_ordered)
            
            return g.fig
            
        except Exception as e:
            logger.error(f"绘制热图失败: {e}")
            # 简化热图备选方案
            try:
                plt.figure(figsize=figsize)
                im = plt.imshow(plot_data_ordered.values, cmap='RdBu_r', aspect='auto')
                plt.colorbar(im, label='Z-score')
                plt.xticks(range(len(plot_data_ordered.columns)), 
                          plot_data_ordered.columns, rotation=45)
                plt.yticks(range(len(plot_data_ordered.index)), 
                          plot_data_ordered.index)
                # 移除标题设置
                plt.tight_layout()
                return plt.gcf()
            except Exception as e2:
                logger.error(f"简化热图也失败: {e2}")
                return None
    
    def convert_ensembl_to_symbols(self, ensembl_ids, batch_size=1000):
        """转换ENSEMBL ID为基因符号"""
        if not HAS_MYGENE:
            logger.warning("mygene未安装，将直接使用ENSEMBL ID")
            return [gene_id.split('.')[0] for gene_id in ensembl_ids]
        
        try:
            mg = mygene.MyGeneInfo()
            clean_ids = [gene_id.split('.')[0] for gene_id in ensembl_ids]
            
            logger.info(f"开始转换 {len(clean_ids)} 个ENSEMBL ID...")
            
            all_converted = {}
            
            # 批量转换
            if len(clean_ids) <= batch_size:
                results = mg.querymany(clean_ids, 
                                     scopes='ensembl.gene', 
                                     fields='symbol', 
                                     species='human',
                                     returnall=True)
                
                for result in results['out']:
                    if 'symbol' in result and 'query' in result:
                        all_converted[result['query']] = result['symbol']
            else:
                # 分批处理
                for i in range(0, len(clean_ids), batch_size):
                    batch = clean_ids[i:i+batch_size]
                    results = mg.querymany(batch, 
                                         scopes='ensembl.gene', 
                                         fields='symbol', 
                                         species='human',
                                         returnall=True)
                    
                    for result in results['out']:
                        if 'symbol' in result and 'query' in result:
                            all_converted[result['query']] = result['symbol']
            
            # 返回转换结果
            symbols = []
            for ensembl_id in clean_ids:
                if ensembl_id in all_converted:
                    symbols.append(all_converted[ensembl_id])
                else:
                    symbols.append(ensembl_id)
            
            success_rate = len(all_converted) / len(clean_ids) * 100
            logger.info(f"基因ID转换完成: {len(all_converted)}/{len(clean_ids)} ({success_rate:.1f}%)")
            
            return symbols
            
        except Exception as e:
            logger.error(f"基因ID转换失败: {e}")
            return [gene_id.split('.')[0] for gene_id in ensembl_ids]
    
    def run_enrichment_analysis(self, gene_list, organism='human', gene_sets=None):
        """运行GO和KEGG富集分析"""
        logger.info(f"开始GO/KEGG富集分析，输入基因数: {len(gene_list)}")
        
        if not HAS_GSEAPY:
            logger.error("gseapy未安装，无法运行富集分析")
            logger.error("请安装: pip install gseapy")
            return None
            
        if len(gene_list) == 0:
            logger.error("基因列表为空")
            return None
        
        # 自动选择基因集
        if gene_sets is None:
            try:
                logger.info("获取可用的基因集...")
                available_libs = gp.get_library_name()
                logger.info(f"总共有 {len(available_libs)} 个可用的基因集")
                
                go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
                kegg_sets = [lib for lib in available_libs if 'KEGG' in lib and 'Human' in lib]
                
                logger.info(f"找到 {len(go_sets)} 个GO生物过程基因集")
                logger.info(f"找到 {len(kegg_sets)} 个KEGG基因集")
                
                gene_sets = []
                if go_sets:
                    selected_go = go_sets[0]  # 选择第一个（通常是最新的）
                    gene_sets.append(selected_go)
                    logger.info(f"选择GO基因集: {selected_go}")
                if kegg_sets:
                    selected_kegg = kegg_sets[0]
                    gene_sets.append(selected_kegg)
                    logger.info(f"选择KEGG基因集: {selected_kegg}")
                
                if not gene_sets:
                    logger.error("找不到合适的GO/KEGG基因集")
                    logger.info("可用基因集示例:")
                    for i, lib in enumerate(available_libs[:10]):
                        logger.info(f"  {i+1}. {lib}")
                    return None
                    
            except Exception as e:
                logger.error(f"获取基因集失败: {e}")
                return None
        
        try:
            # 转换基因ID
            logger.info("转换ENSEMBL ID为基因符号...")
            gene_symbols = self.convert_ensembl_to_symbols(gene_list)
            
            # 过滤有效的基因符号
            valid_symbols = []
            for symbol in gene_symbols:
                if (symbol is not None and pd.notna(symbol) and symbol != '' and 
                    not symbol.startswith('ENSG')):  # 过滤掉未转换的ENSEMBL ID
                    valid_symbols.append(symbol)
            
            logger.info(f"转换结果: {len(gene_list)} -> {len(gene_symbols)} -> {len(valid_symbols)} 个有效基因符号")
            
            if len(valid_symbols) == 0:
                logger.error("没有有效的基因符号用于富集分析")
                return None
            
            if len(valid_symbols) < 5:
                logger.warning(f"基因数量较少({len(valid_symbols)})，富集分析结果可能不可靠")
            
            logger.info(f"基因符号示例: {valid_symbols[:10]}")
            logger.info(f"使用基因集: {gene_sets}")
            
            # 运行富集分析
            logger.info("调用gseapy.enrichr进行富集分析...")
            enr = gp.enrichr(gene_list=valid_symbols,
                           gene_sets=gene_sets,
                           organism='Human',
                           outdir=None,
                           cutoff=0.05)
            
            logger.info("富集分析调用完成，处理结果...")
            
            # 检查结果并正确解析
            if enr and hasattr(enr, 'results'):
                logger.info(f"获得富集分析结果，结果类型: {type(enr.results)}")
                
                # 标准化结果格式
                final_results = {}
                
                if isinstance(enr.results, dict):
                    logger.info(f"结果是字典格式，包含 {len(enr.results)} 个基因集")
                    # 如果是字典格式（每个基因集一个DataFrame）
                    for gene_set_name, result_df in enr.results.items():
                        logger.info(f"处理基因集: {gene_set_name}, 类型: {type(result_df)}")
                        
                        if isinstance(result_df, pd.DataFrame):
                            logger.info(f"  DataFrame形状: {result_df.shape}")
                            if not result_df.empty:
                                # 验证必需列是否存在
                                required_cols = ['Term', 'P-value', 'Adjusted P-value']
                                available_cols = list(result_df.columns)
                                logger.info(f"  可用列: {available_cols}")
                                
                                missing_cols = [col for col in required_cols if col not in available_cols]
                                if missing_cols:
                                    logger.warning(f"  缺少必需列: {missing_cols}")
                                else:
                                    final_results[gene_set_name] = result_df.copy()
                                    logger.info(f"  ✅ 基因集 {gene_set_name}: {len(result_df)} 个显著富集条目")
                            else:
                                logger.info(f"  基因集 {gene_set_name}: DataFrame为空")
                        else:
                            logger.warning(f"  基因集 {gene_set_name}: 不是DataFrame类型")
                    
                elif isinstance(enr.results, pd.DataFrame) and not enr.results.empty:
                    logger.info(f"结果是单一DataFrame格式，形状: {enr.results.shape}")
                    
                    # 按基因集分组
                    if 'Gene_set' in enr.results.columns:
                        for gene_set_name in enr.results['Gene_set'].unique():
                            gene_set_df = enr.results[enr.results['Gene_set'] == gene_set_name].copy()
                            if not gene_set_df.empty:
                                final_results[gene_set_name] = gene_set_df
                                logger.info(f"基因集 {gene_set_name}: {len(gene_set_df)} 个显著富集条目")
                    else:
                        # 如果没有Gene_set列，使用第一个基因集名称
                        first_gene_set = gene_sets[0] if gene_sets else 'Combined_Results'
                        final_results[first_gene_set] = enr.results.copy()
                        logger.info(f"基因集 {first_gene_set}: {len(enr.results)} 个显著富集条目")
                else:
                    logger.warning("富集分析结果为空或格式不正确")
                    logger.info(f"结果类型: {type(enr.results)}")
                    if hasattr(enr.results, 'shape'):
                        logger.info(f"结果形状: {enr.results.shape}")
                
                if final_results:
                    logger.info(f"✅ 富集分析成功完成，获得 {len(final_results)} 个基因集的结果")
                    
                    # 创建结果包装器
                    class EnrichmentResults:
                        def __init__(self, results_dict):
                            self.results = results_dict
                    
                    return EnrichmentResults(final_results)
                else:
                    logger.warning("❌ 没有有效的富集分析结果")
                    return None
            else:
                logger.warning("❌ 富集分析没有返回有效结果")
                logger.info(f"enr类型: {type(enr)}")
                if enr:
                    logger.info(f"enr属性: {dir(enr)}")
                return None
                
        except Exception as e:
            logger.error(f"❌ 富集分析失败: {e}")
            logger.error(f"错误类型: {type(e).__name__}")
            import traceback
            logger.error(f"详细错误信息:\n{traceback.format_exc()}")
            logger.error(f"可能的原因: 网络问题、基因ID格式问题或基因集名称错误")
            return None

    def debug_available_gene_sets(self):
        """🆕 调试：查看所有可用的基因集，重点关注WikiPathways"""
        if not HAS_GSEAPY:
            logger.error("gseapy未安装")
            return
        
        try:
            logger.info("正在获取所有可用基因集...")
            available_libs = gp.get_library_name()
            
            # 按类别分组
            categories = {
                'WikiPathways': [lib for lib in available_libs if 'WikiPathway' in lib],
                'KEGG': [lib for lib in available_libs if 'KEGG' in lib],
                'GO': [lib for lib in available_libs if 'GO_' in lib],
                'MSigDB': [lib for lib in available_libs if 'MSigDB' in lib or 'Hallmark' in lib],
                'Reactome': [lib for lib in available_libs if 'Reactome' in lib],
                'BioPlanet': [lib for lib in available_libs if 'BioPlanet' in lib],
            }
            
            logger.info("="*60)
            logger.info("可用基因集详情")
            logger.info("="*60)
            
            for category, gene_sets in categories.items():
                if gene_sets:
                    logger.info(f"{category} ({len(gene_sets)} 个):")
                    for gene_set in sorted(gene_sets):
                        logger.info(f"  - {gene_set}")
                else:
                    logger.info(f"{category}: 无")
            
            # 🎯 查找最合适的WikiPathways基因集
            wiki_sets = categories.get('WikiPathways', [])
            if wiki_sets:
                recommended_wiki = None
                for pattern in ['2024_Human', '2023_Human', '2022_Human', '2021_Human', '2020_Human']:
                    matching = [lib for lib in wiki_sets if pattern in lib]
                    if matching:
                        recommended_wiki = matching[0]
                        break
                
                if recommended_wiki:
                    logger.info(f"🎯 推荐的WikiPathways基因集: {recommended_wiki}")
                else:
                    logger.info(f"🎯 推荐的WikiPathways基因集: {wiki_sets[0]} (最新可用)")
            else:
                logger.warning("❌ 没有找到WikiPathways基因集")
            
            logger.info("="*60)
            
        except Exception as e:
            logger.error(f"获取基因集列表失败: {e}")

    def _get_available_pathways(self):
        """🆕 获取所有可用的通路数据库，优先WikiPathways"""
        try:
            available_libs = gp.get_library_name()
            logger.info(f"总共有 {len(available_libs)} 个可用的基因集")
            
            # 查找各类基因集
            wiki_sets = [lib for lib in available_libs if 'WikiPathway' in lib and 'Human' in lib]
            kegg_sets = [lib for lib in available_libs if 'KEGG' in lib and 'Human' in lib]
            hallmark_sets = [lib for lib in available_libs if 'Hallmark' in lib]
            go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
            
            logger.info(f"WikiPathways基因集: {wiki_sets}")
            logger.info(f"KEGG基因集: {kegg_sets}")
            logger.info(f"Hallmark基因集: {hallmark_sets}")
            logger.info(f"GO基因集: {go_sets[:3]}...")  # 只显示前3个
            
            # 🎯 优先选择WikiPathways
            selected = []
            
            # 1. 优先WikiPathways - 选择最新版本
            if wiki_sets:
                # 按版本排序，选择最新的
                wiki_sorted = sorted(wiki_sets, reverse=True)  # 按名称排序，新版本在前
                selected.append(wiki_sorted[0])
                logger.info(f"✅ 选择WikiPathways: {wiki_sorted[0]}")
            
            # 2. 添加KEGG
            if kegg_sets and len(selected) < 3:
                kegg_sorted = sorted(kegg_sets, reverse=True)
                selected.append(kegg_sorted[0])
                logger.info(f"✅ 选择KEGG: {kegg_sorted[0]}")
            
            # 3. 添加Hallmark（通常较小，运行快）
            if hallmark_sets and len(selected) < 3:
                selected.append(hallmark_sets[0])
                logger.info(f"✅ 选择Hallmark: {hallmark_sets[0]}")
            
            if not selected:
                logger.error("❌ 没有找到合适的基因集")
                return None
                
            return selected
            
        except Exception as e:
            logger.error(f"获取可用通路失败: {e}")
            return None

    def _standardize_gsea_results(self, gs_res):
        """🆕 标准化GSEA结果字段名，处理不同版本的gseapy + 修复数据类型"""
        if not gs_res or not hasattr(gs_res, 'res2d'):
            return None
        
        result_df = gs_res.res2d.copy()
        if result_df.empty:
            return None
        
        # 打印原始列名以便调试
        logger.info(f"GSEA原始结果列名: {list(result_df.columns)}")
        
        # 字段名映射表 - 处理不同版本的gseapy
        field_mapping = {
            # FDR相关字段
            'fdr': ['fdr', 'FDR q-val', 'FDR_q_val', 'qval', 'FDR', 'padj'],
            # NES相关字段  
            'nes': ['nes', 'NES', 'Normalized Enrichment Score'],
            # ES相关字段
            'es': ['es', 'ES', 'Enrichment Score'],
            # 基因集名称
            'term': ['term', 'Term', 'pathway', 'Pathway', 'gene_set', 'Name'],
            # p值 - 🆕 添加更多可能的p值字段
            'pval': ['pval', 'P-value', 'pvalue', 'P_value', 'NOM p-val', 'nom_pval'],
            # 基因数量
            'size': ['size', 'Size', 'gene_number', 'genes_num', 'Tag %', 'Gene %']
        }
        
        # 执行字段标准化
        for standard_name, possible_names in field_mapping.items():
            found = False
            for possible_name in possible_names:
                if possible_name in result_df.columns:
                    if possible_name != standard_name:
                        result_df[standard_name] = result_df[possible_name]
                        logger.info(f"字段映射: {possible_name} -> {standard_name}")
                    found = True
                    break
            
            if not found:
                logger.warning(f"⚠️ 未找到字段: {standard_name} (尝试过: {possible_names})")
                # 为缺失字段设置默认值
                if standard_name == 'fdr':
                    result_df['fdr'] = 1.0
                elif standard_name == 'nes':
                    result_df['nes'] = 0.0
                elif standard_name == 'pval':
                    result_df['pval'] = 1.0
        
        # 🔧 关键修复：数据类型转换
        numeric_fields = ['fdr', 'nes', 'es', 'pval']
        
        for field in numeric_fields:
            if field in result_df.columns:
                try:
                    # 先处理特殊值
                    original_values = result_df[field].copy()
                    
                    # 将常见的非数值标记转换为NaN
                    result_df[field] = result_df[field].replace({
                        '---': np.nan,
                        '--': np.nan,
                        'nan': np.nan,
                        'NaN': np.nan,
                        'NA': np.nan,
                        'null': np.nan,
                        '': np.nan
                    })
                    
                    # 转换为数值类型
                    result_df[field] = pd.to_numeric(result_df[field], errors='coerce')
                    
                    # 对于关键字段，用合理默认值填充NaN
                    if field == 'fdr':
                        result_df[field] = result_df[field].fillna(1.0)  # FDR默认为1（不显著）
                    elif field == 'nes':
                        result_df[field] = result_df[field].fillna(0.0)  # NES默认为0（无富集）
                    elif field == 'pval':
                        result_df[field] = result_df[field].fillna(1.0)  # p值默认为1（不显著）
                    
                    # 记录转换结果
                    converted_count = (~result_df[field].isna()).sum()
                    total_count = len(result_df)
                    logger.info(f"🔧 {field}字段转换: {converted_count}/{total_count} 个数值成功转换")
                    
                except Exception as e:
                    logger.error(f"❌ {field}字段数据类型转换失败: {e}")
                    # 设置安全默认值
                    if field == 'fdr':
                        result_df[field] = 1.0
                    elif field == 'nes':
                        result_df[field] = 0.0
                    elif field == 'pval':
                        result_df[field] = 1.0
        
        # 验证数据类型
        for field in ['fdr', 'nes']:
            if field in result_df.columns:
                dtype = result_df[field].dtype
                if not pd.api.types.is_numeric_dtype(dtype):
                    logger.warning(f"⚠️ {field}字段仍然不是数值类型: {dtype}")
                else:
                    logger.info(f"✅ {field}字段数据类型: {dtype}")
        numeric_fields = ['fdr', 'nes', 'es', 'pval']
        for field in numeric_fields:
            if field in result_df.columns:
                # 处理非数值标记
                result_df[field] = result_df[field].replace(['---', '--', 'nan', 'NaN'], np.nan)
                # 转换为数值
                result_df[field] = pd.to_numeric(result_df[field], errors='coerce')
                # 填充默认值
                if field == 'fdr':
                    result_df[field] = result_df[field].fillna(1.0)
                elif field in ['nes', 'es']:
                    result_df[field] = result_df[field].fillna(0.0)
                elif field == 'pval':
                    result_df[field] = result_df[field].fillna(1.0)
        
        return result_df

    def run_gsea_analysis(self, de_result, gene_sets=None, **kwargs):
        """🔧 修复版GSEA分析方法 - 确保WikiPathways优先使用和字段名兼容性"""
        
        # 提取关键参数
        ranking_method = kwargs.get('ranking_method', self.config['gsea'].get('ranking_method', 'wald_stat'))
        permutation_num = kwargs.get('permutation_num', self.config['gsea'].get('permutation_num', 5000))
        random_seed = kwargs.get('random_seed', self.config['gsea'].get('random_seed', 42))
        min_size = kwargs.get('min_size', self.config['gsea'].get('min_size', 15))
        max_size = kwargs.get('max_size', self.config['gsea'].get('max_size', 500))
        
        if not HAS_GSEAPY:
            logger.error("gseapy未安装，无法运行GSEA分析")
            return None
        
        if de_result is None or de_result.empty:
            logger.error("没有差异表达结果")
            return None
        
        try:
            logger.info("开始修复版GSEA分析（WikiPathways优先）...")
            logger.info(f"参数: ranking_method={ranking_method}, permutation_num={permutation_num}, seed={random_seed}")
            
            # 1. 严格的数据清理
            ranked_data = de_result.copy()
            
            # 移除无穷大和NaN值
            ranked_data = ranked_data.replace([np.inf, -np.inf], np.nan)
            initial_count = len(ranked_data)
            
            # 根据ranking method选择必要的列
            required_cols = ['log2FoldChange']
            if ranking_method == 'wald_stat':
                required_cols.append('stat')
            elif ranking_method == 'signed_pvalue':
                required_cols.append('pvalue')
            elif ranking_method == 'signal_to_noise':
                required_cols.extend(['pvalue', 'lfcSE'])
            
            ranked_data = ranked_data.dropna(subset=required_cols)
            logger.info(f"数据清理: {initial_count} -> {len(ranked_data)} 个基因")
            
            if len(ranked_data) == 0:
                logger.error("数据清理后没有有效基因")
                return None
            
            # 2. 计算ranking metric - 使用更标准的方法
            if ranking_method == 'wald_stat' and 'stat' in ranked_data.columns:
                # 优先使用DESeq2的Wald统计量（最接近标准GSEA）
                metric = ranked_data['stat'].copy()
                logger.info("使用Wald统计量作为ranking metric")
                
            elif ranking_method == 'log2fc':
                metric = ranked_data['log2FoldChange'].copy()
                logger.info("使用log2FoldChange作为ranking metric")
                
            elif ranking_method == 'signed_pvalue':
                # 更保守的p值处理
                pval_clipped = ranked_data['pvalue'].clip(lower=1e-300, upper=1-1e-10)
                metric = -np.log10(pval_clipped) * np.sign(ranked_data['log2FoldChange'])
                logger.info("使用带符号p值作为ranking metric")
                
            elif ranking_method == 'signal_to_noise':
                # Signal-to-noise ratio
                if 'lfcSE' in ranked_data.columns:
                    se_safe = ranked_data['lfcSE'].clip(lower=1e-6)  # 避免除零
                    metric = ranked_data['log2FoldChange'] / se_safe
                    logger.info("使用信噪比作为ranking metric")
                else:
                    logger.warning("缺少lfcSE列，回退到Wald统计量")
                    metric = ranked_data.get('stat', ranked_data['log2FoldChange'])
            else:
                logger.warning(f"未知的ranking方法: {ranking_method}，使用log2FC")
                metric = ranked_data['log2FoldChange']
            
            # 处理重复值问题 - 添加微小随机扰动
            unique_count = len(metric.unique())
            if unique_count < len(metric) * 0.95:  # 如果重复值超过5%
                logger.warning(f"检测到大量重复值 ({len(metric)-unique_count}/{len(metric)})，添加微小扰动")
                np.random.seed(random_seed)
                metric = metric + np.random.normal(0, abs(metric.std()) * 1e-10, len(metric))
            
            # 检查metric的有效性
            if metric.isna().any():
                logger.warning(f"Ranking metric中有 {metric.isna().sum()} 个NaN值，将被移除")
                valid_idx = ~metric.isna()
                ranked_data = ranked_data[valid_idx]
                metric = metric[valid_idx]
            
            logger.info(f"Ranking metric统计: min={metric.min():.3f}, max={metric.max():.3f}, mean={metric.mean():.3f}")
            
            # 3. 基因ID转换 - 更谨慎的处理
            logger.info("转换基因ID...")
            logger.info("⏳ 正在批量转换基因ID，可能需要几分钟...")
            
            try:
                gene_symbols = self.convert_ensembl_to_symbols(ranked_data.index.tolist(), batch_size=500)  # 减小批量大小
            except Exception as e:
                logger.error(f"基因ID转换失败: {e}")
                logger.info("使用原始基因ID继续分析...")
                gene_symbols = [gene_id.split('.')[0] for gene_id in ranked_data.index]
            
            # 构建有效的基因-得分映射
            gene_score_pairs = []
            conversion_stats = {'total': 0, 'converted': 0, 'duplicates': 0}
            
            for gene_id, symbol, score in zip(ranked_data.index, gene_symbols, metric):
                conversion_stats['total'] += 1
                
                # 严格的基因符号验证
                if (symbol and 
                    pd.notna(symbol) and 
                    isinstance(symbol, str) and
                    symbol.strip() != '' and 
                    not symbol.startswith('ENSG') and 
                    not symbol.startswith('ENS') and
                    len(symbol) > 1):  # 基因符号至少2个字符
                    
                    gene_score_pairs.append({'gene': symbol.strip(), 'score': float(score)})
                    conversion_stats['converted'] += 1
            
            logger.info(f"基因ID转换: {conversion_stats['converted']}/{conversion_stats['total']} "
                       f"({conversion_stats['converted']/conversion_stats['total']*100:.1f}%)")
            
            if len(gene_score_pairs) == 0:
                logger.error("没有有效基因符号用于GSEA")
                return None
            
            # 4. 处理重复基因 - 保留绝对值最大的得分
            gene_scores = {}
            for pair in gene_score_pairs:
                gene = pair['gene']
                score = pair['score']
                
                if gene in gene_scores:
                    # 保留绝对值更大的得分
                    if abs(score) > abs(gene_scores[gene]):
                        gene_scores[gene] = score
                    conversion_stats['duplicates'] += 1
                else:
                    gene_scores[gene] = score
            
            logger.info(f"重复基因处理: 移除 {conversion_stats['duplicates']} 个重复，"
                       f"保留 {len(gene_scores)} 个唯一基因")
            
            # 5. 创建最终的ranking DataFrame
            rnk_df = pd.DataFrame([
                {'gene_symbol': gene, 'ranking_metric': score}
                for gene, score in gene_scores.items()
            ])
            
            # 排序 - 高得分在前
            rnk_df = rnk_df.sort_values('ranking_metric', ascending=False).reset_index(drop=True)
            
            logger.info(f"最终基因列表: {len(rnk_df)} 个基因")
            logger.info(f"得分范围: {rnk_df['ranking_metric'].min():.3f} 到 {rnk_df['ranking_metric'].max():.3f}")
            
            # 显示top/bottom基因
            logger.info("Top 10基因:")
            for _, row in rnk_df.head(10).iterrows():
                logger.info(f"  {row['gene_symbol']}: {row['ranking_metric']:.3f}")
            
            logger.info("Bottom 10基因:")  
            for _, row in rnk_df.tail(10).iterrows():
                logger.info(f"  {row['gene_symbol']}: {row['ranking_metric']:.3f}")
            
            # 6. 🎯 基因集选择 - 确保WikiPathways优先
            if gene_sets is None:
                gene_sets = self._get_available_pathways()  # 使用新方法确保WikiPathways优先
            
            if not gene_sets:
                logger.error("没有可用的基因集")
                return None
                
            logger.info(f"使用基因集: {gene_sets}")
            
            # 7. 运行GSEA - 🔧 修复字段名问题
            gsea_results = {}
            
            for i, gene_set in enumerate(gene_sets, 1):
                try:
                    logger.info(f"运行GSEA ({i}/{len(gene_sets)}): {gene_set}")
                    logger.info("⏳ 正在运行GSEA，请耐心等待...")
                    
                    # 关键：确保参数一致性
                    gs_res = gp.prerank(
                        rnk=rnk_df,                        # 使用我们准备的数据
                        gene_sets=gene_set,                # 基因集
                        processes=1,                       # 单进程确保一致性
                        permutation_num=permutation_num,   # 置换次数
                        outdir=None,                       # 不输出到文件
                        max_size=max_size,                 # 最大基因集大小
                        min_size=min_size,                 # 最小基因集大小  
                        weighted_score_type=1,             # 标准GSEA权重
                        ascending=False,                   # 降序排列
                        seed=random_seed,                  # 固定随机种子
                        verbose=False                      # 减少输出
                    )
                    
                    if gs_res and hasattr(gs_res, 'res2d') and not gs_res.res2d.empty:
                        # 🔧 关键修复：使用标准化方法处理结果
                        result_df = self._standardize_gsea_results(gs_res)
                        
                        if result_df is not None and not result_df.empty:
                            # 添加分析参数信息
                            result_df['ranking_method'] = ranking_method
                            result_df['permutation_num'] = permutation_num
                            result_df['random_seed'] = random_seed
                            result_df['total_genes_input'] = len(rnk_df)
                            
                            gsea_results[gene_set] = result_df
                            
                            # 使用标准化后的字段名统计
                            sig_pos = result_df[(result_df['fdr'] < 0.25) & (result_df['nes'] > 0)]
                            sig_neg = result_df[(result_df['fdr'] < 0.25) & (result_df['nes'] < 0)]
                            
                            logger.info(f"✅ {gene_set}: 总计 {len(result_df)} 个通路")
                            logger.info(f"   显著正向富集: {len(sig_pos)} 个")
                            logger.info(f"   显著负向富集: {len(sig_neg)} 个")
                            
                            # 显示top结果
                            if len(sig_pos) > 0:
                                top_pos = sig_pos.nsmallest(3, 'fdr')
                                logger.info("   Top正向:")
                                for _, row in top_pos.iterrows():
                                    logger.info(f"     {row['term']}: NES={row['nes']:.2f}, FDR={row['fdr']:.3f}")
                            
                            if len(sig_neg) > 0:
                                top_neg = sig_neg.nsmallest(3, 'fdr')
                                logger.info("   Top负向:")
                                for _, row in top_neg.iterrows():
                                    logger.info(f"     {row['term']}: NES={row['nes']:.2f}, FDR={row['fdr']:.3f}")
                        else:
                            logger.warning(f"❌ {gene_set}: 结果标准化失败")
                    else:
                        logger.warning(f"❌ {gene_set}: 没有结果")
                        
                except Exception as e:
                    logger.error(f"❌ GSEA失败 {gene_set}: {e}")
                    # 打印详细错误信息以便调试
                    import traceback
                    logger.error(f"详细错误:\n{traceback.format_exc()}")
                    logger.info("继续处理下一个基因集...")
                    continue
            
            if gsea_results:
                logger.info(f"✅ GSEA分析完成，获得 {len(gsea_results)} 个基因集的结果")
                
                # 创建结果包装器
                class FixedGSEAResults:
                    def __init__(self, results_dict, ranked_genes, analysis_params):
                        self.results = results_dict
                        self.ranked_genes = ranked_genes
                        self.parameters = analysis_params
                        
                    def save_ranked_genes(self, filename):
                        """保存ranking文件，便于与其他工具比较"""
                        try:
                            self.ranked_genes.to_csv(filename, sep='\t', index=False, header=False)
                            logger.info(f"Ranking文件已保存: {filename}")
                        except Exception as e:
                            logger.error(f"保存ranking文件失败: {e}")
                
                analysis_params = {
                    'ranking_method': ranking_method,
                    'permutation_num': permutation_num,
                    'random_seed': random_seed,
                    'min_size': min_size,
                    'max_size': max_size,
                    'total_genes': len(rnk_df),
                    'gene_sets': gene_sets,
                    'conversion_rate': conversion_stats['converted']/conversion_stats['total']
                }
                
                return FixedGSEAResults(gsea_results, rnk_df, analysis_params)
            else:
                logger.warning("❌ 没有GSEA结果")
                return None
                
        except Exception as e:
            logger.error(f"GSEA分析失败: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None
    
    def _select_reproducible_gene_sets(self):
        """选择稳定可重现的基因集"""
        try:
            available_libs = gp.get_library_name()
            
            # 优先选择版本明确、更新频率稳定的基因集
            priority_order = [
                'WikiPathway_2023_Human',     # 🆕 WikiPathways优先
                'WikiPathway_2021_Human',     
                'GO_Biological_Process_2021', 
                'KEGG_2021_Human',
                'MSigDB_Hallmark_2020',
                'WikiPathway_2019_Human',     # 备选版本
                'GO_Biological_Process_2018',
                'KEGG_2019_Human'
            ]
            
            selected = []
            for gene_set in priority_order:
                if gene_set in available_libs:
                    selected.append(gene_set)
                    logger.info(f"选择基因集: {gene_set}")
                    if len(selected) >= 3:  # 增加到3个主要基因集
                        break
            
            if not selected:
                # 备选方案 - 动态查找
                wiki_sets = [lib for lib in available_libs if 'WikiPathway' in lib and 'Human' in lib]
                go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
                kegg_sets = [lib for lib in available_libs if 'KEGG' in lib and ('Human' in lib or '2021' in lib)]
                
                if wiki_sets:
                    selected.append(wiki_sets[0])  # 优先WikiPathways
                if go_sets and len(selected) < 3:
                    selected.append(go_sets[0])
                if kegg_sets and len(selected) < 3:
                    selected.append(kegg_sets[0])
            
            if not selected:
                raise ValueError("无法找到合适的基因集")
            
            return selected
            
        except Exception as e:
            logger.error(f"基因集选择失败: {e}")
            # 默认基因集
            return ['GO_Biological_Process_2021', 'KEGG_2019_Human']
    
    def plot_enrichment(self, enrichment_result, gene_set_name, top_n=20, figsize=(12, 8)):
        """绘制富集分析结果"""
        logger.info(f"开始绘制富集分析图: {gene_set_name}")
        
        if enrichment_result is None:
            logger.error("没有富集分析结果用于绘图")
            return None
            
        if not hasattr(enrichment_result, 'results'):
            logger.error("富集分析结果缺少results属性")
            return None
            
        if gene_set_name not in enrichment_result.results:
            logger.error(f"找不到基因集: {gene_set_name}")
            logger.info(f"可用的基因集: {list(enrichment_result.results.keys())}")
            return None
            
        try:
            result_df = enrichment_result.results[gene_set_name].copy()
            logger.info(f"基因集 {gene_set_name} 数据形状: {result_df.shape}")
            
            if result_df.empty:
                logger.warning(f"基因集 {gene_set_name} 没有富集条目")
                return None
            
            # 检查必需的列
            required_columns = ['Term']
            available_columns = list(result_df.columns)
            logger.info(f"可用列: {available_columns}")
            
            score_column = None
            
            # 查找可用的得分列
            possible_score_columns = ['Combined Score', 'Score', 'Enrichment Score', '-log10(P-value)']
            for col in possible_score_columns:
                if col in result_df.columns:
                    score_column = col
                    logger.info(f"使用得分列: {score_column}")
                    break
            
            if score_column is None:
                logger.error(f"富集结果缺少得分列，可用列: {available_columns}")
                return None
            
            missing_columns = [col for col in required_columns if col not in result_df.columns]
            if missing_columns:
                logger.error(f"富集结果缺少必需的列: {missing_columns}")
                return None
                
            # 选择top_n个条目并排序
            original_length = len(result_df)
            if len(result_df) > top_n:
                result_df = result_df.head(top_n)
                logger.info(f"从 {original_length} 个条目中选择前 {top_n} 个")
            
            # 按得分排序
            result_df = result_df.sort_values(score_column, ascending=True)
            
            plt.figure(figsize=figsize)
            
            # 创建条形图
            y_pos = range(len(result_df))
            bars = plt.barh(y_pos, result_df[score_column])
            
            # 根据p值设置颜色（如果有的话）
            if 'P-value' in result_df.columns or 'Adjusted P-value' in result_df.columns:
                p_col = 'Adjusted P-value' if 'Adjusted P-value' in result_df.columns else 'P-value'
                p_values = pd.to_numeric(result_df[p_col], errors='coerce').fillna(1.0)
                
                if len(p_values) > 0 and p_values.min() < p_values.max():
                    colors = plt.cm.viridis_r((p_values - p_values.min()) / (p_values.max() - p_values.min()))
                    for bar, color in zip(bars, colors):
                        bar.set_color(color)
                        
                    # 添加colorbar
                    try:
                        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis_r, 
                                                 norm=plt.Normalize(vmin=p_values.min(), vmax=p_values.max()))
                        sm.set_array([])
                        cbar = plt.colorbar(sm, ax=plt.gca())
                        cbar.set_label(p_col)
                    except Exception as e:
                        logger.warning(f"添加颜色条失败: {e}")
            
            # 设置标签
            term_labels = []
            for term in result_df['Term']:
                if isinstance(term, str):
                    # 截断过长的标签
                    if len(term) > 60:
                        term_labels.append(term[:57] + '...')
                    else:
                        term_labels.append(term)
                else:
                    term_labels.append(str(term))
            
            plt.yticks(y_pos, term_labels)
            plt.xlabel(score_column)
            plt.title(f'{gene_set_name} Enrichment Analysis (Top {len(result_df)})')
            plt.tight_layout()
            
            logger.info(f"成功生成富集分析图: {gene_set_name}")
            return plt.gcf()
            
        except Exception as e:
            logger.error(f"绘制富集分析图失败: {e}")
            import traceback
            logger.error(f"详细错误:\n{traceback.format_exc()}")
            return None
    
    def analyze_comparison(self, control_group, treatment_group, data_type='gene'):
        """进行组间比较分析 - 增强版"""
        logger.info(f"开始分析 {treatment_group} vs {control_group} ({data_type})")
        
        # 选择数据
        if data_type == 'gene' and self.gene_data is not None:
            data = self.gene_data
        elif data_type == 'mirna' and self.mirna_data is not None:
            data = self.mirna_data
        else:
            logger.error(f"没有可用的{data_type}数据")
            return None
            
        # 创建样本信息 - 传递数据类型
        sample_info = self.create_sample_info(data, data_type)
        
        # 筛选相关样本
        relevant_samples = sample_info[
            sample_info['treatment'].isin([control_group, treatment_group])
        ]
        
        if len(relevant_samples) == 0:
            logger.error(f"没有找到相关样本")
            logger.error(f"期望的分组: {control_group}, {treatment_group}")
            logger.info(f"实际可用的分组: {sample_info['treatment'].unique().tolist()}")
            
            # 尝试提供修复建议
            available_groups = sample_info['treatment'].unique().tolist()
            if 'Unknown' in available_groups:
                logger.warning("部分样本未能正确分组，请检查样本命名或配置文件")
                logger.info("样本详情:")
                for group in available_groups:
                    group_samples = sample_info[sample_info['treatment'] == group]['clean_sample'].tolist()
                    logger.info(f"  {group}: {group_samples}")
            return None
            
        data_subset = data[relevant_samples.index]
        
        # 运行差异表达分析
        de_result = self.run_differential_expression(data_subset, relevant_samples, 
                                                   [control_group, treatment_group], data_type)
        
        if de_result is not None:
            # 保存结果
            comparison_name = f"{treatment_group}_vs_{control_group}_{data_type}"
            self.de_results[comparison_name] = de_result
            
            # 筛选显著基因
            significant_genes = self.filter_significant_genes(de_result)
            
            if significant_genes is not None and len(significant_genes) > 0:
                logger.info(f"发现 {len(significant_genes)} 个显著差异表达的{data_type}")
            else:
                logger.warning(f"没有显著差异表达的{data_type}")
                significant_genes = pd.DataFrame()
                
            return {
                'de_result': de_result,
                'significant_genes': significant_genes,
                'sample_info': relevant_samples,
                'data_subset': data_subset
            }
        
        return None
    
    def run_all_enrichment_analysis(self, run_enrichment=True, run_gsea=True):
        """对所有有显著基因的比较进行富集分析"""
        if not run_enrichment and not run_gsea:
            return
            
        if not HAS_GSEAPY:
            logger.warning("gseapy未安装，跳过富集分析")
            return
            
        gene_comparisons = [comp for comp in self.de_results.keys() if comp.endswith('_gene')]
        
        if not gene_comparisons:
            logger.warning("没有基因差异表达结果")
            return
        
        for comparison_name in gene_comparisons:
            logger.info(f"处理比较: {comparison_name}")
            de_result = self.de_results[comparison_name]
            
            # 富集分析
            if run_enrichment:
                significant_genes = self.filter_significant_genes(de_result)
                
                if significant_genes is not None and len(significant_genes) > 0:
                    logger.info(f"对 {comparison_name} 进行富集分析")
                    
                    gene_list = significant_genes.index.tolist()
                    enrichment_result = self.run_enrichment_analysis(gene_list)
                    
                    if enrichment_result:
                        self.enrichment_results[comparison_name] = enrichment_result
                        logger.info(f"✅ 完成 {comparison_name} 富集分析")
                    else:
                        logger.warning(f"❌ {comparison_name} 富集分析失败")
                else:
                    logger.info(f"{comparison_name} 没有显著差异基因")
            
            # GSEA分析
            if run_gsea:
                logger.info(f"对 {comparison_name} 进行GSEA分析")
                
                gsea_result = self.run_gsea_analysis(de_result)
                
                if gsea_result:
                    self.gsea_results[comparison_name] = gsea_result
                    logger.info(f"✅ 完成 {comparison_name} GSEA分析")
                else:
                    logger.warning(f"❌ {comparison_name} GSEA分析失败")
    
    def save_results(self, output_dir='results'):
        """保存分析结果 - 改进文件命名规则"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 保存差异表达结果 - 使用简化命名
        logger.info("保存差异表达结果...")
        for comparison, result in self.de_results.items():
            try:
                # 提取数据类型
                if comparison.endswith('_gene'):
                    data_type = 'Gene'
                elif comparison.endswith('_mirna'):
                    data_type = 'miRNA'
                else:
                    data_type = 'Unknown'
                
                # 简化文件名：Gene_differential_expression_results.csv
                filepath = output_path / f"{data_type}_differential_expression_results.csv"
                result.to_csv(filepath)
                logger.info(f"保存: {filepath.name}")
            except Exception as e:
                logger.error(f"保存失败 {comparison}: {e}")
        
        # 保存富集分析结果 - 改进命名规则
        if self.enrichment_results:
            logger.info("保存富集分析结果...")
            for comparison, enrichment in self.enrichment_results.items():
                if enrichment and hasattr(enrichment, 'results'):
                    # 提取数据类型
                    if comparison.endswith('_gene'):
                        data_type = 'Gene'
                    elif comparison.endswith('_mirna'):
                        data_type = 'miRNA'
                    else:
                        data_type = 'Unknown'
                    
                    for gene_set, result_df in enrichment.results.items():
                        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
                            try:
                                # 清理基因集名称
                                safe_gene_set = str(gene_set).replace(' ', '_').replace('/', '_').replace('-', '_')
                                # 新命名格式：Gene_KEGG_2021_Human_enrichment.csv
                                filename = f"{data_type}_{safe_gene_set}_enrichment.csv"
                                filepath = output_path / filename
                                result_df.to_csv(filepath, index=False)
                                logger.info(f"保存富集结果: {filename}")
                            except Exception as e:
                                logger.error(f"保存富集结果失败: {e}")
        
        # 保存GSEA结果 - 改进命名规则
        if self.gsea_results:
            logger.info("保存GSEA结果...")
            for comparison, gsea in self.gsea_results.items():
                if gsea and hasattr(gsea, 'results'):
                    # 提取数据类型
                    if comparison.endswith('_gene'):
                        data_type = 'Gene'
                    elif comparison.endswith('_mirna'):
                        data_type = 'miRNA'
                    else:
                        data_type = 'Unknown'
                    
                    for gene_set, result_df in gsea.results.items():
                        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
                            try:
                                # 清理基因集名称
                                safe_gene_set = str(gene_set).replace(' ', '_').replace('/', '_').replace('-', '_')
                                # 新命名格式：Gene_KEGG_2021_Human_GSEA.csv
                                filename = f"{data_type}_{safe_gene_set}_GSEA.csv"
                                filepath = output_path / filename
                                result_df.to_csv(filepath, index=False)
                                logger.info(f"保存GSEA结果: {filename}")
                            except Exception as e:
                                logger.error(f"保存GSEA结果失败: {e}")
                
                # 保存ranked基因列表 - 简化命名
                if hasattr(gsea, 'ranked_genes'):
                    try:
                        # 提取数据类型
                        if comparison.endswith('_gene'):
                            data_type = 'Gene'
                        elif comparison.endswith('_mirna'):
                            data_type = 'miRNA'
                        else:
                            data_type = 'Unknown'
                        
                        # 简化命名：Gene_ranked_genes.rnk
                        filename = f"{data_type}_ranked_genes.rnk"
                        filepath = output_path / filename
                        gsea.save_ranked_genes(filepath)
                        logger.info(f"保存ranked基因: {filename}")
                    except Exception as e:
                        logger.error(f"保存ranked基因失败: {e}")
        
        logger.info(f"所有结果已保存到 {output_path.absolute()}")
    
    def cleanup(self):
        """清理临时文件"""
        try:
            import shutil
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
                logger.info(f"清理临时目录: {self.temp_dir}")
        except Exception as e:
            logger.warning(f"清理临时文件失败: {e}")


def create_sample_config():
    """创建示例配置文件"""
    config = {
        'sample_groups': {
            'Control': ['CT_', 'CT.'],       
            'Treatment': ['BAI_', 'BAI.']    
        },
        'analysis': {
            'min_count': 10,                 
            'min_samples': 2,                
            'padj_thresh': 0.05,             
            'log2fc_thresh': 1.0,            
            'normalization_method': 'deseq2' 
        },
        'deseq2': {
            'fit_type': 'parametric',        
            'test': 'Wald',                  
            'shrink_lfc': True,              
            'alpha': 0.05,                   
            'independent_filtering': True    
        },
        'gsea': {
            'permutation_num': 10000,
            'min_size': 15,
            'max_size': 500,
            'ranking_method': 'wald_stat',
            'random_seed': 42
        }
    }
    
    with open('config.yaml', 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    
    print("示例配置文件已保存为 config.yaml")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='修复版混合Python+R版RNA-seq差异表达分析')
    parser.add_argument('--gene-file', help='基因表达数据文件路径')
    parser.add_argument('--mirna-file', help='miRNA表达数据文件路径') 
    parser.add_argument('--output-dir', default='results', help='输出目录')
    parser.add_argument('--config', help='配置文件路径')
    parser.add_argument('--control', default='Control', help='对照组名称')
    parser.add_argument('--treatment', nargs='+', 
                       default=['Treatment'],
                       help='处理组名称')
    parser.add_argument('--skip-enrichment', action='store_true',
                       help='跳过GO/KEGG富集分析')
    parser.add_argument('--skip-gsea', action='store_true',
                       help='跳过GSEA分析')
    parser.add_argument('--enrichment-only', action='store_true',
                       help='仅进行富集分析')
    parser.add_argument('--padj-thresh', type=float, default=None,
                       help='校正p值阈值（默认0.05）')
    parser.add_argument('--log2fc-thresh', type=float, default=None,
                       help='log2 fold change阈值（默认1.0，条件: |log2FC| >= thresh）')
    parser.add_argument('--min-count', type=int, default=None,
                       help='最小表达量阈值（默认10）')
    parser.add_argument('--min-samples', type=int, default=None,
                       help='最小样本数阈值（默认2）')
    parser.add_argument('--create-config', action='store_true',
                       help='创建示例配置文件')
    # 添加快速模式参数
    parser.add_argument('--quick-mode', action='store_true',
                       help='快速模式：跳过GSEA，减少置换次数')
    parser.add_argument('--check-data', action='store_true',
                       help='仅检查数据格式和样本分组，不进行分析')
    parser.add_argument('--check-deps', action='store_true',
                       help='检查依赖包安装状态')
    
    # GSEA专用参数
    parser.add_argument('--gsea-ranking-method', 
                       choices=['wald_stat', 'log2fc', 'signed_pvalue', 'signal_to_noise'],
                       default=None,
                       help='GSEA排序方法（默认wald_stat）')
    parser.add_argument('--gsea-permutation-num', type=int, default=None,
                       help='GSEA置换次数（默认10000）')
    parser.add_argument('--gsea-random-seed', type=int, default=None,
                       help='GSEA随机种子（默认42）')
    
    # 在argparse部分添加新参数
    parser.add_argument('--run-ppi', action='store_true',
                    help='运行PPI网络分析')
    parser.add_argument('--run-survival', action='store_true',
                    help='运行生存分析')
    parser.add_argument('--run-integration', action='store_true',
                    help='运行整合分析（包含PPI和生存分析）')
    parser.add_argument('--cancer-type', default='BRCA',
                    help='癌症类型（用于生存分析）')
    parser.add_argument('--string-score', type=int, default=400,
                    help='STRING数据库置信度阈值')

    args = parser.parse_args()
    
    # 检查依赖
    if args.check_deps:
        print("检查依赖包状态:")
        print(f"  - R: {'可用' if HAS_R else '不可用'}")
        if HAS_R:
            print(f"    路径: {R_PATH}")
            
            # 创建临时分析器检查DESeq2
            print("  - DESeq2 (R包): 检查中...")
            try:
                temp_analyzer = HybridRNASeqAnalyzer()
                deseq2_available = temp_analyzer.check_deseq2_installation()
                print(f"  - DESeq2 (R包): {'可用' if deseq2_available else '不可用'}")
                temp_analyzer.cleanup()
            except Exception as e:
                print(f"  - DESeq2 (R包): 检查失败 ({e})")
        else:
            print(f"    错误: {R_ERROR}")
            print("  - DESeq2 (R包): 不可用 (R未安装)")
            
        print(f"  - gseapy: {'已安装' if HAS_GSEAPY else '未安装 (需要进行GO/KEGG分析)'}")
        if not HAS_GSEAPY:
            print("    安装命令: pip install gseapy")
        else:
            try:
                import gseapy as gp_test
                print(f"    版本: {gp_test.__version__}")
            except:
                print("    版本: 未知")
                
        print(f"  - sklearn: {'已安装' if HAS_SKLEARN else '未安装'}")
        print(f"  - mygene: {'已安装' if HAS_MYGENE else '未安装 (需要基因ID转换)'}")
        if not HAS_MYGENE:
            print("    安装命令: pip install mygene")
        
        # 🆕 检查可用基因集
        if HAS_GSEAPY:
            print("\n检查可用基因集...")
            try:
                temp_analyzer = HybridRNASeqAnalyzer()
                temp_analyzer.debug_available_gene_sets()
                temp_analyzer.cleanup()
            except Exception as e:
                print(f"检查基因集失败: {e}")
        
        if not HAS_R:
            print("\n安装R的步骤:")
            print("1. 下载并安装R: https://cran.r-project.org/")
            print("2. 在R中安装DESeq2:")
            print("   install.packages('BiocManager')")
            print("   BiocManager::install('DESeq2')")
            print("3. 确保R在系统PATH中")
        
        print(f"\n即使R不可用，程序仍可使用Python备用方法运行")
        print(f"   运行命令: python rnaseq_analyzer.py --gene-file your_data.txt")
        return
    
    # 创建配置文件
    if args.create_config:
        create_sample_config()
        return
    
    # 检查数据格式
    if args.check_data:
        logger.info("="*60)
        logger.info("数据格式检查模式")
        logger.info("="*60)
        
        if not args.gene_file and not args.mirna_file:
            logger.error("请提供基因表达数据或miRNA表达数据文件")
            return
        
        # 创建临时分析器
        temp_analyzer = HybridRNASeqAnalyzer(args.config)
        
        try:
            # 加载数据
            temp_analyzer.load_data(args.gene_file, args.mirna_file)
            
            # 检查基因数据
            if temp_analyzer.gene_data is not None:
                logger.info(f"✅ 基因数据加载成功: {temp_analyzer.gene_data.shape}")
                logger.info("基因数据样本分组测试:")
                gene_sample_info = temp_analyzer.create_sample_info(temp_analyzer.gene_data, 'gene')
                
                groups = gene_sample_info['treatment'].value_counts()
                if len(groups) >= 2 and 'Unknown' not in groups.index:
                    logger.info("✅ 基因数据样本分组成功")
                else:
                    logger.warning("⚠️ 基因数据样本分组可能有问题")
                    if 'Unknown' in groups.index:
                        unknown_samples = gene_sample_info[gene_sample_info['treatment'] == 'Unknown'].index.tolist()
                        logger.warning(f"未分组样本: {unknown_samples}")
            
            # 检查miRNA数据
            if temp_analyzer.mirna_data is not None:
                logger.info(f"✅ miRNA数据加载成功: {temp_analyzer.mirna_data.shape}")
                logger.info("miRNA数据样本分组测试:")
                mirna_sample_info = temp_analyzer.create_sample_info(temp_analyzer.mirna_data, 'mirna')
                
                groups = mirna_sample_info['treatment'].value_counts()
                if len(groups) >= 2 and 'Unknown' not in groups.index:
                    logger.info("✅ miRNA数据样本分组成功")
                else:
                    logger.warning("⚠️ miRNA数据样本分组可能有问题")
                    if 'Unknown' in groups.index:
                        unknown_samples = mirna_sample_info[mirna_sample_info['treatment'] == 'Unknown'].index.tolist()
                        logger.warning(f"未分组样本: {unknown_samples}")
            
            logger.info("="*60)
            logger.info("数据检查完成")
            logger.info("="*60)
            
            if temp_analyzer.gene_data is None and temp_analyzer.mirna_data is None:
                logger.error("❌ 没有成功加载任何数据")
            else:
                logger.info("✅ 数据格式检查通过，可以开始分析")
                logger.info("运行完整分析命令:")
                cmd_parts = [f"python {sys.argv[0]}"]
                if args.gene_file:
                    cmd_parts.append(f"--gene-file {args.gene_file}")
                if args.mirna_file:
                    cmd_parts.append(f"--mirna-file {args.mirna_file}")
                if args.config:
                    cmd_parts.append(f"--config {args.config}")
                cmd_parts.extend(["--control Control", "--treatment Treatment", "--output-dir results"])
                logger.info(" ".join(cmd_parts))
        
        finally:
            temp_analyzer.cleanup()
        return
    
    # 显示分析方法信息  
    if HAS_R:
        logger.info("将使用R+DESeq2进行基因差异表达分析")
    else:
        logger.warning("R不可用，将使用Python备用方法")
        logger.warning("建议安装R和DESeq2以获得更准确的结果")
    
    # 创建输出目录
    Path(args.output_dir).mkdir(exist_ok=True)
    
    # 创建分析器
    analyzer = HybridRNASeqAnalyzer(args.config)
    
    # 快速模式设置
    if args.quick_mode:
        logger.info("🚀 启用快速模式")
        analyzer.config['gsea']['permutation_num'] = 1000  # 大幅减少置换次数
        args.skip_gsea = True  # 跳过GSEA
        logger.info("快速模式设置: 跳过GSEA，减少计算时间")
    
    try:
        # 调整配置参数
        if args.min_count is not None:
            analyzer.config['analysis']['min_count'] = args.min_count
        if args.min_samples is not None:
            analyzer.config['analysis']['min_samples'] = args.min_samples
        if args.padj_thresh is not None:
            analyzer.config['analysis']['padj_thresh'] = args.padj_thresh
        if args.log2fc_thresh is not None:
            analyzer.config['analysis']['log2fc_thresh'] = args.log2fc_thresh
            
        # GSEA专用参数
        if args.gsea_ranking_method is not None:
            analyzer.config['gsea']['ranking_method'] = args.gsea_ranking_method
        if args.gsea_permutation_num is not None:
            analyzer.config['gsea']['permutation_num'] = args.gsea_permutation_num
        if args.gsea_random_seed is not None:
            analyzer.config['gsea']['random_seed'] = args.gsea_random_seed
        
        logger.info(f"分析参数: min_count={analyzer.config['analysis']['min_count']}, "
                   f"min_samples={analyzer.config['analysis']['min_samples']}, "
                   f"padj<{analyzer.config['analysis']['padj_thresh']}, "
                   f"|log2FC|>={analyzer.config['analysis']['log2fc_thresh']}")
        
        logger.info(f"GSEA参数: ranking_method={analyzer.config['gsea']['ranking_method']}, "
                   f"permutation_num={analyzer.config['gsea']['permutation_num']}, "
                   f"random_seed={analyzer.config['gsea']['random_seed']}")
        
        # 预估分析时间
        if not args.skip_gsea and not args.enrichment_only:
            estimated_time = "5-15分钟"
            if analyzer.config['gsea']['permutation_num'] > 5000:
                estimated_time = "15-30分钟"
            logger.info(f"⏰ 预估总分析时间: {estimated_time} (包含GSEA)")
        else:
            logger.info("⏰ 预估总分析时间: 2-5分钟 (不含GSEA)")
        
        if args.enrichment_only:
            logger.info("仅进行富集分析模式")
            # 加载已有结果
            result_files = list(Path(args.output_dir).glob('*_differential_expression_results.csv'))
            if not result_files:
                logger.error("没有找到差异表达结果文件")
                sys.exit(1)
            
            for result_file in result_files:
                comparison_name = result_file.stem.replace('_differential_expression_results', '')
                result_df = pd.read_csv(result_file, index_col=0)
                analyzer.de_results[comparison_name] = result_df
                logger.info(f"加载结果: {comparison_name}")
            
            # 运行富集分析
            analyzer.run_all_enrichment_analysis(
                run_enrichment=not args.skip_enrichment,
                run_gsea=not args.skip_gsea
            )
            
            # 保存结果
            analyzer.save_results(args.output_dir)
        else:
            # 正常分析流程
            # 加载数据
            analyzer.load_data(args.gene_file, args.mirna_file)
            
            if analyzer.gene_data is None and analyzer.mirna_data is None:
                logger.error("请提供基因表达数据或miRNA表达数据")
                sys.exit(1)
            
            logger.info("="*60)
            logger.info("第一阶段: 差异表达分析")
            logger.info("="*60)
            
            # 差异表达分析
            results = {}
            all_comparisons = []
            
            for treatment in args.treatment:
                if analyzer.gene_data is not None:
                    all_comparisons.append((treatment, args.control, 'gene'))
                if analyzer.mirna_data is not None:
                    all_comparisons.append((treatment, args.control, 'mirna'))
            
            logger.info(f"计划进行 {len(all_comparisons)} 个比较分析")
            
            # 执行分析
            for i, (treatment, control, data_type) in enumerate(all_comparisons, 1):
                logger.info(f"进度 {i}/{len(all_comparisons)}: {treatment} vs {control} ({data_type})")
                
                result = analyzer.analyze_comparison(control, treatment, data_type)
                if result:
                    comparison_key = f"{treatment}_vs_{control}_{data_type}"
                    results[comparison_key] = result
                    logger.info(f"完成 {comparison_key}")
                else:
                    logger.warning(f"失败 {treatment}_vs_{control}_{data_type}")
            
            logger.info("="*60)
            logger.info("第二阶段: 生成可视化图表")
            logger.info("="*60)
            
            # 生成图表
            for comparison_key, result in results.items():
                parts = comparison_key.split('_')
                treatment, control, data_type = parts[0], parts[2], parts[3]
                
                # 火山图
                logger.info(f"生成火山图: {comparison_key}")
                fig = analyzer.plot_volcano(result['de_result'])  # 移除标题参数
                if fig:
                    # Modified line for Volcano Plot
                    filename = "mrna_volcano_plot.png" if data_type == 'gene' else "mirna_volcano_plot.png"
                    fig.savefig(f"{args.output_dir}/{filename}", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                
                # MA图
                logger.info(f"生成MA图: {comparison_key}")
                fig = analyzer.plot_ma(result['de_result'])  # 移除标题参数
                if fig:
                    # Modified line for MA Plot
                    filename = "mrna_ma_plot.png" if data_type == 'gene' else "mirna_ma_plot.png"
                    fig.savefig(f"{args.output_dir}/{filename}", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                
                # 热图
                logger.info(f"生成热图: {comparison_key}")
                fig = analyzer.plot_heatmap(result['de_result'], data_type)  # 移除标题参数
                if fig:
                    # Modified line for Heatmap
                    filename = "mrna_heatmap.png" if data_type == 'gene' else "mirna_heatmap.png"
                    fig.savefig(f"{args.output_dir}/{filename}", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
            
            logger.info("="*60)
            logger.info("第三阶段: 功能富集分析")
            logger.info("="*60)
            
            # 统计结果
            gene_comparisons = len([k for k in results.keys() if k.endswith('_gene')])
            mirna_comparisons = len([k for k in results.keys() if k.endswith('_mirna')])
            
            total_significant = 0
            for result in results.values():
                significant = analyzer.filter_significant_genes(result['de_result'])
                if significant is not None:
                    total_significant += len(significant)
            
            logger.info(f"分析总结: 基因比较={gene_comparisons}, miRNA比较={mirna_comparisons}, "
                       f"总显著基因/miRNA={total_significant}")
            
            # 功能分析
            if not args.skip_enrichment or not args.skip_gsea:
                analyzer.run_all_enrichment_analysis(
                    run_enrichment=not args.skip_enrichment,
                    run_gsea=not args.skip_gsea
                )
                
                # 绘制富集分析图
                if analyzer.enrichment_results:
                    logger.info("生成富集分析图表...")
                    for comparison, enrichment in analyzer.enrichment_results.items():
                        if enrichment and hasattr(enrichment, 'results'):
                            for gene_set in enrichment.results.keys():
                                fig = analyzer.plot_enrichment(enrichment, gene_set)
                                if fig:
                                    safe_gene_set = gene_set.replace(' ', '_').replace('/', '_')
                                    fig.savefig(f"{args.output_dir}/{comparison}_{safe_gene_set}_enrichment.png", 
                                               dpi=300, bbox_inches='tight')
                                    plt.close()
            
            logger.info("="*60)
            logger.info("第四阶段: 保存所有结果")
            logger.info("="*60)
            
            # 保存结果
            analyzer.save_results(args.output_dir)
        
        logger.info("🎉 所有分析完成！")
        logger.info(f"结果保存在: {args.output_dir}/")
        
        if HAS_R:
            logger.info("✅ 使用了R+DESeq2，结果与纯R分析一致！")
        else:
            logger.warning("⚠️ 使用了Python备用方法，建议安装R+DESeq2获得更准确结果")
            
        # 显示GSEA一致性信息
        if not args.skip_gsea and analyzer.gsea_results:
            logger.info("="*60)
            logger.info("GSEA分析一致性信息")
            logger.info("="*60)
            for comparison, gsea_result in analyzer.gsea_results.items():
                if hasattr(gsea_result, 'parameters'):
                    params = gsea_result.parameters
                    logger.info(f"{comparison}:")
                    logger.info(f"  排序方法: {params.get('ranking_method', 'unknown')}")
                    logger.info(f"  置换次数: {params.get('permutation_num', 'unknown')}")
                    logger.info(f"  随机种子: {params.get('random_seed', 'unknown')}")
                    logger.info(f"  基因转换率: {params.get('conversion_rate', 0)*100:.1f}%")
                    logger.info(f"  🎯 使用基因集: {params.get('gene_sets', [])}")
                    logger.info(f"  📝 排序文件: {comparison}_ranked_genes.rnk")
            
            logger.info("📝 提示: 使用生成的.rnk文件与其他工具结果比较可验证一致性")
            
            # 检查是否成功使用WikiPathways
            wikipathway_found = False
            for comparison, gsea_result in analyzer.gsea_results.items():
                if hasattr(gsea_result, 'parameters'):
                    gene_sets = gsea_result.parameters.get('gene_sets', [])
                    wiki_sets = [gs for gs in gene_sets if 'WikiPathway' in gs]
                    if wiki_sets:
                        wikipathway_found = True
                        logger.info(f"🎯 成功使用WikiPathways数据库: {wiki_sets}")
                        break
            
            if not wikipathway_found:
                logger.warning("⚠️ 未能使用WikiPathways数据库")
                logger.info("可能原因: WikiPathways在当前gseapy版本中名称已更改")
                logger.info("建议运行: python rnaseq_analyzer.py --check-deps 查看可用基因集")
        # 在保存结果之后，添加生物标志物分析
        if HAS_BIOMARKER_MODULES and analyzer.de_results:
            
            # 选择要分析的差异表达结果（通常是基因数据）
            gene_comparisons = [k for k in analyzer.de_results.keys() if k.endswith('_gene')]
            
            if gene_comparisons and (args.run_ppi or args.run_survival or args.run_integration):
                # 使用第一个基因比较结果
                comparison = gene_comparisons[0]
                de_file = f"{args.output_dir}/{comparison}_differential_expression_results.csv"
                
                logger.info("="*60)
                logger.info("第五阶段: 生物标志物分析")
                logger.info("="*60)
                
                if args.run_integration:
                    # 运行完整的整合分析
                    logger.info("运行整合分析（PPI + 生存分析）...")
                    integration_results = run_integrated_analysis(
                        de_file=de_file,
                        cancer_type=args.cancer_type,
                        output_dir=f"{args.output_dir}/biomarker_integration",
                        padj_threshold=args.padj_thresh or 0.05,
                        log2fc_threshold=args.log2fc_thresh or 1.0,
                        string_score=args.string_score,
                        top_biomarkers=20
                    )
                    logger.info("✅ 整合分析完成")
                    
                else:
                    # 单独运行各个分析
                    if args.run_ppi:
                        logger.info("运行PPI网络分析...")
                        ppi_results = analyze_de_genes_ppi(
                            de_results_file=de_file,
                            output_dir=f"{args.output_dir}/ppi_analysis",
                            padj_threshold=args.padj_thresh or 0.05,
                            log2fc_threshold=args.log2fc_thresh or 1.0,
                            string_score=args.string_score
                        )
                        logger.info("✅ PPI分析完成")
                    
                    if args.run_survival:
                        logger.info("运行生存分析...")
                        # 获取显著基因列表
                        de_df = pd.read_csv(de_file, index_col=0)
                        significant = de_df[
                            (de_df['padj'] < (args.padj_thresh or 0.05)) & 
                            (abs(de_df['log2FoldChange']) >= (args.log2fc_thresh or 1.0))
                        ]
                        gene_list = significant.head(50).index.tolist()
                        
                        survival_results = analyze_gene_survival(
                            gene_list=gene_list,
                            cancer_type=args.cancer_type,
                            output_dir=f"{args.output_dir}/survival_analysis"
                        )
                        logger.info("✅ 生存分析完成")
    finally:
        # 清理临时文件
        analyzer.cleanup()


if __name__ == "__main__":
    main()
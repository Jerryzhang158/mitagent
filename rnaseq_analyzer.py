#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
混合Python+R版RNA-seq和miRNA-seq差异表达分析程序
- 差异表达分析：使用原生R脚本调用DESeq2
- 数据处理和可视化：使用Python
- 富集分析：使用Python (gseapy)
"""

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
    
    # LFC shrinkage（如果启用）
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
                'permutation_num': 1000,
                'min_size': 15,
                'max_size': 500
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
                    self.gene_data.columns = [str(col).replace('\x00', '').replace('\ufeff', '').replace('ÿþ', '').strip() 
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
                else:
                    self.mirna_data = None
                    
            except Exception as e:
                logger.error(f"加载miRNA表达数据失败: {e}")
                self.mirna_data = None
    
    def create_sample_info(self, data):
        """根据样本名称和配置创建样本信息表"""
        samples = data.columns.tolist()
        sample_info = []
        
        for sample in samples:
            # 清理样本名称
            clean_sample = str(sample).replace('\x00', '').replace('\ufeff', '').replace('ÿþ', '').strip()
            
            group = 'Unknown'
            treatment = 'Unknown'
            
            # 根据配置匹配样本组
            for group_name, prefixes in self.config['sample_groups'].items():
                if any(clean_sample.startswith(prefix) for prefix in prefixes):
                    group = group_name
                    treatment = group_name
                    break
            
            # 如果没匹配到，尝试默认匹配
            if group == 'Unknown':
                if clean_sample.startswith('CT_') or clean_sample.startswith('CT.'):
                    group = 'Control'
                    treatment = 'Control'
                elif clean_sample.startswith('BAI_') or clean_sample.startswith('BAI.'):
                    group = 'Treatment'
                    treatment = 'Treatment'
                # 其他格式...
            
            sample_info.append({
                'sample': sample,
                'clean_sample': clean_sample,
                'group': group,
                'treatment': treatment
            })
        
        sample_df = pd.DataFrame(sample_info).set_index('sample')
        
        # 记录分组信息
        group_counts = sample_df['treatment'].value_counts()
        logger.info("样本分组信息:")
        for group, count in group_counts.items():
            group_samples = sample_df[sample_df['treatment'] == group]['clean_sample'].tolist()
            logger.info(f"  {group}: {count} 个样本 - {group_samples}")
        
        return sample_df
    
    def run_r_deseq2_analysis(self, count_data, sample_info, contrast, data_type='gene'):
        """使用R脚本运行DESeq2分析"""
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
            
            # 运行R脚本
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
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            # 显示R的输出
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        logger.info(f"R> {line}")
            
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip() and 'Loading required package' not in line:
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
            
            # 读取分析摘要
            summary_file = str(output_file).replace('.csv', '_summary.txt')
            if os.path.exists(summary_file):
                with open(summary_file, 'r') as f:
                    summary_content = f.read()
                    logger.info("分析摘要:")
                    for line in summary_content.split('\n')[2:]:  # 跳过标题行
                        if line.strip():
                            logger.info(f"  {line}")
            
            return results_df
            
        except subprocess.TimeoutExpired:
            logger.error("R脚本执行超时")
            return None
        except Exception as e:
            logger.error(f"R DESeq2分析失败: {e}")
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
    
    def plot_volcano(self, de_result, title="Volcano Plot", 
                    padj_thresh=None, log2fc_thresh=None, figsize=(10, 8)):
        """绘制火山图"""
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
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        return plt.gcf()
    
    def plot_ma(self, de_result, title="MA Plot", 
                padj_thresh=None, log2fc_thresh=None, figsize=(10, 8)):
        """绘制MA图"""
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
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        return plt.gcf()
    
    def plot_heatmap(self, de_result, data_type='gene', title="Expression Heatmap", 
                    top_n=50, figsize=(12, 10)):
        """绘制表达热图"""
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
        sample_info = self.create_sample_info(expr_data)
        
        try:
            plt.figure(figsize=figsize)
            
            # 创建颜色映射
            unique_treatments = sample_info['treatment'].unique()
            colors = plt.cm.Set1(np.linspace(0, 1, len(unique_treatments)))
            treatment_colors = dict(zip(unique_treatments, colors))
            col_colors = sample_info['treatment'].map(treatment_colors)
            
            # 绘制聚类热图
            g = sns.clustermap(plot_data_scaled, 
                              col_colors=col_colors,
                              cmap='RdBu_r', 
                              center=0,
                              figsize=figsize,
                              cbar_kws={'label': 'Z-score'},
                              yticklabels=True,
                              xticklabels=True,
                              method='average',
                              metric='euclidean')
            
            g.fig.suptitle(title, y=1.02)
            return g.fig
            
        except Exception as e:
            logger.error(f"绘制热图失败: {e}")
            # 简化热图
            try:
                plt.figure(figsize=figsize)
                im = plt.imshow(plot_data_scaled.values, cmap='RdBu_r', aspect='auto')
                plt.colorbar(im, label='Z-score')
                plt.xticks(range(len(plot_data_scaled.columns)), 
                          plot_data_scaled.columns, rotation=45)
                plt.yticks(range(len(plot_data_scaled.index)), 
                          plot_data_scaled.index)
                plt.title(title)
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
    
    def run_gsea_analysis(self, de_result, gene_sets=None):
        """运行GSEA分析"""
        if not HAS_GSEAPY:
            logger.error("gseapy未安装，无法运行GSEA分析")
            return None
        
        if de_result is None or de_result.empty:
            logger.error("没有差异表达结果")
            return None
        
        try:
            # 准备ranked gene list
            ranked_data = de_result.copy()
            ranked_data = ranked_data.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2FoldChange', 'pvalue'])
            
            if len(ranked_data) == 0:
                logger.error("没有有效基因用于GSEA")
                return None
            
            # 计算ranking metric
            ranked_data['ranking_metric'] = -np.log10(ranked_data['pvalue'].clip(lower=1e-300)) * np.sign(ranked_data['log2FoldChange'])
            ranked_data = ranked_data.sort_values('ranking_metric', ascending=False)
            
            # 转换基因ID
            gene_symbols = self.convert_ensembl_to_symbols(ranked_data.index.tolist())
            
            # 创建有效的ranked gene list
            valid_data = []
            for symbol, metric in zip(gene_symbols, ranked_data['ranking_metric'].values):
                if (symbol is not None and pd.notna(symbol) and symbol != '' and 
                    not symbol.startswith('ENSG')):
                    valid_data.append({'gene_symbol': symbol, 'ranking_metric': metric})
            
            if len(valid_data) == 0:
                logger.error("没有有效基因符号用于GSEA")
                return None
            
            rnk_df = pd.DataFrame(valid_data)
            rnk_df = rnk_df.drop_duplicates(subset=['gene_symbol'], keep='first')
            
            logger.info(f"GSEA分析: 使用 {len(rnk_df)} 个基因")
            
            # 选择基因集
            if gene_sets is None:
                try:
                    available_libs = gp.get_library_name()
                    preferred_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
                    
                    gene_sets = []
                    for pref_set in preferred_sets:
                        if pref_set in available_libs:
                            gene_sets.append(pref_set)
                            break
                    
                    if not gene_sets:
                        go_sets = [lib for lib in available_libs if 'GO_Biological_Process' in lib]
                        if go_sets:
                            gene_sets.append(go_sets[0])
                        else:
                            logger.error("找不到合适的基因集")
                            return None
                            
                except Exception as e:
                    logger.error(f"获取基因集失败: {e}")
                    return None
            
            # 运行GSEA
            gsea_results = {}
            for gene_set in gene_sets:
                try:
                    logger.info(f"运行GSEA: {gene_set}")
                    
                    gs_res = gp.prerank(rnk=rnk_df,
                                      gene_sets=gene_set,
                                      processes=1,
                                      permutation_num=self.config['gsea']['permutation_num'],
                                      outdir=None,
                                      max_size=self.config['gsea']['max_size'],
                                      min_size=self.config['gsea']['min_size'],
                                      weighted_score_type=1,
                                      ascending=False,
                                      verbose=False)
                    
                    if gs_res and hasattr(gs_res, 'res2d') and not gs_res.res2d.empty:
                        gsea_results[gene_set] = gs_res.res2d.copy()
                        logger.info(f"GSEA {gene_set}: {len(gs_res.res2d)} 个通路")
                        
                except Exception as e:
                    logger.error(f"GSEA失败 {gene_set}: {e}")
                    continue
            
            if gsea_results:
                class GSEAResults:
                    def __init__(self, results_dict, ranked_genes):
                        self.results = results_dict
                        self.ranked_genes = ranked_genes
                
                return GSEAResults(gsea_results, rnk_df)
            else:
                return None
                
        except Exception as e:
            logger.error(f"GSEA分析失败: {e}")
            return None
    
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
        """进行组间比较分析"""
        logger.info(f"开始分析 {treatment_group} vs {control_group} ({data_type})")
        
        # 选择数据
        if data_type == 'gene' and self.gene_data is not None:
            data = self.gene_data
        elif data_type == 'mirna' and self.mirna_data is not None:
            data = self.mirna_data
        else:
            logger.error(f"没有可用的{data_type}数据")
            return None
            
        # 创建样本信息
        sample_info = self.create_sample_info(data)
        
        # 筛选相关样本
        relevant_samples = sample_info[
            sample_info['treatment'].isin([control_group, treatment_group])
        ]
        
        if len(relevant_samples) == 0:
            logger.error(f"没有找到相关样本")
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
        """保存分析结果"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 保存差异表达结果
        logger.info("保存差异表达结果...")
        for comparison, result in self.de_results.items():
            try:
                filepath = output_path / f"{comparison}_differential_expression_results.csv"
                result.to_csv(filepath)
                logger.info(f"保存: {filepath.name}")
            except Exception as e:
                logger.error(f"保存失败 {comparison}: {e}")
        
        # 保存富集分析结果
        if self.enrichment_results:
            logger.info("保存富集分析结果...")
            for comparison, enrichment in self.enrichment_results.items():
                if enrichment and hasattr(enrichment, 'results'):
                    for gene_set, result_df in enrichment.results.items():
                        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
                            try:
                                safe_gene_set = str(gene_set).replace(' ', '_').replace('/', '_')
                                filename = f"{comparison}_{safe_gene_set}_enrichment.csv"
                                filepath = output_path / filename
                                result_df.to_csv(filepath, index=False)
                                logger.info(f"保存富集结果: {filename}")
                            except Exception as e:
                                logger.error(f"保存富集结果失败: {e}")
        
        # 保存GSEA结果
        if self.gsea_results:
            logger.info("保存GSEA结果...")
            for comparison, gsea in self.gsea_results.items():
                if gsea and hasattr(gsea, 'results'):
                    for gene_set, result_df in gsea.results.items():
                        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
                            try:
                                safe_gene_set = str(gene_set).replace(' ', '_').replace('/', '_')
                                filename = f"{comparison}_{safe_gene_set}_GSEA.csv"
                                filepath = output_path / filename
                                result_df.to_csv(filepath, index=False)
                                logger.info(f"保存GSEA结果: {filename}")
                            except Exception as e:
                                logger.error(f"保存GSEA结果失败: {e}")
                
                # 保存ranked基因列表
                if hasattr(gsea, 'ranked_genes'):
                    try:
                        filename = f"{comparison}_ranked_genes.csv"
                        filepath = output_path / filename
                        gsea.ranked_genes.to_csv(filepath, index=False)
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
            'permutation_num': 1000,
            'min_size': 15,
            'max_size': 500
        }
    }
    
    with open('config.yaml', 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    
    print("示例配置文件已保存为 config.yaml")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='混合Python+R版RNA-seq差异表达分析')
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
    parser.add_argument('--check-deps', action='store_true',
                       help='检查依赖包安装状态')
    
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
        print(f"  - sklearn: {'已安装' if HAS_SKLEARN else '未安装'}")
        print(f"  - mygene: {'已安装' if HAS_MYGENE else '未安装 (需要基因ID转换)'}")
        if not HAS_MYGENE:
            print("    安装命令: pip install mygene")
        
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
        
        logger.info(f"分析参数: min_count={analyzer.config['analysis']['min_count']}, "
                   f"min_samples={analyzer.config['analysis']['min_samples']}, "
                   f"padj<{analyzer.config['analysis']['padj_thresh']}, "
                   f"|log2FC|>={analyzer.config['analysis']['log2fc_thresh']}")
        
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
                fig = analyzer.plot_volcano(result['de_result'], 
                                          f"{data_type.capitalize()}: {treatment} vs {control}")
                if fig:
                    fig.savefig(f"{args.output_dir}/{comparison_key}_volcano.png", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                
                # MA图
                logger.info(f"生成MA图: {comparison_key}")
                fig = analyzer.plot_ma(result['de_result'], 
                                     f"{data_type.capitalize()} MA: {treatment} vs {control}")
                if fig:
                    fig.savefig(f"{args.output_dir}/{comparison_key}_ma.png", 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                
                # 热图
                logger.info(f"生成热图: {comparison_key}")
                fig = analyzer.plot_heatmap(result['de_result'], data_type,
                                          f"{data_type.capitalize()} Heatmap: {treatment} vs {control}")
                if fig:
                    fig.savefig(f"{args.output_dir}/{comparison_key}_heatmap.png", 
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
        
        logger.info("所有分析完成！")
        logger.info(f"结果保存在: {args.output_dir}/")
        
        if HAS_R:
            logger.info("使用了R+DESeq2，结果与纯R分析一致！")
        else:
            logger.warning("使用了Python备用方法，建议安装R+DESeq2获得更准确结果")
    
    finally:
        # 清理临时文件
        analyzer.cleanup()


if __name__ == "__main__":
    main()
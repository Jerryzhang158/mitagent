"""
通用工具函数 - miRNA Research Pipeline
包含日志、文件处理、数据验证等通用功能
"""

import os
import pandas as pd
import csv
import requests
from datetime import datetime
from typing import Optional, List, Dict, Any, Tuple

class Logger:
    """日志管理器"""
    
    def __init__(self, log_file: str, console_output: bool = True):
        self.log_file = log_file
        self.console_output = console_output
        
        # 确保日志目录存在
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    def log(self, message: str, level: str = "INFO"):
        """记录日志"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        formatted_message = f"[{timestamp}] {level}: {message}"
        
        # 写入文件
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(formatted_message + "\n")
        
        # 控制台输出
        if self.console_output:
            print(formatted_message)
    
    def info(self, message: str):
        """信息日志"""
        self.log(message, "INFO")
    
    def warning(self, message: str):
        """警告日志"""
        self.log(message, "WARNING")
    
    def error(self, message: str):
        """错误日志"""
        self.log(message, "ERROR")
    
    def debug(self, message: str):
        """调试日志"""
        self.log(message, "DEBUG")

class FileUtils:
    """文件处理工具"""
    
    @staticmethod
    def smart_read_csv(file_path: str, encoding: str = 'utf-8') -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        智能读取CSV文件，自动检测分隔符和编码
        
        Returns:
            tuple: (DataFrame, metadata_dict)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        metadata = {
            'file_size': os.path.getsize(file_path),
            'delimiter': None,
            'encoding': encoding,
            'columns': [],
            'shape': (0, 0)
        }
        
        try:
            # 读取文件头部样本
            with open(file_path, 'r', encoding=encoding, errors='ignore') as f:
                first_line = f.readline().strip()
                sample = f.read(2048)
            
            # 检测分隔符
            if '\t' in first_line and first_line.count('\t') > first_line.count(','):
                delimiter = '\t'
            elif ',' in first_line:
                delimiter = ','
            else:
                # 使用csv.Sniffer自动检测
                sniffer = csv.Sniffer()
                delimiter = sniffer.sniff(sample).delimiter
            
            metadata['delimiter'] = delimiter
            
            # 读取数据
            df = pd.read_csv(file_path, delimiter=delimiter, encoding=encoding, low_memory=False)
            
            metadata['columns'] = list(df.columns)
            metadata['shape'] = df.shape
            
            return df, metadata
            
        except Exception as e:
            raise Exception(f"Error reading CSV file {file_path}: {e}")
    
    @staticmethod
    def load_gene_list(file_path: str) -> Optional[set]:
        """加载基因列表文件"""
        if not file_path or not os.path.exists(file_path):
            return None
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                genes = set([line.strip() for line in f.readlines() if line.strip()])
            return genes
        except Exception as e:
            raise Exception(f"Error loading gene list from {file_path}: {e}")
    
    @staticmethod
    def save_excel_with_summary(df: pd.DataFrame, file_path: str, 
                               summary_data: Dict[str, any] = None,
                               validation_stats: Dict[str, int] = None):
        """保存DataFrame到Excel，包含摘要信息"""
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
            # 主数据
            df.to_excel(writer, sheet_name='MTI_Results', index=False)
            
            # 摘要数据 - 修复标量值问题
            if summary_data:
                # 将标量值转换为列表，以便创建DataFrame
                summary_items = []
                for key, value in summary_data.items():
                    summary_items.append({'Metric': key, 'Value': value})
                
                if summary_items:  # 确保有数据
                    summary_df = pd.DataFrame(summary_items)
                    summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # 验证统计
            if validation_stats:
                validation_df = pd.DataFrame(list(validation_stats.items()), 
                                           columns=['Validation_Type', 'Count'])
                validation_df.to_excel(writer, sheet_name='Validation_Summary', index=False)

class DataValidator:
    """数据验证工具"""
    
    @staticmethod
    def validate_database_columns(df: pd.DataFrame, required_columns: List[str], 
                                 database_name: str) -> bool:
        """验证数据库文件是否包含必需的列"""
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing required columns in {database_name}: {missing_columns}")
        
        return True
    
    @staticmethod
    def validate_mirtarbase_format(df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """验证miRTarBase文件格式"""
        required_columns = ['miRNA', 'Target Gene']
        missing_cols = [col for col in required_columns if col not in df.columns]
        
        return len(missing_cols) == 0, missing_cols
    
    @staticmethod
    def apply_column_mapping(df: pd.DataFrame, column_mapping: Dict[str, str]) -> pd.DataFrame:
        """应用列名映射"""
        df_copy = df.copy()
        
        for old_col, new_col in column_mapping.items():
            if old_col in df_copy.columns:
                df_copy.rename(columns={old_col: new_col}, inplace=True)
        
        return df_copy

class TextProcessor:
    """文本处理工具"""
    
    @staticmethod
    def extract_mirna_core(mirna_name: str) -> Optional[str]:
        """提取miRNA核心名称用于搜索"""
        if not mirna_name:
            return None
        
        import re
        
        # 提取miR-XXX部分
        match = re.search(r"hsa-(miR-\d+[a-z]*)", mirna_name)
        if match:
            return match.group(1)
        
        # 如果没有hsa-前缀
        if mirna_name.startswith('hsa-'):
            return mirna_name[4:]
        
        return mirna_name
    
    @staticmethod
    def safe_get_text(data_dict: Dict, key: str, max_length: int = 500) -> str:
        """安全获取文本数据，限制长度"""
        try:
            value = data_dict.get(key, '')
            if isinstance(value, str):
                return value[:max_length] + "..." if len(value) > max_length else value
            elif isinstance(value, (list, tuple)):
                return "; ".join(str(v) for v in value)[:max_length]
            else:
                return str(value)[:max_length] if value else ''
        except Exception:
            return ''
    
    @staticmethod
    def safe_get_list(data_dict: Dict, key: str, max_items: int = 5, 
                     max_length: int = 300) -> str:
        """安全获取列表数据，转换为分号分隔的字符串"""
        try:
            value = data_dict.get(key, [])
            if isinstance(value, (list, tuple)):
                limited_items = value[:max_items]
                result = "; ".join(str(item) for item in limited_items if item)
                return result[:max_length] + "..." if len(result) > max_length else result
            elif isinstance(value, str):
                return value[:max_length] + "..." if len(value) > max_length else value
            else:
                return str(value)[:max_length] if value else ''
        except Exception:
            return ''

class ResultsIntegrator:
    """结果整合工具"""
    
    @staticmethod
    def merge_mti_results(results_list: List[Dict], unique_key_func=None) -> List[Dict]:
        if not unique_key_func:
            unique_key_func = lambda x: (x['Gene'], x['miRNA'])
        
        unique_mtis = {}
        
        for mti in results_list:
            key = unique_key_func(mti)
            
            if key not in unique_mtis:
                unique_mtis[key] = mti.copy()
            else:
                existing = unique_mtis[key]
                
                # 更新方向信息
                if existing['Direction'] != mti['Direction']:
                    existing['Direction'] = 'Bidirectional'
                
                # 合并Database_Sources（追加而不是覆盖）
                existing_sources = set(existing.get('Database_Sources', '').split(','))
                new_sources = set(mti.get('Database_Sources', '').split(','))
                merged_sources = existing_sources | new_sources
                merged_sources.discard('')
                existing['Database_Sources'] = ','.join(sorted(merged_sources))
                
                # 数据库数量用合并后的source数重新计算
                existing['Databases'] = len(merged_sources)
                
                # Priority取更高的
                priority_order = {'High': 3, 'Medium': 2, 'Low': 1}
                if priority_order.get(mti.get('Priority', 'Low'), 1) > priority_order.get(existing.get('Priority', 'Low'), 1):
                    existing['Priority'] = mti['Priority']
                
                # miRTarBase验证信息：有就保留
                if mti.get('miRTarBase_Validation'):
                    existing['miRTarBase_Validation'] = True
                    existing['miRTarBase_Support_Type'] = mti.get('miRTarBase_Support_Type')
                    existing['Validation_Strength'] = mti.get('Validation_Strength')
        
        return list(unique_mtis.values())
    
    @staticmethod
    def create_summary_stats(mti_results: List[Dict]) -> Dict[str, int]:
        """创建MTI结果的统计摘要"""
        if not mti_results:
            return {}
        
        stats = {
            'Total MTIs': len(mti_results),
            'In 3 databases': sum(1 for m in mti_results if m.get('Databases') == 3),
            'In 2 databases': sum(1 for m in mti_results if m.get('Databases') == 2),
            'In 1 database': sum(1 for m in mti_results if m.get('Databases') == 1),
            'Gene→miRNA': sum(1 for m in mti_results if m.get('Direction') == 'Gene→miRNA'),
            'miRNA→Gene': sum(1 for m in mti_results if m.get('Direction') == 'miRNA→Gene'),
            'Bidirectional': sum(1 for m in mti_results if m.get('Direction') == 'Bidirectional'),
            'With Experimental Validation': sum(1 for m in mti_results if m.get('miRTarBase_Validation', False))
        }
        
        return stats

def create_timestamp() -> str:
    """创建时间戳字符串"""
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def check_ollama_connection(base_url: str = "http://localhost:11434") -> Tuple[bool, List[str]]:
    """检查Ollama连接和可用模型"""
    import requests
    
    try:
        response = requests.get(f"{base_url}/api/tags", timeout=5)
        if response.status_code == 200:
            models = response.json().get('models', [])
            model_names = [m['name'] for m in models]
            return True, model_names
        else:
            return False, []
    except Exception:
        return False, []
"""
miRTarBase验证重新识别工具
用于更新现有结果CSV文件中的miRTarBase验证信息
"""

import pandas as pd
import argparse
import os
from typing import Dict, Optional


class MiRTarBaseValidator:
    """miRTarBase验证器 - 简化版"""
    
    def __init__(self, mirtarbase_file: str):
        """初始化验证器"""
        self.mirtarbase_df = None
        self.load_mirtarbase(mirtarbase_file)
    
    def load_mirtarbase(self, mirtarbase_file: str):
        """加载miRTarBase数据"""
        if not os.path.exists(mirtarbase_file):
            raise FileNotFoundError(f"miRTarBase file not found: {mirtarbase_file}")
        
        print(f"Loading miRTarBase from: {mirtarbase_file}")
        
        # 智能读取CSV/TSV
        for sep in ['\t', ',', '|']:
            try:
                df = pd.read_csv(mirtarbase_file, sep=sep, low_memory=False)
                if len(df.columns) > 1:
                    break
            except:
                continue
        
        print(f"Loaded {len(df)} records with {len(df.columns)} columns")
        print(f"Columns: {list(df.columns)}")
        
        # 检查必需列
        required_cols = ['miRNA', 'Target Gene']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"⚠️  Warning: Missing columns {missing_cols}")
            print("Available columns:", list(df.columns))
            raise ValueError(f"Required columns not found: {missing_cols}")
        
        self.mirtarbase_df = df
        print(f"✅ miRTarBase loaded successfully")
        
        # 统计信息
        if 'Support Type' in df.columns:
            print("\nSupport Type distribution:")
            print(df['Support Type'].value_counts().head(10))
    
    def validate_pair(self, mirna: str, gene: str) -> Dict[str, any]:
        """验证单个miRNA-gene对"""
        if self.mirtarbase_df is None:
            return {
                'miRTarBase_Validation': False,
                'Validation_Strength': 'No miRTarBase data'
            }
        
        # 查找匹配记录
        matches = self.mirtarbase_df[
            (self.mirtarbase_df['miRNA'] == mirna) & 
            (self.mirtarbase_df['Target Gene'] == gene)
        ]
        
        if matches.empty:
            return {
                'miRTarBase_Validation': False,
                'Validation_Strength': 'Not Validated'
            }
        
        # 收集所有Support Type
        all_support_types = matches['Support Type'].dropna().unique().tolist()
        
        # 收集所有Experiments
        all_experiments = []
        if 'Experiments' in matches.columns:
            for exp_str in matches['Experiments'].dropna().unique():
                if pd.notna(exp_str):
                    # 分割多个实验
                    all_experiments.extend([e.strip() for e in str(exp_str).replace(';', ',').split(',')])
        
        all_experiments = list(set(all_experiments))
        
        # 定义强验证标准
        strong_experiment_keywords = [
            'luciferase reporter assay',
            'luciferase assay', 
            'reporter assay',
            'qrt-pcr',
            'qpcr',
            'western blot',
            'western',
            'immunoblot'
        ]
        
        # 检查强实验证据
        has_strong_experiment = False
        if all_experiments:
            experiments_lower = [exp.lower() for exp in all_experiments]
            has_strong_experiment = any(
                any(keyword in exp for keyword in strong_experiment_keywords)
                for exp in experiments_lower
            )
        
        # 检查Functional MTI
        has_functional_mti = any(
            'functional mti' in st.lower() and 'weak' not in st.lower()
            for st in all_support_types
        )
        
        # 检查Non-Functional
        has_nonfunctional = any(
            'non-functional' in st.lower()
            for st in all_support_types
        )
        
        # 决定验证强度
        if has_functional_mti or has_strong_experiment:
            validation_strength = 'Strong'
            has_validation = True
        elif not has_nonfunctional:
            validation_strength = 'Moderate'
            has_validation = True
        elif has_nonfunctional and len(all_support_types) > 1:
            validation_strength = 'Weak'
            has_validation = True
        else:
            validation_strength = 'Not Validated'
            has_validation = False
        
        return {
            'miRTarBase_Validation': has_validation,
            'Validation_Strength': validation_strength,
            'num_records': len(matches),
            'all_support_types': '; '.join(all_support_types[:3]),  # 最多显示3个
            'all_experiments': '; '.join(all_experiments[:5])  # 最多显示5个
        }


def re_validate_csv(input_csv: str, mirtarbase_file: str, output_csv: Optional[str] = None):
    """重新验证CSV文件中的所有MTI"""
    
    # 读取输入CSV
    print(f"\nReading input CSV: {input_csv}")
    df = pd.read_csv(input_csv)
    print(f"Found {len(df)} MTI records")
    
    # 检查必需列
    if 'miRNA' not in df.columns or 'Target_Gene' not in df.columns:
        raise ValueError("Input CSV must contain 'miRNA' and 'Target_Gene' columns")
    
    # 初始化验证器
    validator = MiRTarBaseValidator(mirtarbase_file)
    
    # 重新验证所有记录
    print("\nRe-validating all MTI pairs...")
    
    new_validation = []
    new_strength = []
    num_records = []
    all_support = []
    all_experiments = []
    
    for idx, row in df.iterrows():
        mirna = row['miRNA']
        gene = row['Target_Gene']
        
        result = validator.validate_pair(mirna, gene)
        
        new_validation.append(result['miRTarBase_Validation'])
        new_strength.append(result['Validation_Strength'])
        num_records.append(result.get('num_records', 0))
        all_support.append(result.get('all_support_types', 'N/A'))
        all_experiments.append(result.get('all_experiments', 'N/A'))
        
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(df)} records...")
    
    # 更新DataFrame
    df['miRTarBase_Validation'] = new_validation
    df['Validation_Strength'] = new_strength
    df['miRTarBase_Num_Records'] = num_records
    df['miRTarBase_Support_Types'] = all_support
    df['miRTarBase_Experiments'] = all_experiments
    
    # 统计
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    print(f"Total MTIs: {len(df)}")
    print(f"\nValidation Status:")
    print(f"  ✅ Validated: {sum(new_validation)}")
    print(f"  ❌ Not Validated: {sum(~pd.Series(new_validation))}")
    print(f"\nValidation Strength Distribution:")
    for strength, count in pd.Series(new_strength).value_counts().items():
        print(f"  {strength}: {count}")
    
    # 保存结果
    if output_csv is None:
        base, ext = os.path.splitext(input_csv)
        output_csv = f"{base}_revalidated{ext}"
    
    df.to_csv(output_csv, index=False)
    print(f"\n✅ Results saved to: {output_csv}")
    
    return df


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Re-validate miRTarBase information in existing results CSV'
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input CSV file with MTI results')
    parser.add_argument('-m', '--mirtarbase', required=True,
                       help='miRTarBase database file')
    parser.add_argument('-o', '--output', 
                       help='Output CSV file (default: input_revalidated.csv)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("miRTarBase Re-validation Tool")
    print("="*60)
    
    # 执行重新验证
    re_validate_csv(args.input, args.mirtarbase, args.output)
    
    print("\n🎉 Re-validation completed successfully!")


if __name__ == "__main__":
    main()
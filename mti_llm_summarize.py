import os
import sys
import json
import pandas as pd
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import io
import re
from dataclasses import dataclass
from langchain_community.llms import Ollama
from langchain_core.callbacks.manager import CallbackManager
from langchain_core.callbacks.streaming_stdout import StreamingStdOutCallbackHandler

# 设置UTF-8编码
if __name__ == '__main__':
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8')

def extract_mirna_simplified(mirna_full):
    """提取简化的miRNA名称，与pipeline文件命名一致"""
    # 匹配 hsa-miR-数字[字母] 格式
    match = re.search(r"hsa-(miR-\d+[a-z]*)", mirna_full)
    if match:
        return match.group(1)
    
    # 如果没匹配到，尝试移除hsa-前缀
    if mirna_full.startswith('hsa-'):
        return mirna_full[4:]
    
    return mirna_full

@dataclass
class MTISummary:
    """MTI总结结果数据类"""
    mirna: str
    gene: str
    relationship_summary: str
    functional_summary: str
    evidence_strength: str
    key_findings: List[str]
    contradictions: List[str]
    clinical_relevance: str
    research_gaps: List[str]
    overall_assessment: str

class MTILLMSummarizer:
    """使用LLM总结MTI关系的工具"""
    
    def __init__(self, model_name="qwen3.5:9b", temperature=0):
        """初始化LLM总结器"""
        self.model_name = model_name
        self.temperature = temperature
        self.llm = self._setup_llm()
        
    def _setup_llm(self):
        """设置LLM模型"""
        try:
            callback_manager = CallbackManager([StreamingStdOutCallbackHandler()])
            llm = Ollama(
                model=self.model_name,
                callback_manager=callback_manager,
                temperature= 0,
                top_p=0.9,
                num_ctx=4096,  # 增加上下文窗口以处理多个abstracts
                repeat_penalty=1.1
            )
            return llm
        except Exception as e:
            print(f"Error setting up LLM: {e}")
            return None
    
    def read_abstracts(self, file_path: str) -> List[Tuple[str, str]]:
        """原始的abstracts读取函数，保持向后兼容"""
        return self.improved_read_abstracts(file_path)
    
    def improved_read_abstracts(self, file_path: str) -> List[Tuple[str, str]]:
        """改进的abstracts读取函数，专门处理标准格式的abstracts"""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                text = f.read()
            
            abstracts = []
            
            # 按双换行符分割文本块
            blocks = [block.strip() for block in text.split('\n\n') if block.strip()]
            
            current_title = None
            current_content = []
            
            for block in blocks:
                lines = block.split('\n')
                
                # 检查是否是标题（方括号格式或单独一行）
                first_line = lines[0].strip()
                
                # 处理方括号标题格式：[Title here]
                if first_line.startswith('[') and first_line.endswith(']'):
                    # 保存之前的abstract
                    if current_title and current_content:
                        content_text = '\n'.join(current_content).strip()
                        if content_text:
                            abstracts.append((current_title, content_text))
                    
                    # 开始新的abstract
                    current_title = first_line[1:-1].strip()  # 移除方括号
                    current_content = []
                    
                    # 处理标题后的内容
                    if len(lines) > 1:
                        remaining_lines = lines[1:]
                        # 跳过"Abstract"标签
                        if remaining_lines and remaining_lines[0].strip().lower() == 'abstract':
                            remaining_lines = remaining_lines[1:]
                        current_content.extend(remaining_lines)
                
                # 检查是否是"Abstract"标签开始的内容块
                elif first_line.lower() == 'abstract':
                    if len(lines) > 1:
                        current_content.extend(lines[1:])
                    else:
                        current_content.append(block)
                
                # 其他情况作为内容添加
                else:
                    if current_title:
                        current_content.append(block)
                    else:
                        # 如果没有明确标题，使用第一行作为标题
                        if len(lines) >= 2:
                            current_title = first_line
                            current_content = lines[1:]
                        else:
                            # 单行内容，使用通用标题
                            current_title = "Research Article"
                            current_content = [block]
            
            # 处理最后一个abstract
            if current_title and current_content:
                content_text = '\n'.join(current_content).strip()
                if content_text:
                    abstracts.append((current_title, content_text))
            
            # 清理和验证abstracts
            cleaned_abstracts = []
            for title, content in abstracts:
                # 移除多余的"Abstract"标签
                if content.lower().startswith('abstract'):
                    content = content[8:].strip()
                
                # 确保内容不为空且有意义
                if len(content) > 50:  # 至少50个字符
                    cleaned_abstracts.append((title, content))
            
            return cleaned_abstracts
            
        except Exception as e:
            print(f"Error reading abstracts from {file_path}: {e}")
            return []
    
    def extract_mti_info_from_filename(self, filename: str) -> Tuple[str, str]:
        """从文件名提取基因和miRNA信息"""
        try:
            # 移除文件扩展名
            basename = os.path.splitext(filename)[0]
            
            # 处理不同的文件名格式
            if '_' in basename:
                # 格式: CYP1B1_miR-27b.txt 或 gene_mirna.txt
                parts = basename.split('_')
                if len(parts) >= 2:
                    gene = parts[0]
                    mirna = '_'.join(parts[1:])  # 处理可能有多个下划线的情况
                    return gene, mirna
            
            # 如果不能从文件名解析，返回空值
            return "", ""
            
        except Exception as e:
            print(f"Error extracting MTI info from filename {filename}: {e}")
            return "", ""
    
    def create_mti_summary_prompt(self, mirna: str, gene: str, abstracts: List[Tuple[str, str]], 
                                  validation_scores: Dict, functions: List[str]) -> str:
        """创建原始MTI总结的prompt（保持向后兼容）"""
        
        # 限制abstracts数量以避免过长prompt
        max_abstracts = 10
        if len(abstracts) > max_abstracts:
            abstracts = abstracts[:max_abstracts]
            abstracts_note = f"\nNote: Showing top {max_abstracts} abstracts from {len(abstracts)} available."
        else:
            abstracts_note = ""
        
        # 准备abstracts文本
        abstracts_text = ""
        for i, (title, content) in enumerate(abstracts, 1):
            # 截断过长的abstract
            content_preview = content[:500] + "..." if len(content) > 500 else content
            abstracts_text += f"\nAbstract {i}:\nTitle: {title}\nContent: {content_preview}\n"
        
        # 准备评分信息
        scores_text = f"""
Validation Scores:
- Traditional Score: {validation_scores.get('traditional_score', 'N/A')}%
- Evidence-Based Score: {validation_scores.get('evidence_based_score', 'N/A')}%
- Final Score: {validation_scores.get('final_score', 'N/A')}%
- Confidence Level: {validation_scores.get('confidence_level', 'N/A')}
- Negative Patterns Found: {validation_scores.get('negative_patterns_count', 0)}
"""
        
        # 准备功能评分
        if functions and 'function_scores' in validation_scores:
            scores_text += "\nFunction-specific Scores:\n"
            for func in functions:
                score = validation_scores['function_scores'].get(func, 'N/A')
                scores_text += f"- {func}: {score}%\n"
        
        prompt = f"""You are an expert in miRNA-target gene interaction research. Please provide a comprehensive summary of the relationship between {mirna} and {gene} based on the following abstracts and validation scores.

{scores_text}{abstracts_note}

Research Abstracts:
{abstracts_text}

Please analyze and provide a structured summary in the following JSON format:
{{
    "relationship_summary": "A comprehensive summary of the {mirna}-{gene} interaction, including direct targeting evidence, regulatory mechanisms, and validation methods used across studies",
    
    "functional_summary": "Summary of biological functions and pathways affected by this interaction, with emphasis on: {', '.join(functions) if functions else 'general biological processes'}",
    
    "evidence_strength": "Assessment of overall evidence quality (Strong/Moderate/Weak) with explanation based on the validation scores and literature consistency",
    
    "key_findings": [
        "List 3-5 most important findings about this interaction from the abstracts"
    ],
    
    "contradictions": [
        "List any contradictory findings or inconsistencies between studies, if any"
    ],
    
    "clinical_relevance": "Summary of potential clinical implications, disease associations, or therapeutic relevance mentioned in the studies",
    
    "research_gaps": [
        "List 2-3 areas where more research is needed based on current evidence"
    ],
    
    "overall_assessment": "A final assessment integrating the computational validation scores with literature evidence, providing a confidence rating for this MTI"
}}

Focus on:
1. Consistency between computational predictions (scores) and experimental evidence
2. Strength of experimental validation methods used
3. Biological plausibility and functional coherence
4. Any concerns raised by negative patterns in the validation

Provide your response in the specified JSON format."""
        
        return prompt
    
    def create_mti_functional_prompt(self, mirna: str, gene: str, abstracts: List[Tuple[str, str]], 
                                    validation_scores: Dict, functions: List[str]) -> str:
        """创建专门用于基因和miRNA功能/通路分析的prompt"""
        
        # 准备abstracts文本
        abstracts_text = ""
        for i, (title, content) in enumerate(abstracts, 1):
            abstracts_text += f"\nAbstract {i}:\nTitle: {title}\nContent: {content}\n"
            abstracts_text += "-" * 80 + "\n"
        
        # 准备评分信息（简化版）
        scores_text = f"""
Computational Validation Scores:
- Final Score: {validation_scores.get('final_score', 'N/A')}%
- Confidence Level: {validation_scores.get('confidence_level', 'N/A')}
"""
        
        prompt = f"""You are a systems biology expert specializing in miRNA-gene regulatory networks, cellular pathways, and functional genomics. Your task is to extract and analyze the biological functions and pathway information for {gene} and {mirna} from the provided research abstracts.

{scores_text}

Research Literature:
{abstracts_text}

Please analyze the abstracts and provide a comprehensive functional assessment in this EXACT JSON format:

{{
    "gene_analysis": {{
        "gene_name": "{gene}",
        "primary_functions": [
            "List 3-5 key biological functions of {gene} mentioned in the literature"
        ],
        "cellular_processes": [
            "List cellular processes that {gene} is involved in (e.g., cell cycle, apoptosis, metabolism)"
        ],
        "signaling_pathways": [
            "List specific signaling pathways where {gene} plays a role (e.g., MAPK, PI3K/AKT, Wnt)"
        ],
        "tissue_expression": "Tissues/cell types where {gene} is predominantly expressed",
        "disease_associations": [
            "Diseases or pathological conditions associated with {gene} dysfunction"
        ]
    }},
    
    "mirna_analysis": {{
        "mirna_name": "{mirna}",
        "regulatory_functions": [
            "List 3-5 key regulatory functions of {mirna} mentioned in the literature"
        ],
        "target_pathways": [
            "List biological pathways that {mirna} regulates"
        ],
        "biological_processes": [
            "List biological processes controlled by {mirna} (focus on: {', '.join(functions) if functions else 'general processes'})"
        ],
        "expression_patterns": "Conditions or contexts where {mirna} expression changes",
        "disease_relevance": [
            "Diseases where {mirna} dysregulation is implicated"
        ]
    }},
    
    "interaction_analysis": {{
        "regulatory_relationship": "How {mirna} regulates {gene} (direct targeting, indirect effects, feedback loops)",
        "functional_consequences": [
            "What happens when {mirna} targets {gene} - downstream effects"
        ],
        "pathway_integration": [
            "Common pathways where both {gene} and {mirna} are involved"
        ],
        "biological_significance": "Overall biological importance of this {mirna}-{gene} interaction"
    }},
    
    "pathway_analysis": {{
        "shared_pathways": [
            "Pathways where both {gene} and {mirna} are active"
        ],
        "upstream_regulators": [
            "Factors that regulate {gene} expression or {mirna} expression"
        ],
        "downstream_targets": [
            "Genes/processes affected by {gene} or {mirna} activity"
        ],
        "crosstalk_networks": [
            "Other miRNAs or genes that interact with this pair"
        ]
    }},
    
    "functional_summary": {{
        "gene_role": "Concise summary of {gene}'s main cellular role",
        "mirna_role": "Concise summary of {mirna}'s main regulatory role", 
        "interaction_impact": "How this interaction affects cellular function",
        "therapeutic_potential": "Potential for targeting this interaction therapeutically",
        "research_priority": "HIGH/MEDIUM/LOW - importance for further functional studies"
    }},
    
    "evidence_quality": {{
        "experimental_methods": [
            "List experimental approaches used to study these functions"
        ],
        "validation_level": "WELL-VALIDATED/PARTIALLY-VALIDATED/PREDICTED",
        "literature_consistency": "Are findings consistent across studies?",
        "functional_gaps": [
            "Aspects of function/pathway involvement that need more research"
        ]
    }}
}}

ANALYSIS GUIDELINES:
1. **Gene Functions**: Focus on enzymatic activity, protein interactions, transcriptional roles
2. **miRNA Functions**: Emphasize regulatory mechanisms, target specificity, expression control
3. **Pathway Integration**: Look for canonical pathways (KEGG, Reactome style names)
4. **Functional Validation**: Distinguish between computational predictions and experimental evidence
5. **Biological Context**: Consider tissue specificity, developmental stages, disease states

SPECIFIC FOCUS AREAS:
- Metabolic pathways and enzyme functions
- Cell signaling cascades and regulatory networks  
- Disease mechanisms and therapeutic targets
- Developmental processes and tissue homeostasis
- Stress responses and environmental adaptations

Extract information conservatively - only include functions and pathways that are explicitly mentioned or strongly implied in the abstracts. If information is not available, state "Not mentioned in provided literature" rather than speculating."""

        return prompt
    
    def parse_llm_response(self, response: str) -> Dict:
        """解析LLM响应"""
        try:
            # 查找JSON部分
            response = response.strip()
            start_idx = response.find('{')
            end_idx = response.rfind('}') + 1
            
            if start_idx != -1 and end_idx > start_idx:
                json_str = response[start_idx:end_idx]
                return json.loads(json_str)
            else:
                # 如果没有找到JSON，返回原始响应
                return {"raw_response": response, "parse_error": "No JSON found in response"}
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON: {e}")
            return {"raw_response": response, "parse_error": str(e)}
    
    def summarize_mti(self, mirna: str, gene: str, abstracts_file: str, 
                      validation_scores: Dict, functions: List[str] = None) -> Dict:
        """为单个MTI生成总结（原始方法，保持向后兼容）"""
        print(f"\nSummarizing {mirna} → {gene} interaction...")
        
        # 读取abstracts
        abstracts = self.read_abstracts(abstracts_file)
        if not abstracts:
            print(f"No abstracts found for {mirna}-{gene}")
            return {
                "error": "No abstracts found",
                "mirna": mirna,
                "gene": gene
            }
        
        print(f"Found {len(abstracts)} abstracts")
        
        # 创建prompt
        prompt = self.create_mti_summary_prompt(mirna, gene, abstracts, validation_scores, functions)
        
        # 调用LLM
        try:
            response = self.llm.invoke(prompt)
            result = self.parse_llm_response(response)
            
            # 添加元数据
            result['mirna'] = mirna
            result['gene'] = gene
            result['abstract_count'] = len(abstracts)
            result['validation_scores'] = validation_scores
            
            return result
            
        except Exception as e:
            print(f"Error generating summary: {e}")
            return {
                "error": str(e),
                "mirna": mirna,
                "gene": gene
            }
    
    def enhanced_summarize_mti_functions(self, mirna: str, gene: str, abstracts_file: str, 
                                       validation_scores: Dict, functions: List[str] = None) -> Dict:
        """增强版MTI功能分析总结"""
        print(f"\nAnalyzing functions and pathways for {mirna} → {gene} interaction...")
        
        # 读取abstracts（使用改进的解析器）
        abstracts = self.improved_read_abstracts(abstracts_file)
        if not abstracts:
            print(f"No abstracts found for {mirna}-{gene}")
            return {
                "error": "No abstracts found",
                "mirna": mirna,
                "gene": gene
            }
        
        print(f"Found {len(abstracts)} abstracts for functional analysis")
        
        # 创建功能分析prompt
        prompt = self.create_mti_functional_prompt(mirna, gene, abstracts, validation_scores, functions)
        
        # 调用LLM
        try:
            response = self.llm.invoke(prompt)
            result = self.parse_llm_response(response)
            
            # 添加元数据
            result['mirna'] = mirna
            result['gene'] = gene
            result['abstract_count'] = len(abstracts)
            result['analysis_type'] = 'functional_pathway_analysis'
            result['validation_scores'] = validation_scores
            
            return result
            
        except Exception as e:
            print(f"Error generating functional analysis: {e}")
            return {
                "error": str(e),
                "mirna": mirna,
                "gene": gene
            }
    
    def batch_summarize(self, validation_results_file: str, abstracts_dir: str, 
                        functions: List[str] = None, output_dir: str = "mti_summaries") -> Dict:
        """批量处理MTI总结（修复文件名匹配问题）"""
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 读取验证结果
        try:
            validation_df = pd.read_csv(validation_results_file)
            print(f"Loaded validation results: {len(validation_df)} entries")
        except Exception as e:
            print(f"Error reading validation results: {e}")
            return {}
        
        all_summaries = {}
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        print(f"Processing {len(validation_df)} MTIs for standard analysis...")
        
        # 处理每个MTI
        for _, row in validation_df.iterrows():
            mirna_full = row['miRNA']
            gene = row['Gene'] if 'Gene' in row else row.get('Target_Gene', '')
            
            # 关键修复：提取简化的miRNA名称
            mirna_simplified = extract_mirna_simplified(mirna_full)
            
            print(f"\nProcessing: {mirna_full} -> {mirna_simplified}, Gene: {gene}")
            
            # 构建abstracts文件路径（修复后的逻辑）
            possible_files = [
                # 优先尝试pipeline生成的标准格式
                os.path.join(abstracts_dir, f"{gene}_{mirna_simplified}.txt"),
                
                # 备用格式
                os.path.join(abstracts_dir, f"{mirna_simplified}_{gene}.txt"),
                
                # 原有的复杂格式（向后兼容）
                os.path.join(abstracts_dir, "direct_mti", f"{gene}_{mirna_full}_direct.txt"),
                os.path.join(abstracts_dir, "functional", f"{gene}_{mirna_full}_functional_combined.txt"),
                os.path.join(abstracts_dir, f"{gene}_{mirna_full}.txt"),
                os.path.join(abstracts_dir, f"{mirna_full}_{gene}.txt")
            ]
            
            abstracts_file = None
            for file_path in possible_files:
                print(f"  Checking: {os.path.basename(file_path)}")
                if os.path.exists(file_path):
                    abstracts_file = file_path
                    print(f"  ✅ Found: {os.path.basename(file_path)}")
                    break
            
            if not abstracts_file:
                print(f"  ❌ No abstracts file found for {mirna_full}-{gene}")
                print(f"     Tried simplified: {mirna_simplified}")
                continue
            
            # 准备验证分数
            validation_scores = {
                'traditional_score': row.get('Traditional_Score', row.get('Overall_Score', 0)),
                'evidence_based_score': row.get('Evidence_Based_Score', 0),
                'final_score': row.get('Final_Score', row.get('Functional_Prediction_Score', 0)),
                'confidence_level': row.get('Confidence_Level', row.get('Prediction_Level', 'Unknown')),
                'negative_patterns_count': row.get('Negative_Patterns_Count', 0)
            }
            
            # 添加功能分数
            if functions:
                function_scores = {}
                for func in functions:
                    clean_func = func.replace(' ', '_').replace('/', '_')
                    score_col = f'Function_Score_{clean_func}'
                    if score_col in row:
                        function_scores[func] = row[score_col]
                validation_scores['function_scores'] = function_scores
            
            # 生成总结（使用简化名称）
            print(f"  Generating summary for {mirna_simplified} -> {gene}")
            summary = self.summarize_mti(mirna_simplified, gene, abstracts_file, validation_scores, functions)
            
            # 同时生成功能分析
            print(f"  Generating functional analysis for {mirna_simplified} -> {gene}")
            functional_analysis = self.enhanced_summarize_mti_functions(
                mirna_simplified, gene, abstracts_file, validation_scores, functions
            )
            
            # 保存单个总结 - 使用简化名称作为键
            mti_key = f"{mirna_simplified}_{gene}"
            all_summaries[mti_key] = {
                "standard_summary": summary,
                "functional_analysis": functional_analysis
            }
            
            # 保存单个MTI的详细报告
            single_report_file = os.path.join(output_dir, f"{mti_key}_summary_{timestamp}.json")
            with open(single_report_file, 'w', encoding='utf-8') as f:
                json.dump(all_summaries[mti_key], f, indent=2, ensure_ascii=False)
            
            print(f"  ✅ Completed analysis for {mti_key}")
        
        # 保存所有总结
        all_summaries_file = os.path.join(output_dir, f"all_mti_summaries_{timestamp}.json")
        with open(all_summaries_file, 'w', encoding='utf-8') as f:
            json.dump(all_summaries, f, indent=2, ensure_ascii=False)
        
        print(f"\n🎉 Generated {len(all_summaries)} MTI summaries")
        print(f"Results saved to: {output_dir}")
        
        # 生成综合报告
        self._generate_comprehensive_report(all_summaries, output_dir, timestamp)
        
        return all_summaries
    
    def batch_functional_analysis(self, validation_results_file: str, abstracts_dir: str, 
                                functions: List[str] = None, output_dir: str = "mti_functional_summaries") -> Dict:
        """批量功能分析处理（修复文件名匹配）"""
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 读取验证结果
        try:
            validation_df = pd.read_csv(validation_results_file)
            print(f"Loaded validation results: {len(validation_df)} entries")
        except Exception as e:
            print(f"Error reading validation results: {e}")
            return {}
        
        all_functional_analyses = {}
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        print(f"Processing {len(validation_df)} MTIs for functional analysis...")
        
        # 处理每个MTI
        for _, row in validation_df.iterrows():
            mirna_full = row['miRNA'] 
            gene = row['Gene'] if 'Gene' in row else row.get('Target_Gene', '')
            
            # 关键修复：提取简化的miRNA名称
            mirna_simplified = extract_mirna_simplified(mirna_full)
            
            print(f"\nProcessing functional analysis: {mirna_full} -> {mirna_simplified}, Gene: {gene}")
            
            # 构建abstracts文件路径（修复后的逻辑）
            possible_filenames = [
                f"{gene}_{mirna_simplified}.txt",  # 标准格式
                f"{mirna_simplified}_{gene}.txt",  # 备用格式
                f"{gene}_{mirna_full}.txt",        # 完整名称格式（向后兼容）
                f"{mirna_full}_{gene}.txt",
                f"{gene}_{mirna_simplified}_direct.txt",
                f"{gene}_{mirna_simplified}_functional_combined.txt"
            ]
            
            abstracts_file = None
            for filename in possible_filenames:
                # 直接在主目录查找
                file_path = os.path.join(abstracts_dir, filename)
                print(f"  Checking: {filename}")
                if os.path.exists(file_path):
                    abstracts_file = file_path
                    print(f"  ✅ Found: {filename}")
                    break
                
                # 也检查子目录
                for subdir in ['direct_mti', 'functional', '']:
                    if subdir:
                        sub_path = os.path.join(abstracts_dir, subdir, filename)
                        if os.path.exists(sub_path):
                            abstracts_file = sub_path
                            print(f"  ✅ Found in {subdir}/: {filename}")
                            break
                
                if abstracts_file:
                    break
            
            if not abstracts_file:
                print(f"  ❌ No abstracts file found for {mirna_full}-{gene}")
                print(f"     Tried simplified: {mirna_simplified}")
                continue
            
            # 准备验证分数
            validation_scores = {
                'final_score': row.get('Final_Score', row.get('Overall_Score', 0)),
                'confidence_level': row.get('Confidence_Level', 'Unknown')
            }
            
            # 生成功能分析
            print(f"  Generating functional analysis for {mirna_simplified} -> {gene}")
            functional_analysis = self.enhanced_summarize_mti_functions(
                mirna_simplified, gene, abstracts_file, validation_scores, functions
            )
            
            # 保存单个分析 - 使用简化名称作为键
            mti_key = f"{mirna_simplified}_{gene}"
            all_functional_analyses[mti_key] = functional_analysis
            
            # 保存单个MTI的功能分析报告
            single_report_file = os.path.join(output_dir, f"{mti_key}_functional_analysis_{timestamp}.json")
            with open(single_report_file, 'w', encoding='utf-8') as f:
                json.dump(functional_analysis, f, indent=2, ensure_ascii=False)
            
            print(f"  ✅ Completed functional analysis for {mti_key}")
        
        # 保存所有功能分析
        all_analyses_file = os.path.join(output_dir, f"all_functional_analyses_{timestamp}.json")
        with open(all_analyses_file, 'w', encoding='utf-8') as f:
            json.dump(all_functional_analyses, f, indent=2, ensure_ascii=False)
        
        print(f"\n🎉 Generated {len(all_functional_analyses)} functional analyses")
        print(f"Results saved to: {output_dir}")
        
        # 生成功能分析报告
        self._generate_functional_report(all_functional_analyses, output_dir, timestamp)
        
        return all_functional_analyses
    
    def _generate_comprehensive_report(self, all_summaries: Dict, output_dir: str, timestamp: str):
        """生成综合报告"""
        report_file = os.path.join(output_dir, f"comprehensive_mti_report_{timestamp}.txt")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*100 + "\n")
            f.write("COMPREHENSIVE MTI ANALYSIS REPORT\n")
            f.write("="*100 + "\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total MTIs analyzed: {len(all_summaries)}\n\n")
            
            # 处理可能的两层结构（新版本有standard_summary和functional_analysis）
            strong_evidence = []
            moderate_evidence = []
            weak_evidence = []
            
            for mti_key, summary_data in all_summaries.items():
                # 兼容新旧格式
                if isinstance(summary_data, dict) and 'standard_summary' in summary_data:
                    summary = summary_data['standard_summary']
                else:
                    summary = summary_data
                
                if 'error' in summary:
                    continue
                
                evidence_strength = summary.get('evidence_strength', '').lower()
                if 'strong' in evidence_strength:
                    strong_evidence.append((mti_key, summary))
                elif 'moderate' in evidence_strength:
                    moderate_evidence.append((mti_key, summary))
                else:
                    weak_evidence.append((mti_key, summary))
            
            # 写入强证据MTIs
            if strong_evidence:
                f.write("="*80 + "\n")
                f.write("HIGH CONFIDENCE MTIs (Strong Evidence)\n")
                f.write("="*80 + "\n\n")
                
                for mti_key, summary in strong_evidence:
                    mirna = summary.get('mirna', '')
                    gene = summary.get('gene', '')
                    f.write(f"\n{mirna} → {gene}\n")
                    f.write("-"*60 + "\n")
                    f.write(f"Relationship: {summary.get('relationship_summary', 'N/A')}\n\n")
                    f.write(f"Clinical Relevance: {summary.get('clinical_relevance', 'N/A')}\n\n")
                    f.write(f"Overall Assessment: {summary.get('overall_assessment', 'N/A')}\n")
                    f.write("-"*60 + "\n")
            
            # 写入中等证据MTIs
            if moderate_evidence:
                f.write("\n" + "="*80 + "\n")
                f.write("MODERATE CONFIDENCE MTIs\n")
                f.write("="*80 + "\n\n")
                
                for mti_key, summary in moderate_evidence[:10]:  # 只显示前10个
                    mirna = summary.get('mirna', '')
                    gene = summary.get('gene', '')
                    f.write(f"\n{mirna} → {gene}\n")
                    key_findings = summary.get('key_findings', [])
                    if key_findings:
                        f.write(f"Key Findings: {'; '.join(key_findings[:2])}\n")
            
            # 统计信息
            f.write("\n" + "="*80 + "\n")
            f.write("SUMMARY STATISTICS\n")
            f.write("="*80 + "\n")
            f.write(f"Strong Evidence MTIs: {len(strong_evidence)}\n")
            f.write(f"Moderate Evidence MTIs: {len(moderate_evidence)}\n")
            f.write(f"Weak Evidence MTIs: {len(weak_evidence)}\n")
            f.write(f"Failed to analyze: {len(all_summaries) - len(strong_evidence) - len(moderate_evidence) - len(weak_evidence)}\n")
        
        print(f"Comprehensive report saved to: {report_file}")
    
    def _generate_functional_report(self, all_analyses: Dict, output_dir: str, timestamp: str):
        """生成功能分析综合报告"""
        report_file = os.path.join(output_dir, f"functional_analysis_report_{timestamp}.txt")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*100 + "\n")
            f.write("miRNA-GENE FUNCTIONAL AND PATHWAY ANALYSIS REPORT\n")
            f.write("="*100 + "\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total MTI pairs analyzed: {len(all_analyses)}\n\n")
            
            # 按研究优先级分类
            high_priority = []
            medium_priority = []
            low_priority = []
            
            for mti_key, analysis in all_analyses.items():
                if 'error' in analysis:
                    continue
                    
                priority = analysis.get('functional_summary', {}).get('research_priority', 'LOW').upper()
                if priority == 'HIGH':
                    high_priority.append((mti_key, analysis))
                elif priority == 'MEDIUM':
                    medium_priority.append((mti_key, analysis))
                else:
                    low_priority.append((mti_key, analysis))
            
            # 高优先级MTIs
            if high_priority:
                f.write("="*80 + "\n")
                f.write("HIGH PRIORITY MTIs FOR FUNCTIONAL STUDIES\n") 
                f.write("="*80 + "\n\n")
                
                for mti_key, analysis in high_priority:
                    mirna = analysis.get('mirna', '')
                    gene = analysis.get('gene', '')
                    f.write(f"\n{mirna} → {gene}\n")
                    f.write("-"*60 + "\n")
                    
                    gene_role = analysis.get('functional_summary', {}).get('gene_role', 'N/A')
                    mirna_role = analysis.get('functional_summary', {}).get('mirna_role', 'N/A') 
                    interaction_impact = analysis.get('functional_summary', {}).get('interaction_impact', 'N/A')
                    
                    f.write(f"Gene Function: {gene_role}\n")
                    f.write(f"miRNA Function: {mirna_role}\n")
                    f.write(f"Interaction Impact: {interaction_impact}\n")
                    
                    # 共同通路
                    shared_pathways = analysis.get('pathway_analysis', {}).get('shared_pathways', [])
                    if shared_pathways:
                        f.write(f"Key Pathways: {', '.join(shared_pathways[:3])}\n")
                    
                    f.write("-"*60 + "\n")
            
            # 通路统计
            pathway_counts = {}
            for analysis in all_analyses.values():
                if 'error' in analysis:
                    continue
                pathways = analysis.get('pathway_analysis', {}).get('shared_pathways', [])
                for pathway in pathways:
                    pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1
            
            f.write("\n" + "="*80 + "\n")
            f.write("TOP BIOLOGICAL PATHWAYS\n")
            f.write("="*80 + "\n")
            
            sorted_pathways = sorted(pathway_counts.items(), key=lambda x: x[1], reverse=True)
            for pathway, count in sorted_pathways[:10]:
                f.write(f"{pathway}: {count} MTI pairs\n")
            
            # 统计信息
            f.write("\n" + "="*80 + "\n")
            f.write("SUMMARY STATISTICS\n")
            f.write("="*80 + "\n")
            f.write(f"High Priority MTIs: {len(high_priority)}\n")
            f.write(f"Medium Priority MTIs: {len(medium_priority)}\n") 
            f.write(f"Low Priority MTIs: {len(low_priority)}\n")
            f.write(f"Most Common Pathways: {len(pathway_counts)} unique pathways identified\n")
        
        print(f"Functional analysis report saved to: {report_file}")


def main():
    """主函数，用于独立运行测试"""
    import argparse
    
    parser = argparse.ArgumentParser(description='MTI LLM Summarizer - Fixed Version')
    parser.add_argument('-v', '--validation-results', required=True, 
                        help='Path to validation results CSV file')
    parser.add_argument('-a', '--abstracts-dir', required=True,
                        help='Directory containing abstracts files')
    parser.add_argument('-f', '--functions', 
                        help='Comma-separated list of functions to analyze')
    parser.add_argument('-o', '--output-dir', default='mti_summaries',
                        help='Output directory for summaries')
    parser.add_argument('--mode', choices=['standard', 'functional', 'both'], default='both',
                        help='Analysis mode: standard (original), functional (pathway focus), or both')
    
    args = parser.parse_args()
    
    # 解析功能列表
    functions = None
    if args.functions:
        functions = [f.strip() for f in args.functions.split(',')]
    
    # 创建总结器并运行
    print("="*80)
    print("MTI LLM SUMMARIZER - FIXED VERSION")
    print("="*80)
    print(f"Mode: {args.mode}")
    print(f"Validation file: {args.validation_results}")
    print(f"Abstracts directory: {args.abstracts_dir}")
    print(f"Output directory: {args.output_dir}")
    if functions:
        print(f"Functions: {', '.join(functions)}")
    print("="*80)
    
    summarizer = MTILLMSummarizer()
    
    if not summarizer.llm:
        print("❌ ERROR: Could not initialize LLM. Please check Ollama is running.")
        return
    
    if args.mode == 'standard':
        summarizer.batch_summarize(
            validation_results_file=args.validation_results,
            abstracts_dir=args.abstracts_dir,
            functions=functions,
            output_dir=args.output_dir
        )
    elif args.mode == 'functional':
        summarizer.batch_functional_analysis(
            validation_results_file=args.validation_results,
            abstracts_dir=args.abstracts_dir,
            functions=functions,
            output_dir=args.output_dir
        )
    else:  # both
        print("Running standard analysis...")
        summarizer.batch_summarize(
            validation_results_file=args.validation_results,
            abstracts_dir=args.abstracts_dir,
            functions=functions,
            output_dir=args.output_dir
        )
        
        print("\nRunning functional analysis...")
        functional_output_dir = args.output_dir + "_functional"
        summarizer.batch_functional_analysis(
            validation_results_file=args.validation_results,
            abstracts_dir=args.abstracts_dir,
            functions=functions,
            output_dir=functional_output_dir
        )


if __name__ == "__main__":
    main()

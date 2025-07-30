import re
import json
import sys
import os
import argparse
import numpy as np
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import nltk
from nltk.tokenize import sent_tokenize
from typing import List, Dict, Tuple, Optional
import pandas as pd
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')

# 下载必要的NLTK数据
try:
    nltk.data.find('tokenizers/punkt')
except LookupError:
    print("下载NLTK punkt数据...")
    nltk.download('punkt', quiet=True)

@dataclass
class MTIValidation:
    """MTI验证结果数据类"""
    mirna: str
    target_gene: str
    overall_score: float
    gene_relation_score: float
    function_scores: Dict[str, float]
    evidence_count: int
    confidence_level: str
    key_evidence: List[str]
    detailed_scores: Dict[str, float]
    function_evidence: Dict[str, List[str]]  # 每个功能的证据
    relation_evidence: List[str]  # MTI关系的证据
    experimental_evidence: List[str]  # 实验证据

class MTIValidator:
    """基于BERT的miRNA-Target Gene Interaction验证器"""
    
    def __init__(self, model_name='NeuML/pubmedbert-base-embeddings'):
        """
        初始化MTI验证器
        Args:
            model_name: Sentence-BERT模型名称
        """
        print(f"初始化MTI验证器，加载模型: {model_name}...")
        self.model = SentenceTransformer(model_name)
        
        # 评分权重配置
        self.weights = {
            'direct_relation': 0.35,      # 直接关系证据
            'functional_relevance': 0.25,  # 功能相关性
            'co_occurrence': 0.15,        # 共现频率
            'semantic_coherence': 0.15,   # 语义连贯性
            'experimental_evidence': 0.10  # 实验证据
        }
        
        # 关系类型关键词
        self.relation_keywords = {
            'direct_targeting': {
                'keywords': ['target', 'targets', 'targeting', 'bind', 'binds', 'binding', 
                           'direct', 'directly', '3\'UTR', '3\' UTR', 'seed sequence'],
                'weight': 1.0
            },
            'regulation': {
                'keywords': ['regulate', 'regulates', 'regulation', 'modulate', 'modulates',
                           'control', 'controls', 'affect', 'affects'],
                'weight': 0.8
            },
            'inhibition': {
                'keywords': ['inhibit', 'inhibits', 'suppress', 'suppresses', 'repress',
                           'represses', 'downregulate', 'downregulates', 'decrease'],
                'weight': 0.9
            },
            'activation': {
                'keywords': ['activate', 'activates', 'upregulate', 'upregulates', 
                           'increase', 'promote', 'promotes', 'enhance', 'enhances'],
                'weight': 0.8
            }
        }
        
        # 实验证据关键词
        self.experimental_keywords = [
            'luciferase', 'reporter', 'transfection', 'overexpression', 'knockdown',
            'western blot', 'qPCR', 'RT-PCR', 'immunoprecipitation', 'pull-down',
            'CLIP', 'RIP', 'AGO', 'validation', 'confirmed', 'demonstrated'
        ]
    
    def parse_abstracts(self, text: str) -> List[str]:
        """解析多个摘要文本"""
        # 多种分割模式
        abstracts = []
        
        # 尝试按空行分割
        blocks = re.split(r'\n\s*\n', text)
        
        for block in blocks:
            block = block.strip()
            if len(block) > 50:  # 过滤太短的文本
                abstracts.append(block)
        
        # 如果分割效果不好，尝试其他模式
        if len(abstracts) <= 1:
            # 尝试按Abstract关键词分割
            abstracts = re.split(r'(?i)(?:^|\n)(?:abstract|title)[:：\s]*\n?', text)
            abstracts = [a.strip() for a in abstracts if len(a.strip()) > 50]
        
        return abstracts if abstracts else [text]
    
    def validate_mti(self, abstracts: List[str], mirna: str, target_gene: str, 
                    functions: List[str] = None) -> MTIValidation:
        """
        验证单个MTI的有效性
        
        Args:
            abstracts: 摘要文本列表
            mirna: miRNA名称
            target_gene: 目标基因
            functions: 要评估的生物学功能列表
            
        Returns:
            MTIValidation对象
        """
        if functions is None:
            functions = []
        
        # 初始化评分组件
        scores = {
            'direct_relation': 0,
            'functional_relevance': 0,
            'co_occurrence': 0,
            'semantic_coherence': 0,
            'experimental_evidence': 0
        }
        
        function_scores = {func: 0 for func in functions}
        function_evidence = {func: [] for func in functions}
        key_evidence = []
        relation_evidence = []
        experimental_evidence = []
        evidence_count = 0
        
        # 分析每个摘要
        for abstract in abstracts:
            # 1. 直接关系分析
            relation_score, rel_evidence = self._analyze_direct_relation(
                abstract, mirna, target_gene
            )
            scores['direct_relation'] = max(scores['direct_relation'], relation_score)
            if rel_evidence:
                relation_evidence.extend(rel_evidence[:3])
                key_evidence.extend(rel_evidence[:2])
            
            # 2. 功能相关性分析
            for func in functions:
                func_score, func_evidence = self._analyze_functional_relevance_with_evidence(
                    abstract, mirna, target_gene, func
                )
                function_scores[func] = max(function_scores[func], func_score)
                if func_evidence:
                    function_evidence[func].extend(func_evidence[:2])
            
            # 3. 共现分析
            co_score = self._analyze_cooccurrence(abstract, mirna, target_gene)
            scores['co_occurrence'] = max(scores['co_occurrence'], co_score)
            
            # 4. 语义连贯性分析
            semantic_score = self._analyze_semantic_coherence(
                abstract, mirna, target_gene
            )
            scores['semantic_coherence'] = max(scores['semantic_coherence'], semantic_score)
            
            # 5. 实验证据分析
            exp_score, exp_evidence = self._analyze_experimental_evidence_with_details(
                abstract, mirna, target_gene
            )
            scores['experimental_evidence'] = max(scores['experimental_evidence'], exp_score)
            if exp_evidence:
                experimental_evidence.extend(exp_evidence[:2])
            
            if relation_score > 0.5 or exp_score > 0.5:
                evidence_count += 1
        
        # 计算功能相关性总分
        if function_scores:
            scores['functional_relevance'] = np.mean(list(function_scores.values()))
        
        # 计算基因关系得分
        gene_relation_score = (
            scores['direct_relation'] * 0.5 +
            scores['co_occurrence'] * 0.2 +
            scores['semantic_coherence'] * 0.2 +
            scores['experimental_evidence'] * 0.1
        )
        
        # 计算总分
        overall_score = sum(scores[k] * self.weights[k] for k in scores.keys())
        
        # 确定置信度等级
        confidence_level = self._determine_confidence_level(
            overall_score, gene_relation_score, evidence_count
        )
        
        # 去重和限制证据数量
        key_evidence = list(dict.fromkeys(key_evidence))[:5]
        relation_evidence = list(dict.fromkeys(relation_evidence))[:5]
        experimental_evidence = list(dict.fromkeys(experimental_evidence))[:5]
        
        for func in function_evidence:
            function_evidence[func] = list(dict.fromkeys(function_evidence[func]))[:5]
        
        return MTIValidation(
            mirna=mirna,
            target_gene=target_gene,
            overall_score=overall_score * 100,
            gene_relation_score=gene_relation_score * 100,
            function_scores={k: v * 100 for k, v in function_scores.items()},
            evidence_count=evidence_count,
            confidence_level=confidence_level,
            key_evidence=key_evidence,
            detailed_scores={k: v * 100 for k, v in scores.items()},
            function_evidence=function_evidence,
            relation_evidence=relation_evidence,
            experimental_evidence=experimental_evidence
        )
    
    def _analyze_direct_relation(self, text: str, mirna: str, gene: str) -> Tuple[float, List[str]]:
        """分析直接靶向关系"""
        sentences = sent_tokenize(text)
        max_score = 0
        evidence_sentences = []
        
        for sent in sentences:
            sent_lower = sent.lower()
            mirna_lower = mirna.lower()
            gene_lower = gene.lower()
            
            # 检查是否同时包含miRNA和基因
            if mirna_lower in sent_lower and gene_lower in sent_lower:
                score = 0.3  # 基础共现分
                
                # 检查关系类型
                for rel_type, rel_info in self.relation_keywords.items():
                    if any(kw in sent_lower for kw in rel_info['keywords']):
                        score = max(score, 0.5 * rel_info['weight'])
                        
                        # 直接证据加分
                        if 'direct' in sent_lower or '3\'utr' in sent_lower or 'binding site' in sent_lower:
                            score = min(1.0, score + 0.3)
                        
                        evidence_sentences.append(sent)
                        break
                
                max_score = max(max_score, score)
        
        # 使用语义相似度增强
        if sentences:
            query = f"{mirna} targets {gene}"
            query_embedding = self.model.encode([query])
            sent_embeddings = self.model.encode(sentences)
            
            similarities = cosine_similarity(query_embedding, sent_embeddings)[0]
            semantic_score = np.max(similarities)
            
            # 结合语义得分
            max_score = max_score * 0.5 + semantic_score * 0.5
        
        return max_score, evidence_sentences
    
    def _analyze_functional_relevance(self, text: str, mirna: str, gene: str, function: str) -> float:
        """分析功能相关性"""
        # 分割复合功能（如 "immune ferroptosis" 分为 "immune" 和 "ferroptosis"）
        function_parts = function.split()
        
        # 构建功能相关查询
        queries = [
            f"{mirna} {function}",
            f"{gene} {function}",
            f"{mirna} regulates {gene} in {function}",
            f"{mirna} {gene} {function}"
        ]
        
        # 为复合功能添加额外查询
        if len(function_parts) > 1:
            for part in function_parts:
                queries.extend([
                    f"{mirna} {part}",
                    f"{gene} {part}"
                ])
        
        # 编码查询和文本
        query_embeddings = self.model.encode(queries)
        sentences = sent_tokenize(text)
        
        if not sentences:
            return 0
        
        sent_embeddings = self.model.encode(sentences)
        
        # 计算最大相似度
        max_similarities = []
        for q_emb in query_embeddings:
            similarities = cosine_similarity([q_emb], sent_embeddings)[0]
            max_similarities.append(np.max(similarities))
        
        # 检查关键词匹配
        text_lower = text.lower()
        keyword_score = 0
        
        # 对于复合功能，检查所有部分
        if len(function_parts) > 1:
            # 检查完整功能短语
            if function.lower() in text_lower:
                keyword_score = 1.0
            else:
                # 检查各个部分
                parts_found = sum(1 for part in function_parts if part.lower() in text_lower)
                keyword_score = parts_found / len(function_parts) * 0.8
        else:
            # 单一功能词
            if function.lower() in text_lower:
                keyword_score = 0.5
                if mirna.lower() in text_lower or gene.lower() in text_lower:
                    keyword_score = 0.8
                if mirna.lower() in text_lower and gene.lower() in text_lower:
                    keyword_score = 1.0
        
        # 组合得分
        semantic_score = np.max(max_similarities)
        return semantic_score * 0.7 + keyword_score * 0.3
    
    def _analyze_functional_relevance_with_evidence(self, text: str, mirna: str, gene: str, function: str) -> Tuple[float, List[str]]:
        """分析功能相关性并返回证据"""
        sentences = sent_tokenize(text)
        evidence_sentences = []
        
        # 分割复合功能
        function_parts = function.split()
        
        # 构建功能相关查询
        queries = [
            f"{mirna} {function}",
            f"{gene} {function}",
            f"{mirna} regulates {gene} in {function}",
            f"{mirna} {gene} {function}"
        ]
        
        # 为复合功能添加额外查询
        if len(function_parts) > 1:
            for part in function_parts:
                queries.extend([
                    f"{mirna} {part}",
                    f"{gene} {part}"
                ])
        
        # 编码查询和句子
        if not sentences:
            return 0, []
        
        query_embeddings = self.model.encode(queries)
        sent_embeddings = self.model.encode(sentences)
        
        # 找出最相关的句子
        best_sentences = []
        for q_emb in query_embeddings:
            similarities = cosine_similarity([q_emb], sent_embeddings)[0]
            top_indices = np.argsort(similarities)[-3:][::-1]  # Top 3
            
            for idx in top_indices:
                if similarities[idx] > 0.5:  # 相似度阈值
                    sent = sentences[idx]
                    sent_lower = sent.lower()
                    
                    # 额外检查关键词
                    if (function.lower() in sent_lower or 
                        any(part.lower() in sent_lower for part in function_parts)):
                        
                        # 优先选择同时包含miRNA或基因的句子
                        priority = 0
                        if mirna.lower() in sent_lower:
                            priority += 1
                        if gene.lower() in sent_lower:
                            priority += 1
                        
                        best_sentences.append((similarities[idx] + priority * 0.1, sent))
        
        # 去重并排序
        unique_sentences = {}
        for score, sent in best_sentences:
            if sent not in unique_sentences or score > unique_sentences[sent]:
                unique_sentences[sent] = score
        
        # 选择最好的证据
        sorted_evidence = sorted(unique_sentences.items(), key=lambda x: x[1], reverse=True)
        evidence_sentences = [sent for sent, _ in sorted_evidence[:3]]
        
        # 计算得分
        score = self._analyze_functional_relevance(text, mirna, gene, function)
        
        return score, evidence_sentences
    
    def _analyze_experimental_evidence(self, text: str, mirna: str, gene: str) -> float:
        """分析实验证据"""
        text_lower = text.lower()
        score = 0
        
        # 检查实验方法关键词
        exp_methods_found = sum(1 for kw in self.experimental_keywords 
                              if kw in text_lower)
        
        # 检查是否在同一句子中提到实验验证
        sentences = sent_tokenize(text)
        for sent in sentences:
            sent_lower = sent.lower()
            if (mirna.lower() in sent_lower or gene.lower() in sent_lower):
                if any(kw in sent_lower for kw in self.experimental_keywords):
                    score = max(score, 0.7)
                    if mirna.lower() in sent_lower and gene.lower() in sent_lower:
                        score = 1.0
                        break
        
        # 结合整体实验方法提及
        method_score = min(1.0, exp_methods_found / 3)
        
        return score * 0.7 + method_score * 0.3
    
    def _analyze_experimental_evidence_with_details(self, text: str, mirna: str, gene: str) -> Tuple[float, List[str]]:
        """分析实验证据并返回具体句子"""
        sentences = sent_tokenize(text)
        evidence_sentences = []
        max_score = 0
        
        for sent in sentences:
            sent_lower = sent.lower()
            
            # 检查实验方法
            exp_methods = [method for method in self.experimental_keywords 
                          if method in sent_lower]
            
            if exp_methods:
                # 检查是否提到miRNA或基因
                mentions_mirna = mirna.lower() in sent_lower
                mentions_gene = gene.lower() in sent_lower
                
                if mentions_mirna or mentions_gene:
                    score = 0.5
                    
                    # 高价值实验方法
                    high_value_methods = ['luciferase', 'reporter', 'transfection', 
                                        'overexpression', 'knockdown', 'chip', 'clip']
                    if any(method in exp_methods for method in high_value_methods):
                        score = 0.8
                    
                    # 同时提到两者
                    if mentions_mirna and mentions_gene:
                        score = min(1.0, score + 0.2)
                    
                    # 包含验证词汇
                    if any(word in sent_lower for word in ['confirmed', 'validated', 
                                                           'demonstrated', 'showed']):
                        score = min(1.0, score + 0.1)
                    
                    if score > 0.6:
                        evidence_sentences.append(sent)
                    
                    max_score = max(max_score, score)
        
        return max_score, evidence_sentences
        """分析功能相关性"""
        # 分割复合功能（如 "immune ferroptosis" 分为 "immune" 和 "ferroptosis"）
        function_parts = function.split()
        
        # 构建功能相关查询
        queries = [
            f"{mirna} {function}",
            f"{gene} {function}",
            f"{mirna} regulates {gene} in {function}",
            f"{mirna} {gene} {function}"
        ]
        
        # 为复合功能添加额外查询
        if len(function_parts) > 1:
            for part in function_parts:
                queries.extend([
                    f"{mirna} {part}",
                    f"{gene} {part}"
                ])
        
        # 编码查询和文本
        query_embeddings = self.model.encode(queries)
        sentences = sent_tokenize(text)
        
        if not sentences:
            return 0
        
        sent_embeddings = self.model.encode(sentences)
        
        # 计算最大相似度
        max_similarities = []
        for q_emb in query_embeddings:
            similarities = cosine_similarity([q_emb], sent_embeddings)[0]
            max_similarities.append(np.max(similarities))
        
        # 检查关键词匹配
        text_lower = text.lower()
        keyword_score = 0
        
        # 对于复合功能，检查所有部分
        if len(function_parts) > 1:
            # 检查完整功能短语
            if function.lower() in text_lower:
                keyword_score = 1.0
            else:
                # 检查各个部分
                parts_found = sum(1 for part in function_parts if part.lower() in text_lower)
                keyword_score = parts_found / len(function_parts) * 0.8
        else:
            # 单一功能词
            if function.lower() in text_lower:
                keyword_score = 0.5
                if mirna.lower() in text_lower or gene.lower() in text_lower:
                    keyword_score = 0.8
                if mirna.lower() in text_lower and gene.lower() in text_lower:
                    keyword_score = 1.0
        
        # 组合得分
        semantic_score = np.max(max_similarities)
        return semantic_score * 0.6 + keyword_score * 0.4
    
    def _analyze_cooccurrence(self, text: str, mirna: str, gene: str) -> float:
        """分析共现频率"""
        text_lower = text.lower()
        mirna_lower = mirna.lower()
        gene_lower = gene.lower()
        
        # 计算出现次数
        mirna_count = text_lower.count(mirna_lower)
        gene_count = len(re.findall(r'\b' + re.escape(gene_lower) + r'\b', text_lower))
        
        # 计算句子级共现
        sentences = sent_tokenize(text)
        co_occur_sents = sum(1 for sent in sentences 
                           if mirna_lower in sent.lower() and gene_lower in sent.lower())
        
        # 归一化得分
        occurrence_score = min(1.0, (mirna_count + gene_count) / 20)
        co_occur_score = min(1.0, co_occur_sents / 3)
        
        return occurrence_score * 0.4 + co_occur_score * 0.6
    
    def _analyze_semantic_coherence(self, text: str, mirna: str, gene: str) -> float:
        """分析语义连贯性"""
        # 构建上下文查询
        context_queries = [
            f"{mirna} regulation mechanism",
            f"{gene} expression regulation",
            f"microRNA target interaction",
            f"{mirna} {gene} pathway"
        ]
        
        # 编码
        query_embeddings = self.model.encode(context_queries)
        text_embedding = self.model.encode([text])
        
        # 计算相似度
        similarities = cosine_similarity(query_embeddings, text_embedding)
        
        return np.mean(similarities)
    
    def _analyze_experimental_evidence(self, text: str, mirna: str, gene: str) -> float:
        """分析实验证据"""
        text_lower = text.lower()
        score = 0
        
        # 检查实验方法关键词
        exp_methods_found = sum(1 for kw in self.experimental_keywords 
                              if kw in text_lower)
        
        # 检查是否在同一句子中提到实验验证
        sentences = sent_tokenize(text)
        for sent in sentences:
            sent_lower = sent.lower()
            if (mirna.lower() in sent_lower or gene.lower() in sent_lower):
                if any(kw in sent_lower for kw in self.experimental_keywords):
                    score = max(score, 0.7)
                    if mirna.lower() in sent_lower and gene.lower() in sent_lower:
                        score = 1.0
                        break
        
        # 结合整体实验方法提及
        method_score = min(1.0, exp_methods_found / 3)
        
        return score * 0.7 + method_score * 0.3
    
    def _determine_confidence_level(self, overall_score: float, gene_score: float, 
                                  evidence_count: int) -> str:
        """确定置信度等级"""
        if overall_score > 0.8 and gene_score > 0.7 and evidence_count >= 2:
            return "High Confidence"
        elif overall_score > 0.6 and gene_score > 0.5:
            return "Medium Confidence"
        elif overall_score > 0.4 or gene_score > 0.4:
            return "Low Confidence"
        else:
            return "Insufficient Evidence"
    
    def normalize_mirna_name(self, mirna: str) -> str:
        """标准化miRNA名称格式"""
        # 转换为标准格式 miR-XXX
        mirna_lower = mirna.lower()
        if mirna_lower.startswith('mir-'):
            return 'miR-' + mirna[4:]
        elif mirna_lower.startswith('microrna-'):
            return 'miR-' + mirna[9:]
        elif mirna_lower.startswith('let-'):
            return mirna.lower()  # let-7系列保持小写
        return mirna
    
    def validate_mti_batch(self, abstracts_dict: Dict[str, List[str]], 
                         functions: List[str] = None) -> pd.DataFrame:
        """
        批量验证MTI
        
        Args:
            abstracts_dict: {(mirna, gene): [abstracts]}格式的字典
            functions: 要评估的功能列表
            
        Returns:
            包含所有MTI验证结果的DataFrame
        """
        results = []
        
        for (mirna, gene), abstracts in abstracts_dict.items():
            # 标准化miRNA名称
            mirna = self.normalize_mirna_name(mirna)
            
            print(f"验证 {mirna} -> {gene}...")
            validation = self.validate_mti(abstracts, mirna, gene, functions)
            
            # 构建结果行
            row = {
                'miRNA': mirna,
                'Target_Gene': gene,
                'Overall_Score': validation.overall_score,
                'Gene_Relation_Score': validation.gene_relation_score,
                'Evidence_Count': validation.evidence_count,
                'Confidence_Level': validation.confidence_level,
                'Valid': validation.overall_score >= 50,  # 50分以上认为有效
                'Abstract_Count': len(abstracts)
            }
            
            # 添加每个功能的独立得分
            for func, score in validation.function_scores.items():
                # 清理功能名称，替换空格为下划线
                clean_func = func.replace(' ', '_').replace('/', '_')
                row[f'Function_Score_{clean_func}'] = score
            
            # 添加详细得分
            row['Score_Direct_Relation'] = validation.detailed_scores.get('direct_relation', 0)
            row['Score_Co_occurrence'] = validation.detailed_scores.get('co_occurrence', 0)
            row['Score_Semantic_Coherence'] = validation.detailed_scores.get('semantic_coherence', 0)
            row['Score_Experimental_Evidence'] = validation.detailed_scores.get('experimental_evidence', 0)
            
            # 添加关键证据（最多3条，用 ||| 分隔）
            evidence_text = ' ||| '.join(validation.key_evidence[:3]) if validation.key_evidence else ''
            row['Key_Evidence'] = evidence_text
            
            # 添加证据摘要
            row['Evidence_Summary'] = self._generate_evidence_summary(validation)
            
            results.append(row)
        
        return pd.DataFrame(results)
    
    def _generate_evidence_summary(self, validation: MTIValidation) -> str:
        """生成证据摘要"""
        summary_parts = []
        
        # 直接关系证据
        if validation.detailed_scores.get('direct_relation', 0) > 70:
            summary_parts.append("Strong direct targeting evidence")
        elif validation.detailed_scores.get('direct_relation', 0) > 40:
            summary_parts.append("Moderate targeting evidence")
        
        # 实验证据
        if validation.detailed_scores.get('experimental_evidence', 0) > 70:
            summary_parts.append("Strong experimental validation")
        elif validation.detailed_scores.get('experimental_evidence', 0) > 40:
            summary_parts.append("Some experimental evidence")
        
        # 功能相关
        high_func_scores = [func for func, score in validation.function_scores.items() if score > 70]
        if high_func_scores:
            summary_parts.append(f"High relevance to: {', '.join(high_func_scores)}")
        
        return '; '.join(summary_parts) if summary_parts else "Limited evidence"

def main():
    parser = argparse.ArgumentParser(description='MTI验证工具 - 基于BERT的miRNA-靶基因相互作用验证')
    parser.add_argument('-i', '--input', help='输入文件路径（单个MTI模式需要）')
    parser.add_argument('-o', '--output', default='mti_validation_results.csv', help='输出CSV文件路径')
    parser.add_argument('-m', '--mirna', help='miRNA名称（单个MTI模式需要）')
    parser.add_argument('-g', '--gene', help='目标基因名称（单个MTI模式需要）')
    parser.add_argument('-f', '--functions', default='', help='功能列表，逗号分隔 (例如: "immune,ferroptosis,steroid regulation")')
    parser.add_argument('--batch', action='store_true', help='批量处理模式')
    parser.add_argument('--batch-file', help='批量处理输入文件(CSV格式: mirna,gene,file_path)')
    
    args = parser.parse_args()
    
    # 验证参数
    if args.batch:
        if not args.batch_file:
            parser.error("批量模式需要提供 --batch-file 参数")
    else:
        if not all([args.input, args.mirna, args.gene]):
            parser.error("单个MTI模式需要提供 -i/--input, -m/--mirna, -g/--gene 参数")
    
    # 初始化验证器
    validator = MTIValidator()
    
    if args.batch and args.batch_file:
        # 批量处理模式
        print("批量处理模式...")
        batch_df = pd.read_csv(args.batch_file)
        abstracts_dict = {}
        
        for _, row in batch_df.iterrows():
            mirna = row['mirna']
            gene = row['gene']
            file_path = row['file_path']
            
            with open(file_path, 'r', encoding='utf-8') as f:
                text = f.read()
            abstracts = validator.parse_abstracts(text)
            abstracts_dict[(mirna, gene)] = abstracts
        
        # 执行批量验证
        functions = [f.strip() for f in args.functions.split(',') if f.strip()]
        results_df = validator.validate_mti_batch(abstracts_dict, functions)
        
        # 保存结果
        results_df.to_csv(args.output, index=False)
        print(f"批量验证完成！结果保存至: {args.output}")
        
        # 打印摘要
        print("\n验证摘要:")
        print(f"总MTI数: {len(results_df)}")
        print(f"有效MTI数: {results_df['Valid'].sum()}")
        print(f"平均得分: {results_df['Overall_Score'].mean():.2f}")
        
        # 生成详细报告
        detailed_report = args.output.replace('.csv', '_detailed_report.txt')
        with open(detailed_report, 'w', encoding='utf-8') as f:
            f.write("="*120 + "\n")
            f.write("MTI批量验证详细报告\n")
            f.write("="*120 + "\n\n")
            
            f.write(f"分析日期: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总MTI数量: {len(results_df)}\n")
            f.write(f"有效MTI数量: {results_df['Valid'].sum()} ({results_df['Valid'].sum()/len(results_df)*100:.1f}%)\n")
            f.write(f"平均总体得分: {results_df['Overall_Score'].mean():.2f}/100\n")
            f.write(f"平均基因关系得分: {results_df['Gene_Relation_Score'].mean():.2f}/100\n\n")
            
            # 按得分排序
            sorted_df = results_df.sort_values('Overall_Score', ascending=False)
            
            # 置信度分布
            f.write("置信度等级分布:\n")
            f.write("-"*120 + "\n")
            conf_dist = sorted_df['Confidence_Level'].value_counts()
            for level, count in conf_dist.items():
                f.write(f"  {level}: {count} ({count/len(sorted_df)*100:.1f}%)\n")
            f.write("\n")
            
            # 高置信度MTI
            high_conf = sorted_df[sorted_df['Confidence_Level'] == 'High Confidence']
            if not high_conf.empty:
                f.write("高置信度MTI:\n")
                f.write("-"*100 + "\n")
                for _, row in high_conf.iterrows():
                    f.write(f"{row['miRNA']} -> {row['Target_Gene']}: ")
                    f.write(f"总分={row['Overall_Score']:.1f}, 基因关系分={row['Gene_Relation_Score']:.1f}\n")
                    if pd.notna(row.get('Key_Evidence', '')):
                        f.write(f"  关键证据: {row['Key_Evidence'][:200]}...\n")
                f.write("\n")
            
            # 功能分析
            if functions:
                f.write("功能相关性分析:\n")
                f.write("-"*100 + "\n")
                for func in functions:
                    clean_func = func.replace(' ', '_').replace('/', '_')
                    col_name = f'Function_Score_{clean_func}'
                    if col_name in sorted_df.columns:
                        high_func = sorted_df[sorted_df[col_name] > 70]
                        if not high_func.empty:
                            f.write(f"\n与'{func}'高度相关的MTI:\n")
                            for _, row in high_func.head(5).iterrows():
                                f.write(f"  {row['miRNA']} -> {row['Target_Gene']}: {row[col_name]:.1f}分\n")
                
                f.write("\n")
            
            # 添加更多详细信息
            f.write("="*120 + "\n")
            f.write("所有MTI详细信息:\n")
            f.write("="*120 + "\n\n")
            
            for idx, row in sorted_df.iterrows():
                f.write(f"\n{idx+1}. {row['miRNA']} → {row['Target_Gene']}\n")
                f.write("-"*80 + "\n")
                f.write(f"总体得分: {row['Overall_Score']:.2f}/100\n")
                f.write(f"基因关系得分: {row['Gene_Relation_Score']:.2f}/100\n")
                f.write(f"置信度: {row['Confidence_Level']}\n")
                f.write(f"分析摘要数: {row['Abstract_Count']}\n")
                
                # 功能得分
                if functions:
                    f.write("\n功能得分:\n")
                    for func in functions:
                        clean_func = func.replace(' ', '_').replace('/', '_')
                        col_name = f'Function_Score_{clean_func}'
                        if col_name in row:
                            f.write(f"  {func}: {row[col_name]:.2f}/100\n")
                
                # 关键证据
                if pd.notna(row.get('Key_Evidence', '')) and row['Key_Evidence']:
                    f.write("\n关键证据:\n")
                    evidences = row['Key_Evidence'].split(' ||| ')
                    for i, evidence in enumerate(evidences, 1):
                        f.write(f"  {i}. {evidence}\n")
                
                f.write("\n")
        
        print(f"详细报告已保存至: {detailed_report}")
        
    else:
        # 单个MTI验证模式
        with open(args.input, 'r', encoding='utf-8') as f:
            text = f.read()
        
        # 解析摘要
        abstracts = validator.parse_abstracts(text)
        print(f"解析到 {len(abstracts)} 个摘要")
        
        # 解析功能
        functions = [f.strip() for f in args.functions.split(',') if f.strip()]
        
        # 执行验证
        print(f"\n验证 {args.mirna} -> {args.gene}")
        validation = validator.validate_mti(abstracts, args.mirna, args.gene, functions)
        
        # 生成详细报告
        report_file = args.output.replace('.csv', '.txt')
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write(f"MTI验证报告: {validation.mirna} -> {validation.target_gene}\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"总体得分: {validation.overall_score:.2f}/100\n")
            f.write(f"基因关系得分: {validation.gene_relation_score:.2f}/100\n")
            f.write(f"置信度等级: {validation.confidence_level}\n")
            f.write(f"证据数量: {validation.evidence_count}\n")
            f.write(f"验证结果: {'有效' if validation.overall_score >= 50 else '证据不足'}\n\n")
            
            f.write("-"*80 + "\n")
            f.write("详细得分:\n")
            f.write("-"*80 + "\n")
            for score_type, score in validation.detailed_scores.items():
                f.write(f"{score_type}: {score:.2f}/100\n")
            
            if validation.function_scores:
                f.write("\n" + "-"*80 + "\n")
                f.write("功能相关性得分:\n")
                f.write("-"*80 + "\n")
                for func, score in validation.function_scores.items():
                    f.write(f"{func}: {score:.2f}/100\n")
            
            if validation.key_evidence:
                f.write("\n" + "-"*80 + "\n")
                f.write("关键证据:\n")
                f.write("-"*80 + "\n")
                for i, evidence in enumerate(validation.key_evidence, 1):
                    f.write(f"{i}. {evidence}\n\n")
        
        print(f"\n详细报告已保存至: {report_file}")
        
        # 保存简化的CSV结果
        result_df = pd.DataFrame([{
            'miRNA': validation.mirna,
            'Target_Gene': validation.target_gene,
            'Overall_Score': validation.overall_score,
            'Gene_Relation_Score': validation.gene_relation_score,
            'Confidence_Level': validation.confidence_level,
            'Valid': validation.overall_score >= 50
        }])
        result_df.to_csv(args.output, index=False)
        print(f"CSV结果已保存至: {args.output}") == validator.validate_mti(abstracts, args.mirna, args.gene, functions)
        
        # 生成详细报告
        report_file = args.output.replace('.csv', '.txt')
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write(f"MTI验证报告: {validation.mirna} -> {validation.target_gene}\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"总体得分: {validation.overall_score:.2f}/100\n")
            f.write(f"基因关系得分: {validation.gene_relation_score:.2f}/100\n")
            f.write(f"置信度等级: {validation.confidence_level}\n")
            f.write(f"证据数量: {validation.evidence_count}\n")
            f.write(f"验证结果: {'有效' if validation.overall_score >= 50 else '证据不足'}\n\n")
            
            f.write("-"*80 + "\n")
            f.write("详细得分:\n")
            f.write("-"*80 + "\n")
            for score_type, score in validation.detailed_scores.items():
                f.write(f"{score_type}: {score:.2f}/100\n")
            
            if validation.function_scores:
                f.write("\n" + "-"*80 + "\n")
                f.write("功能相关性得分:\n")
                f.write("-"*80 + "\n")
                for func, score in validation.function_scores.items():
                    f.write(f"{func}: {score:.2f}/100\n")
            
            if validation.key_evidence:
                f.write("\n" + "-"*80 + "\n")
                f.write("关键证据:\n")
                f.write("-"*80 + "\n")
                for i, evidence in enumerate(validation.key_evidence, 1):
                    f.write(f"{i}. {evidence}\n\n")
        
        print(f"\n详细报告已保存至: {report_file}")
        
        # 保存简化的CSV结果
        result_df = pd.DataFrame([{
            'miRNA': validation.mirna,
            'Target_Gene': validation.target_gene,
            'Overall_Score': validation.overall_score,
            'Gene_Relation_Score': validation.gene_relation_score,
            'Confidence_Level': validation.confidence_level,
            'Valid': validation.overall_score >= 50
        }])
        result_df.to_csv(args.output, index=False)
        print(f"CSV结果已保存至: {args.output}")

if __name__ == "__main__":
    main()
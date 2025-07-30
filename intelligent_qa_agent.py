"""
智能问答Agent - miRNA Research Pipeline
双层架构：问题路由 + 内容检索回答
修复版：支持完整数据分析，解决最高分等问题
"""

import os
import json
import pandas as pd
import glob
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from enum import Enum
import re

from config import PipelineConfig
from utils import Logger, TextProcessor
from utils import check_ollama_connection

class QuestionType(Enum):
    """问题类型枚举"""
    MTI_DISCOVERY = "mti_discovery"           # MTI发现相关
    EXPERIMENTAL_VALIDATION = "validation"    # 实验验证相关
    LITERATURE_EVIDENCE = "literature"       # 文献证据相关
    AI_ANALYSIS = "ai_analysis"              # AI分析相关
    NETWORK_ANALYSIS = "network"             # 网络分析相关
    PIPELINE_STATUS = "pipeline_status"      # 流程状态相关
    COMPARATIVE = "comparative"              # 比较分析相关
    FUNCTIONAL = "functional"                # 功能分析相关
    GENERAL_SUMMARY = "general_summary"      # 综合总结相关
    UNKNOWN = "unknown"                      # 未知类型

@dataclass
class ContentSource:
    """内容源定义"""
    file_path: str
    file_type: str
    description: str
    priority: int = 1

@dataclass
class QuestionContext:
    """问题上下文"""
    question: str
    question_type: QuestionType
    content_sources: List[ContentSource]
    confidence: float
    routing_reasoning: str
    needs_complete_data: bool = False  # 新增：是否需要完整数据分析

class QuestionRoutingAgent:
    """问题路由Agent - 第一层"""
    
    def __init__(self, results_dir: str, logger: Optional[Logger] = None):
        self.results_dir = results_dir
        self.logger = logger
        self.config = PipelineConfig()
        
        # 问题模式匹配规则
        self.question_patterns = {
            QuestionType.MTI_DISCOVERY: [
                r"(多少|how many|number|count).*mti",
                r"发现.*mirna.*基因|discover.*mirna.*gene",
                r"(数量|count|number).*interaction",
                r"哪些.*mirna.*target|which.*mirna.*target",
                r"基因.*mirna.*关系|gene.*mirna.*relationship",
                r"strongest.*possibility|最强.*可能|highest.*score",
                r"top.*mti|best.*mti|最好.*mti",
                r"ferroptosis.*mti|immune.*mti|相关.*mti",
                r"tell me about.*mti|关于.*mti",
                r"which.*mti.*strongest|哪个.*mti.*最强"
            ],
            QuestionType.EXPERIMENTAL_VALIDATION: [
                r"实验.*验证|experimental.*validation",
                r"mirtarbase.*结果|mirtarbase.*result",
                r"(支持|support).*证据|evidence.*support",
                r"验证.*强度|validation.*strength",
                r"实验.*证明|experimental.*proof",
                r"validated.*mti|verified.*interaction",
                r"strongest.*validation|最强.*验证"
            ],
            QuestionType.LITERATURE_EVIDENCE: [
                r"文献.*支持|literature.*evidence|literature.*support",
                r"pubmed.*结果|pubmed.*result",
                r"(多少|how many).*文章|article.*count",
                r"研究.*报告|research.*report",
                r"文献.*挖掘|literature.*mining",
                r"paper.*found|找到.*论文"
            ],
            QuestionType.AI_ANALYSIS: [
                r"bert.*分析|bert.*result|bert.*analysis",
                r"llm.*总结|llm.*analysis|llm.*summary",
                r"ai.*评估|ai.*assessment|ai.*analysis",
                r"机器学习.*结果|machine.*learning.*result",
                r"智能.*分析|intelligent.*analysis",
                r"score|打分|评分",
                r"confidence|置信度"
            ],
            QuestionType.NETWORK_ANALYSIS: [
                r"网络.*分析|network.*analysis",
                r"cytoscape.*结果|cytoscape.*result",
                r"节点.*边|node.*edge",
                r"调控.*网络|regulatory.*network",
                r"网络.*可视化|network.*visualization",
                r"graph|图.*分析"
            ],
            QuestionType.PIPELINE_STATUS: [
                r"流程.*状态|pipeline.*status",
                r"执行.*结果|execution.*result",
                r"错误.*警告|error.*warning",
                r"日志.*信息|log.*information",
                r"运行.*情况|running.*status",
                r"完成.*情况|completion.*status"
            ],
            QuestionType.COMPARATIVE: [
                r"(比较|compare|对比).*数据库|database.*comparison",
                r"(哪个|which).*更好|better|best",
                r"差异.*分析|difference.*analysis|differential.*analysis",
                r"相似.*性|similarity",
                r"(排名|rank|top|highest|strongest|best)",
                r"versus|vs|相比|对比"
            ],
            QuestionType.FUNCTIONAL: [
                r"功能.*分析|functional.*analysis",
                r"通路.*pathway|pathway.*analysis",
                r"生物.*过程|biological.*process",
                r"细胞.*功能|cellular.*function",
                r"分子.*机制|molecular.*mechanism",
                r"ferroptosis|immune|apoptosis|proliferation"
            ],
            QuestionType.GENERAL_SUMMARY: [
                r"总结.*结果|summarize.*result|summary",
                r"整体.*分析|overall.*analysis",
                r"主要.*发现|key.*finding|main.*finding",
                r"结论.*建议|conclusion.*recommendation",
                r"完整.*报告|complete.*report",
                r"tell me about|关于|about",
                r"what.*found|发现了.*什么"
            ]
        }
        
        # 新增：需要完整数据分析的关键词
        self.complete_data_keywords = [
            'highest', 'lowest', 'maximum', 'minimum', 'top', 'best', 'worst',
            'strongest', 'weakest', 'most', 'least', 'rank', 'ranking',
            '最高', '最低', '最大', '最小', '最强', '最弱', '排名', '排行',
            'which.*highest', 'which.*best', 'which.*strongest'
        ]
    
    def identify_question_type(self, question: str) -> Tuple[QuestionType, float, str, bool]:
        """识别问题类型 - 修复版，返回是否需要完整数据分析"""
        question_lower = question.lower()
        
        # 检测是否需要完整数据分析
        needs_complete_data = any(re.search(keyword, question_lower) for keyword in self.complete_data_keywords)
        
        # 计算每种类型的匹配分数
        type_scores = {}
        matched_patterns = {}
        
        for q_type, patterns in self.question_patterns.items():
            score = 0
            matches = []
            
            for pattern in patterns:
                if re.search(pattern, question_lower):
                    score += 1
                    matches.append(pattern)
            
            if score > 0:
                type_scores[q_type] = score / len(patterns)  # 归一化分数
                matched_patterns[q_type] = matches
        
        if not type_scores:
            # Fallback路由
            q_type, confidence, reasoning = self._fallback_routing(question)
            return q_type, confidence, reasoning, needs_complete_data
        
        # 找到最高分的类型
        best_type = max(type_scores.keys(), key=lambda k: type_scores[k])
        confidence = type_scores[best_type]
        
        reasoning = f"Matched {len(matched_patterns[best_type])} patterns for {best_type.value}"
        if needs_complete_data:
            reasoning += " [Requires complete data analysis]"
        
        return best_type, confidence, reasoning, needs_complete_data
    
    def map_type_to_sources(self, question_type: QuestionType) -> List[ContentSource]:
        """将问题类型映射到内容源 - 增强版，优先加载详细数据"""
        sources = []
        
        if question_type == QuestionType.MTI_DISCOVERY:
            # 优先级1：包含详细评分的集成验证结果
            sources.append(
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Integrated validation results with Overall_Score (PRIORITY FOR RANKING)",
                    priority=1
                )
            )
            # 优先级2：BERT验证结果（包含评分）
            sources.append(
                ContentSource(
                    file_path=self._find_file("bert_validation", "mti_validation_results.csv"),
                    file_type="csv", 
                    description="BERT validation results with confidence scores",
                    priority=2
                )
            )
            # 优先级3：MTI选择结果
            sources.append(
                ContentSource(
                    file_path=self._find_file("mirna_selection", "mti_selection_results_*.xlsx"),
                    file_type="excel",
                    description="MTI selection results with database sources",
                    priority=3
                )
            )
        
        elif question_type == QuestionType.EXPERIMENTAL_VALIDATION:
            # 优先加载包含验证信息的详细结果
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Experimental validation data with miRTarBase info",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("mirna_selection", "mti_selection_results_*.xlsx"),
                    file_type="excel",
                    description="MTI selection with validation details",
                    priority=2
                )
            ])
        
        elif question_type == QuestionType.LITERATURE_EVIDENCE:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("pubmed_articles", "mining_summary.csv"),
                    file_type="csv",
                    description="Literature mining summary with article counts",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Integrated results with literature information",
                    priority=2
                ),
                ContentSource(
                    file_path=self._find_directory("pubmed_articles"),
                    file_type="text_directory",
                    description="PubMed abstracts collection",
                    priority=3
                )
            ])
        
        elif question_type == QuestionType.AI_ANALYSIS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Complete AI analysis results with Overall_Score (PRIORITY FOR RANKING)",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("bert_validation", "mti_validation_results.csv"),
                    file_type="csv",
                    description="BERT analysis results (secondary)",
                    priority=2
                ),
                ContentSource(
                    file_path=self._find_directory("llm_summaries"),
                    file_type="json_directory",
                    description="LLM analysis summaries",
                    priority=3
                )
            ])
        
        elif question_type == QuestionType.NETWORK_ANALYSIS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("cytoscape_network", "nodes.csv"),
                    file_type="csv",
                    description="Network nodes with properties",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("cytoscape_network", "edges.csv"),
                    file_type="csv",
                    description="Network edges with scores", 
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="MTI data for network context",
                    priority=2
                )
            ])
        
        elif question_type == QuestionType.PIPELINE_STATUS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file(".", "pipeline_log.txt"),
                    file_type="text",
                    description="Pipeline execution log",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file(".", "pipeline_report.txt"),
                    file_type="text",
                    description="Pipeline summary report",
                    priority=2
                )
            ])
        
        elif question_type in [QuestionType.COMPARATIVE, QuestionType.FUNCTIONAL, QuestionType.GENERAL_SUMMARY]:
            # 这些类型需要详细的多数据源，优先加载包含完整信息的文件
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Complete analysis results with all metrics",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("bert_validation", "mti_validation_results.csv"),
                    file_type="csv",
                    description="BERT validation with detailed scores",
                    priority=2
                ),
                ContentSource(
                    file_path=self._find_file(".", "pipeline_report.txt"),
                    file_type="text",
                    description="Pipeline summary report",
                    priority=3
                ),
                ContentSource(
                    file_path=self._find_file("mirna_selection", "mti_selection_results_*.xlsx"),
                    file_type="excel",
                    description="MTI selection with database details",
                    priority=4
                ),
                ContentSource(
                    file_path=self._find_directory("llm_summaries"),
                    file_type="json_directory",
                    description="Detailed LLM analyses",
                    priority=5
                )
            ])
        
        # 过滤并添加fallback机制
        valid_sources = []
        fallback_sources = []
        
        for source in sources:
            if source.file_path and os.path.exists(source.file_path):
                valid_sources.append(source)
            else:
                fallback_sources.append(source)
        
        # 如果没有找到特定的有效源，尝试查找任何包含MTI数据的文件
        if not valid_sources or len(valid_sources) < 2:
            fallback_files = [
                ("bert_validation/integrated_validation_results.csv", "csv", "Integrated validation results"),
                ("bert_validation/mti_validation_results.csv", "csv", "BERT validation results"),
                ("mirna_selection/mti_selection_results_*.xlsx", "excel", "MTI selection results"),
                ("pipeline_log.txt", "text", "Pipeline execution log"),
                ("pipeline_report.txt", "text", "Pipeline summary report")
            ]
            
            for filename, file_type, description in fallback_files:
                found_file = self._find_file_advanced(filename)
                if found_file and found_file not in [s.file_path for s in valid_sources]:
                    valid_sources.append(ContentSource(
                        file_path=found_file,
                        file_type=file_type,
                        description=f"Fallback: {description}",
                        priority=10
                    ))
        
        # 按优先级排序
        valid_sources.sort(key=lambda x: x.priority)
        
        return valid_sources
    
    def _find_file(self, subdir: str, pattern: str) -> Optional[str]:
        """查找文件 - 增强版，支持多种路径格式"""
        # 尝试多种路径组合
        search_paths = []
        
        # 标准路径
        if subdir != ".":
            search_paths.append(os.path.join(self.results_dir, subdir, pattern))
        
        # 直接在results_dir中搜索
        search_paths.append(os.path.join(self.results_dir, pattern))
        
        # 如果subdir是"."，直接在results_dir中搜索
        if subdir == ".":
            search_paths.append(os.path.join(self.results_dir, pattern))
        
        # 尝试不同的文件模式
        if "*" not in pattern:
            # 如果pattern不包含通配符，添加一些变体
            base_name = pattern
            if "." in base_name:
                name_part, ext_part = base_name.rsplit(".", 1)
                search_paths.extend([
                    os.path.join(self.results_dir, subdir, f"{name_part}*.{ext_part}"),
                    os.path.join(self.results_dir, f"{name_part}*.{ext_part}")
                ])
        
        # 搜索所有可能的路径
        for search_path in search_paths:
            import glob
            matches = glob.glob(search_path)
            if matches:
                # 返回最新的文件（按修改时间排序）
                latest_file = max(matches, key=lambda x: os.path.getmtime(x))
                return latest_file
        
        return None
    
    def _find_file_advanced(self, filename: str) -> Optional[str]:
        """高级文件查找"""
        # 递归搜索整个results_dir
        for root, dirs, files in os.walk(self.results_dir):
            for file in files:
                if filename.replace("*", "") in file:
                    return os.path.join(root, file)
        return None
    
    def _find_directory(self, subdir: str) -> Optional[str]:
        """查找目录 - 增强版"""
        # 尝试多种路径
        possible_paths = [
            os.path.join(self.results_dir, subdir),
            os.path.join(self.results_dir, subdir.replace("_", "-")),  # 处理下划线/连字符变体
            os.path.join(self.results_dir, subdir.lower()),
            os.path.join(self.results_dir, subdir.upper())
        ]
        
        for path in possible_paths:
            if os.path.exists(path) and os.path.isdir(path):
                return path
        
        return None
    
    def route_question(self, question: str) -> QuestionContext:
        """路由问题到相应的内容源 - 修复版"""
        # 识别问题类型（现在返回4个值）
        q_type, confidence, reasoning, needs_complete_data = self.identify_question_type(question)
        
        # 如果识别失败，尝试基于关键词进行fallback路由
        if q_type == QuestionType.UNKNOWN or confidence < 0.3:
            q_type, confidence, reasoning = self._fallback_routing(question)
            # 重新检查是否需要完整数据
            needs_complete_data = any(re.search(keyword, question.lower()) for keyword in self.complete_data_keywords)
        
        # 映射到内容源
        sources = self.map_type_to_sources(q_type)
        
        if self.logger:
            self.logger.info(f"Question routing - Type: {q_type.value}, Confidence: {confidence:.2f}")
            self.logger.info(f"Needs complete data: {needs_complete_data}")
            self.logger.info(f"Found {len(sources)} relevant content sources")
        
        return QuestionContext(
            question=question,
            question_type=q_type,
            content_sources=sources,
            confidence=confidence,
            routing_reasoning=reasoning,
            needs_complete_data=needs_complete_data
        )
    
    def _fallback_routing(self, question: str) -> Tuple[QuestionType, float, str]:
        """Fallback路由机制 - 基于关键词的简单匹配"""
        question_lower = question.lower()
        
        # 关键词映射
        keyword_mapping = {
            QuestionType.MTI_DISCOVERY: [
                'mti', 'interaction', 'mirna', 'gene', 'target', 'discover', 'found',
                '发现', '相互作用', '基因', '目标', 'ferroptosis', 'immune', '多少', 'how many',
                'strongest', 'best', 'top', '最强', '最好'
            ],
            QuestionType.EXPERIMENTAL_VALIDATION: [
                'validation', 'experimental', 'mirtarbase', 'verify', 'evidence',
                '验证', '实验', '证据', '支持'
            ],
            QuestionType.LITERATURE_EVIDENCE: [
                'literature', 'pubmed', 'paper', 'article', 'reference',
                '文献', '论文', '文章', '参考'
            ],
            QuestionType.AI_ANALYSIS: [
                'bert', 'llm', 'ai', 'score', 'analysis', 'confidence',
                '分析', '评分', '置信度', 'ai'
            ],
            QuestionType.NETWORK_ANALYSIS: [
                'network', 'cytoscape', 'node', 'edge', 'graph',
                '网络', '节点', '边', '图'
            ],
            QuestionType.FUNCTIONAL: [
                'function', 'pathway', 'biological', 'cellular', 'ferroptosis', 'immune',
                '功能', '通路', '生物', '细胞', '铁死亡', '免疫'
            ]
        }
        
        best_type = QuestionType.MTI_DISCOVERY  # 默认类型
        best_score = 0
        best_matches = []
        
        # 计算关键词匹配分数
        for q_type, keywords in keyword_mapping.items():
            matches = [kw for kw in keywords if kw in question_lower]
            score = len(matches) / len(keywords)  # 归一化分数
            
            if score > best_score:
                best_score = score
                best_type = q_type
                best_matches = matches
        
        # 至少有一些匹配才返回，否则默认为MTI_DISCOVERY
        if best_score > 0:
            confidence = min(0.6, best_score * 2)  # 限制fallback的confidence
            reasoning = f"Fallback routing based on keywords: {best_matches}"
        else:
            confidence = 0.3  # 最低置信度
            reasoning = "Default routing to MTI_DISCOVERY"
        
        return best_type, confidence, reasoning

class ContentRetrievalAgent:
    """内容检索和回答Agent - 第二层"""
    
    def __init__(self, logger: Optional[Logger] = None):
        self.logger = logger
        self.config = PipelineConfig()
        self.llm = self._initialize_llm()
        
        # 内容缓存
        self.content_cache = {}
    
    def _initialize_llm(self):
        """初始化LLM - 修复版，直接初始化而不依赖其他模块"""
        try:
            is_connected, models = check_ollama_connection()
            if not is_connected:
                if self.logger:
                    self.logger.error("Ollama not available for QA Agent")
                print("⚠️ Warning: Ollama not connected. Please start Ollama: ollama serve")
                return None
            
            if self.logger:
                self.logger.info(f"Ollama connected. Available models: {models}")
            
            # 直接初始化LLM，不依赖其他模块
            try:
                from langchain_community.llms import Ollama
                
                # 尝试使用可用的模型
                preferred_models = ['llama3.1', 'llama3', 'llama2', 'mistral', 'codellama']
                selected_model = None
                
                for model in preferred_models:
                    if any(model in available_model.lower() for available_model in models):
                        selected_model = model
                        break
                
                if not selected_model and models:
                    # 如果没有首选模型，使用第一个可用模型
                    selected_model = models[0].split(':')[0]  # 移除tag部分
                
                if selected_model:
                    llm = Ollama(model=selected_model, temperature=0.1)
                    if self.logger:
                        self.logger.info(f"Initialized LLM with model: {selected_model}")
                    print(f"✅ LLM initialized with model: {selected_model}")
                    return llm
                else:
                    if self.logger:
                        self.logger.error("No suitable LLM models found")
                    print("❌ No suitable LLM models found. Please install: ollama pull llama3.1")
                    return None
                    
            except ImportError:
                # Fallback: 尝试使用其他模块的LLM
                try:
                    from mti_llm_summarize import MTILLMSummarizer
                    summarizer = MTILLMSummarizer()
                    if summarizer.llm:
                        if self.logger:
                            self.logger.info("Using LLM from MTILLMSummarizer")
                        print("✅ Using LLM from existing summarizer")
                        return summarizer.llm
                except ImportError:
                    pass
                
                if self.logger:
                    self.logger.error("Cannot import LangChain Ollama. Please install: pip install langchain-community")
                print("❌ Cannot import LangChain. Please install: pip install langchain-community")
                return None
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to initialize LLM for QA: {e}")
            print(f"❌ LLM initialization failed: {e}")
            return None
    
    def _extract_complete_scores(self, df: pd.DataFrame) -> str:
        """提取完整的分数分析 - 修复版，优先Overall_Score"""
        try:
            # 查找所有分数列
            all_score_columns = [col for col in df.columns if any(keyword in col.lower() 
                                for keyword in ['score', 'confidence', 'overall', 'rating', 'strength'])]
            
            if not all_score_columns:
                return "No score columns found in the data."
            
            # 按优先级排序：Overall_Score优先
            overall_columns = [col for col in all_score_columns if 'overall' in col.lower()]
            other_columns = [col for col in all_score_columns if 'overall' not in col.lower()]
            prioritized_columns = overall_columns + other_columns
            
            analysis_text = f"\n=== 🎯 COMPLETE DATA ANALYSIS ({len(df)} total rows) ===\n"
            analysis_text += f"Found {len(prioritized_columns)} score columns\n"
            
            for i, col in enumerate(prioritized_columns):
                if pd.api.types.is_numeric_dtype(df[col]):
                    try:
                        valid_data = df[col].dropna()
                        if valid_data.empty:
                            continue
                            
                        max_value = valid_data.max()
                        max_idx = valid_data.idxmax()
                        max_row = df.loc[max_idx]
                        
                        # 特别标记Overall_Score
                        if 'overall' in col.lower():
                            marker = "🏆 [MAIN SCORE - USE THIS FOR RANKING]"
                        else:
                            marker = f"📊 [ADDITIONAL SCORE {i+1}]"
                        
                        analysis_text += f"\n{marker} '{col}':\n"
                        analysis_text += f"   📈 MAXIMUM VALUE: {max_value:.6f}\n"
                        
                        # 找miRNA和基因列
                        mirna_col = None
                        gene_col = None
                        for c in df.columns:
                            if 'mirna' in c.lower() and not mirna_col:
                                mirna_col = c
                            elif any(keyword in c.lower() for keyword in ['gene', 'target']) and not gene_col:
                                gene_col = c
                        
                        if mirna_col and gene_col:
                            analysis_text += f"   🎯 WINNING PAIR: {max_row[mirna_col]} ➜ {max_row[gene_col]}\n"
                        
                        # 显示Top 5（减少输出）
                        top_5 = df.nlargest(5, col)
                        analysis_text += f"\n   🏅 TOP 5 for {col}:\n"
                        for j, (idx, row) in enumerate(top_5.iterrows(), 1):
                            if mirna_col and gene_col:
                                analysis_text += f"      {j}. {row[mirna_col]} ➜ {row[gene_col]} | {row[col]:.6f}\n"
                            else:
                                analysis_text += f"      {j}. Row {idx} | {row[col]:.6f}\n"
                        
                        analysis_text += f"   📊 Stats: Mean={valid_data.mean():.4f}, Total={len(valid_data)} rows\n"
                        
                        if 'overall' in col.lower():
                            analysis_text += f"   ⭐ FOR 'HIGHEST SCORE' QUESTIONS: Use this {col} result!\n"
                        
                        analysis_text += "\n" + "-"*60 + "\n"
                        
                    except Exception as e:
                        continue
            
            # 添加明确指导
            analysis_text += "\n🔍 IMPORTANT FOR RANKING QUESTIONS:\n"
            if overall_columns:
                analysis_text += f"✅ Use '{overall_columns[0]}' (marked with 🏆) for main ranking\n"
                analysis_text += "✅ This column represents the comprehensive MTI score\n"
            else:
                analysis_text += "⚠️ No Overall_Score found, using other score columns\n"
            
            return analysis_text
            
        except Exception as e:
            return f"❌ Error: {str(e)}"
    
    def load_content(self, sources: List[ContentSource]) -> Dict[str, Any]:
        """加载相关内容 - 修复版，支持完整数据分析"""
        loaded_content = {}
        
        for source in sorted(sources, key=lambda x: x.priority):
            if source.file_path in self.content_cache:
                loaded_content[source.description] = self.content_cache[source.file_path]
                continue
            
            try:
                if source.file_type == "csv":
                    df = pd.read_csv(source.file_path)
                    content = {
                        'type': 'dataframe',
                        'shape': df.shape,
                        'columns': list(df.columns),
                        'sample_data': df.head().to_dict('records') if len(df) > 0 else [],
                        'summary_stats': self._get_dataframe_summary(df),
                        'complete_score_analysis': self._extract_complete_scores(df)  # 新增
                    }
                    
                elif source.file_type == "excel":
                    df = pd.read_excel(source.file_path, sheet_name=None)  # 读取所有sheet
                    content = {
                        'type': 'excel',
                        'sheets': list(df.keys()),
                        'data': {}
                    }
                    for sheet_name, sheet_df in df.items():
                        content['data'][sheet_name] = {
                            'shape': sheet_df.shape,
                            'columns': list(sheet_df.columns),
                            'sample_data': sheet_df.head().to_dict('records') if len(sheet_df) > 0 else [],
                            'complete_score_analysis': self._extract_complete_scores(sheet_df)  # 新增
                        }
                
                elif source.file_type == "text":
                    with open(source.file_path, 'r', encoding='utf-8') as f:
                        text_content = f.read()
                    content = {
                        'type': 'text',
                        'length': len(text_content),
                        'content': text_content[:5000]  # 限制长度避免token过多
                    }
                
                elif source.file_type == "json_directory":
                    content = self._load_json_directory(source.file_path)
                
                elif source.file_type == "text_directory":
                    content = self._load_text_directory_summary(source.file_path)
                
                else:
                    continue
                
                # 缓存内容
                self.content_cache[source.file_path] = content
                loaded_content[source.description] = content
                
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"Failed to load {source.file_path}: {e}")
                continue
        
        return loaded_content
    
    def _get_dataframe_summary(self, df: pd.DataFrame) -> Dict[str, Any]:
        """获取DataFrame摘要统计"""
        summary = {
            'row_count': len(df),
            'column_count': len(df.columns)
        }
        
        # 数值列统计
        numeric_cols = df.select_dtypes(include=['number']).columns
        if len(numeric_cols) > 0:
            summary['numeric_summary'] = df[numeric_cols].describe().to_dict()
        
        # 分类列统计
        categorical_cols = df.select_dtypes(include=['object']).columns
        if len(categorical_cols) > 0:
            summary['categorical_summary'] = {}
            for col in categorical_cols:
                summary['categorical_summary'][col] = df[col].value_counts().head().to_dict()
        
        return summary
    
    def _load_json_directory(self, dir_path: str) -> Dict[str, Any]:
        """加载JSON目录内容"""
        content = {
            'type': 'json_directory',
            'files': [],
            'summary': {}
        }
        
        json_files = glob.glob(os.path.join(dir_path, "*.json"))
        
        for json_file in json_files[:10]:  # 限制文件数量
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                filename = os.path.basename(json_file)
                content['files'].append({
                    'filename': filename,
                    'keys': list(data.keys()) if isinstance(data, dict) else 'non-dict',
                    'sample': str(data)[:500]  # 样本内容
                })
            except Exception as e:
                continue
        
        content['summary']['total_files'] = len(json_files)
        content['summary']['loaded_files'] = len(content['files'])
        
        return content
    
    def _load_text_directory_summary(self, dir_path: str) -> Dict[str, Any]:
        """加载文本目录摘要"""
        content = {
            'type': 'text_directory',
            'summary': {}
        }
        
        txt_files = glob.glob(os.path.join(dir_path, "*.txt"))
        
        total_size = sum(os.path.getsize(f) for f in txt_files if os.path.exists(f))
        
        content['summary'] = {
            'total_files': len(txt_files),
            'total_size_mb': total_size / (1024 * 1024),
            'avg_file_size_kb': (total_size / len(txt_files) / 1024) if txt_files else 0,
            'sample_files': [os.path.basename(f) for f in txt_files[:5]]
        }
        
        return content
    
    def generate_answer(self, question: str, content: Dict[str, Any], 
                       question_type: QuestionType, needs_complete_data: bool = False) -> str:
        """生成智能回答 - 修复版，支持完整数据分析"""
        print(f"🤖 Generating answer for question: '{question}'")
        print(f"🎯 Question type: {question_type.value}")
        print(f"📊 Content sources: {len(content)}")
        print(f"🔍 Needs complete data analysis: {needs_complete_data}")
        
        if not self.llm:
            error_msg = "Sorry, LLM is not available for answering questions. Please check Ollama connection and model installation."
            print(f"❌ {error_msg}")
            return error_msg
        
        # 构建提示
        prompt = self._build_answer_prompt(question, content, question_type, needs_complete_data)
        
        # 调试：显示prompt的开头部分
        print(f"📝 Prompt preview: {prompt[:200]}...")
        
        try:
            print("🔄 Calling LLM...")
            # 使用LLM生成回答
            response = self.llm.invoke(prompt)
            
            print("✅ LLM response received")
            
            if hasattr(response, 'content'):
                answer = response.content
            else:
                answer = str(response)
            
            # 验证回答不为空
            if not answer or answer.strip() == "":
                answer = f"I received your question '{question}' but generated an empty response. Please try rephrasing your question or check the LLM configuration."
            
            print(f"📄 Answer length: {len(answer)} characters")
            return answer
                
        except Exception as e:
            error_msg = f"I apologize, but I encountered an error while generating the answer: {str(e)}"
            print(f"❌ LLM error: {e}")
            if self.logger:
                self.logger.error(f"Failed to generate answer for question '{question}': {e}")
            return error_msg
    
    def _build_answer_prompt(self, question: str, content: Dict[str, Any], 
                           question_type: QuestionType, needs_complete_data: bool = False) -> str:
        """构建回答提示 - 修复版，明确指导使用正确的分数列"""
        print(f"🔨 Building prompt for question: '{question}'")
        
        prompt_parts = [
            "You are an intelligent research assistant specializing in miRNA analysis.",
            "Based on the provided pipeline analysis results, answer the user's question accurately.",
            "",
            "=== USER QUESTION ===",
            f"Question: {question}",
            f"Question Type: {question_type.value}",
            f"Requires Complete Data Analysis: {needs_complete_data}",
            "",
            "=== AVAILABLE DATA ==="
        ]
        
        if not content:
            prompt_parts.append("No data sources available.")
        else:
            for source_name, source_content in content.items():
                prompt_parts.append(f"\n--- {source_name} ---")
                
                if source_content['type'] == 'dataframe':
                    prompt_parts.append(f"Shape: {source_content['shape']}")
                    prompt_parts.append(f"Columns: {', '.join(source_content['columns'])}")
                    
                    # 优先显示完整分析结果
                    if 'complete_score_analysis' in source_content and source_content['complete_score_analysis']:
                        prompt_parts.append(source_content['complete_score_analysis'])
                    
                    # 样本数据作为补充
                    if source_content['sample_data']:
                        prompt_parts.append("Sample rows (for reference only):")
                        for i, row in enumerate(source_content['sample_data'][:2]):
                            prompt_parts.append(f"  {i+1}. {row}")
                
                elif source_content['type'] == 'excel':
                    prompt_parts.append(f"Excel sheets: {', '.join(source_content['sheets'])}")
                    for sheet_name, sheet_data in source_content['data'].items():
                        prompt_parts.append(f"\nSheet '{sheet_name}' Analysis:")
                        if 'complete_score_analysis' in sheet_data and sheet_data['complete_score_analysis']:
                            prompt_parts.append(sheet_data['complete_score_analysis'])
                
                elif source_content['type'] == 'text':
                    prompt_parts.append(f"Text content: {source_content['content'][:800]}")
                
                elif source_content['type'] == 'json_directory':
                    summary = source_content['summary']
                    prompt_parts.append(f"JSON directory: {summary['total_files']} files")
                
                elif source_content['type'] == 'text_directory':
                    summary = source_content['summary']
                    prompt_parts.append(f"Text directory: {summary['total_files']} files")
        
        # 关键指令
        prompt_parts.extend([
            "",
            "=== CRITICAL INSTRUCTIONS ===",
            f"Answer this question: '{question}'",
            ""
        ])
        
        if needs_complete_data:
            prompt_parts.extend([
                "🚨 CRITICAL for highest/lowest score questions:",
                "1. Use the section marked '🏆 [MAIN SCORE - USE THIS FOR RANKING]'",
                "2. Look for Overall_Score columns - these are the primary ranking criteria",
                "3. Use the 'WINNING PAIR' information from the main score section",
                "4. Ignore experimental evidence or function scores for overall ranking",
                "5. The complete analysis shows ALL data, not just samples",
                "6. When you see '⭐ FOR HIGHEST SCORE QUESTIONS: Use this result!', follow that guidance",
                ""
            ])
        
        prompt_parts.extend([
            "ANSWER REQUIREMENTS:",
            "- Use the complete data analysis above (not sample data)",
            "- For 'highest score' questions, use Overall_Score marked with 🏆",
            "- Provide exact miRNA-gene pair and score value",
            "- Be specific and accurate",
            "- When multiple score columns exist, prioritize the one marked as MAIN SCORE",
            "",
            "=== YOUR ANSWER ==="
        ])
        
        final_prompt = "\n".join(prompt_parts)
        
        # 调试信息
        print(f"📏 Prompt length: {len(final_prompt)} characters")
        print(f"📋 Contains complete data analysis: {needs_complete_data}")
        
        return final_prompt

class IntelligentQAAgent:
    """智能问答Agent - 主控制器"""
    
    def __init__(self, results_dir: str, logger: Optional[Logger] = None):
        self.results_dir = results_dir
        self.logger = logger
        
        # 初始化两层Agent
        self.routing_agent = QuestionRoutingAgent(results_dir, logger)
        self.retrieval_agent = ContentRetrievalAgent(logger)
        
        # 对话历史
        self.conversation_history = []
    
    def ask_question(self, question: str) -> Dict[str, Any]:
        """处理用户问题 - 修复版，支持完整数据分析"""
        if self.logger:
            self.logger.info(f"Processing question: {question}")
        
        # 第一层：问题路由
        context = self.routing_agent.route_question(question)
        
        # 即使没有找到内容源，也尝试提供基本回答
        if not context.content_sources:
            # 尝试从日志文件获取基本信息
            basic_answer = self._generate_basic_answer(question, context.question_type)
            
            return {
                'question': question,
                'answer': basic_answer,
                'confidence': context.confidence,
                'sources_used': [],
                'question_type': context.question_type.value,
                'note': 'Answer based on limited available information',
                'complete_data_analysis': context.needs_complete_data
            }
        
        # 第二层：内容加载和回答生成
        content = self.retrieval_agent.load_content(context.content_sources)
        
        if not content:
            # 如果内容加载失败，也尝试提供基本回答
            basic_answer = self._generate_basic_answer(question, context.question_type)
            
            return {
                'question': question,
                'answer': basic_answer,
                'confidence': context.confidence,
                'sources_used': [s.file_path for s in context.content_sources],
                'question_type': context.question_type.value,
                'note': 'Content loading failed, answer based on limited information',
                'complete_data_analysis': context.needs_complete_data
            }
        
        # 生成回答（传递完整数据需求标志）
        answer = self.retrieval_agent.generate_answer(
            question, content, context.question_type, context.needs_complete_data
        )
        
        # 保存到对话历史
        qa_pair = {
            'question': question,
            'answer': answer,
            'question_type': context.question_type.value,
            'confidence': context.confidence,
            'sources_used': [s.file_path for s in context.content_sources],
            'timestamp': pd.Timestamp.now().isoformat(),
            'complete_data_analysis': context.needs_complete_data
        }
        
        self.conversation_history.append(qa_pair)
        
        return qa_pair
    
    def _generate_basic_answer(self, question: str, question_type: QuestionType) -> str:
        """生成基本回答（当无法获取详细内容时）"""
        basic_answers = {
            QuestionType.MTI_DISCOVERY: 
                "I understand you're asking about miRNA-target interactions (MTIs). To provide accurate numbers and details, I need access to the pipeline analysis results. Please ensure your pipeline has completed successfully and the results files are available.",
            
            QuestionType.EXPERIMENTAL_VALIDATION:
                "You're asking about experimental validation. The pipeline typically uses miRTarBase to provide experimental evidence for MTI relationships. Please check if your analysis included miRTarBase data and completed successfully.",
            
            QuestionType.LITERATURE_EVIDENCE:
                "You're inquiring about literature support. The pipeline mines PubMed for relevant articles supporting each MTI. Please ensure the literature mining step completed successfully in your analysis.",
            
            QuestionType.FUNCTIONAL:
                "I see you're asking about functional analysis. This relates to the biological pathways and functions associated with the miRNAs and genes in your analysis. Please check if the functional analysis components completed successfully.",
            
            QuestionType.AI_ANALYSIS:
                "You're asking about AI analysis results (BERT or LLM analysis). These provide confidence scores and detailed insights into the MTI relationships. Please ensure these analysis steps completed successfully.",
            
            QuestionType.NETWORK_ANALYSIS:
                "You're inquiring about network analysis results. This involves the construction and analysis of miRNA-gene regulatory networks. Please check if the network generation step completed successfully.",
            
            QuestionType.COMPARATIVE:
                "You're asking for comparative analysis. This typically involves comparing results across different databases or conditions. Please ensure your pipeline analysis completed with multiple data sources.",
            
            QuestionType.GENERAL_SUMMARY:
                "You're looking for a summary of results. I'd be happy to provide this once I can access your pipeline analysis results. Please ensure the analysis completed successfully and all result files are available.",
            
            QuestionType.PIPELINE_STATUS:
                "You're asking about pipeline execution status. Please check the pipeline_log.txt file in your results directory for detailed execution information.",
            
            QuestionType.UNKNOWN:
                "I'm trying to understand your question better. Could you please rephrase it or provide more specific details? For example, you could ask about MTI numbers, experimental validation, literature support, or analysis results."
        }
        
        base_answer = basic_answers.get(question_type, basic_answers[QuestionType.UNKNOWN])
        
        # 添加一些通用的帮助信息
        help_text = "\n\n💡 Helpful tips:\n"
        help_text += "• Make sure your pipeline analysis completed successfully\n"
        help_text += "• Check that result files exist in your pipeline output directory\n"
        help_text += "• Try rephrasing your question with specific terms like 'MTI', 'miRNA', 'validation', etc.\n"
        help_text += "• You can ask questions like: 'How many MTIs were found?', 'What is the experimental validation rate?', 'Summarize the main findings'"
        
        return base_answer + help_text
    
    def interactive_session(self):
        """启动交互式问答会话"""
        print("\n" + "="*60)
        print("🤖 Intelligent miRNA Research Assistant")
        print("="*60)
        print("Ask me anything about your pipeline results!")
        print("Type 'quit', 'exit', or 'bye' to end the session.")
        print("Type 'help' to see example questions.")
        print("-"*60)
        
        while True:
            try:
                question = input("\n🔬 Your Question: ").strip()
                
                if question.lower() in ['quit', 'exit', 'bye', 'q']:
                    print("\n👋 Thank you for using the miRNA Research Assistant!")
                    break
                
                if question.lower() == 'help':
                    self._show_example_questions()
                    continue
                
                if not question:
                    print("⚠️ Please enter a question.")
                    continue
                
                print("\n🤖 Processing your question...")
                result = self.ask_question(question)
                
                print(f"\n📊 Answer (Type: {result['question_type']}, Confidence: {result['confidence']:.1%}):")
                if result.get('complete_data_analysis'):
                    print("🔍 Used complete data analysis")
                print("-" * 40)
                print(result['answer'])
                
                if result['sources_used']:
                    print(f"\n📁 Sources used: {len(result['sources_used'])} files")
                
            except KeyboardInterrupt:
                print("\n\n👋 Session interrupted. Goodbye!")
                break
            except Exception as e:
                print(f"\n❌ An error occurred: {e}")
                if self.logger:
                    self.logger.error(f"Interactive session error: {e}")
    
    def _show_example_questions(self):
        """显示示例问题"""
        examples = [
            "How many MTIs were discovered in total?",
            "Which MTI pair has the highest overall score?",  # 改进的问法
            "What are the top 5 MTIs by Overall_Score?",  # 新增
            "Which miRNAs have the strongest experimental validation?",
            "What is the literature support for the top MTIs?",
            "What did the AI analysis reveal about ferroptosis pathways?",
            "Show me the network analysis results",
            "What errors or warnings occurred during the pipeline?",
            "Compare the results from different databases",
            "Summarize the main findings of this analysis",
            "Which miRNA has the highest confidence score?",  # 新增
            "What are the strongest MTI interactions found?"  # 新增
        ]
        
        print("\n💡 Example Questions:")
        print("-" * 30)
        for i, example in enumerate(examples, 1):
            print(f"{i:2d}. {example}")
        
        print("\n🎯 Tips for best results:")
        print("• For ranking questions, use terms like 'highest', 'top', 'strongest'")
        print("• Ask about 'Overall_Score' for comprehensive MTI rankings")
        print("• Be specific about what type of score you want to know about")
    
    def save_conversation_history(self, output_file: Optional[str] = None):
        """保存对话历史"""
        if not self.conversation_history:
            return
        
        if output_file is None:
            output_file = os.path.join(self.results_dir, "qa_conversation_history.json")
        
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(self.conversation_history, f, ensure_ascii=False, indent=2)
            
            if self.logger:
                self.logger.info(f"Conversation history saved to {output_file}")
        
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to save conversation history: {e}")
    
    def get_agent_status(self) -> Dict[str, Any]:
        """获取Agent状态"""
        return {
            'results_dir': self.results_dir,
            'llm_available': self.retrieval_agent.llm is not None,
            'conversation_count': len(self.conversation_history),
            'content_cache_size': len(self.retrieval_agent.content_cache),
            'supported_question_types': [qt.value for qt in QuestionType]
        }
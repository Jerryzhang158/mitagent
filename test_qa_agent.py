"""
智能问答Agent - miRNA Research Pipeline
双层架构：问题路由 + 内容检索回答
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

class QuestionRoutingAgent:
    """问题路由Agent - 第一层"""
    
    def __init__(self, results_dir: str, logger: Optional[Logger] = None):
        self.results_dir = results_dir
        self.logger = logger
        self.config = PipelineConfig()
        
        # 问题模式匹配规则
        self.question_patterns = {
            QuestionType.MTI_DISCOVERY: [
                r"(多少|how many).*mti",
                r"发现.*mirna.*基因",
                r"(数量|count|number).*interaction",
                r"哪些.*mirna.*target",
                r"基因.*mirna.*关系"
            ],
            QuestionType.EXPERIMENTAL_VALIDATION: [
                r"实验.*验证|experimental.*validation",
                r"mirtarbase.*结果",
                r"(支持|support).*证据",
                r"验证.*强度|validation.*strength",
                r"实验.*证明"
            ],
            QuestionType.LITERATURE_EVIDENCE: [
                r"文献.*支持|literature.*evidence",
                r"pubmed.*结果",
                r"(多少|how many).*文章|article",
                r"研究.*报告",
                r"文献.*挖掘"
            ],
            QuestionType.AI_ANALYSIS: [
                r"bert.*分析|bert.*result",
                r"llm.*总结|llm.*analysis",
                r"ai.*评估|ai.*assessment",
                r"机器学习.*结果",
                r"智能.*分析"
            ],
            QuestionType.NETWORK_ANALYSIS: [
                r"网络.*分析|network.*analysis",
                r"cytoscape.*结果",
                r"节点.*边|node.*edge",
                r"调控.*网络",
                r"网络.*可视化"
            ],
            QuestionType.PIPELINE_STATUS: [
                r"流程.*状态|pipeline.*status",
                r"执行.*结果",
                r"错误.*警告|error.*warning",
                r"日志.*信息",
                r"运行.*情况"
            ],
            QuestionType.COMPARATIVE: [
                r"(比较|compare|对比).*数据库",
                r"(哪个|which).*更好|better",
                r"差异.*分析|difference.*analysis",
                r"相似.*性|similarity",
                r"(排名|rank|top)"
            ],
            QuestionType.FUNCTIONAL: [
                r"功能.*分析|functional.*analysis",
                r"通路.*pathway",
                r"生物.*过程|biological.*process",
                r"细胞.*功能",
                r"分子.*机制"
            ],
            QuestionType.GENERAL_SUMMARY: [
                r"总结.*结果|summarize.*result",
                r"整体.*分析|overall.*analysis",
                r"主要.*发现|key.*finding",
                r"结论.*建议|conclusion.*recommendation",
                r"完整.*报告"
            ]
        }
    
    def identify_question_type(self, question: str) -> Tuple[QuestionType, float, str]:
        """识别问题类型"""
        question_lower = question.lower()
        
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
            return QuestionType.UNKNOWN, 0.0, "No patterns matched"
        
        # 找到最高分的类型
        best_type = max(type_scores.keys(), key=lambda k: type_scores[k])
        confidence = type_scores[best_type]
        
        reasoning = f"Matched {len(matched_patterns[best_type])} patterns for {best_type.value}"
        
        return best_type, confidence, reasoning
    
    def map_type_to_sources(self, question_type: QuestionType) -> List[ContentSource]:
        """将问题类型映射到内容源"""
        sources = []
        
        if question_type == QuestionType.MTI_DISCOVERY:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("mirna_selection", "mti_selection_results_*.xlsx"),
                    file_type="excel",
                    description="MTI selection results",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Integrated validation results",
                    priority=2
                )
            ])
        
        elif question_type == QuestionType.EXPERIMENTAL_VALIDATION:
            sources.append(
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Experimental validation data",
                    priority=1
                )
            )
        
        elif question_type == QuestionType.LITERATURE_EVIDENCE:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("pubmed_articles", "mining_summary.csv"),
                    file_type="csv",
                    description="Literature mining summary",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_directory("pubmed_articles"),
                    file_type="text_directory",
                    description="PubMed abstracts collection",
                    priority=2
                )
            ])
        
        elif question_type == QuestionType.AI_ANALYSIS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="BERT analysis results",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_directory("llm_summaries"),
                    file_type="json_directory",
                    description="LLM analysis summaries",
                    priority=1
                )
            ])
        
        elif question_type == QuestionType.NETWORK_ANALYSIS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("cytoscape_network", "nodes.csv"),
                    file_type="csv",
                    description="Network nodes",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file("cytoscape_network", "edges.csv"),
                    file_type="csv",
                    description="Network edges", 
                    priority=1
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
            # 这些类型需要多个数据源
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Complete analysis results",
                    priority=1
                ),
                ContentSource(
                    file_path=self._find_file(".", "pipeline_report.txt"),
                    file_type="text",
                    description="Pipeline summary report",
                    priority=2
                ),
                ContentSource(
                    file_path=self._find_directory("llm_summaries"),
                    file_type="json_directory",
                    description="Detailed LLM analyses",
                    priority=3
                )
            ])
        
        # 过滤掉不存在的文件
        valid_sources = [s for s in sources if s.file_path and os.path.exists(s.file_path)]
        
        return valid_sources
    
    def _find_file(self, subdir: str, pattern: str) -> Optional[str]:
        """查找文件"""
        search_path = os.path.join(self.results_dir, subdir, pattern)
        matches = glob.glob(search_path)
        
        if matches:
            return matches[0]  # 返回第一个匹配的文件
        
        return None
    
    def _find_directory(self, subdir: str) -> Optional[str]:
        """查找目录"""
        dir_path = os.path.join(self.results_dir, subdir)
        if os.path.exists(dir_path) and os.path.isdir(dir_path):
            return dir_path
        return None
    
    def route_question(self, question: str) -> QuestionContext:
        """路由问题到相应的内容源"""
        # 识别问题类型
        q_type, confidence, reasoning = self.identify_question_type(question)
        
        # 映射到内容源
        sources = self.map_type_to_sources(q_type)
        
        if self.logger:
            self.logger.info(f"Question routing - Type: {q_type.value}, Confidence: {confidence:.2f}")
            self.logger.info(f"Found {len(sources)} relevant content sources")
        
        return QuestionContext(
            question=question,
            question_type=q_type,
            content_sources=sources,
            confidence=confidence,
            routing_reasoning=reasoning
        )

class ContentRetrievalAgent:
    """内容检索和回答Agent - 第二层"""
    
    def __init__(self, logger: Optional[Logger] = None):
        self.logger = logger
        self.config = PipelineConfig()
        self.llm = self._initialize_llm()
        
        # 内容缓存
        self.content_cache = {}
    
    def _initialize_llm(self):
        """初始化LLM"""
        try:
            is_connected, models = check_ollama_connection()
            if not is_connected:
                if self.logger:
                    self.logger.error("Ollama not available for QA Agent")
                return None
            
            # 使用与其他模块相同的LLM初始化方式
            from mti_llm_summarize import MTILLMSummarizer
            summarizer = MTILLMSummarizer()
            return summarizer.llm
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to initialize LLM for QA: {e}")
            return None
    
    def load_content(self, sources: List[ContentSource]) -> Dict[str, Any]:
        """加载相关内容"""
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
                        'summary_stats': self._get_dataframe_summary(df)
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
                            'sample_data': sheet_df.head().to_dict('records') if len(sheet_df) > 0 else []
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
                       question_type: QuestionType) -> str:
        """生成智能回答"""
        if not self.llm:
            return "Sorry, LLM is not available for answering questions."
        
        # 构建提示
        prompt = self._build_answer_prompt(question, content, question_type)
        
        try:
            # 使用LLM生成回答
            response = self.llm.invoke(prompt)
            
            if hasattr(response, 'content'):
                return response.content
            else:
                return str(response)
                
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to generate answer: {e}")
            return f"I apologize, but I encountered an error while generating the answer: {str(e)}"
    
    def _build_answer_prompt(self, question: str, content: Dict[str, Any], 
                           question_type: QuestionType) -> str:
        """构建回答提示"""
        prompt_parts = [
            "You are an intelligent research assistant specializing in miRNA analysis.",
            "Based on the provided pipeline analysis results, please answer the user's question accurately and comprehensively.",
            "",
            f"User Question: {question}",
            f"Question Type: {question_type.value}",
            "",
            "Available Data Sources:"
        ]
        
        # 添加内容摘要
        for source_name, source_content in content.items():
            prompt_parts.append(f"\n--- {source_name} ---")
            
            if source_content['type'] == 'dataframe':
                prompt_parts.append(f"Data shape: {source_content['shape']}")
                prompt_parts.append(f"Columns: {', '.join(source_content['columns'])}")
                if source_content['sample_data']:
                    prompt_parts.append("Sample data:")
                    for i, row in enumerate(source_content['sample_data'][:3]):
                        prompt_parts.append(f"  Row {i+1}: {row}")
                if 'summary_stats' in source_content:
                    prompt_parts.append(f"Summary statistics: {source_content['summary_stats']}")
            
            elif source_content['type'] == 'excel':
                prompt_parts.append(f"Excel sheets: {', '.join(source_content['sheets'])}")
                for sheet_name, sheet_data in source_content['data'].items():
                    prompt_parts.append(f"  {sheet_name}: {sheet_data['shape']} shape, columns: {', '.join(sheet_data['columns'])}")
            
            elif source_content['type'] == 'text':
                prompt_parts.append(f"Text length: {source_content['length']} characters")
                prompt_parts.append(f"Content preview: {source_content['content'][:1000]}")
            
            elif source_content['type'] == 'json_directory':
                prompt_parts.append(f"JSON directory with {source_content['summary']['total_files']} files")
                for file_info in source_content['files'][:3]:
                    prompt_parts.append(f"  {file_info['filename']}: {file_info['sample'][:200]}")
            
            elif source_content['type'] == 'text_directory':
                summary = source_content['summary']
                prompt_parts.append(f"Text directory: {summary['total_files']} files, {summary['total_size_mb']:.1f} MB total")
        
        # 添加回答指导
        prompt_parts.extend([
            "",
            "Instructions:",
            "1. Answer the question directly and accurately based on the provided data",
            "2. Use specific numbers and facts from the analysis results",
            "3. If the question cannot be fully answered with the available data, explain what information is missing",
            "4. Provide context and interpretation to help the user understand the significance",
            "5. Format your response clearly with appropriate structure",
            "6. If relevant, suggest follow-up questions or additional analyses",
            "",
            "Answer:"
        ])
        
        return "\n".join(prompt_parts)

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
        """处理用户问题"""
        if self.logger:
            self.logger.info(f"Processing question: {question}")
        
        # 第一层：问题路由
        context = self.routing_agent.route_question(question)
        
        if not context.content_sources:
            return {
                'question': question,
                'answer': "I couldn't find relevant data sources to answer your question. Please make sure the pipeline has completed successfully and try rephrasing your question.",
                'confidence': 0.0,
                'sources_used': [],
                'question_type': context.question_type.value
            }
        
        # 第二层：内容加载和回答生成
        content = self.retrieval_agent.load_content(context.content_sources)
        
        if not content:
            return {
                'question': question,
                'answer': "I found relevant data sources but couldn't load the content. Please check if the pipeline results are accessible.",
                'confidence': context.confidence,
                'sources_used': [s.file_path for s in context.content_sources],
                'question_type': context.question_type.value
            }
        
        # 生成回答
        answer = self.retrieval_agent.generate_answer(question, content, context.question_type)
        
        # 保存到对话历史
        qa_pair = {
            'question': question,
            'answer': answer,
            'question_type': context.question_type.value,
            'confidence': context.confidence,
            'sources_used': [s.file_path for s in context.content_sources],
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        self.conversation_history.append(qa_pair)
        
        return qa_pair
    
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
            "Which miRNAs have the strongest experimental validation?",
            "What is the literature support for the top MTIs?",
            "What did the AI analysis reveal about ferroptosis pathways?",
            "Show me the network analysis results",
            "What errors or warnings occurred during the pipeline?",
            "Compare the results from different databases",
            "Summarize the main findings of this analysis"
        ]
        
        print("\n💡 Example Questions:")
        print("-" * 30)
        for i, example in enumerate(examples, 1):
            print(f"{i}. {example}")
    
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
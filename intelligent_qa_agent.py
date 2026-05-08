"""
Optimized Intelligent Q&A Agent - miRNA Research Pipeline
Major optimizations:
1. LLM-driven question classification
2. Smart token management to avoid input length limits
3. Layered data loading strategy
4. Optimized MTI discovery logic
"""

import os
import json
import pandas as pd
import glob
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from enum import Enum
import re
import numpy as np

from config import PipelineConfig
from utils import Logger, TextProcessor
from utils import check_ollama_connection

# Token counter (for estimating LLM input length)
def count_tokens(text: str, model: str = "gpt-3.5-turbo") -> int:
    """Estimate token count for text"""
    try:
        import tiktoken
        encoding = tiktoken.encoding_for_model(model)
        return len(encoding.encode(str(text)))
    except:
        # Fallback: rough estimation (1 token ≈ 4 characters)
        return len(str(text)) // 4

class QuestionType(Enum):
    """Question type enumeration"""
    MTI_DISCOVERY = "mti_discovery"           
    EXPERIMENTAL_VALIDATION = "validation"    
    LITERATURE_EVIDENCE = "literature"       
    AI_ANALYSIS = "ai_analysis"              
    NETWORK_ANALYSIS = "network"             
    PIPELINE_STATUS = "pipeline_status"      
    COMPARATIVE = "comparative"              
    FUNCTIONAL = "functional"                
    GENERAL_SUMMARY = "general_summary"      
    UNKNOWN = "unknown"                      

@dataclass
class ContentSource:
    """Content source definition"""
    file_path: str
    file_type: str
    description: str
    priority: int = 1
    max_tokens: int = 2000  # New: max token limit per source

@dataclass
class QuestionContext:
    """Question context"""
    question: str
    question_type: QuestionType
    content_sources: List[ContentSource]
    confidence: float
    reasoning: str
    needs_complete_data: bool = False
    complexity_level: int = 1  # New: complexity level (1-3)

class OptimizedQuestionRoutingAgent:
    """Optimized Question Routing Agent - Uses LLM for classification"""
    
    def __init__(self, results_dir: str, logger: Optional[Logger] = None):
        self.results_dir = results_dir
        self.logger = logger
        self.config = PipelineConfig()
        self.llm = self._initialize_classification_llm()
        
        # Token limit configuration
        self.max_total_tokens = 6000  # Total input limit
        self.classification_tokens = 1000  # Reserved tokens for classification
        self.content_tokens = self.max_total_tokens - self.classification_tokens
        
        # Question type descriptions (for LLM classification)
        self.type_descriptions = {
            QuestionType.MTI_DISCOVERY: "Questions about discovering miRNA-target interactions, finding MTI pairs, counting interactions, asking about strongest/highest score MTIs",
            QuestionType.EXPERIMENTAL_VALIDATION: "Questions about experimental validation, miRTarBase results, validation strength, experimental evidence",
            QuestionType.LITERATURE_EVIDENCE: "Questions about literature support, PubMed results, research papers, publication evidence",
            QuestionType.AI_ANALYSIS: "Questions about AI analysis results, BERT scores, LLM analysis, confidence scores, AI assessments",
            QuestionType.NETWORK_ANALYSIS: "Questions about network analysis, Cytoscape results, nodes, edges, network visualization",
            QuestionType.PIPELINE_STATUS: "Questions about pipeline execution, logs, errors, warnings, completion status",
            QuestionType.COMPARATIVE: "Questions comparing databases, asking which is better, differential analysis, rankings",
            QuestionType.FUNCTIONAL: "Questions about biological functions, pathways, cellular processes, molecular mechanisms",
            QuestionType.GENERAL_SUMMARY: "Questions asking for summaries, overall findings, conclusions, comprehensive reports"
        }
    
    def _initialize_classification_llm(self):
        """Initialize LLM for classification"""
        try:
            is_connected, models = check_ollama_connection()
            if not is_connected:
                if self.logger:
                    self.logger.warning("Ollama not available for classification")
                return None
            
            from langchain_community.llms import Ollama
            
            # Choose smaller model for classification (better speed)
            preferred_models = ['qwen3.5:9b', 'llama3:8b', 'qwen3:8b', 'mistral:7b']
            selected_model = None
            
            for model in preferred_models:
                if any(model.split(':')[0] in available_model.lower() for available_model in models):
                    selected_model = model
                    break
            
            if not selected_model and models:
                selected_model = models[0]
            
            if selected_model:
                llm = Ollama(model=selected_model, temperature=0.0)  # Low temperature for consistent classification
                if self.logger:
                    self.logger.info(f"Classification LLM initialized: {selected_model}")
                return llm
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to initialize classification LLM: {e}")
        
        return None
    
    def classify_question_with_llm(self, question: str) -> Tuple[QuestionType, float, str, bool, int]:
        """Use LLM for question classification"""
        if not self.llm:
            # Fallback to original regex method
            return self._fallback_classification(question)
        
        # Build classification prompt
        classification_prompt = self._build_classification_prompt(question)
        
        try:
            response = self.llm.invoke(classification_prompt)
            result = self._parse_classification_response(response, question)
            return result
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"LLM classification failed: {e}")
            return self._fallback_classification(question)
    
    def _build_classification_prompt(self, question: str) -> str:
        """Build classification prompt"""
        prompt_parts = [
            "You are a miRNA research question classifier. Classify the user's question into one of the predefined types.",
            "",
            "QUESTION TYPES AND DESCRIPTIONS:",
        ]
        
        for q_type, description in self.type_descriptions.items():
            prompt_parts.append(f"- {q_type.value}: {description}")
        
        prompt_parts.extend([
            "",
            "ANALYSIS CRITERIA:",
            "- needs_complete_data: True if question asks for 'highest', 'lowest', 'top', 'best', 'strongest', 'ranking'",
            "- complexity_level: 1=simple facts, 2=analysis needed, 3=comprehensive research",
            "",
            f"USER QUESTION: {question}",
            "",
            "Respond with EXACTLY this JSON format:",
            "{",
            '  "question_type": "type_name",',
            '  "confidence": 0.95,',
            '  "reasoning": "Brief explanation",',
            '  "needs_complete_data": true/false,',
            '  "complexity_level": 1-3',
            "}"
        ])
        
        return "\n".join(prompt_parts)
    
    def _parse_classification_response(self, response: str, question: str) -> Tuple[QuestionType, float, str, bool, int]:
        """Parse LLM classification response"""
        try:
            # Clean response
            if hasattr(response, 'content'):
                response_text = response.content
            else:
                response_text = str(response)
            
            # Extract JSON
            json_start = response_text.find('{')
            json_end = response_text.rfind('}') + 1
            if json_start >= 0 and json_end > json_start:
                json_text = response_text[json_start:json_end]
                result = json.loads(json_text)
                
                # Validate and convert result
                q_type_str = result.get('question_type', 'unknown')
                try:
                    q_type = QuestionType(q_type_str)
                except ValueError:
                    q_type = QuestionType.MTI_DISCOVERY  # Default type
                
                confidence = float(result.get('confidence', 0.5))
                reasoning = result.get('reasoning', 'LLM classification')
                needs_complete_data = result.get('needs_complete_data', False)
                complexity_level = int(result.get('complexity_level', 1))
                
                return q_type, confidence, reasoning, needs_complete_data, complexity_level
                
        except Exception as e:
            if self.logger:
                self.logger.warning(f"Failed to parse LLM classification: {e}")
        
        # Fallback
        return self._fallback_classification(question)
    
    def _fallback_classification(self, question: str) -> Tuple[QuestionType, float, str, bool, int]:
        """Backup classification method (keyword-based)"""
        question_lower = question.lower()
        
        # Detect if complete data needed
        complete_data_keywords = [
            'highest', 'lowest', 'maximum', 'minimum', 'top', 'best', 'worst',
            'strongest', 'weakest', 'most', 'least', 'rank', 'ranking'
        ]
        needs_complete_data = any(keyword in question_lower for keyword in complete_data_keywords)
        
        # Simple keyword-based classification
        if any(kw in question_lower for kw in ['mti', 'interaction', 'target', 'mirna', 'gene']):
            return QuestionType.MTI_DISCOVERY, 0.7, "Keyword-based classification", needs_complete_data, 2
        elif any(kw in question_lower for kw in ['validation', 'experimental', 'mirtarbase']):
            return QuestionType.EXPERIMENTAL_VALIDATION, 0.6, "Keyword-based classification", needs_complete_data, 2
        else:
            return QuestionType.MTI_DISCOVERY, 0.5, "Default classification", needs_complete_data, 1
    
    def route_question(self, question: str) -> QuestionContext:
        """Route question to appropriate content sources"""
        # Use LLM for classification
        q_type, confidence, reasoning, needs_complete_data, complexity_level = self.classify_question_with_llm(question)
        
        # Adjust content source strategy based on complexity
        sources = self._map_type_to_optimized_sources(q_type, complexity_level, needs_complete_data)
        
        if self.logger:
            self.logger.info(f"Question classified as {q_type.value} (confidence: {confidence:.2f}, complexity: {complexity_level})")
        
        return QuestionContext(
            question=question,
            question_type=q_type,
            content_sources=sources,
            confidence=confidence,
            reasoning=reasoning,
            needs_complete_data=needs_complete_data,
            complexity_level=complexity_level
        )
    
    def _map_type_to_optimized_sources(self, question_type: QuestionType, 
                                      complexity_level: int, needs_complete_data: bool) -> List[ContentSource]:
        """Optimized content source mapping - adjust based on complexity and data needs"""
        sources = []
        
        # Adjust token allocation based on complexity
        if complexity_level == 1:  # Simple questions
            max_tokens_per_source = 1500
        elif complexity_level == 2:  # Medium complexity
            max_tokens_per_source = 2000
        else:  # High complexity
            max_tokens_per_source = 2500
        
        if question_type == QuestionType.MTI_DISCOVERY:
            # Optimized MTI discovery logic
            if needs_complete_data:
                # When complete data needed, prioritize files with Overall_Score
                sources.append(ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Complete MTI analysis with Overall_Score (PRIMARY FOR RANKING)",
                    priority=1,
                    max_tokens=max_tokens_per_source
                ))
            else:
                # When complete data not needed, load summary info
                sources.append(ContentSource(
                    file_path=self._find_file("bert_validation", "mti_validation_results.csv"),
                    file_type="csv",
                    description="MTI validation summary",
                    priority=1,
                    max_tokens=max_tokens_per_source // 2
                ))
            
            # Add additional sources based on complexity
            if complexity_level >= 2:
                sources.append(ContentSource(
                    file_path=self._find_file("mirna_selection", "mti_selection_results_*.xlsx"),
                    file_type="excel",
                    description="MTI selection details",
                    priority=2,
                    max_tokens=max_tokens_per_source // 2
                ))
        
        elif question_type == QuestionType.EXPERIMENTAL_VALIDATION:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="Experimental validation data",
                    priority=1,
                    max_tokens=max_tokens_per_source
                )
            ])
        
        elif question_type == QuestionType.AI_ANALYSIS:
            sources.extend([
                ContentSource(
                    file_path=self._find_file("bert_validation", "integrated_validation_results.csv"),
                    file_type="csv",
                    description="AI analysis results with comprehensive scores",
                    priority=1,
                    max_tokens=max_tokens_per_source
                )
            ])
            
            if complexity_level >= 3:
                sources.append(ContentSource(
                    file_path=self._find_directory("llm_summaries"),
                    file_type="json_directory",
                    description="Detailed LLM analysis summaries",
                    priority=2,
                    max_tokens=max_tokens_per_source // 2
                ))
        
        # Other type mappings...
        
        # Filter valid sources
        valid_sources = [s for s in sources if s.file_path and os.path.exists(s.file_path)]
        return sorted(valid_sources, key=lambda x: x.priority)
    
    def _find_file(self, subdir: str, pattern: str) -> Optional[str]:
        """Find file - reuse original logic"""
        search_paths = []
        
        if subdir != ".":
            search_paths.append(os.path.join(self.results_dir, subdir, pattern))
        search_paths.append(os.path.join(self.results_dir, pattern))
        
        for search_path in search_paths:
            matches = glob.glob(search_path)
            if matches:
                return max(matches, key=lambda x: os.path.getmtime(x))
        return None
    
    def _find_directory(self, subdir: str) -> Optional[str]:
        """Find directory"""
        possible_paths = [
            os.path.join(self.results_dir, subdir),
            os.path.join(self.results_dir, subdir.lower()),
            os.path.join(self.results_dir, subdir.upper())
        ]
        
        for path in possible_paths:
            if os.path.exists(path) and os.path.isdir(path):
                return path
        return None

class OptimizedContentRetrievalAgent:
    """Optimized Content Retrieval Agent - Smart token management"""
    
    def __init__(self, logger: Optional[Logger] = None):
        self.logger = logger
        self.config = PipelineConfig()
        self.llm = self._initialize_llm()
        self.content_cache = {}
        
        # Token management configuration
        self.max_context_tokens = 6000
        self.answer_tokens_reserve = 1500
        self.available_content_tokens = self.max_context_tokens - self.answer_tokens_reserve
    
    def _initialize_llm(self):
        """Initialize LLM - reuse original logic"""
        try:
            is_connected, models = check_ollama_connection()
            if not is_connected:
                return None
            
            from langchain_community.llms import Ollama
            preferred_models = ['qwen3.5:9b', 'llama3.1', 'llama3', 'mistral']
            
            for model in preferred_models:
                if any(model.split(':')[0] in available_model.lower() for available_model in models):
                    return Ollama(model=model, temperature=0.1)
            
            if models:
                return Ollama(model=models[0], temperature=0.1)
                
        except Exception as e:
            if self.logger:
                self.logger.error(f"Failed to initialize LLM: {e}")
        
        return None
    
    def load_content_optimized(self, sources: List[ContentSource], 
                             needs_complete_data: bool = False) -> Dict[str, Any]:
        """Optimized content loading - smart token management"""
        loaded_content = {}
        total_tokens_used = 0
        
        for source in sorted(sources, key=lambda x: x.priority):
            if total_tokens_used >= self.available_content_tokens:
                break
            
            remaining_tokens = min(source.max_tokens, 
                                 self.available_content_tokens - total_tokens_used)
            
            try:
                content = self._load_single_source_optimized(source, remaining_tokens, needs_complete_data)
                if content:
                    loaded_content[source.description] = content
                    # Estimate actual token usage
                    content_tokens = self._estimate_content_tokens(content)
                    total_tokens_used += content_tokens
                    
            except Exception as e:
                if self.logger:
                    self.logger.warning(f"Failed to load {source.file_path}: {e}")
                continue
        
        if self.logger:
            self.logger.info(f"Loaded content using ~{total_tokens_used} tokens")
        
        return loaded_content
    
    def _load_single_source_optimized(self, source: ContentSource, 
                                    max_tokens: int, needs_complete_data: bool) -> Optional[Dict[str, Any]]:
        """Optimized single source loading"""
        if source.file_path in self.content_cache:
            cached_content = self.content_cache[source.file_path]
            return self._truncate_content(cached_content, max_tokens, needs_complete_data)
        
        try:
            if source.file_type == "csv":
                df = pd.read_csv(source.file_path)
                content = self._process_csv_optimized(df, max_tokens, needs_complete_data)
                
            elif source.file_type == "excel":
                df_dict = pd.read_excel(source.file_path, sheet_name=None)
                content = self._process_excel_optimized(df_dict, max_tokens, needs_complete_data)
                
            elif source.file_type == "text":
                with open(source.file_path, 'r', encoding='utf-8') as f:
                    text_content = f.read()
                content = self._process_text_optimized(text_content, max_tokens)
                
            else:
                return None
            
            # Cache complete content
            self.content_cache[source.file_path] = content
            return content
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error loading {source.file_path}: {e}")
            return None
    
    def _process_csv_optimized(self, df: pd.DataFrame, 
                             max_tokens: int, needs_complete_data: bool) -> Dict[str, Any]:
        """Optimized CSV processing"""
        content = {
            'type': 'dataframe',
            'shape': df.shape,
            'columns': list(df.columns)
        }
        
        if needs_complete_data:
            # When complete data needed, provide full score analysis
            content['complete_score_analysis'] = self._extract_complete_scores_optimized(df, max_tokens)
            content['sample_data'] = df.head(3).to_dict('records')  # Reduce sample data
        else:
            # When complete data not needed, only provide summary
            content['summary_stats'] = self._get_dataframe_summary_optimized(df)
            content['sample_data'] = df.head(5).to_dict('records')
        
        return content
    
    def _extract_complete_scores_optimized(self, df: pd.DataFrame, max_tokens: int) -> str:
        """Optimized complete score extraction - token control"""
        try:
            score_columns = [col for col in df.columns if any(keyword in col.lower() 
                            for keyword in ['score', 'confidence', 'overall', 'rating'])]
            
            if not score_columns:
                return "No score columns found."
            
            # Prioritize Overall_Score
            overall_columns = [col for col in score_columns if 'overall' in col.lower()]
            other_columns = [col for col in score_columns if 'overall' not in col.lower()]
            prioritized_columns = overall_columns + other_columns
            
            analysis_parts = [f"\n=== COMPLETE DATA ANALYSIS ({len(df)} total rows) ==="]
            
            for i, col in enumerate(prioritized_columns[:3]):  # Limit number of columns
                if pd.api.types.is_numeric_dtype(df[col]):
                    valid_data = df[col].dropna()
                    if valid_data.empty:
                        continue
                    
                    max_value = valid_data.max()
                    max_idx = valid_data.idxmax()
                    max_row = df.loc[max_idx]
                    
                    marker = "[MAIN SCORE]" if 'overall' in col.lower() else f"[SCORE {i+1}]"
                    
                    analysis_parts.append(f"\n{marker} '{col}':")
                    analysis_parts.append(f"   MAXIMUM: {max_value:.6f}")
                    
                    # Find miRNA and gene columns
                    mirna_col = next((c for c in df.columns if 'mirna' in c.lower()), None)
                    gene_col = next((c for c in df.columns if any(kw in c.lower() for kw in ['gene', 'target'])), None)
                    
                    if mirna_col and gene_col:
                        analysis_parts.append(f"   TOP PAIR: {max_row[mirna_col]} -> {max_row[gene_col]}")
                    
                    # Top 3 instead of 5 to save tokens
                    top_3 = df.nlargest(3, col)
                    analysis_parts.append(f"\n   TOP 3:")
                    for j, (idx, row) in enumerate(top_3.iterrows(), 1):
                        if mirna_col and gene_col:
                            analysis_parts.append(f"      {j}. {row[mirna_col]} -> {row[gene_col]} | {row[col]:.6f}")
                        else:
                            analysis_parts.append(f"      {j}. Row {idx} | {row[col]:.6f}")
                    
                    # Check token limit
                    current_text = "\n".join(analysis_parts)
                    if count_tokens(current_text) > max_tokens:
                        break
            
            if overall_columns:
                analysis_parts.append(f"\nUSE '{overall_columns[0]}' FOR RANKING QUESTIONS")
            
            return "\n".join(analysis_parts)
            
        except Exception as e:
            return f"Error in score analysis: {str(e)}"
    
    def _get_dataframe_summary_optimized(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Optimized DataFrame summary"""
        summary = {
            'row_count': len(df),
            'column_count': len(df.columns),
            'key_columns': [col for col in df.columns if any(kw in col.lower() 
                          for kw in ['mirna', 'gene', 'target', 'score', 'confidence'])]
        }
        
        # Only include most important statistical info
        numeric_cols = df.select_dtypes(include=['number']).columns
        if len(numeric_cols) > 0:
            # Only select first 3 numeric columns for basic stats
            important_cols = [col for col in numeric_cols if any(kw in col.lower() 
                            for kw in ['score', 'confidence', 'overall'])][:3]
            if important_cols:
                summary['key_stats'] = df[important_cols].describe().loc[['mean', 'max', 'min']].to_dict()
        
        return summary
    
    def _process_excel_optimized(self, df_dict: Dict[str, pd.DataFrame], 
                               max_tokens: int, needs_complete_data: bool) -> Dict[str, Any]:
        """Optimized Excel processing"""
        content = {
            'type': 'excel',
            'sheets': list(df_dict.keys()),
            'data': {}
        }
        
        tokens_per_sheet = max_tokens // len(df_dict) if len(df_dict) > 0 else max_tokens
        
        for sheet_name, sheet_df in df_dict.items():
            if needs_complete_data:
                content['data'][sheet_name] = {
                    'complete_score_analysis': self._extract_complete_scores_optimized(sheet_df, tokens_per_sheet)
                }
            else:
                content['data'][sheet_name] = {
                    'summary': self._get_dataframe_summary_optimized(sheet_df),
                    'sample_data': sheet_df.head(3).to_dict('records')
                }
        
        return content
    
    def _process_text_optimized(self, text_content: str, max_tokens: int) -> Dict[str, Any]:
        """Optimized text processing"""
        # Truncate text based on token limit
        max_chars = max_tokens * 4  # Rough estimation
        
        content = {
            'type': 'text',
            'total_length': len(text_content),
            'content': text_content[:max_chars] if len(text_content) > max_chars else text_content
        }
        
        if len(text_content) > max_chars:
            content['truncated'] = True
            content['truncated_at'] = f"{max_chars} characters (~{max_tokens} tokens)"
        
        return content
    
    def _truncate_content(self, content: Dict[str, Any], 
                         max_tokens: int, needs_complete_data: bool) -> Dict[str, Any]:
        """Truncate cached content to fit token limit"""
        # Simplified version: return original content
        return content
    
    def _estimate_content_tokens(self, content: Dict[str, Any]) -> int:
        """Estimate token count for content"""
        content_str = json.dumps(content, ensure_ascii=False)
        return count_tokens(content_str)
    
    def generate_answer_optimized(self, question: str, content: Dict[str, Any], 
                                question_context: QuestionContext) -> str:
        """Optimized answer generation - token control"""
        if not self.llm:
            return "LLM is not available. Please check Ollama connection."
        
        # Build optimized prompt
        prompt = self._build_optimized_prompt(question, content, question_context)
        
        # Check prompt length
        prompt_tokens = count_tokens(prompt)
        if prompt_tokens > self.max_context_tokens - self.answer_tokens_reserve:
            prompt = self._truncate_prompt(prompt, self.max_context_tokens - self.answer_tokens_reserve)
        
        try:
            response = self.llm.invoke(prompt)
            answer = response.content if hasattr(response, 'content') else str(response)
            
            if not answer.strip():
                return f"Received your question '{question}' but generated an empty response. Please try rephrasing your question."
            
            return answer
            
        except Exception as e:
            return f"Error generating answer: {str(e)}"
    
    def _build_optimized_prompt(self, question: str, content: Dict[str, Any], 
                              context: QuestionContext) -> str:
        """Build optimized prompt"""
        prompt_parts = [
            "You are an intelligent miRNA research assistant. Answer user questions accurately based on the provided pipeline analysis results.",
            "",
            f"Question: {question}",
            f"Type: {context.question_type.value}",
            f"Complexity: {context.complexity_level}/3",
            f"Needs Complete Data: {context.needs_complete_data}",
            "",
            "=== AVAILABLE DATA ==="
        ]
        
        # Adjust content display based on complete data needs
        for source_name, source_content in content.items():
            prompt_parts.append(f"\n--- {source_name} ---")
            
            if source_content.get('type') == 'dataframe':
                prompt_parts.append(f"Data shape: {source_content['shape']}")
                
                if context.needs_complete_data and 'complete_score_analysis' in source_content:
                    # Complete data analysis priority
                    prompt_parts.append(source_content['complete_score_analysis'])
                else:
                    # Summary info
                    if 'summary_stats' in source_content:
                        summary = source_content['summary_stats']
                        prompt_parts.append(f"Key statistics: {summary}")
                    
                    if source_content.get('sample_data'):
                        prompt_parts.append("Sample data:")
                        for i, row in enumerate(source_content['sample_data'][:2]):
                            prompt_parts.append(f"  {i+1}. {row}")
        
        # Key instructions
        prompt_parts.extend([
            "",
            "=== KEY REQUIREMENTS ===",
            f"Answer this question: '{question}'"
        ])
        
        if context.needs_complete_data:
            prompt_parts.extend([
                "CRITICAL for ranking questions:",
                "1. Use the '[MAIN SCORE]' Overall_Score column for rankings",
                "2. Provide specific miRNA-gene pairs and score values",
                "3. Use complete data analysis results (not sample data)"
            ])
        
        prompt_parts.extend([
            "",
            "Requirements: Be accurate, specific, and based on provided data",
            "",
            "=== YOUR ANSWER ==="
        ])
        
        return "\n".join(prompt_parts)
    
    def _truncate_prompt(self, prompt: str, max_tokens: int) -> str:
        """Truncate overly long prompts"""
        current_tokens = count_tokens(prompt)
        if current_tokens <= max_tokens:
            return prompt
        
        # Simple truncation strategy: keep beginning and end, truncate middle
        lines = prompt.split('\n')
        keep_start = len(lines) // 4
        keep_end = len(lines) // 4
        
        truncated_lines = (lines[:keep_start] + 
                          ["\n[... content truncated to fit token limit ...]\n"] + 
                          lines[-keep_end:])
        
        return "\n".join(truncated_lines)

class OptimizedIntelligentQAAgent:
    """Optimized Intelligent Q&A Agent main controller"""
    
    def __init__(self, results_dir: str, logger: Optional[Logger] = None):
        self.results_dir = results_dir
        self.logger = logger
        
        # Initialize optimized agents
        self.routing_agent = OptimizedQuestionRoutingAgent(results_dir, logger)
        self.retrieval_agent = OptimizedContentRetrievalAgent(logger)
        
        self.conversation_history = []
    
    def ask_question(self, question: str) -> Dict[str, Any]:
        """Process user question - optimized version"""
        if self.logger:
            self.logger.info(f"Processing optimized question: {question}")
        
        # First layer: LLM-driven question routing
        context = self.routing_agent.route_question(question)
        
        # Second layer: Smart content loading
        content = self.retrieval_agent.load_content_optimized(
            context.content_sources, 
            context.needs_complete_data
        )
        
        if not content:
            return {
                'question': question,
                'answer': "Sorry, unable to load relevant data files. Please check if the pipeline completed successfully and result files exist.",
                'confidence': 0.0,
                'question_type': context.question_type.value,
                'token_usage': 'minimal'
            }
        
        # Third layer: Optimized answer generation
        answer = self.retrieval_agent.generate_answer_optimized(question, content, context)
        
        result = {
            'question': question,
            'answer': answer,
            'confidence': context.confidence,
            'question_type': context.question_type.value,
            'complexity_level': context.complexity_level,
            'needs_complete_data': context.needs_complete_data,
            'sources_used': len(context.content_sources),
            'llm_classification': True,
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        self.conversation_history.append(result)
        return result
    
    def interactive_session(self):
        """Start interactive Q&A session"""
        print("\n" + "="*60)
        print("Intelligent miRNA Research Assistant")
        print("="*60)
        print("Ask me anything about your pipeline results!")
        print("Type 'quit', 'exit', or 'bye' to end the session.")
        print("Type 'help' to see example questions.")
        print("-"*60)
        
        while True:
            try:
                question = input("\nYour Question: ").strip()
                
                if question.lower() in ['quit', 'exit', 'bye', 'q']:
                    print("\nThank you for using the miRNA Research Assistant!")
                    break
                
                if question.lower() == 'help':
                    self._show_example_questions()
                    continue
                
                if not question:
                    print("Please enter a question.")
                    continue
                
                print("\nProcessing your question...")
                result = self.ask_question(question)
                
                print(f"\nAnswer (Type: {result['question_type']}, Confidence: {result['confidence']:.1%}):")
                if result.get('needs_complete_data'):
                    print("Used complete data analysis")
                print("-" * 40)
                print(result['answer'])
                
                if result['sources_used']:
                    print(f"\nSources used: {result['sources_used']} files")
                
            except KeyboardInterrupt:
                print("\n\nSession interrupted. Goodbye!")
                break
            except Exception as e:
                print(f"\nAn error occurred: {e}")
                if self.logger:
                    self.logger.error(f"Interactive session error: {e}")
    
    def _show_example_questions(self):
        """Show example questions"""
        examples = [
            "How many MTIs were discovered in total?",
            "Which MTI pair has the highest overall score?",
            "What are the top 5 MTIs by Overall_Score?",
            "Which miRNAs have the strongest experimental validation?",
            "What is the literature support for the top MTIs?",
            "What did the AI analysis reveal about ferroptosis pathways?",
            "Show me the network analysis results",
            "What errors or warnings occurred during the pipeline?",
            "Compare the results from different databases",
            "Summarize the main findings of this analysis",
            "Which miRNA has the highest confidence score?",
            "What are the strongest MTI interactions found?"
        ]
        
        print("\nExample Questions:")
        print("-" * 30)
        for i, example in enumerate(examples, 1):
            print(f"{i:2d}. {example}")
        
        print("\nTips for best results:")
        print("- For ranking questions, use terms like 'highest', 'top', 'strongest'")
        print("- Ask about 'Overall_Score' for comprehensive MTI rankings")
        print("- Be specific about what type of score you want to know about")
    
    def save_conversation_history(self, output_file: Optional[str] = None):
        """Save conversation history"""
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
    
    def get_optimization_status(self) -> Dict[str, Any]:
        """Get optimization status"""
        return {
            'llm_classification_available': self.routing_agent.llm is not None,
            'answer_generation_available': self.retrieval_agent.llm is not None,
            'max_context_tokens': self.retrieval_agent.max_context_tokens,
            'content_cache_size': len(self.retrieval_agent.content_cache),
            'conversation_count': len(self.conversation_history)
        }
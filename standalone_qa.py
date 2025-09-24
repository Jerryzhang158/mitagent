#!/usr/bin/env python3
"""
Optimized Standalone Q&A Assistant - miRNA Research Pipeline
Direct intelligent Q&A using completed pipeline results
Supports LLM-driven question classification and smart token management
"""

import os
import argparse
import sys
import time
import json
from typing import Optional, Dict, Any

def find_latest_pipeline_results(base_dir: str = "pipeline_results") -> Optional[str]:
    """Auto-find the latest pipeline results directory"""
    if not os.path.exists(base_dir):
        return None
    
    # Look for all run_* directories
    run_dirs = []
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and item.startswith('run_'):
            # Check if it contains pipeline result files
            if has_pipeline_files(item_path):
                run_dirs.append(item_path)
    
    if not run_dirs:
        return None
    
    # Return the latest (sorted by modification time)
    latest_dir = max(run_dirs, key=lambda x: os.path.getmtime(x))
    return latest_dir

def has_pipeline_files(results_dir: str) -> bool:
    """Check if directory contains pipeline result files"""
    required_indicators = [
        "pipeline_log.txt",
        "mirna_selection",
        "pubmed_articles",
        "bert_validation"
    ]
    
    found_count = 0
    for indicator in required_indicators:
        if os.path.exists(os.path.join(results_dir, indicator)):
            found_count += 1
    
    return found_count >= 2  # Need at least 2 indicator files/directories

def validate_pipeline_results(results_dir: str) -> dict:
    """Validate pipeline results completeness"""
    validation = {
        'valid': True,
        'issues': [],
        'available_components': [],
        'missing_components': [],
        'optimization_ready': False
    }
    
    # Check each component's results
    components = {
        'MTI Selection': ['mirna_selection'],
        'Literature Mining': ['pubmed_articles/mining_summary.csv'],
        'BERT Validation': ['bert_validation/mti_validation_results.csv', 
                           'bert_validation/integrated_validation_results.csv'],
        'LLM Analysis': ['llm_summaries'],
        'Network Analysis': ['cytoscape_network'],
        'Pipeline Logs': ['pipeline_log.txt', 'pipeline_report.txt']
    }
    
    for component_name, files in components.items():
        component_exists = False
        for file_path in files:
            full_path = os.path.join(results_dir, file_path)
            if os.path.exists(full_path):
                component_exists = True
                break
        
        if component_exists:
            validation['available_components'].append(component_name)
        else:
            validation['missing_components'].append(component_name)
    
    # Need at least basic components
    essential_components = ['MTI Selection', 'Pipeline Logs']
    missing_essential = [c for c in essential_components if c in validation['missing_components']]
    
    if missing_essential:
        validation['valid'] = False
        validation['issues'].append(f"Missing essential components: {', '.join(missing_essential)}")
    
    # Check if ready for optimized Q&A (needs detailed data files)
    optimization_files = [
        'bert_validation/integrated_validation_results.csv',
        'bert_validation/mti_validation_results.csv'
    ]
    
    optimization_count = sum(1 for f in optimization_files if os.path.exists(os.path.join(results_dir, f)))
    validation['optimization_ready'] = optimization_count >= 1
    
    if not validation['optimization_ready']:
        validation['issues'].append("Missing key files for optimized Q&A (integrated_validation_results.csv)")
    
    return validation

def show_pipeline_summary(results_dir: str) -> dict:
    """Display pipeline results summary"""
    print("\n" + "="*60)
    print("PIPELINE RESULTS SUMMARY")
    print("="*60)
    
    # Validate results completeness
    validation = validate_pipeline_results(results_dir)
    
    print(f"Results Directory: {os.path.basename(results_dir)}")
    print(f"Full Path: {results_dir}")
    print(f"Status: {'Valid' if validation['valid'] else 'Issues Found'}")
    print(f"Optimization Ready: {'Yes' if validation['optimization_ready'] else 'Limited'}")
    
    print(f"\nAvailable Components ({len(validation['available_components'])}):")
    for component in validation['available_components']:
        print(f"  ✓ {component}")
    
    if validation['missing_components']:
        print(f"\nMissing Components ({len(validation['missing_components'])}):")
        for component in validation['missing_components']:
            print(f"  ✗ {component}")
    
    if validation['issues']:
        print(f"\nIssues Found:")
        for issue in validation['issues']:
            print(f"  ! {issue}")
    
    # Try to read some statistics
    stats = get_pipeline_statistics(results_dir)
    if stats:
        print(f"\nKey Statistics:")
        for stat_name, stat_value in stats.items():
            print(f"  • {stat_name}: {stat_value}")
    
    return validation

def get_pipeline_statistics(results_dir: str) -> dict:
    """Get pipeline statistics"""
    stats = {}
    
    try:
        # Extract stats from log file
        log_file = os.path.join(results_dir, "pipeline_log.txt")
        if os.path.exists(log_file):
            with open(log_file, 'r', encoding='utf-8') as f:
                log_content = f.read()
            stats.update(extract_stats_from_log(log_content))
        
        # Extract stats from CSV files
        csv_stats = extract_stats_from_csv(results_dir)
        stats.update(csv_stats)
        
    except Exception as e:
        print(f"  ! Could not read statistics: {e}")
    
    return stats

def extract_stats_from_log(log_content: str) -> dict:
    """Extract key statistics from log"""
    stats = {}
    import re
    
    # Extract MTI count
    mti_match = re.search(r'Found (\d+) MTIs', log_content)
    if mti_match:
        stats['Total MTIs'] = mti_match.group(1)
    
    # Extract articles count
    articles_match = re.search(r'Total articles found: (\d+)', log_content)
    if articles_match:
        stats['Articles Found'] = articles_match.group(1)
    
    # Extract validation count
    validation_match = re.search(r'MTIs validated: (\d+)', log_content)
    if validation_match:
        stats['MTIs Validated'] = validation_match.group(1)
    
    return stats

def extract_stats_from_csv(results_dir: str) -> dict:
    """Extract statistics from CSV files"""
    stats = {}
    
    try:
        import pandas as pd
        
        # Try to read integrated validation results
        integrated_file = os.path.join(results_dir, 'bert_validation', 'integrated_validation_results.csv')
        if os.path.exists(integrated_file):
            df = pd.read_csv(integrated_file)
            stats['Total Analyzed MTIs'] = str(len(df))
            
            # Look for score columns
            score_cols = [col for col in df.columns if 'score' in col.lower() or 'confidence' in col.lower()]
            if score_cols:
                stats['Score Columns Available'] = str(len(score_cols))
                
                # Look for Overall_Score
                overall_cols = [col for col in score_cols if 'overall' in col.lower()]
                if overall_cols:
                    max_score = df[overall_cols[0]].max()
                    stats['Highest Overall Score'] = f"{max_score:.4f}"
        
        # Try to read network files
        nodes_file = os.path.join(results_dir, 'cytoscape_network', 'nodes.csv')
        if os.path.exists(nodes_file):
            nodes_df = pd.read_csv(nodes_file)
            stats['Network Nodes'] = str(len(nodes_df))
        
        edges_file = os.path.join(results_dir, 'cytoscape_network', 'edges.csv')
        if os.path.exists(edges_file):
            edges_df = pd.read_csv(edges_file)
            stats['Network Edges'] = str(len(edges_df))
            
    except ImportError:
        stats['Note'] = "pandas not available for detailed stats"
    except Exception as e:
        stats['CSV Stats Error'] = str(e)
    
    return stats

def check_dependencies() -> dict:
    """Check dependency availability"""
    deps = {
        'pandas': False,
        'ollama_available': False,
        'langchain': False,
        'tiktoken': False
    }
    
    try:
        import pandas
        deps['pandas'] = True
    except ImportError:
        pass
    
    try:
        from utils import check_ollama_connection
        is_connected, models = check_ollama_connection()
        deps['ollama_available'] = is_connected
        deps['ollama_models'] = models if is_connected else []
    except:
        deps['ollama_available'] = False
        deps['ollama_models'] = []
    
    try:
        import langchain_community
        deps['langchain'] = True
    except ImportError:
        pass
    
    try:
        import tiktoken
        deps['tiktoken'] = True
    except ImportError:
        pass
    
    return deps

def show_optimization_status(deps: dict):
    """Show optimization features status"""
    print("\nOPTIMIZATION STATUS")
    print("="*30)
    
    print(f"Data Processing: {'✓' if deps['pandas'] else '✗'} pandas")
    print(f"LLM Classification: {'✓' if deps['ollama_available'] and deps['langchain'] else '✗'} Ollama + LangChain")
    print(f"Token Management: {'✓' if deps['tiktoken'] else '✗'} tiktoken")
    
    if deps['ollama_available'] and deps['ollama_models']:
        print(f"Available Models: {', '.join(deps['ollama_models'][:3])}")
        if len(deps['ollama_models']) > 3:
            print(f"    ... and {len(deps['ollama_models']) - 3} more")
    
    # Calculate optimization score
    optimization_score = sum([
        deps['pandas'],
        deps['ollama_available'] and deps['langchain'],
        deps['tiktoken']
    ])
    
    if optimization_score == 3:
        level = "Full Optimization"
        desc = "All advanced features available"
    elif optimization_score == 2:
        level = "Partial Optimization"
        desc = "Most features available, some limitations"
    elif optimization_score == 1:
        level = "Basic Mode"
        desc = "Limited functionality, basic Q&A only"
    else:
        level = "Minimal Mode"
        desc = "Severe limitations, may not work properly"
    
    print(f"\nOptimization Level: {level}")
    print(f"   {desc}")
    
    return optimization_score

def install_recommendations(deps: dict):
    """Show installation recommendations"""
    missing = []
    
    if not deps['pandas']:
        missing.append("pip install pandas")
    if not deps['langchain']:
        missing.append("pip install langchain-community")
    if not deps['tiktoken']:
        missing.append("pip install tiktoken")
    if not deps['ollama_available']:
        missing.append("# Install and start Ollama:\n# Download from https://ollama.ai\n# ollama serve\n# ollama pull llama3.1")
    
    if missing:
        print("\nTO ENABLE FULL OPTIMIZATION:")
        print("-" * 35)
        for cmd in missing:
            if cmd.startswith('#'):
                print(cmd)
            else:
                print(f"  {cmd}")

def start_qa_session(results_dir: str, single_question: Optional[str] = None, force_basic: bool = False):
    """Start Q&A session - supports optimized and basic versions"""
    
    # Check dependencies
    deps = check_dependencies()
    optimization_score = sum([deps['pandas'], deps['ollama_available'] and deps['langchain'], deps['tiktoken']])
    
    try:
        # Try to use optimized version
        if optimization_score >= 2 and not force_basic:
            print("\nLoading Optimized Q&A Assistant...")
            success = start_optimized_qa_session(results_dir, single_question, deps)
            if success:
                return True
            else:
                print("! Optimized version failed, falling back to basic version...")
        
        # Fallback to basic version
        print("\nLoading Basic Q&A Assistant...")
        return start_basic_qa_session(results_dir, single_question)
        
    except Exception as e:
        print(f"✗ Error starting Q&A session: {e}")
        import traceback
        print("Full error traceback:")
        print(traceback.format_exc())
        return False

def start_optimized_qa_session(results_dir: str, single_question: Optional[str], deps: dict) -> bool:
    """Start optimized Q&A session"""
    try:
        # Dynamic import of optimized agent
        if 'intelligent_qa_agent' in sys.modules:
            # Reload module to get latest version
            import importlib
            importlib.reload(sys.modules['intelligent_qa_agent'])
        
        # Try to import optimized version
        try:
            from intelligent_qa_agent import OptimizedIntelligentQAAgent
            agent_class = OptimizedIntelligentQAAgent
            print("✓ Using OptimizedIntelligentQAAgent")
        except ImportError:
            # If no optimized version, try original
            from intelligent_qa_agent import IntelligentQAAgent
            agent_class = IntelligentQAAgent
            print("! Using IntelligentQAAgent (optimization features limited)")
        
        from utils import Logger
        
        # Create logger
        log_file = os.path.join(results_dir, "qa_debug.log")
        logger = Logger(log_file, console_output=False)  # Disable console output for cleaner experience
        
        # Create agent
        qa_agent = agent_class(results_dir, logger)
        
        # Check status
        if hasattr(qa_agent, 'get_optimization_status'):
            status = qa_agent.get_optimization_status()
            print(f"Optimization Status:")
            print(f"  LLM Classification: {'✓' if status.get('llm_classification_available', False) else '✗'}")
            print(f"  Answer Generation: {'✓' if status.get('answer_generation_available', False) else '✗'}")
            print(f"  Max Context Tokens: {status.get('max_context_tokens', 'N/A')}")
        else:
            status = qa_agent.get_agent_status()
            print(f"Agent Status:")
            print(f"  LLM Available: {'✓' if status['llm_available'] else '✗'}")
        
        # Execute Q&A
        if single_question:
            return execute_single_question(qa_agent, single_question, log_file)
        else:
            return execute_interactive_session(qa_agent, log_file)
            
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False
    except Exception as e:
        print(f"✗ Optimization session error: {e}")
        return False

def start_basic_qa_session(results_dir: str, single_question: Optional[str]) -> bool:
    """Start basic Q&A session"""
    try:
        # Basic version implementation - direct file reading and simple answers
        print("Basic Q&A mode - limited functionality")
        
        if single_question:
            answer = generate_basic_answer(results_dir, single_question)
            print(f"\nQuestion: {single_question}")
            print("-" * 60)
            print(f"Answer: {answer}")
            return True
        else:
            return execute_basic_interactive(results_dir)
            
    except Exception as e:
        print(f"✗ Basic session error: {e}")
        return False

def generate_basic_answer(results_dir: str, question: str) -> str:
    """Generate basic answer (no LLM dependency)"""
    question_lower = question.lower()
    
    # Simple keyword matching
    try:
        if any(kw in question_lower for kw in ['how many', 'count', 'number']):
            return analyze_counts(results_dir)
        elif any(kw in question_lower for kw in ['highest', 'best', 'top', 'strongest']):
            return analyze_top_results(results_dir)
        elif any(kw in question_lower for kw in ['summary', 'summarize']):
            return generate_summary(results_dir)
        else:
            return get_general_info(results_dir)
    except Exception as e:
        return f"Unable to analyze data files: {str(e)}"

def analyze_counts(results_dir: str) -> str:
    """Analyze count information"""
    try:
        import pandas as pd
        
        info = []
        
        # Check integrated validation results
        integrated_file = os.path.join(results_dir, 'bert_validation', 'integrated_validation_results.csv')
        if os.path.exists(integrated_file):
            df = pd.read_csv(integrated_file)
            info.append(f"Total analyzed MTI interactions: {len(df)}")
            
        # Check network files
        nodes_file = os.path.join(results_dir, 'cytoscape_network', 'nodes.csv')
        if os.path.exists(nodes_file):
            nodes_df = pd.read_csv(nodes_file)
            info.append(f"Network contains {len(nodes_df)} nodes")
        
        edges_file = os.path.join(results_dir, 'cytoscape_network', 'edges.csv')
        if os.path.exists(edges_file):
            edges_df = pd.read_csv(edges_file)
            info.append(f"Network contains {len(edges_df)} edges")
            
        return "\n".join(info) if info else "No count information found"
        
    except ImportError:
        return "Need pandas to analyze count information"
    except Exception as e:
        return f"Count analysis error: {str(e)}"

def analyze_top_results(results_dir: str) -> str:
    """Analyze highest score results"""
    try:
        import pandas as pd
        
        integrated_file = os.path.join(results_dir, 'bert_validation', 'integrated_validation_results.csv')
        if not os.path.exists(integrated_file):
            return "Integrated validation results file not found"
        
        df = pd.read_csv(integrated_file)
        
        # Look for score columns
        score_cols = [col for col in df.columns if 'score' in col.lower()]
        overall_cols = [col for col in score_cols if 'overall' in col.lower()]
        
        if overall_cols:
            main_score_col = overall_cols[0]
        elif score_cols:
            main_score_col = score_cols[0]
        else:
            return "No score columns found"
        
        # Get highest score
        max_score = df[main_score_col].max()
        max_idx = df[main_score_col].idxmax()
        max_row = df.loc[max_idx]
        
        # Find miRNA and gene columns
        mirna_col = None
        gene_col = None
        for col in df.columns:
            if 'mirna' in col.lower() and not mirna_col:
                mirna_col = col
            elif any(kw in col.lower() for kw in ['gene', 'target']) and not gene_col:
                gene_col = col
        
        result = f"Highest score: {max_score:.6f} (from {main_score_col})\n"
        
        if mirna_col and gene_col:
            result += f"Best MTI pair: {max_row[mirna_col]} -> {max_row[gene_col]}\n"
        
        # Top 3
        top_3 = df.nlargest(3, main_score_col)
        result += "\nTop 3 MTIs:\n"
        for i, (idx, row) in enumerate(top_3.iterrows(), 1):
            if mirna_col and gene_col:
                result += f"{i}. {row[mirna_col]} -> {row[gene_col]}: {row[main_score_col]:.6f}\n"
            else:
                result += f"{i}. Row {idx}: {row[main_score_col]:.6f}\n"
        
        return result
        
    except ImportError:
        return "Need pandas to analyze highest score results"
    except Exception as e:
        return f"Highest score analysis error: {str(e)}"

def generate_summary(results_dir: str) -> str:
    """Generate summary"""
    stats = get_pipeline_statistics(results_dir)
    validation = validate_pipeline_results(results_dir)
    
    summary = "=== miRNA Pipeline Results Summary ===\n\n"
    summary += f"Available components: {', '.join(validation['available_components'])}\n"
    
    if stats:
        summary += "\nMain statistics:\n"
        for key, value in stats.items():
            summary += f"- {key}: {value}\n"
    
    if validation['issues']:
        summary += f"\nNotes: {len(validation['issues'])} issues found\n"
    
    return summary

def get_general_info(results_dir: str) -> str:
    """Get general information"""
    return f"""
This is a miRNA pipeline analysis result.

Available question types:
- Count questions: "How many MTIs were found?"
- Ranking questions: "Which MTI has the highest score?"  
- Summary questions: "Summarize the results"

Results directory: {results_dir}

For more detailed analysis, install full dependencies:
- pip install pandas langchain-community tiktoken
- Install and start Ollama
"""

def execute_single_question(qa_agent, question: str, log_file: str) -> bool:
    """Execute single question"""
    print(f"\nProcessing question: {question}")
    
    result = qa_agent.ask_question(question)
    
    print("="*60)
    print("RESULT:")
    print("="*60)
    print(f"Question Type: {result.get('question_type', 'unknown')}")
    print(f"Confidence: {result.get('confidence', 0):.1%}")
    
    if result.get('needs_complete_data'):
        print("Used complete data analysis")
    
    print("\nAnswer:")
    print("-" * 40)
    print(result['answer'])
    
    if result.get('sources_used'):
        count = result['sources_used'] if isinstance(result['sources_used'], int) else len(result['sources_used'])
        print(f"\nData sources used: {count}")
    
    # Save history
    if hasattr(qa_agent, 'save_conversation_history'):
        qa_agent.save_conversation_history()
    
    print(f"\nDebug log: {log_file}")
    return True

def execute_interactive_session(qa_agent, log_file: str) -> bool:
    """Execute interactive session"""
    if hasattr(qa_agent, 'interactive_session'):
        qa_agent.interactive_session()
    else:
        # Simple interaction implementation
        print("\n" + "="*60)
        print("Interactive Q&A Session")
        print("="*60)
        print("Type 'quit' to exit")
        
        while True:
            try:
                question = input("\nYour Question: ").strip()
                if question.lower() in ['quit', 'exit', 'q']:
                    break
                
                if not question:
                    continue
                
                result = qa_agent.ask_question(question)
                print(f"\nAnswer ({result.get('confidence', 0):.1%} confidence):")
                print("-" * 40)
                print(result['answer'])
                
            except KeyboardInterrupt:
                print("\n\nSession interrupted!")
                break
    
    # Save history
    if hasattr(qa_agent, 'save_conversation_history'):
        qa_agent.save_conversation_history()
        print(f"\nConversation saved")
    
    print(f"Debug log: {log_file}")
    return True

def execute_basic_interactive(results_dir: str) -> bool:
    """Execute basic interaction"""
    print("\n" + "="*60)
    print("Basic Interactive Q&A")
    print("="*60)
    print("Limited functionality - install dependencies for full features")
    print("Type 'quit' to exit")
    
    while True:
        try:
            question = input("\nQuestion: ").strip()
            if question.lower() in ['quit', 'exit', 'q']:
                break
            
            if not question:
                continue
            
            answer = generate_basic_answer(results_dir, question)
            print(f"\nAnswer:")
            print("-" * 40)
            print(answer)
            
        except KeyboardInterrupt:
            print("\n\nSession ended!")
            break
    
    return True

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Optimized Intelligent Q&A Assistant for miRNA Pipeline Results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-find latest results and start interactive Q&A
  python standalone_qa.py
  
  # Use specific results directory
  python standalone_qa.py -d pipeline_results/run_20240101_120000
  
  # Ask a single question
  python standalone_qa.py --ask "What is the highest scoring MTI?"
  
  # Check optimization status
  python standalone_qa.py --check-deps
  
  # Force basic mode (no LLM)
  python standalone_qa.py --basic-mode
        """
    )
    
    parser.add_argument('-d', '--results-dir', 
                       help='Path to pipeline results directory')
    
    parser.add_argument('--ask', type=str,
                       help='Ask a single question (non-interactive mode)')
    
    parser.add_argument('--list-results', action='store_true',
                       help='List all available pipeline result directories')
    
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate the results without starting Q&A')
    
    parser.add_argument('--check-deps', action='store_true',
                       help='Check optimization dependencies')
    
    parser.add_argument('--basic-mode', action='store_true',
                       help='Force basic mode (no LLM optimization)')
    
    parser.add_argument('--base-dir', default='pipeline_results',
                       help='Base directory to search for pipeline results')
    
    args = parser.parse_args()
    
    print("="*60)
    print("miRNA Pipeline - Optimized Q&A Assistant")
    print("="*60)
    
    # Check dependencies
    if args.check_deps:
        deps = check_dependencies()
        optimization_score = show_optimization_status(deps)
        install_recommendations(deps)
        return 0 if optimization_score >= 2 else 1
    
    # List available result directories
    if args.list_results:
        print(f"Searching for pipeline results in: {args.base_dir}")
        
        if not os.path.exists(args.base_dir):
            print(f"✗ Base directory not found: {args.base_dir}")
            return 1
        
        found_dirs = []
        for item in os.listdir(args.base_dir):
            item_path = os.path.join(args.base_dir, item)
            if os.path.isdir(item_path) and has_pipeline_files(item_path):
                found_dirs.append(item_path)
        
        if found_dirs:
            print(f"\nFound {len(found_dirs)} pipeline result directories:")
            for i, result_dir in enumerate(sorted(found_dirs), 1):
                print(f"  {i}. {os.path.basename(result_dir)}")
                print(f"     Path: {result_dir}")
                mtime = os.path.getmtime(result_dir)
                print(f"     Modified: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mtime))}")
        else:
            print("✗ No pipeline result directories found")
        
        return 0
    
    # Determine results directory
    if args.results_dir:
        results_dir = args.results_dir
        if not os.path.exists(results_dir):
            print(f"✗ Results directory not found: {results_dir}")
            return 1
        print(f"Using specified directory: {results_dir}")
    else:
        results_dir = find_latest_pipeline_results(args.base_dir)
        if not results_dir:
            print(f"✗ No pipeline results found in {args.base_dir}")
            print("Run pipeline first or use --list-results")
            return 1
        print(f"Auto-detected: {os.path.basename(results_dir)}")
    
    # Show results summary
    validation = show_pipeline_summary(results_dir)
    
    # Show optimization status
    deps = check_dependencies()
    optimization_score = show_optimization_status(deps)
    
    if optimization_score < 2:
        install_recommendations(deps)
    
    # Validation mode
    if args.validate_only:
        success = validation['valid'] and validation['optimization_ready']
        print(f"\nValidation complete: {'Fully ready' if success else 'Limited functionality'}")
        return 0 if success else 1
    
    # Start Q&A
    print(f"\nStarting Q&A Assistant...")
    success = start_qa_session(results_dir, args.ask, args.basic_mode)
    
    if success:
        print(f"\nQ&A session completed!")
        return 0
    else:
        print(f"\nQ&A session failed")
        return 1

if __name__ == "__main__":
    exit(main())
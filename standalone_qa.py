#!/usr/bin/env python3
"""
独立问答助手 - miRNA Research Pipeline
直接使用已完成的pipeline结果进行智能问答
"""

import os
import argparse
import sys
from typing import Optional

def find_latest_pipeline_results(base_dir: str = "pipeline_results") -> Optional[str]:
    """自动找到最新的pipeline结果目录"""
    if not os.path.exists(base_dir):
        return None
    
    # 查找所有run_*目录
    run_dirs = []
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and item.startswith('run_'):
            # 检查是否包含pipeline结果文件
            if has_pipeline_files(item_path):
                run_dirs.append(item_path)
    
    if not run_dirs:
        return None
    
    # 返回最新的（按修改时间排序）
    latest_dir = max(run_dirs, key=lambda x: os.path.getmtime(x))
    return latest_dir

def has_pipeline_files(results_dir: str) -> bool:
    """检查目录是否包含pipeline结果文件"""
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
    
    return found_count >= 2  # 至少要有2个指标文件/目录

def validate_pipeline_results(results_dir: str) -> dict:
    """验证pipeline结果的完整性"""
    validation = {
        'valid': True,
        'issues': [],
        'available_components': [],
        'missing_components': []
    }
    
    # 检查各个组件的结果
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
    
    # 至少需要有基本组件
    essential_components = ['MTI Selection', 'Pipeline Logs']
    missing_essential = [c for c in essential_components if c in validation['missing_components']]
    
    if missing_essential:
        validation['valid'] = False
        validation['issues'].append(f"Missing essential components: {', '.join(missing_essential)}")
    
    return validation

def show_pipeline_summary(results_dir: str):
    """显示pipeline结果摘要"""
    print("\n" + "="*60)
    print("📊 PIPELINE RESULTS SUMMARY")
    print("="*60)
    
    # 验证结果完整性
    validation = validate_pipeline_results(results_dir)
    
    print(f"📁 Results Directory: {os.path.basename(results_dir)}")
    print(f"📍 Full Path: {results_dir}")
    print(f"✅ Status: {'Valid' if validation['valid'] else 'Issues Found'}")
    
    print(f"\n🔬 Available Components ({len(validation['available_components'])}):")
    for component in validation['available_components']:
        print(f"  ✅ {component}")
    
    if validation['missing_components']:
        print(f"\n⚠️ Missing Components ({len(validation['missing_components'])}):")
        for component in validation['missing_components']:
            print(f"  ❌ {component}")
    
    if validation['issues']:
        print(f"\n🚨 Issues Found:")
        for issue in validation['issues']:
            print(f"  ⚠️ {issue}")
    
    # 尝试读取一些统计信息
    try:
        log_file = os.path.join(results_dir, "pipeline_log.txt")
        if os.path.exists(log_file):
            with open(log_file, 'r', encoding='utf-8') as f:
                log_content = f.read()
                
            # 提取一些关键统计信息
            stats = extract_stats_from_log(log_content)
            if stats:
                print(f"\n📈 Key Statistics:")
                for stat_name, stat_value in stats.items():
                    print(f"  📊 {stat_name}: {stat_value}")
    
    except Exception as e:
        print(f"  ⚠️ Could not read pipeline statistics: {e}")
    
    return validation['valid']

def extract_stats_from_log(log_content: str) -> dict:
    """从日志中提取关键统计信息"""
    stats = {}
    
    import re
    
    # 提取MTI数量
    mti_match = re.search(r'Found (\d+) MTIs', log_content)
    if mti_match:
        stats['Total MTIs'] = mti_match.group(1)
    
    # 提取文献数量
    articles_match = re.search(r'Total articles found: (\d+)', log_content)
    if articles_match:
        stats['Articles Found'] = articles_match.group(1)
    
    # 提取验证数量
    validation_match = re.search(r'MTIs validated: (\d+)', log_content)
    if validation_match:
        stats['MTIs Validated'] = validation_match.group(1)
    
    # 提取网络节点数
    nodes_match = re.search(r'Generated (\d+) nodes', log_content)
    if nodes_match:
        stats['Network Nodes'] = nodes_match.group(1)
    
    return stats

def start_qa_session(results_dir: str, single_question: Optional[str] = None):
    """启动问答会话"""
    try:
        # 导入QA Agent
        from intelligent_qa_agent import IntelligentQAAgent
        from utils import Logger
        
        # 创建临时日志用于调试
        log_file = os.path.join(results_dir, "qa_debug.log")
        logger = Logger(log_file, console_output=True)  # 启用控制台输出
        
        # 创建QA Agent
        print("\n🤖 Initializing Intelligent Q&A Assistant...")
        qa_agent = IntelligentQAAgent(results_dir, logger)
        
        # 检查Agent状态
        status = qa_agent.get_agent_status()
        print(f"📊 QA Agent Status:")
        print(f"  🧠 LLM Available: {status['llm_available']}")
        print(f"  💾 Content Cache: {status['content_cache_size']} entries")
        print(f"  🎯 Supported Questions: {len(status['supported_question_types'])} types")
        
        if not status['llm_available']:
            print("\n⚠️ Warning: LLM is not available. Responses may be limited.")
            print("💡 To enable full functionality:")
            print("  1. Start Ollama: ollama serve")
            print("  2. Install model: ollama pull llama3.1")
            print("  3. Or install langchain-community: pip install langchain-community")
        
        if single_question:
            # 单次问题模式
            print(f"\n🔬 Processing your question: {single_question}")
            
            # 添加调试信息
            print("="*60)
            print("🔍 DEBUG INFORMATION:")
            print("="*60)
            
            result = qa_agent.ask_question(single_question)
            
            print("="*60)
            print("📊 FINAL RESULT:")
            print("="*60)
            
            print(f"📊 Answer (Type: {result['question_type']}, Confidence: {result['confidence']:.1%}):")
            print("-" * 60)
            print(result['answer'])
            
            if result['sources_used']:
                print(f"\n📁 Data sources used: {len(result['sources_used'])}")
                for i, source in enumerate(result['sources_used'][:3], 1):
                    print(f"  {i}. {os.path.basename(source)}")
                if len(result['sources_used']) > 3:
                    print(f"  ... and {len(result['sources_used']) - 3} more")
            
            # 保存单次问答历史
            qa_agent.save_conversation_history()
            
            # 提示查看调试日志
            print(f"\n🔍 Debug log saved to: {log_file}")
        
        else:
            # 交互式问答模式
            qa_agent.interactive_session()
            
            # 保存对话历史
            qa_agent.save_conversation_history()
            
            # 显示对话统计
            if qa_agent.conversation_history:
                print(f"\n📋 Session Summary:")
                print(f"  💬 Questions Asked: {len(qa_agent.conversation_history)}")
                print(f"  💾 Conversation Saved: qa_conversation_history.json")
                print(f"  🔍 Debug Log: {log_file}")
    
    except ImportError as e:
        print(f"❌ Error: Required modules not found - {e}")
        print("💡 Make sure all pipeline files are in the same directory")
        return False
    
    except Exception as e:
        print(f"❌ Error starting Q&A session: {e}")
        import traceback
        print("🔍 Full error traceback:")
        print(traceback.format_exc())
        return False
    
    return True

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Intelligent Q&A Assistant for miRNA Pipeline Results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-find latest results and start interactive Q&A
  python standalone_qa.py
  
  # Use specific results directory
  python standalone_qa.py -d pipeline_results/run_20240101_120000
  
  # Ask a single question
  python standalone_qa.py --ask "How many MTIs were discovered?"
  
  # List available result directories
  python standalone_qa.py --list-results
        """
    )
    
    parser.add_argument('-d', '--results-dir', 
                       help='Path to pipeline results directory')
    
    parser.add_argument('--ask', '--ask-question', type=str,
                       help='Ask a single question (non-interactive mode)')
    
    parser.add_argument('--list-results', action='store_true',
                       help='List all available pipeline result directories')
    
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate the results without starting Q&A')
    
    parser.add_argument('--base-dir', default='pipeline_results',
                       help='Base directory to search for pipeline results')
    
    args = parser.parse_args()
    
    print("="*60)
    print("🤖 miRNA Pipeline - Standalone Q&A Assistant")
    print("="*60)
    
    # 列出可用结果目录
    if args.list_results:
        print(f"🔍 Searching for pipeline results in: {args.base_dir}")
        
        if not os.path.exists(args.base_dir):
            print(f"❌ Base directory not found: {args.base_dir}")
            return 1
        
        found_dirs = []
        for item in os.listdir(args.base_dir):
            item_path = os.path.join(args.base_dir, item)
            if os.path.isdir(item_path) and has_pipeline_files(item_path):
                found_dirs.append(item_path)
        
        if found_dirs:
            print(f"\n📁 Found {len(found_dirs)} pipeline result directories:")
            for i, result_dir in enumerate(sorted(found_dirs), 1):
                print(f"  {i}. {os.path.basename(result_dir)}")
                print(f"     Path: {result_dir}")
                # 显示修改时间
                import time
                mtime = os.path.getmtime(result_dir)
                print(f"     Modified: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mtime))}")
        else:
            print("❌ No pipeline result directories found")
            print("💡 Make sure you have run the pipeline at least once")
        
        return 0
    
    # 确定要使用的结果目录
    if args.results_dir:
        results_dir = args.results_dir
        if not os.path.exists(results_dir):
            print(f"❌ Error: Results directory not found: {results_dir}")
            return 1
        print(f"📁 Using specified results directory: {results_dir}")
    else:
        # 自动查找最新结果
        results_dir = find_latest_pipeline_results(args.base_dir)
        if not results_dir:
            print(f"❌ Error: No pipeline results found in {args.base_dir}")
            print("💡 Available options:")
            print("  1. Run the pipeline first to generate results")
            print("  2. Specify results directory with -d option")
            print("  3. Use --list-results to see available directories")
            return 1
        print(f"🔍 Auto-detected latest results: {os.path.basename(results_dir)}")
    
    # 显示结果摘要和验证
    is_valid = show_pipeline_summary(results_dir)
    
    if not is_valid:
        print("\n❌ Pipeline results validation failed!")
        print("💡 The Q&A assistant may not work properly with incomplete results.")
        
        user_input = input("\nDo you want to continue anyway? (y/N): ").strip().lower()
        if user_input not in ['y', 'yes']:
            print("👋 Exiting. Please use complete pipeline results.")
            return 1
    
    # 如果只是验证模式，到此结束
    if args.validate_only:
        print(f"\n✅ Validation complete. Results are {'valid' if is_valid else 'incomplete but usable'}.")
        return 0 if is_valid else 1
    
    # 启动问答会话
    print(f"\n🚀 Starting Q&A Assistant...")
    success = start_qa_session(results_dir, args.ask)
    
    if success:
        print(f"\n✅ Q&A session completed successfully!")
        return 0
    else:
        print(f"\n❌ Q&A session failed. Please check the error messages above.")
        return 1

if __name__ == "__main__":
    exit(main())
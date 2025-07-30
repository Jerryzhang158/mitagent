#!/usr/bin/env python3
"""
QA Agent诊断工具
帮助诊断智能问答Agent的问题
"""

import os
import sys
import glob
from typing import List, Dict

def check_results_directory(results_dir: str) -> Dict[str, any]:
    """检查结果目录的完整性"""
    diagnosis = {
        'directory_exists': os.path.exists(results_dir),
        'files_found': {},
        'directories_found': {},
        'issues': [],
        'recommendations': []
    }
    
    if not diagnosis['directory_exists']:
        diagnosis['issues'].append(f"Results directory does not exist: {results_dir}")
        diagnosis['recommendations'].append("Make sure you specify the correct pipeline results directory")
        return diagnosis
    
    print(f"📁 Checking directory: {results_dir}")
    
    # 检查预期的文件
    expected_files = {
        'pipeline_log.txt': 'Pipeline execution log',
        'pipeline_report.txt': 'Pipeline summary report',
        'mirna_selection/mti_selection_results_*.xlsx': 'MTI selection results',
        'bert_validation/mti_validation_results.csv': 'BERT validation results',
        'bert_validation/integrated_validation_results.csv': 'Integrated validation results',
        'pubmed_articles/mining_summary.csv': 'Literature mining summary',
        'llm_summaries/': 'LLM analysis results (directory)',
        'cytoscape_network/nodes.csv': 'Network nodes',
        'cytoscape_network/edges.csv': 'Network edges'
    }
    
    # 检查预期的目录
    expected_dirs = [
        'mirna_selection',
        'pubmed_articles', 
        'bert_validation',
        'llm_summaries',
        'cytoscape_network'
    ]
    
    # 检查文件
    for file_pattern, description in expected_files.items():
        if '/' in file_pattern and not file_pattern.endswith('/'):
            # 这是一个具体文件
            file_path = os.path.join(results_dir, file_pattern)
            if '*' in file_path:
                matches = glob.glob(file_path)
                found = len(matches) > 0
                if found:
                    diagnosis['files_found'][file_pattern] = matches[0]
                    print(f"  ✅ {description}: {os.path.basename(matches[0])}")
                else:
                    print(f"  ❌ {description}: Not found")
            else:
                found = os.path.exists(file_path)
                if found:
                    diagnosis['files_found'][file_pattern] = file_path
                    print(f"  ✅ {description}: Found")
                else:
                    print(f"  ❌ {description}: Not found")
        elif file_pattern.endswith('/'):
            # 这是一个目录
            dir_path = os.path.join(results_dir, file_pattern.rstrip('/'))
            found = os.path.exists(dir_path) and os.path.isdir(dir_path)
            if found:
                diagnosis['directories_found'][file_pattern] = dir_path
                # 计算目录中的文件数量
                files_count = len([f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))])
                print(f"  ✅ {description}: Found ({files_count} files)")
            else:
                print(f"  ❌ {description}: Not found")
        else:
            # 根目录中的文件
            file_path = os.path.join(results_dir, file_pattern)
            found = os.path.exists(file_path)
            if found:
                diagnosis['files_found'][file_pattern] = file_path
                print(f"  ✅ {description}: Found")
            else:
                print(f"  ❌ {description}: Not found")
    
    # 分析问题
    essential_files = ['pipeline_log.txt', 'mirna_selection/mti_selection_results_*.xlsx']
    missing_essential = []
    
    for essential in essential_files:
        if essential not in diagnosis['files_found']:
            if '*' in essential:
                # 检查通配符匹配
                pattern_path = os.path.join(results_dir, essential)
                if not glob.glob(pattern_path):
                    missing_essential.append(essential)
            else:
                missing_essential.append(essential)
    
    if missing_essential:
        diagnosis['issues'].append(f"Missing essential files: {missing_essential}")
        diagnosis['recommendations'].append("The pipeline may not have completed successfully")
    
    # 计算完整性分数
    total_expected = len(expected_files)
    found_count = len(diagnosis['files_found']) + len(diagnosis['directories_found'])
    diagnosis['completeness'] = found_count / total_expected * 100
    
    return diagnosis

def test_question_routing(results_dir: str):
    """测试问题路由功能"""
    print("\n🧪 Testing question routing...")
    
    try:
        from intelligent_qa_agent import QuestionRoutingAgent
        from utils import Logger
        
        # 创建路由Agent
        routing_agent = QuestionRoutingAgent(results_dir)
        
        # 测试问题
        test_questions = [
            "Tell me about which mti has the strongest possibility",
            "ferroptosis相关的MTI有多少个",
            "How many MTIs were discovered?",
            "What is the experimental validation rate?",
            "Show me the literature support",
            "Summarize the main findings"
        ]
        
        print("Testing question routing for sample questions:")
        for question in test_questions:
            try:
                context = routing_agent.route_question(question)
                print(f"  📝 '{question}'")
                print(f"    → Type: {context.question_type.value}")
                print(f"    → Confidence: {context.confidence:.1%}")
                print(f"    → Sources found: {len(context.content_sources)}")
                if context.content_sources:
                    for source in context.content_sources[:2]:  # 显示前2个源
                        status = "✅" if os.path.exists(source.file_path) else "❌"
                        print(f"      {status} {source.description}")
                print()
            except Exception as e:
                print(f"    ❌ Error: {e}")
    
    except ImportError as e:
        print(f"❌ Cannot import QA Agent modules: {e}")
        print("Make sure all pipeline files are in the current directory")

def test_content_loading(results_dir: str):
    """测试内容加载功能"""
    print("\n🧪 Testing content loading...")
    
    try:
        from intelligent_qa_agent import IntelligentQAAgent
        
        qa_agent = IntelligentQAAgent(results_dir)
        
        # 测试加载一些基本内容
        basic_question = "How many MTIs were found?"
        context = qa_agent.routing_agent.route_question(basic_question)
        
        if context.content_sources:
            print(f"Attempting to load {len(context.content_sources)} content sources:")
            content = qa_agent.retrieval_agent.load_content(context.content_sources)
            
            if content:
                print(f"✅ Successfully loaded {len(content)} content sources:")
                for name, source_content in content.items():
                    print(f"  📊 {name}: {source_content['type']}")
            else:
                print("❌ Failed to load any content")
        else:
            print("❌ No content sources found for the test question")
    
    except Exception as e:
        print(f"❌ Content loading test failed: {e}")

def check_llm_availability():
    """检查LLM可用性"""
    print("\n🧪 Checking LLM availability...")
    
    try:
        from utils import check_ollama_connection
        is_connected, models = check_ollama_connection()
        
        if is_connected:
            print("✅ Ollama connection: OK")
            print(f"  Available models: {models}")
            
            # 检查推荐的模型
            recommended_models = ['llama3.1', 'llama3', 'mistral', 'codellama']
            available_recommended = [m for m in models if any(rec in m.lower() for rec in recommended_models)]
            
            if available_recommended:
                print(f"  ✅ Recommended models available: {available_recommended}")
            else:
                print(f"  ⚠️ No recommended models found. Consider installing: ollama pull llama3.1")
        else:
            print("❌ Ollama connection: Failed")
            print("  💡 Try: ollama serve")
            print("  💡 Install model: ollama pull llama3.1")
    
    except Exception as e:
        print(f"❌ LLM availability check failed: {e}")

def provide_recommendations(diagnosis: Dict[str, any]):
    """提供修复建议"""
    print("\n💡 RECOMMENDATIONS:")
    print("-" * 40)
    
    completeness = diagnosis.get('completeness', 0)
    
    if completeness >= 80:
        print("✅ Your pipeline results look good!")
        print("  The QA Agent should work properly.")
    elif completeness >= 50:
        print("⚠️ Your pipeline results are partially complete.")
        print("  The QA Agent will work but with limited functionality.")
        print("  Consider re-running the pipeline to completion.")
    else:
        print("❌ Your pipeline results appear incomplete.")
        print("  The QA Agent may not work properly.")
        print("  Please check if the pipeline completed successfully.")
    
    # 具体建议
    recommendations = diagnosis.get('recommendations', [])
    for i, rec in enumerate(recommendations, 1):
        print(f"  {i}. {rec}")
    
    # 通用建议
    if completeness < 80:
        print("\n📋 General troubleshooting steps:")
        print("  1. Check pipeline_log.txt for error messages")
        print("  2. Ensure the pipeline ran to completion") 
        print("  3. Verify all input files were valid")
        print("  4. Re-run the pipeline if necessary")
    
    print("\n🎯 To test the QA Agent with your current results:")
    print("  python ask_assistant.py")
    print("  python standalone_qa.py")

def main():
    """主函数"""
    print("="*60)
    print("🔍 QA AGENT DIAGNOSTIC TOOL")
    print("="*60)
    
    # 获取结果目录
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        # 自动查找结果目录
        possible_dirs = []
        
        # 检查常见位置
        base_dirs = ["pipeline_results", "results", "output", "."]
        for base_dir in base_dirs:
            if os.path.exists(base_dir):
                if base_dir == ".":
                    possible_dirs.append(base_dir)
                else:
                    # 查找run_*目录
                    items = [os.path.join(base_dir, item) for item in os.listdir(base_dir)]
                    run_dirs = [item for item in items if os.path.isdir(item) and 'run_' in os.path.basename(item)]
                    possible_dirs.extend(run_dirs)
        
        if possible_dirs:
            # 使用最新的目录
            results_dir = max(possible_dirs, key=lambda x: os.path.getmtime(x))
            print(f"🔍 Auto-detected results directory: {results_dir}")
        else:
            print("❌ No results directory found!")
            print("Usage: python diagnose_qa.py [results_directory]")
            return 1
    
    # 运行诊断
    print(f"\n📋 Diagnosing results in: {results_dir}")
    
    # 1. 检查结果目录
    diagnosis = check_results_directory(results_dir)
    
    # 2. 测试问题路由
    test_question_routing(results_dir)
    
    # 3. 测试内容加载
    test_content_loading(results_dir)
    
    # 4. 检查LLM可用性
    check_llm_availability()
    
    # 5. 提供建议
    provide_recommendations(diagnosis)
    
    print("\n" + "="*60)
    print(f"🏁 DIAGNOSIS COMPLETE")
    print(f"📊 Results completeness: {diagnosis.get('completeness', 0):.1f}%")
    
    if diagnosis.get('completeness', 0) >= 50:
        print("✅ QA Agent should be functional")
        return 0
    else:
        print("⚠️ QA Agent may have limited functionality")
        return 1

if __name__ == "__main__":
    exit(main())
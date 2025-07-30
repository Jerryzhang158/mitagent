#!/usr/bin/env python3
"""
简单问答助手 - 最简单的使用方式
自动找到最新的pipeline结果并启动问答
"""

import os
import sys

def quick_start():
    """快速启动问答助手"""
    print("🤖 Quick Start - miRNA Research Assistant")
    print("="*50)
    
    # 自动查找结果目录
    base_dirs = [
        "pipeline_results",
        ".",
        "results", 
        "output"
    ]
    
    results_dir = None
    for base_dir in base_dirs:
        if os.path.exists(base_dir):
            # 查找最新的run_*目录
            try:
                items = [os.path.join(base_dir, item) for item in os.listdir(base_dir)]
                run_dirs = [item for item in items if os.path.isdir(item) and 'run_' in os.path.basename(item)]
                
                if run_dirs:
                    results_dir = max(run_dirs, key=lambda x: os.path.getmtime(x))
                    break
                elif base_dir == "." and any(f.endswith('.txt') or f.endswith('.csv') for f in os.listdir(base_dir)):
                    # 当前目录就是结果目录
                    results_dir = base_dir
                    break
            except:
                continue
    
    if not results_dir:
        print("❌ No pipeline results found!")
        print("\n💡 Please make sure you have:")
        print("  1. Run the miRNA pipeline at least once")
        print("  2. Or place this script in your results directory")
        print("  3. Or use the full version: python standalone_qa.py")
        return 1
    
    print(f"📁 Found results: {results_dir}")
    
    try:
        from intelligent_qa_agent import IntelligentQAAgent
        from utils import Logger
        
        # 创建调试日志
        log_file = os.path.join(results_dir, "qa_quick_debug.log")
        logger = Logger(log_file, console_output=True)
        
        print("🚀 Starting interactive Q&A session...")
        print("="*50)
        
        # 创建并启动QA Agent
        qa_agent = IntelligentQAAgent(results_dir, logger)
        
        # 检查状态
        status = qa_agent.get_agent_status()
        print(f"🧠 LLM Available: {status['llm_available']}")
        
        if not status['llm_available']:
            print("⚠️ LLM not available. Please check Ollama setup.")
            print("💡 Try: ollama serve && ollama pull llama3.1")
        
        qa_agent.interactive_session()
        qa_agent.save_conversation_history()
        
        print(f"\n✅ Session completed!")
        print(f"💾 Chat history: qa_conversation_history.json")
        print(f"🔍 Debug log: {log_file}")
        return 0
        
    except ImportError:
        print("❌ Required files not found in current directory")
        print("\n💡 Please make sure these files are present:")
        print("  - intelligent_qa_agent.py")
        print("  - utils.py")
        print("  - config.py")
        print("  - Other pipeline modules")
        return 1
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        print("🔍 Full error traceback:")
        print(traceback.format_exc())
        return 1

if __name__ == "__main__":
    # 如果有命令行参数，直接传给用户作为问题
    if len(sys.argv) > 1:
        question = " ".join(sys.argv[1:])
        print(f"🤖 Quick Question: {question}")
        print("="*50)
        
        try:
            from intelligent_qa_agent import IntelligentQAAgent
            from utils import Logger
            
            # 自动找到结果目录
            results_dir = None
            for base_dir in ["pipeline_results", ".", "results", "output"]:
                if os.path.exists(base_dir):
                    items = [os.path.join(base_dir, item) for item in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, item))]
                    run_dirs = [item for item in items if 'run_' in os.path.basename(item)]
                    if run_dirs:
                        results_dir = max(run_dirs, key=lambda x: os.path.getmtime(x))
                        break
            
            if not results_dir:
                results_dir = "."  # 最后尝试当前目录
            
            print(f"📁 Using results from: {results_dir}")
            
            # 创建调试日志
            log_file = os.path.join(results_dir, "qa_single_debug.log")
            logger = Logger(log_file, console_output=True)
            
            print("🔍 DEBUG MODE ENABLED - Detailed processing information:")
            print("="*60)
            
            qa_agent = IntelligentQAAgent(results_dir, logger)
            result = qa_agent.ask_question(question)
            
            print("="*60)
            print("📊 FINAL ANSWER:")
            print("-" * 40)
            print(result['answer'])
            print(f"\n🎯 Confidence: {result['confidence']:.1%}")
            print(f"🏷️ Question Type: {result['question_type']}")
            
            if result.get('sources_used'):
                print(f"📁 Sources: {len(result['sources_used'])} files")
            
            print(f"\n🔍 Debug log saved to: {log_file}")
            
        except Exception as e:
            print(f"❌ Error: {e}")
            import traceback
            print("🔍 Full error traceback:")
            print(traceback.format_exc())
            exit(1)
    else:
        exit(quick_start())
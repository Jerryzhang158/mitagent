#!/usr/bin/env python3
"""
模块测试脚本 - 验证所有模块是否可以正常导入和工作
"""
import os
def test_imports():
    """测试所有模块导入"""
    print("🧪 Testing module imports...")
    
    modules_status = {}
    
    # 测试必需模块
    required_modules = [
        ('config', 'PipelineConfig'),
        ('utils', 'Logger'),
        ('mirna_matcher', 'MiRNAMatcher'),
        ('refseq_cache', 'RefSeqManager'),
        ('mti_selection', 'MTISelector'),
        ('literature_mining', 'LiteratureMiner'),
        ('main_pipeline', 'MiRNAPipeline')
    ]
    
    for module_name, class_name in required_modules:
        try:
            module = __import__(module_name)
            cls = getattr(module, class_name)
            modules_status[module_name] = {'status': 'OK', 'class': cls}
            print(f"✅ {module_name}.{class_name}")
        except ImportError as e:
            modules_status[module_name] = {'status': 'MISSING', 'error': str(e)}
            print(f"❌ {module_name}: {e}")
        except AttributeError as e:
            modules_status[module_name] = {'status': 'CLASS_MISSING', 'error': str(e)}
            print(f"⚠️ {module_name}: {e}")
    
    # 测试可选模块
    optional_modules = [
        'mti_llm_summarize',
        'mti_cytoscape_network', 
        'mirna_bert'
    ]
    
    print("\n🔧 Testing optional modules...")
    for module_name in optional_modules:
        try:
            __import__(module_name)
            modules_status[module_name] = {'status': 'OK'}
            print(f"✅ {module_name}")
        except ImportError as e:
            modules_status[module_name] = {'status': 'OPTIONAL_MISSING', 'error': str(e)}
            print(f"⚠️ {module_name}: {e} (optional)")
    
    return modules_status

def test_basic_functionality():
    """测试基本功能"""
    print("\n🧪 Testing basic functionality...")
    
    try:
        # 测试配置
        from config import PipelineConfig
        config = PipelineConfig()
        print(f"✅ Config loaded - Default functions: {config.DEFAULT_FUNCTIONS}")
        
        # 测试Logger
        from utils import Logger
        import tempfile
        import os
        
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            logger = Logger(tmp_file.name, console_output=False)
            logger.info("Test message")
            
            # 检查日志文件
            if os.path.exists(tmp_file.name) and os.path.getsize(tmp_file.name) > 0:
                print("✅ Logger working")
            else:
                print("❌ Logger not working")
            
            # 清理临时文件
            os.unlink(tmp_file.name)
        
        # 测试miRNA匹配器
        from mirna_matcher import MiRNAMatcher
        matcher = MiRNAMatcher()
        normalized = matcher.normalize_mirna_name("hsa-miR-21-5p")
        if normalized:
            print(f"✅ miRNA matcher working - normalized: {normalized}")
        else:
            print("❌ miRNA matcher not working")
        
        # 测试RefSeq管理器
        from refseq_cache import RefSeqManager
        manager = RefSeqManager()
        cache_stats = manager.get_cache_info()
        print(f"✅ RefSeq manager working - cache stats: {cache_stats}")
        
        print("✅ Basic functionality tests passed")
        return True
        
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False

def test_file_utils():
    """测试文件工具"""
    print("\n🧪 Testing file utilities...")
    
    try:
        from utils import FileUtils
        import tempfile
        import pandas as pd
        
        # 创建测试数据
        test_data = pd.DataFrame({
            'Gene': ['GENE1', 'GENE2', 'GENE3'],
            'miRNA': ['hsa-miR-21', 'hsa-miR-155', 'hsa-miR-200c'],
            'Score': [85, 92, 78]
        })
        
        # 测试CSV读写
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_file:
            test_data.to_csv(tmp_file.name, index=False)
            
            # 使用FileUtils读取
            df, metadata = FileUtils.smart_read_csv(tmp_file.name)
            
            if len(df) == 3 and 'Gene' in df.columns:
                print(f"✅ FileUtils working - detected delimiter: '{metadata['delimiter']}'")
            else:
                print("❌ FileUtils not working properly")
            
            # 清理
            import os
            os.unlink(tmp_file.name)
        
        return True
        
    except Exception as e:
        print(f"❌ File utils test failed: {e}")
        return False

def test_pipeline_initialization():
    """测试Pipeline初始化"""
    print("\n🧪 Testing pipeline initialization...")
    
    try:
        from main_pipeline import MiRNAPipeline
        import tempfile
        import shutil
        
        # 创建临时输出目录
        with tempfile.TemporaryDirectory() as temp_dir:
            pipeline = MiRNAPipeline(output_dir=temp_dir, min_databases=2)
            
            # 检查结果目录是否创建
            if hasattr(pipeline, 'results_dir') and os.path.exists(pipeline.results_dir):
                print("✅ Pipeline initialization working")
                
                # 检查子目录
                subdirs = ['mirna_selection', 'pubmed_articles', 'bert_validation', 
                          'llm_summaries', 'cytoscape_network']
                all_dirs_exist = all(os.path.exists(os.path.join(pipeline.results_dir, subdir)) 
                                   for subdir in subdirs)
                
                if all_dirs_exist:
                    print("✅ All output subdirectories created")
                else:
                    print("⚠️ Some output subdirectories missing")
                
                return True
            else:
                print("❌ Pipeline initialization failed")
                return False
        
    except Exception as e:
        print(f"❌ Pipeline initialization test failed: {e}")
        return False

def main():
    """主测试函数"""
    print("="*60)
    print("🧬 ENHANCED MIRNA PIPELINE - MODULE TESTS")
    print("="*60)
    
    # 运行所有测试
    tests_passed = 0
    total_tests = 4
    
    # 测试1: 模块导入
    modules_status = test_imports()
    required_modules_ok = all(status['status'] == 'OK' 
                            for module, status in modules_status.items() 
                            if 'mti_' not in module and 'mirna_bert' not in module)
    if required_modules_ok:
        tests_passed += 1
        print("\n✅ Module import test: PASSED")
    else:
        print("\n❌ Module import test: FAILED")
    
    # 测试2: 基本功能
    if test_basic_functionality():
        tests_passed += 1
        print("✅ Basic functionality test: PASSED")
    else:
        print("❌ Basic functionality test: FAILED")
    
    # 测试3: 文件工具
    if test_file_utils():
        tests_passed += 1
        print("✅ File utilities test: PASSED")
    else:
        print("❌ File utilities test: FAILED")
    
    # 测试4: Pipeline初始化
    if test_pipeline_initialization():
        tests_passed += 1
        print("✅ Pipeline initialization test: PASSED")
    else:
        print("❌ Pipeline initialization test: FAILED")
    
    # 总结
    print("\n" + "="*60)
    print(f"📊 TEST SUMMARY: {tests_passed}/{total_tests} tests passed")
    print("="*60)
    
    if tests_passed == total_tests:
        print("🎉 All tests passed! The modular pipeline is ready to use.")
        
        print("\n💡 Next steps:")
        print("1. Prepare your data files (genes.txt, mirnas.txt, database files)")
        print("2. Install optional dependencies if needed:")
        print("   - For BERT: pip install torch transformers")
        print("   - For LLM: Install and start Ollama")
        print("3. Run the pipeline: python main_pipeline.py --help")
        
        return 0
    else:
        print("⚠️ Some tests failed. Please check the error messages above.")
        print("💡 You may still be able to use the pipeline with reduced functionality.")
        return 1

if __name__ == "__main__":
    exit(main())
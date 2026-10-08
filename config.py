"""
配置文件 - miRNA Research Pipeline
包含所有默认配置参数和路径设置
"""

import os

class PipelineConfig:
    """Pipeline配置类"""

    # 默认参数
    DEFAULT_MIN_DATABASES = 2
    DEFAULT_MAX_ARTICLES = 50
    DEFAULT_SCORE_THRESHOLD = 50
    DEFAULT_SIMILARITY_THRESHOLD = 0.8
    DEFAULT_BATCH_SIZE = 200

    # 默认功能列表
    DEFAULT_FUNCTIONS = ["cell proliferation", "apoptosis", "migration", "invasion"]

    # 文件路径配置
    CACHE_DIR = "cache"
    REFSEQ_CACHE_FILE = "refseq_cache.json"

    # API配置
    MYGENE_BASE_URL = "https://mygene.info/v3/query"
    PUBMED_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    PUBMED_BACKEND = os.environ.get("MITAGENT_PUBMED_BACKEND", "local").lower()
    PUBMED_LOCAL_DB = os.environ.get(
        "MITAGENT_PUBMED_LOCAL_DB",
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "pubmed_local_backend",
            "data",
            "pubmed_mirna_baseline_2026",
            "pubmed_mirna.sqlite",
        ),
    )
    OLLAMA_BASE_URL = "http://localhost:11434"

    # 请求延迟配置 (秒)
    API_DELAY = 0.5
    PUBMED_DELAY = 1.0
    MYGENE_DELAY = 0.5

    # 超时配置 (秒)
    DEFAULT_TIMEOUT = 30
    PUBMED_TIMEOUT = 30
    MYGENE_TIMEOUT = 10
    OLLAMA_TIMEOUT = 5

    # 文件编码
    DEFAULT_ENCODING = 'utf-8'

    # 目录结构
    OUTPUT_SUBDIRS = [
        "mirna_selection",
        "pubmed_articles",
        "bert_validation",
        "llm_summaries",
        "cytoscape_network"
    ]

    # miRNA命名模式
    MIRNA_PATTERNS = {
        'version_pattern': r'\.\d+$',
        'core_pattern': r'(mir-\d+[a-z]*(?:-\d+[a-z]*)?)',
        'let_pattern': r'(let-\d+[a-z]*(?:-\d+[a-z]*)?)',
        'human_prefix': 'hsa-'
    }

    # 数据库列映射
    DATABASE_COLUMNS = {
        'targetscan': {
            'mirna': 'miRNA',
            'gene': 'Gene Symbol'
        },
        'mirdb': {
            'mirna': 'miRNA',
            'refseq': 'RefSeq',
            'score': 'Score'
        },
        'mirwalk': {
            'mirna': 'miRNA',
            'gene': 'Genesymbol',
            'probability': 'binding_probability'
        },
        'mirtarbase': {
            'mirna': 'miRNA',
            'gene': 'Target Gene',
            'support_type': 'Support Type',
            'experiments': 'Experiments',
            'pmids': 'References (PMID)'
        }
    }

    # miRTarBase列名映射（用于自动修复）
    MIRTARBASE_COLUMN_MAPPING = {
        'miRNA_name': 'miRNA',
        'mirna': 'miRNA',
        'Target_Gene': 'Target Gene',
        'Gene_Symbol': 'Target Gene',
        'target_gene': 'Target Gene',
        'Target_Gene_Symbol': 'Target Gene',
        'Gene': 'Target Gene'
    }

    # 验证强度配置
    VALIDATION_STRENGTH = {
        'strong_keywords': ['Functional MTI'],
        'weak_keywords': ['Functional MTI (Weak)'],
        'negative_keywords': ['Non-Functional MTI']
    }

    @classmethod
    def get_cache_path(cls):
        """获取缓存文件完整路径"""
        os.makedirs(cls.CACHE_DIR, exist_ok=True)
        return os.path.join(cls.CACHE_DIR, cls.REFSEQ_CACHE_FILE)

    @classmethod
    def create_output_structure(cls, base_dir):
        """创建输出目录结构"""
        os.makedirs(base_dir, exist_ok=True)
        for subdir in cls.OUTPUT_SUBDIRS:
            os.makedirs(os.path.join(base_dir, subdir), exist_ok=True)
        return base_dir

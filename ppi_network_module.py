#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPI网络分析模块
用于RNA-seq分析结果的蛋白质相互作用网络构建和分析

主要功能：
1. STRING数据库API调用
2. PPI网络构建
3. Hub基因识别
4. 网络可视化
5. 功能模块检测

依赖安装：
pip install networkx pyvis requests matplotlib pandas numpy
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
import warnings
import community
warnings.filterwarnings('ignore')

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class PPINetworkAnalyzer:
    """PPI网络分析类"""
    
    def __init__(self, species=9606):
        """
        初始化PPI分析器
        
        Args:
            species: 物种ID (9606=human, 10090=mouse, 10116=rat)
        """
        self.species = species
        self.string_api_url = "https://string-db.org/api"
        self.version = "11.5"
        self.network = None
        self.hub_genes = {}
        self.modules = {}
        
        # 默认参数
        self.default_params = {
            'required_score': 400,  # STRING置信度阈值 (0-1000)
            'network_type': 'functional',  # 功能网络
            'add_nodes': 0,  # 不添加额外节点
            'show_query_node_labels': 0
        }
        
        logger.info(f"PPI分析器初始化完成 (物种: {species})")
    
    def get_string_ids(self, genes: List[str]) -> Dict[str, str]:
        """
        将基因符号转换为STRING ID
        
        Args:
            genes: 基因符号列表
            
        Returns:
            基因符号到STRING ID的映射
        """
        logger.info(f"转换 {len(genes)} 个基因符号为STRING ID...")
        
        # 批量查询，每次最多100个
        batch_size = 100
        all_mappings = {}
        
        for i in range(0, len(genes), batch_size):
            batch = genes[i:i+batch_size]
            
            params = {
                'identifiers': '\r'.join(batch),
                'species': self.species,
                'limit': 1,  # 每个基因只返回最佳匹配
                'echo_query': 1,
                'format': 'json'
            }
            
            try:
                response = requests.post(f"{self.string_api_url}/json/get_string_ids", data=params)
                
                if response.status_code == 200:
                    results = response.json()
                    for item in results:
                        query = item.get('queryItem', '')
                        string_id = item.get('stringId', '')
                        preferred_name = item.get('preferredName', query)
                        
                        if query and string_id:
                            all_mappings[query] = {
                                'string_id': string_id,
                                'preferred_name': preferred_name
                            }
                else:
                    logger.warning(f"STRING API返回错误: {response.status_code}")
                    
            except Exception as e:
                logger.error(f"获取STRING ID失败: {e}")
        
        logger.info(f"成功映射 {len(all_mappings)}/{len(genes)} 个基因")
        return all_mappings
    
    def fetch_ppi_network(self, genes: List[str], score_threshold: int = 400) -> pd.DataFrame:
        """
        从STRING数据库获取PPI网络
        
        Args:
            genes: 基因符号列表
            score_threshold: 置信度阈值 (0-1000)
            
        Returns:
            PPI边的DataFrame
        """
        logger.info(f"从STRING获取PPI网络 (基因数: {len(genes)}, 阈值: {score_threshold})...")
        
        if len(genes) == 0:
            logger.error("基因列表为空")
            return pd.DataFrame()
        
        # 处理基因列表
        clean_genes = []
        for gene in genes[:500]:  # 限制最多500个基因，避免API超时
            if gene and str(gene) != 'nan':
                # 清理基因名
                gene_clean = str(gene).strip().upper()
                if not gene_clean.startswith('ENSG'):  # 跳过Ensembl ID
                    clean_genes.append(gene_clean)
        
        if len(clean_genes) == 0:
            logger.error("没有有效的基因符号")
            return pd.DataFrame()
        
        # 构建请求参数
        params = {
            'identifiers': '%0d'.join(clean_genes),
            'species': self.species,
            'required_score': score_threshold,
            'network_type': 'functional',
            'add_nodes': 0,  # 不添加额外节点
            'format': 'tsv'
        }
        
        try:
            # 获取网络数据
            response = requests.post(f"{self.string_api_url}/tsv/network", data=params)
            
            if response.status_code == 200:
                # 解析TSV数据
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    # 解析header
                    header = lines[0].split('\t')
                    
                    # 解析数据行
                    data = []
                    for line in lines[1:]:
                        values = line.split('\t')
                        if len(values) == len(header):
                            data.append(values)
                    
                    if data:
                        df = pd.DataFrame(data, columns=header)
                        
                        # 转换得分为数值
                        if 'score' in df.columns:
                            df['score'] = pd.to_numeric(df['score'], errors='coerce')
                        
                        logger.info(f"获得 {len(df)} 条PPI相互作用")
                        return df
                    else:
                        logger.warning("没有找到PPI相互作用")
                        return pd.DataFrame()
                else:
                    logger.warning("STRING返回空结果")
                    return pd.DataFrame()
            else:
                logger.error(f"STRING API错误: {response.status_code}")
                logger.error(f"响应: {response.text[:500]}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"获取PPI网络失败: {e}")
            return pd.DataFrame()
    
    def build_network(self, ppi_data: pd.DataFrame) -> nx.Graph:
        """
        构建NetworkX图对象
        
        Args:
            ppi_data: PPI数据DataFrame
            
        Returns:
            NetworkX图对象
        """
        if ppi_data.empty:
            logger.warning("PPI数据为空，返回空网络")
            return nx.Graph()
        
        logger.info("构建网络图...")
        
        G = nx.Graph()
        
        # 添加边和权重
        for _, row in ppi_data.iterrows():
            # 获取节点名称
            node1 = row.get('preferredName_A', row.get('stringId_A', 'Unknown_A'))
            node2 = row.get('preferredName_B', row.get('stringId_B', 'Unknown_B'))
            score = row.get('score', 0)
            
            if node1 and node2:
                G.add_edge(node1, node2, weight=float(score)/1000.0)  # 归一化得分
        
        self.network = G
        logger.info(f"网络构建完成: {G.number_of_nodes()} 个节点, {G.number_of_edges()} 条边")
        
        return G
    
    def calculate_centrality_measures(self, G: nx.Graph = None) -> pd.DataFrame:
        """
        计算网络中心性指标
        
        Args:
            G: NetworkX图对象
            
        Returns:
            包含中心性指标的DataFrame
        """
        if G is None:
            G = self.network
        
        if G is None or G.number_of_nodes() == 0:
            logger.warning("网络为空")
            return pd.DataFrame()
        
        logger.info("计算中心性指标...")
        
        # 计算各种中心性
        metrics = {
            'degree': nx.degree_centrality(G),
            'betweenness': nx.betweenness_centrality(G),
            'closeness': nx.closeness_centrality(G),
            'eigenvector': nx.eigenvector_centrality_numpy(G) if G.number_of_nodes() > 2 else {}
        }
        
        # 转换为DataFrame
        df = pd.DataFrame(metrics).fillna(0)
        
        # 添加度数（连接数）
        df['degree_count'] = pd.Series(dict(G.degree()))
        
        # 计算综合得分
        df['hub_score'] = (
            df['degree'] * 0.3 +
            df['betweenness'] * 0.3 +
            df['closeness'] * 0.2 +
            df['eigenvector'] * 0.2
        )
        
        # 排序
        df = df.sort_values('hub_score', ascending=False)
        
        logger.info(f"计算完成，识别出 {len(df)} 个节点的中心性")
        
        return df
    
    def identify_hub_genes(self, centrality_df: pd.DataFrame, 
                          method: str = 'top_percent',
                          threshold: float = 0.1) -> List[str]:
        """
        识别hub基因
        
        Args:
            centrality_df: 中心性指标DataFrame
            method: 识别方法 ('top_percent', 'top_n', 'threshold')
            threshold: 阈值参数
            
        Returns:
            Hub基因列表
        """
        if centrality_df.empty:
            return []
        
        logger.info(f"识别hub基因 (方法: {method}, 阈值: {threshold})...")
        
        if method == 'top_percent':
            # 选择前x%的基因
            n_hubs = max(1, int(len(centrality_df) * threshold))
            hubs = centrality_df.head(n_hubs).index.tolist()
            
        elif method == 'top_n':
            # 选择前n个基因
            n_hubs = min(int(threshold), len(centrality_df))
            hubs = centrality_df.head(n_hubs).index.tolist()
            
        elif method == 'threshold':
            # 基于hub_score阈值
            hubs = centrality_df[centrality_df['hub_score'] > threshold].index.tolist()
            
        else:
            # 默认选择前10%
            n_hubs = max(1, int(len(centrality_df) * 0.1))
            hubs = centrality_df.head(n_hubs).index.tolist()
        
        self.hub_genes = {
            'genes': hubs,
            'scores': centrality_df.loc[hubs, 'hub_score'].to_dict()
        }
        
        logger.info(f"识别出 {len(hubs)} 个hub基因")
        logger.info(f"Top 5 hubs: {hubs[:5]}")
        
        return hubs
    
    def detect_modules(self, G: nx.Graph = None, resolution: float = 1.0) -> Dict:
        """
        检测网络模块（社区）
        
        Args:
            G: NetworkX图对象
            resolution: 分辨率参数（越高模块越多）
            
        Returns:
            模块检测结果
        """
        if G is None:
            G = self.network
        
        if G is None or G.number_of_nodes() < 3:
            logger.warning("网络节点太少，无法检测模块")
            return {}
        
        logger.info("检测网络模块...")
        
        # 尝试不同的社区检测方法
        communities = None
        
        # 方法1: 尝试使用python-louvain (新版本)
        try:
            import community as community_louvain
            communities = community_louvain.best_partition(G, resolution=resolution)
            logger.info("使用Louvain算法检测社区")
        except ImportError:
            pass
        
        # 方法2: 尝试旧版本python-louvain
        if communities is None:
            try:
                import community
                communities = community.best_partition(G, resolution=resolution)
                logger.info("使用Louvain算法检测社区（旧版本）")
            except ImportError:
                pass
        
        # 方法3: 使用NetworkX内置的贪婪模块度算法
        if communities is None:
            logger.info("使用NetworkX贪婪模块度算法检测社区")
            from networkx.algorithms import community
            communities_generator = community.greedy_modularity_communities(G, resolution=resolution)
            communities = {}
            for i, comm in enumerate(communities_generator):
                for node in comm:
                    communities[node] = i
        
        # 整理模块信息
        modules = {}
        for node, module_id in communities.items():
            if module_id not in modules:
                modules[module_id] = []
            modules[module_id].append(node)
        
        # 计算模块统计
        module_stats = {}
        for module_id, nodes in modules.items():
            subgraph = G.subgraph(nodes)
            module_stats[module_id] = {
                'size': len(nodes),
                'density': nx.density(subgraph) if len(nodes) > 1 else 0,
                'nodes': nodes
            }
        
        self.modules = module_stats
        logger.info(f"检测到 {len(modules)} 个模块")
        
        # 显示最大的几个模块
        sorted_modules = sorted(module_stats.items(), 
                              key=lambda x: x[1]['size'], 
                              reverse=True)
        for i, (mod_id, stats) in enumerate(sorted_modules[:3]):
            logger.info(f"  模块 {mod_id}: {stats['size']} 个节点, 密度: {stats['density']:.3f}")
        
        return module_stats
    
    def visualize_network(self, G: nx.Graph = None, 
                         highlight_nodes: List[str] = None,
                         output_file: str = None,
                         layout: str = 'spring',
                         figsize: Tuple[int, int] = (12, 10)) -> None:
        """
        可视化PPI网络
        
        Args:
            G: NetworkX图对象
            highlight_nodes: 要高亮的节点（如hub基因）
            output_file: 输出文件路径
            layout: 布局算法
            figsize: 图像大小
        """
        if G is None:
            G = self.network
        
        if G is None or G.number_of_nodes() == 0:
            logger.warning("网络为空，无法可视化")
            return
        
        logger.info(f"绘制网络图 (节点: {G.number_of_nodes()}, 边: {G.number_of_edges()})...")
        
        plt.figure(figsize=figsize)
        
        # 选择布局
        if layout == 'spring':
            pos = nx.spring_layout(G, k=1/np.sqrt(G.number_of_nodes()), iterations=50)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        else:
            pos = nx.spring_layout(G)
        
        # 节点大小基于度数
        node_sizes = [300 * G.degree(node) for node in G.nodes()]
        
        # 节点颜色
        if highlight_nodes:
            node_colors = ['red' if node in highlight_nodes else 'lightblue' 
                          for node in G.nodes()]
        else:
            node_colors = 'lightblue'
        
        # 绘制网络
        nx.draw_networkx_nodes(G, pos, 
                              node_size=node_sizes,
                              node_color=node_colors,
                              alpha=0.7)
        
        # 绘制边
        edges = G.edges()
        weights = [G[u][v].get('weight', 1) for u, v in edges]
        nx.draw_networkx_edges(G, pos, 
                              width=[w*2 for w in weights],
                              alpha=0.3)
        
        # 添加标签（只显示重要节点）
        labels = {}
        if G.number_of_nodes() <= 30:
            labels = {node: node for node in G.nodes()}
        else:
            # 只显示度数最高的节点
            degrees = dict(G.degree())
            top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:20]
            labels = {node: node for node, _ in top_nodes}
        
        nx.draw_networkx_labels(G, pos, labels, font_size=8)
        
        plt.title(f"PPI Network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
        plt.axis('off')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"网络图已保存: {output_file}")
        
        plt.show()
    
    def create_interactive_network(self, G: nx.Graph = None,
                                  highlight_nodes: List[str] = None,
                                  output_file: str = 'ppi_network.html') -> None:
        """
        创建交互式网络可视化（HTML）
        
        Args:
            G: NetworkX图对象
            highlight_nodes: 要高亮的节点
            output_file: 输出HTML文件路径
        """
        try:
            from pyvis.network import Network
        except ImportError:
            logger.error("pyvis未安装，无法创建交互式网络")
            logger.info("安装命令: pip install pyvis")
            return
        
        if G is None:
            G = self.network
        
        if G is None or G.number_of_nodes() == 0:
            logger.warning("网络为空")
            return
        
        logger.info("创建交互式网络...")
        
        # 创建pyvis网络
        net = Network(height='750px', width='100%', 
                     bgcolor='#ffffff', font_color='black')
        
        # 设置物理引擎
        net.barnes_hut(gravity=-80000, central_gravity=0.3, 
                       spring_length=100, spring_strength=0.001)
        
        # 添加节点
        for node in G.nodes():
            degree = G.degree(node)
            size = min(50, 10 + degree * 2)  # 基于度数的大小
            
            if highlight_nodes and node in highlight_nodes:
                color = '#ff4444'  # 红色表示hub基因
                title = f"{node} (Hub Gene)\nDegree: {degree}"
            else:
                color = '#97c2fc'  # 浅蓝色
                title = f"{node}\nDegree: {degree}"
            
            net.add_node(node, label=node, title=title, 
                        size=size, color=color)
        
        # 添加边
        for edge in G.edges(data=True):
            weight = edge[2].get('weight', 1)
            net.add_edge(edge[0], edge[1], value=weight)
        
        # 添加控制按钮
        net.show_buttons(filter_=['physics'])
        
        # 保存HTML
        net.save_graph(output_file)
        logger.info(f"交互式网络已保存: {output_file}")
    
    def run_enrichment_on_modules(self, modules: Dict = None, 
                                 organism: str = 'human') -> Dict:
        """
        对每个模块进行功能富集分析
        
        Args:
            modules: 模块字典
            organism: 物种
            
        Returns:
            富集分析结果
        """
        if modules is None:
            modules = self.modules
        
        if not modules:
            logger.warning("没有模块可分析")
            return {}
        
        logger.info("对网络模块进行功能富集分析...")
        
        try:
            import gseapy as gp
        except ImportError:
            logger.error("gseapy未安装，无法进行富集分析")
            return {}
        
        enrichment_results = {}
        
        # 只分析较大的模块
        significant_modules = {k: v for k, v in modules.items() 
                             if v['size'] >= 5}
        
        for module_id, module_info in significant_modules.items():
            genes = module_info['nodes']
            logger.info(f"分析模块 {module_id} ({len(genes)} 个基因)...")
            
            try:
                # GO富集分析
                enr = gp.enrichr(gene_list=genes,
                               gene_sets=['GO_Biological_Process_2021'],
                               organism='Human',
                               outdir=None,
                               cutoff=0.05)
                
                if enr and hasattr(enr, 'results'):
                    enrichment_results[module_id] = enr.results
                    logger.info(f"  模块 {module_id} 富集分析完成")
                    
            except Exception as e:
                logger.error(f"  模块 {module_id} 富集分析失败: {e}")
        
        return enrichment_results
    
    def export_results(self, output_dir: str = 'ppi_results') -> None:
        """
        导出所有分析结果
        
        Args:
            output_dir: 输出目录
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        logger.info(f"导出结果到 {output_path}...")
        
        # 导出网络
        if self.network:
            # 导出边列表
            edge_list = []
            for u, v, data in self.network.edges(data=True):
                edge_list.append({
                    'source': u,
                    'target': v,
                    'weight': data.get('weight', 1)
                })
            
            edge_df = pd.DataFrame(edge_list)
            edge_df.to_csv(output_path / 'network_edges.csv', index=False)
            logger.info("  导出网络边列表")
            
            # 导出节点属性
            if hasattr(self, 'centrality_df'):
                self.centrality_df.to_csv(output_path / 'node_centrality.csv')
                logger.info("  导出节点中心性")
        
        # 导出hub基因
        if self.hub_genes:
            hub_df = pd.DataFrame({
                'gene': self.hub_genes['genes'],
                'hub_score': [self.hub_genes['scores'][g] 
                            for g in self.hub_genes['genes']]
            })
            hub_df.to_csv(output_path / 'hub_genes.csv', index=False)
            logger.info("  导出hub基因")
        
        # 导出模块
        if self.modules:
            with open(output_path / 'modules.json', 'w') as f:
                # 转换为可序列化格式
                modules_export = {}
                for mod_id, mod_info in self.modules.items():
                    modules_export[str(mod_id)] = {
                        'size': mod_info['size'],
                        'density': mod_info['density'],
                        'nodes': mod_info['nodes']
                    }
                json.dump(modules_export, f, indent=2)
            logger.info("  导出网络模块")
        
        logger.info(f"所有结果已导出到 {output_path}")


def analyze_de_genes_ppi(de_results_file: str, 
                         output_dir: str = 'ppi_analysis',
                         padj_threshold: float = 0.05,
                         log2fc_threshold: float = 1.0,
                         string_score: int = 400) -> Dict:
    """
    便捷函数：从差异表达结果直接进行PPI分析
    
    Args:
        de_results_file: 差异表达结果CSV文件路径
        output_dir: 输出目录
        padj_threshold: 校正p值阈值
        log2fc_threshold: log2 fold change阈值
        string_score: STRING置信度阈值
        
    Returns:
        分析结果字典
    """
    logger.info("="*60)
    logger.info("开始PPI网络分析")
    logger.info("="*60)
    
    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # 读取差异表达结果
    logger.info(f"读取差异表达结果: {de_results_file}")
    de_df = pd.read_csv(de_results_file, index_col=0)
    
    # 筛选显著基因
    significant = de_df[
        (de_df['padj'] < padj_threshold) & 
        (abs(de_df['log2FoldChange']) >= log2fc_threshold)
    ]
    
    logger.info(f"显著差异基因: {len(significant)} 个")
    
    if len(significant) == 0:
        logger.error("没有显著差异基因")
        return {}
    
    # 获取基因列表
    gene_list = significant.index.tolist()
    
    # 创建PPI分析器
    ppi = PPINetworkAnalyzer()
    
    # 获取PPI网络
    ppi_data = ppi.fetch_ppi_network(gene_list, string_score)
    
    if ppi_data.empty:
        logger.error("未获得PPI数据")
        return {}
    
    # 构建网络
    G = ppi.build_network(ppi_data)
    
    # 计算中心性
    centrality_df = ppi.calculate_centrality_measures(G)
    
    # 识别hub基因
    hub_genes = ppi.identify_hub_genes(centrality_df, method='top_n', threshold=20)
    
    # 检测模块
    modules = ppi.detect_modules(G)
    
    # 可视化
    ppi.visualize_network(G, 
                         highlight_nodes=hub_genes,
                         output_file=str(output_path / 'ppi_network.png'))
    
    # 创建交互式网络
    ppi.create_interactive_network(G,
                                  highlight_nodes=hub_genes,
                                  output_file=str(output_path / 'ppi_network_interactive.html'))
    
    # 导出结果
    ppi.export_results(str(output_path))
    
    logger.info("="*60)
    logger.info("PPI分析完成")
    logger.info("="*60)
    
    return {
        'network': G,
        'centrality': centrality_df,
        'hub_genes': hub_genes,
        'modules': modules,
        'ppi_analyzer': ppi
    }


if __name__ == "__main__":
    # 测试代码
    import argparse
    
    parser = argparse.ArgumentParser(description='PPI网络分析模块')
    parser.add_argument('--de-file', required=True, help='差异表达结果文件')
    parser.add_argument('--output', default='ppi_results', help='输出目录')
    parser.add_argument('--padj', type=float, default=0.05, help='padj阈值')
    parser.add_argument('--log2fc', type=float, default=1.0, help='log2FC阈值')
    parser.add_argument('--string-score', type=int, default=400, help='STRING置信度')
    
    args = parser.parse_args()
    
    # 运行分析
    results = analyze_de_genes_ppi(
        args.de_file,
        args.output,
        args.padj,
        args.log2fc,
        args.string_score
    )
    
    print(f"分析完成！结果保存在: {args.output}")
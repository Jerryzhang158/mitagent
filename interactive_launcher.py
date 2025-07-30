#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
纯Python版RNA-seq和miRNA-seq分析程序交互式启动器
简化用户操作，提供图形化界面选择，无需R环境
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import subprocess
import sys
from pathlib import Path

class RNASeqLauncher:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("纯Python版 RNA-seq & miRNA-seq 差异表达分析程序")
        self.root.geometry("800x600")
        
        # 变量
        self.gene_file = tk.StringVar()
        self.mirna_file = tk.StringVar()
        self.output_dir = tk.StringVar(value="results")
        self.control_group = tk.StringVar(value="Control")
        
        # 处理组选择变量
        self.treatment_vars = {
            'G_Low': tk.BooleanVar(value=True),
            'G_Medium': tk.BooleanVar(value=True), 
            'G_High': tk.BooleanVar(value=True),
            'D_Low': tk.BooleanVar(value=True),
            'D_Medium': tk.BooleanVar(value=True),
            'D_High': tk.BooleanVar(value=True)
        }
        
        self.setup_ui()
        
    def setup_ui(self):
        """设置用户界面"""
        # 标题
        title_label = tk.Label(self.root, text="纯Python版 RNA-seq & miRNA-seq 差异表达分析", 
                              font=("Arial", 16, "bold"))
        title_label.pack(pady=10)
        
        # 主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # 文件选择区域
        file_frame = ttk.LabelFrame(main_frame, text="数据文件选择", padding="10")
        file_frame.pack(fill=tk.X, pady=5)
        
        # 基因表达文件
        gene_frame = ttk.Frame(file_frame)
        gene_frame.pack(fill=tk.X, pady=2)
        ttk.Label(gene_frame, text="基因表达数据:").pack(side=tk.LEFT)
        ttk.Entry(gene_frame, textvariable=self.gene_file, width=50).pack(side=tk.LEFT, padx=5)
        ttk.Button(gene_frame, text="浏览", 
                  command=lambda: self.browse_file(self.gene_file, "基因表达数据", 
                                                 [("Text files", "*.txt"), ("All files", "*.*")])).pack(side=tk.LEFT)
        
        # miRNA表达文件
        mirna_frame = ttk.Frame(file_frame)
        mirna_frame.pack(fill=tk.X, pady=2)
        ttk.Label(mirna_frame, text="miRNA表达数据:").pack(side=tk.LEFT)
        ttk.Entry(mirna_frame, textvariable=self.mirna_file, width=50).pack(side=tk.LEFT, padx=5)
        ttk.Button(mirna_frame, text="浏览", 
                  command=lambda: self.browse_file(self.mirna_file, "miRNA表达数据", 
                                                 [("CSV files", "*.csv"), ("All files", "*.*")])).pack(side=tk.LEFT)
        
        # 输出目录
        output_frame = ttk.Frame(file_frame)
        output_frame.pack(fill=tk.X, pady=2)
        ttk.Label(output_frame, text="输出目录:").pack(side=tk.LEFT)
        ttk.Entry(output_frame, textvariable=self.output_dir, width=50).pack(side=tk.LEFT, padx=5)
        ttk.Button(output_frame, text="浏览", 
                  command=self.browse_directory).pack(side=tk.LEFT)
        
        # 分析参数区域
        param_frame = ttk.LabelFrame(main_frame, text="分析参数", padding="10")
        param_frame.pack(fill=tk.X, pady=5)
        
        # 对照组设置
        control_frame = ttk.Frame(param_frame)
        control_frame.pack(fill=tk.X, pady=2)
        ttk.Label(control_frame, text="对照组:").pack(side=tk.LEFT)
        control_combo = ttk.Combobox(control_frame, textvariable=self.control_group,
                                   values=['Control', 'G_Low', 'G_Medium', 'G_High', 
                                          'D_Low', 'D_Medium', 'D_High'])
        control_combo.pack(side=tk.LEFT, padx=5)
        
        # 分析选项区域
        option_frame = ttk.LabelFrame(param_frame, text="分析选项", padding="5")
        option_frame.pack(fill=tk.X, pady=5)
        
        self.skip_enrichment = tk.BooleanVar(value=False)
        skip_cb = ttk.Checkbutton(option_frame, text="跳过GO/KEGG富集分析（加快运行速度）", 
                                 variable=self.skip_enrichment)
        skip_cb.pack(anchor='w', pady=2)
        
        # 处理组选择
        treatment_frame = ttk.LabelFrame(param_frame, text="处理组选择", padding="5")
        treatment_frame.pack(fill=tk.X, pady=5)
        
        treatment_grid = ttk.Frame(treatment_frame)
        treatment_grid.pack()
        
        row = 0
        col = 0
        for treatment, var in self.treatment_vars.items():
            cb = ttk.Checkbutton(treatment_grid, text=treatment, variable=var)
            cb.grid(row=row, column=col, sticky='w', padx=10, pady=2)
            col += 1
            if col > 2:
                col = 0
                row += 1
        
        # 快捷按钮
        quick_frame = ttk.Frame(treatment_frame)
        quick_frame.pack(pady=5)
        ttk.Button(quick_frame, text="选择全部", command=self.select_all_treatments).pack(side=tk.LEFT, padx=5)
        ttk.Button(quick_frame, text="清除全部", command=self.clear_all_treatments).pack(side=tk.LEFT, padx=5)
        ttk.Button(quick_frame, text="只选G组", command=self.select_g_treatments).pack(side=tk.LEFT, padx=5)
        ttk.Button(quick_frame, text="只选D组", command=self.select_d_treatments).pack(side=tk.LEFT, padx=5)
        
        # 进度和日志区域
        log_frame = ttk.LabelFrame(main_frame, text="运行日志", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # 文本框和滚动条
        text_frame = ttk.Frame(log_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        self.log_text = tk.Text(text_frame, height=10, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # 控制按钮
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(button_frame, text="检查环境", command=self.check_environment).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="开始分析", command=self.start_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清除日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="退出", command=self.root.quit).pack(side=tk.RIGHT, padx=5)
        
        # 初始日志
        self.log("欢迎使用纯Python版RNA-seq & miRNA-seq差异表达分析程序!")
        self.log("无需R环境，配置简单，功能强大。")
        self.log("请选择数据文件并设置分析参数。")
        
    def browse_file(self, var, title, filetypes):
        """浏览文件"""
        filename = filedialog.askopenfilename(title=f"选择{title}", filetypes=filetypes)
        if filename:
            var.set(filename)
            self.log(f"已选择{title}: {os.path.basename(filename)}")
    
    def browse_directory(self):
        """浏览目录"""
        directory = filedialog.askdirectory(title="选择输出目录")
        if directory:
            self.output_dir.set(directory)
            self.log(f"输出目录设置为: {directory}")
    
    def select_all_treatments(self):
        """选择所有处理组"""
        for var in self.treatment_vars.values():
            var.set(True)
        self.log("已选择所有处理组")
    
    def clear_all_treatments(self):
        """清除所有处理组选择"""
        for var in self.treatment_vars.values():
            var.set(False)
        self.log("已清除所有处理组选择")
    
    def select_g_treatments(self):
        """只选择G组处理"""
        for name, var in self.treatment_vars.items():
            var.set(name.startswith('G_'))
        self.log("已选择G组处理组")
    
    def select_d_treatments(self):
        """只选择D组处理"""
        for name, var in self.treatment_vars.items():
            var.set(name.startswith('D_'))
        self.log("已选择D组处理组")
    
    def log(self, message):
        """添加日志信息"""
        self.log_text.insert(tk.END, f"{message}\n")
        self.log_text.see(tk.END)
        self.root.update()
    
    def clear_log(self):
        """清除日志"""
        self.log_text.delete(1.0, tk.END)
    
    def check_environment(self):
        """检查运行环境"""
        self.log("正在检查运行环境...")
        
        # 检查Python包
        required_packages = ['pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels', 'gseapy']
        missing_packages = []
        
        for package in required_packages:
            try:
                __import__(package)
                self.log(f"✓ {package} 已安装")
            except ImportError:
                missing_packages.append(package)
                self.log(f"✗ {package} 未安装")
        
        # 检查可选包
        optional_packages = ['sklearn']
        for package in optional_packages:
            try:
                __import__(package)
                self.log(f"✓ {package} 已安装 (可选)")
            except ImportError:
                self.log(f"△ {package} 未安装 (可选，用于高级标准化)")
        
        if missing_packages:
            self.log(f"缺失的必需包: {', '.join(missing_packages)}")
            self.log("安装命令: pip install " + " ".join(missing_packages))
            self.log("请参考使用说明安装缺失的包")
        else:
            self.log("✓ 环境检查完成，所有必需依赖都已安装")
    
    def validate_inputs(self):
        """验证输入"""
        if not self.gene_file.get() and not self.mirna_file.get():
            messagebox.showerror("错误", "请至少选择一个数据文件（基因表达或miRNA表达）")
            return False
        
        if self.gene_file.get() and not os.path.exists(self.gene_file.get()):
            messagebox.showerror("错误", "基因表达数据文件不存在")
            return False
            
        if self.mirna_file.get() and not os.path.exists(self.mirna_file.get()):
            messagebox.showerror("错误", "miRNA表达数据文件不存在")
            return False
        
        # 检查是否选择了处理组
        selected_treatments = [name for name, var in self.treatment_vars.items() if var.get()]
        if not selected_treatments:
            messagebox.showerror("错误", "请至少选择一个处理组")
            return False
        
        return True
    
    def start_analysis(self):
        """开始分析"""
        if not self.validate_inputs():
            return
        
        self.log("开始准备分析...")
        
        # 构建命令
        cmd = [sys.executable, "rnaseq_analyzer.py"]
        
        if self.gene_file.get():
            cmd.extend(["--gene-file", self.gene_file.get()])
        
        if self.mirna_file.get():
            cmd.extend(["--mirna-file", self.mirna_file.get()])
        
        cmd.extend(["--output-dir", self.output_dir.get()])
        cmd.extend(["--control", self.control_group.get()])
        
        # 添加选中的处理组
        selected_treatments = [name for name, var in self.treatment_vars.items() if var.get()]
        cmd.extend(["--treatment"] + selected_treatments)
        
        # 添加分析选项
        if self.skip_enrichment.get():
            cmd.extend(["--skip-enrichment"])
        
        self.log(f"执行命令: {' '.join(cmd)}")
        
        # 创建输出目录
        Path(self.output_dir.get()).mkdir(exist_ok=True)
        
        try:
            # 运行分析程序
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                     stderr=subprocess.STDOUT, 
                                     universal_newlines=True, bufsize=1)
            
            # 实时显示输出
            for line in iter(process.stdout.readline, ''):
                self.log(line.strip())
            
            process.wait()
            
            if process.returncode == 0:
                self.log("✓ 分析完成！")
                self.log(f"结果已保存在: {self.output_dir.get()}")
                
                # 询问是否打开结果目录
                if messagebox.askyesno("完成", "分析完成！是否打开结果目录？"):
                    if sys.platform.startswith('darwin'):  # macOS
                        subprocess.call(['open', self.output_dir.get()])
                    elif sys.platform.startswith('win'):   # Windows
                        subprocess.call(['explorer', self.output_dir.get()])
                    else:  # Linux
                        subprocess.call(['xdg-open', self.output_dir.get()])
            else:
                self.log("✗ 分析过程中出现错误")
                
        except Exception as e:
            self.log(f"✗ 启动分析失败: {str(e)}")
            messagebox.showerror("错误", f"启动分析失败: {str(e)}")
    
    def run(self):
        """运行程序"""
        self.root.mainloop()


def main():
    """主函数"""
    # 检查是否存在主程序文件
    if not os.path.exists("rnaseq_analyzer.py"):
        messagebox.showerror("错误", "找不到主程序文件 rnaseq_analyzer.py\n请确保两个文件在同一目录下")
        return
    
    app = RNASeqLauncher()
    app.run()


if __name__ == "__main__":
    main()
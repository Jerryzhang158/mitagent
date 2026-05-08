"""
visualize.py — 可视化模块
- 水平条形图（Top-K 关键词频率）
- 词云图
- 年度趋势折线图
"""

import logging
import os
from collections import Counter
from pathlib import Path

import config

log = logging.getLogger(__name__)

COLORS = [
    "#2563EB", "#16A34A", "#DC2626", "#9333EA",
    "#EA580C", "#0891B2", "#65A30D", "#DB2777",
]


def _ensure_dir(output_dir: str):
    Path(output_dir).mkdir(parents=True, exist_ok=True)


# ── 条形图 ────────────────────────────────────────────────────────────────────

def plot_barplot(
    counter: Counter,
    title: str = "Top Keywords",
    filename: str = "barplot.png",
    output_dir: str = None,
    top_k: int = None,
    color: str = "#2563EB",
):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        log.error("matplotlib 未安装")
        return

    output_dir = output_dir or config.OUTPUT_DIR
    _ensure_dir(output_dir)
    top_k = top_k or config.TOP_K_DISPLAY

    items = counter.most_common(top_k)
    if not items:
        return

    labels = [kw for kw, _ in reversed(items)]
    values = [cnt for _, cnt in reversed(items)]

    fig_height = max(8, top_k * 0.28)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    bars = ax.barh(labels, values, color=color, alpha=0.85, edgecolor="white", linewidth=0.5)

    # 数值标注
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width() + max(values) * 0.005,
            bar.get_y() + bar.get_height() / 2,
            str(val), va="center", ha="left", fontsize=8, color="#555"
        )

    ax.set_xlabel("Frequency (# documents)", fontsize=11)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=9)
    ax.set_xlim(0, max(values) * 1.12)

    plt.tight_layout()
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"条形图已保存: {path}")


def plot_all_barplots(method_counters: dict[str, Counter], borda: Counter, output_dir: str = None):
    output_dir = output_dir or config.OUTPUT_DIR
    for i, (method, counter) in enumerate(method_counters.items()):
        plot_barplot(
            counter,
            title=f"Top {config.TOP_K_DISPLAY} Keywords — {method.upper()}",
            filename=f"barplot_{method}.png",
            output_dir=output_dir,
            color=COLORS[i % len(COLORS)],
        )
    # Borda 融合图
    plot_barplot(
        borda,
        title=f"Top {config.TOP_K_DISPLAY} Keywords — Borda Merged",
        filename="barplot_borda.png",
        output_dir=output_dir,
        color="#374151",
    )


# ── 词云图 ────────────────────────────────────────────────────────────────────

def plot_wordcloud(
    counter: Counter,
    title: str = "Keyword Word Cloud",
    filename: str = "wordcloud.png",
    output_dir: str = None,
):
    try:
        from wordcloud import WordCloud
        import matplotlib.pyplot as plt
    except ImportError:
        log.error("wordcloud 未安装，跳过词云图")
        return

    output_dir = output_dir or config.OUTPUT_DIR
    _ensure_dir(output_dir)

    wc = WordCloud(
        width=1600, height=900,
        background_color="white",
        colormap="Blues",
        max_words=200,
        prefer_horizontal=0.85,
        collocations=False,
    ).generate_from_frequencies(dict(counter))

    fig, ax = plt.subplots(figsize=(16, 9))
    ax.imshow(wc, interpolation="bilinear")
    ax.axis("off")
    ax.set_title(title, fontsize=16, fontweight="bold", pad=10)

    plt.tight_layout()
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"词云图已保存: {path}")


# ── 年度趋势折线图 ─────────────────────────────────────────────────────────────

def plot_trend(
    trend: dict[str, dict[str, int]],
    title: str = "Keyword Trend (by Year)",
    filename: str = "trend.png",
    output_dir: str = None,
):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import numpy as np
    except ImportError:
        log.error("matplotlib 未安装")
        return

    output_dir = output_dir or config.OUTPUT_DIR
    _ensure_dir(output_dir)

    if not trend:
        return

    # 取所有年份并排序
    all_years = sorted({yr for counts in trend.values() for yr in counts})
    if len(all_years) < 2:
        log.warning("年份数据不足，跳过趋势图")
        return

    fig, ax = plt.subplots(figsize=(14, 7))
    cmap = cm.get_cmap("tab20", len(trend))

    for i, (kw, yr_counts) in enumerate(trend.items()):
        y_vals = [yr_counts.get(yr, 0) for yr in all_years]
        ax.plot(all_years, y_vals, marker="o", linewidth=2,
                markersize=5, label=kw, color=cmap(i), alpha=0.85)

    ax.set_xlabel("Year", fontsize=11)
    ax.set_ylabel("Frequency (# documents)", fontsize=11)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, ncol=1)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"趋势图已保存: {path}")

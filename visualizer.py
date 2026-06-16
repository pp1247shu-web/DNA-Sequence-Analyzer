"""
visualizer.py
-------------
Creates publication-quality charts for nucleotide analysis results.
Uses matplotlib with a clean, professional aesthetic suitable for reports.
"""

import os
from pathlib import Path
from datetime import datetime

import matplotlib
matplotlib.use("Agg")           # non-interactive backend (safe for all environments)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# --------------------------------------------------------------------------- #
#  Color palette — follows standard bioinformatics conventions
# --------------------------------------------------------------------------- #
BASE_COLORS = {
    "A": "#2ECC71",   # green
    "T": "#E74C3C",   # red
    "G": "#3498DB",   # blue
    "C": "#F39C12",   # orange
}


def _ensure_dir(directory: str) -> Path:
    p = Path(directory)
    p.mkdir(parents=True, exist_ok=True)
    return p


# --------------------------------------------------------------------------- #
#  Single-sequence bar chart
# --------------------------------------------------------------------------- #

def plot_base_counts(result: dict, output_dir: str = "charts") -> str:
    """
    Generate a bar chart of nucleotide counts for a single sequence.

    Parameters
    ----------
    result     : output dict from analyzer.analyze_sequence()
    output_dir : directory where the PNG will be saved

    Returns
    -------
    str : absolute path to the saved PNG file
    """
    _ensure_dir(output_dir)

    bases  = ["A", "T", "G", "C"]
    counts = [result["base_counts"][b] for b in bases]
    colors = [BASE_COLORS[b] for b in bases]
    freqs  = [result["nucleotide_frequency"][b] for b in bases]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Nucleotide Analysis — {result['name'][:60]}",
        fontsize=14, fontweight="bold", y=1.01
    )

    # --- Left: absolute counts ---
    ax1 = axes[0]
    bars = ax1.bar(bases, counts, color=colors, edgecolor="white",
                   linewidth=1.2, zorder=3, width=0.5)
    ax1.set_title("Nucleotide Counts", fontsize=12)
    ax1.set_xlabel("Nucleotide", fontsize=11)
    ax1.set_ylabel("Count", fontsize=11)
    ax1.set_ylim(0, max(counts) * 1.18)
    ax1.yaxis.grid(True, linestyle="--", alpha=0.6, zorder=0)
    ax1.set_axisbelow(True)

    for bar, count in zip(bars, counts):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.02,
            str(count), ha="center", va="bottom", fontsize=11, fontweight="bold"
        )

    # --- Right: percentage pie chart ---
    ax2 = axes[1]
    wedges, texts, autotexts = ax2.pie(
        freqs, labels=bases, colors=colors,
        autopct="%1.2f%%", startangle=140,
        textprops={"fontsize": 11},
        wedgeprops={"edgecolor": "white", "linewidth": 1.5}
    )
    for at in autotexts:
        at.set_fontsize(10)
    ax2.set_title("Nucleotide Frequency (%)", fontsize=12)

    # GC / AT annotation below pie
    ax2.text(
        0, -1.35,
        f"GC content: {result['gc_content']}%   |   AT content: {result['at_content']}%",
        ha="center", va="center", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#EBF5FB", edgecolor="#AED6F1")
    )

    plt.tight_layout()

    safe_name = "".join(c if c.isalnum() else "_" for c in result["name"])[:40]
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename   = f"{safe_name}_{timestamp}.png"
    filepath   = Path(output_dir) / filename

    plt.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return str(filepath.resolve())


# --------------------------------------------------------------------------- #
#  Multi-sequence comparison chart
# --------------------------------------------------------------------------- #

def plot_gc_comparison(results: list[dict], output_dir: str = "charts") -> str:
    """
    Horizontal bar chart comparing GC% across multiple sequences.
    Ideal for comparing genomes or genes side-by-side.

    Returns
    -------
    str : absolute path to the saved PNG file
    """
    _ensure_dir(output_dir)

    names   = [r["name"][:35] for r in results]
    gc_vals = [r["gc_content"]  for r in results]

    fig, ax = plt.subplots(figsize=(11, max(4, len(results) * 0.7 + 2)))

    bar_colors = [
        "#27AE60" if 40 <= v <= 60 else
        "#E74C3C" if v > 70 else
        "#F39C12"
        for v in gc_vals
    ]

    bars = ax.barh(names, gc_vals, color=bar_colors, edgecolor="white",
                   linewidth=1.0, height=0.55, zorder=3)

    # Typical GC range annotation
    ax.axvline(40, color="#7F8C8D", linestyle="--", linewidth=1, alpha=0.7, label="Typical range (40–60%)")
    ax.axvline(60, color="#7F8C8D", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvspan(40, 60, alpha=0.07, color="#2ECC71", zorder=0)

    ax.set_xlabel("GC Content (%)", fontsize=12)
    ax.set_title("GC Content Comparison Across Sequences", fontsize=13, fontweight="bold")
    ax.set_xlim(0, 105)
    ax.xaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.invert_yaxis()

    for bar, val in zip(bars, gc_vals):
        ax.text(
            bar.get_width() + 1.0, bar.get_y() + bar.get_height() / 2,
            f"{val}%", va="center", fontsize=9, fontweight="bold"
        )

    legend_patch = mpatches.Patch(color="#2ECC71", alpha=0.3, label="Typical range (40–60%)")
    ax.legend(handles=[legend_patch], fontsize=9, loc="lower right")

    plt.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filepath  = Path(output_dir) / f"gc_comparison_{timestamp}.png"
    plt.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return str(filepath.resolve())

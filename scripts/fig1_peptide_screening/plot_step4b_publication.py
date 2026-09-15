#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_step4b_publication.py

Purpose
-------
Generate a publication-style figure showing the physicochemical properties
of final deduplicated candidate peptides.

Input
-----
CSV file, typically:
    pipeline_results/step4b_final_dedup.csv

Output
------
Three files will be saved:
    Fig3_1_physicochemical_properties.png
    Fig3_1_physicochemical_properties.pdf
    Fig3_1_physicochemical_properties.svg

Recommended usage
-----------------
python plot_step4b_publication.py \
    -i pipeline_results/step4b_final_dedup.csv \
    -o pipeline_results
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm


def choose_font():
    """
    Choose a clean sans-serif font available on the system.
    Falls back to DejaVu Sans.
    """
    candidates = [
        "Arial",
        "Helvetica",
        "Liberation Sans",
        "Nimbus Sans",
        "DejaVu Sans"
    ]
    available = {f.name for f in fm.fontManager.ttflist}
    for font_name in candidates:
        if font_name in available:
            return font_name
    return "DejaVu Sans"


FONT_FAMILY = choose_font()

plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "regular",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

# A clean, publication-friendly palette
BAR_COLOR = "#5B8DB8"
MEAN_COLOR = "#E07A5F"
MEDIAN_COLOR = "#3D405B"
TEXT_COLOR = "#222222"


def clean_axis(ax):
    """
    Remove unnecessary decorations for a clean publication style.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.tick_params(axis="both", which="both", direction="out")
    ax.grid(False)
    ax.set_axisbelow(True)


def add_panel_label(ax, label):
    """
    Add panel label such as A, B, C...
    """
    ax.text(
        -0.12, 1.05, label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="right",
        va="top",
        color=TEXT_COLOR
    )


def add_stats_text(ax, mean_value, median_value, loc="upper right"):
    """
    Add mean and median text box inside the panel.
    """
    if loc == "upper right":
        x, y, ha = 0.97, 0.95, "right"
    elif loc == "upper left":
        x, y, ha = 0.03, 0.95, "left"
    else:
        x, y, ha = 0.97, 0.95, "right"

    stat_text = f"Mean = {mean_value:.2f}\nMedian = {median_value:.2f}"
    ax.text(
        x, y, stat_text,
        transform=ax.transAxes,
        ha=ha,
        va="top",
        fontsize=8.5,
        color=TEXT_COLOR,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="none",
            alpha=0.92
        )
    )


def plot_length_panel(ax, data):
    """
    Panel A: discrete bar plot for peptide length.
    """
    values = pd.Series(data).dropna().astype(int)
    counts = values.value_counts().sort_index()

    bars = ax.bar(
        counts.index.astype(str),
        counts.values,
        color=BAR_COLOR,
        edgecolor="white",
        linewidth=0.8,
        width=0.8
    )

    for bar, y in zip(bars, counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y + 0.3,
            str(int(y)),
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=TEXT_COLOR
        )

    ax.set_title("Length", pad=8)
    ax.set_xlabel("Peptide length, aa")
    ax.set_ylabel("Count")
    clean_axis(ax)
    ax.set_ylim(0, max(counts.values) * 1.18)


def plot_hist_panel(ax, data, title, xlabel, bins=10, stat_loc="upper right"):
    """
    Generic histogram panel with mean and median lines.
    """
    values = pd.Series(data).dropna().astype(float)

    ax.hist(
        values,
        bins=bins,
        color=BAR_COLOR,
        edgecolor="white",
        linewidth=0.8
    )

    mean_value = values.mean()
    median_value = values.median()

    ax.axvline(mean_value, color=MEAN_COLOR, linewidth=1.6)
    ax.axvline(median_value, color=MEDIAN_COLOR, linewidth=1.5, linestyle="--")

    ax.set_title(title, pad=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    clean_axis(ax)
    add_stats_text(ax, mean_value, median_value, loc=stat_loc)


def make_figure(df, outdir):
    """
    Create the 6-panel publication-style figure.
    """
    required_cols = [
        "length",
        "molecular_weight",
        "net_charge_ph7",
        "gravy",
        "instability_index",
        "source_protein_count"
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing)
        )

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.0))
    axes = axes.flatten()

    fig.suptitle(
        "Physicochemical properties of final deduplicated candidate peptides",
        fontsize=13,
        fontweight="bold",
        y=0.97
    )

    # A
    plot_length_panel(axes[0], df["length"])
    add_panel_label(axes[0], "")

    # B
    plot_hist_panel(
        axes[1],
        df["molecular_weight"],
        title="Molecular weight",
        xlabel="Molecular weight, Da",
        bins=10,
        stat_loc="upper right"
    )
    add_panel_label(axes[1], "")

    # C
    plot_hist_panel(
        axes[2],
        df["net_charge_ph7"],
        title="Net charge",
        xlabel="Net charge at pH 7.0",
        bins=8,
        stat_loc="upper right"
    )
    add_panel_label(axes[2], "")

    # D
    plot_hist_panel(
        axes[3],
        df["gravy"],
        title="Hydrophobicity",
        xlabel="GRAVY index",
        bins=10,
        stat_loc="upper right"
    )
    add_panel_label(axes[3], "")

    # E
    plot_hist_panel(
        axes[4],
        df["instability_index"],
        title="Predicted stability",
        xlabel="Instability index",
        bins=10,
        stat_loc="upper left"
    )
    add_panel_label(axes[4], "")

    # F
    plot_hist_panel(
        axes[5],
        df["source_protein_count"],
        title="Source protein coverage",
        xlabel="Source protein count",
        bins=10,
        stat_loc="upper right"
    )
    add_panel_label(axes[5], "")

    fig.subplots_adjust(
        left=0.07,
        right=0.98,
        bottom=0.09,
        top=0.88,
        wspace=0.28,
        hspace=0.40
    )

    base = os.path.join(outdir, "Fig3_1_physicochemical_properties")
    png_path = base + ".png"
    pdf_path = base + ".pdf"
    svg_path = base + ".svg"

    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    return png_path, pdf_path, svg_path


def summarize(df):
    """
    Print a short summary to terminal.
    """
    print("\nSummary of final deduplicated candidate peptides")
    print("-" * 55)
    print(f"Number of peptides      : {len(df)}")
    print(f"Length range            : {df['length'].min()}–{df['length'].max()} aa")
    print(
        f"Molecular weight range  : "
        f"{df['molecular_weight'].min():.2f}–{df['molecular_weight'].max():.2f} Da"
    )
    print(
        f"Net charge range        : "
        f"{df['net_charge_ph7'].min():.3f}–{df['net_charge_ph7'].max():.3f}"
    )
    print(
        f"GRAVY range             : "
        f"{df['gravy'].min():.3f}–{df['gravy'].max():.3f}"
    )
    print(
        f"Instability index range : "
        f"{df['instability_index'].min():.2f}–{df['instability_index'].max():.2f}"
    )
    print(
        f"Source protein count    : "
        f"{df['source_protein_count'].min()}–{df['source_protein_count'].max()}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot publication-style physicochemical property figure"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file, e.g. pipeline_results/step4b_final_dedup.csv"
    )
    parser.add_argument(
        "-o", "--outdir",
        default="pipeline_results",
        help="Output directory"
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("\nReading input file...")
    print(f"Input: {args.input}")

    df = pd.read_csv(args.input)
    if df.empty:
        raise ValueError("Input CSV is empty.")

    summarize(df)

    print("\nGenerating figure...")
    png_path, pdf_path, svg_path = make_figure(df, args.outdir)

    print("\nDone.")
    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")
    print(f"SVG: {svg_path}")


if __name__ == "__main__":
    main()

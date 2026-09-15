#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4_select.py
综合评分+ 多样性感知选取，最终输出 40–60 条候选肽。
输出：step4_final_candidates.csv, step4_final_candidates.fasta,分布图
用法：python step4_select.py -i pipeline_results/step3_filtered.csv
"""

import argparse
import os
import sys
import time
from collections import Counter

import pandas as pd
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams["axes.unicode_minus"] = False
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("提示：未找到 matplotlib/seaborn，跳过绘图。")

HYDROPHOBIC = frozenset("VILFMWA")
AROMATIC    = frozenset("FWY")
KMER_K      = 2


def kmer_set(seq, k=KMER_K):
    return frozenset(seq[i:i+k] for i in range(len(seq) - k + 1))


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def bio_score(row):
    """
    生物活性潜力评分（满分 100）。
    基于文献中ACE 抑制肽、屏障保护肽的序列特征经验规则。
    参考：Li et al. (2018) Food Chem; Nongonierma & FitzGerald (2017) J Funct Foods
    """
    score = 50.0
    seq = row["sequence"]
    n= len(seq)
    c   = Counter(seq)

    # 芳香族残基（有利于疏水口袋结合）
    aro = sum(c.get(aa, 0) for aa in AROMATIC)
    if aro >= 1:
        score += 8
    if aro >= 2:
        score += 5

    # 适中疏水性（GRAVY -1.5~0.0最佳）
    gravy = row.get("gravy", np.nan)
    if pd.notna(gravy):
        if -1.5 <= gravy <= 0.0:
            score += 10
        elif gravy < -2.0 or gravy > 0.3:
            score -= 8

    # 适中净电荷（轻微正电荷有利于与靶蛋白结合）
    charge = row.get("net_charge_ph7", np.nan)
    if pd.notna(charge):
        if 0.0 <= charge <= 2.0:
            score += 8
        elif charge > 3.0 or charge < -1.0:
            score -= 8

    # 含Lys 或 Arg（有利于氢键和盐桥）
    if c.get("K", 0) + c.get("R", 0) >= 1:
        score += 5

    # 来源蛋白数量多→跨蛋白保守→ 可能具有结构/功能意义
    src = row.get("source_protein_count", 1)
    if pd.notna(src):
        score += min(10, int(src) * 0.5)

    # 长度适中（6–10 aa 是研究最多的活性短肽范围）
    if 6 <= n <= 10:
        score += 8
    elif n < 5 or n > 12:
        score -= 5

    # 合成评分加权贡献
    synth = row.get("synth_score", 70.0)
    if pd.notna(synth):
        score += (synth - 70.0) * 0.1

    # 不稳定指数（越低越好）
    instab = row.get("instability_index", 30.0)
    if pd.notna(instab):
        if instab < 25:
            score += 5
        elif instab > 35:
            score -= 3

    return round(max(0.0, min(100.0, score)), 2)


def composite_score(row):
    """
    综合评分 = 生物活性潜力(60%) + 合成可行性(25%) + 来源保守性(15%)
    """
    b = row.get("bio_score", 50.0)
    s = row.get("synth_score", 70.0)

    src = row.get("source_protein_count", 1)
    src_score = min(100.0, float(src) * 5.0) if pd.notna(src) else 50.0

    return round(0.60 * b + 0.25 * s + 0.15 * src_score, 2)


def diversity_select(df, target_n, sim_threshold=0.60):
    """
    贪心多样性选取：
    按综合评分从高到低排序，逐条加入候选集，
    若与已选肽段的最大 Jaccard 相似度超过阈值则跳过。
    """
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    sequences= df["sequence"].tolist()
    kmer_sets_list = [kmer_set(seq) for seq in sequences]

    selected_pos= []
    selected_kmers = []

    for i in range(len(df)):
        if len(selected_pos) >= target_n:
            break
        ks = kmer_sets_list[i]
        max_sim = max(jaccard(ks, sk) for sk in selected_kmers) if selected_kmers else 0.0
        if max_sim <= sim_threshold:
            selected_pos.append(i)
            selected_kmers.append(ks)

    # 数量不足时放宽阈值补充
    if len(selected_pos) < target_n:
        relaxed = sim_threshold + 0.15
        print(f"  多样性选取数量不足，放宽相似度阈值至 {relaxed:.2f} 补充...")
        for i in range(len(df)):
            if len(selected_pos) >= target_n:
                break
            if i in selected_pos:
                continue
            ks = kmer_sets_list[i]
            max_sim = max(jaccard(ks, sk) for sk in selected_kmers) if selected_kmers else 0.0
            if max_sim <= relaxed:
                selected_pos.append(i)
                selected_kmers.append(ks)

    return df.iloc[selected_pos].reset_index(drop=True)


def write_fasta(df, filepath):
    with open(filepath, "w") as f:
        for i, row in df.iterrows():
            rank= i + 1
            seq    = row["sequence"]
            mw     = row.get("molecular_weight", "")
            charge = row.get("net_charge_ph7", "")
            gravy  = row.get("gravy", "")
            cscore = row.get("composite_score", "")
            header = (f">candidate_{rank:03d} | seq={seq} | "f"MW={mw} | charge={charge} | "
                      f"GRAVY={gravy} | score={cscore}")
            f.write(header + "\n")
            f.write(seq + "\n")


def plot_distributions(df, outdir):
    """理化性质分布图（英文标签）。"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Final Candidate Peptides — Physicochemical Property Distributions",
                 fontsize=13, fontweight="bold")

    plot_cfg = [
        ("length","Peptide Length (aa)",   "steelblue"),
        ("molecular_weight", "Molecular Weight (Da)",  "darkorange"),
        ("net_charge_ph7",   "Net Charge (pH 7.0)",    "seagreen"),
        ("gravy",            "GRAVY Index",             "mediumpurple"),
    ]

    for ax, (col, xlabel, color) in zip(axes.flatten(), plot_cfg):
        if col not in df.columns:
            ax.set_visible(False)
            continue
        data = df[col].dropna()
        ax.hist(data, bins=min(20, len(data)), color=color,
                edgecolor="white", linewidth=0.6)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"Distribution of {xlabel}", fontsize=10)
        mean_val = data.mean()
        ax.axvline(mean_val, color="crimson", linestyle="--",
                   linewidth=1.2, label=f"Mean={mean_val:.2f}")
        ax.legend(fontsize=8)
        sns.despine(ax=ax)

    plt.tight_layout()
    out_path = os.path.join(outdir, "step4_distributions.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  分布图: {out_path}")


def plot_funnel(funnel_csv, outdir):
    """筛选漏斗图（英文标签，避免中文字体问题）。"""
    if not os.path.exists(funnel_csv):
        return

    label_map = {
        "读取输入（Stage2输出）":"S0: Input (Stage2 output)",
        "长度 5–12 aa":              "S1: Length 5-12 aa",
        "分子量 450.0–1500.0 Da":    "S2: MW 450-1500 Da",
        "净电荷 -1.0–4.0":           "S3: Net charge -1.0~4.0",
        "GRAVY -2.5–0.3":            "S4: GRAVY -2.5~0.3",
        "不稳定指数 < 40.0":         "S5: Instability < 40",
        "不含 Cys（合成复杂）":      "S6: No Cys",
        "Pro 比例 ≤ 0.3":            "S7: Pro ratio <= 30%",
        "连续相同残基 ≤ 3":          "S8: Max same-aa run <= 3",
        "连续强疏水残基 ≤ 4":        "S9: Max hydrophobic run <= 4",
        "无毒性规则模式":            "S10: No toxic pattern",
        "无过敏原规则模式":          "S11: No allergen pattern",
        "合成可行性评分 ≥ 60.0":     "S12: Synth score >= 60",
    }

    df = pd.read_csv(funnel_csv)
    labels = [label_map.get(str(row["筛选条件"]), str(row["筛选条件"]))
              for _, row in df.iterrows()]
    counts = df["保留数量"].tolist()

    fig, ax = plt.subplots(figsize=(11, 7))
    bars = ax.barh(labels[::-1], counts[::-1], color="steelblue", edgecolor="white")
    ax.set_xlabel("Number of Peptides Retained", fontsize=11)
    ax.set_title("Peptide Screening Funnel", fontsize=13, fontweight="bold")
    ax.bar_label(bars, labels=[f"{v:,}" for v in counts[::-1]],
                 padding=4, fontsize=9)
    ax.set_xscale("log")
    sns.despine(ax=ax)
    plt.tight_layout()

    out_path = os.path.join(outdir, "step4_funnel.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  漏斗图: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Step4:综合评分 & 多样性选取")
    parser.add_argument("-i", "--input",required=True,
                        help="step3_filtered.csv")
    parser.add_argument("-o", "--outdir", default="pipeline_results")
    parser.add_argument("--target", type=int, default=50,
                        help="最终候选肽数量（默认 50，建议 40–60）")
    parser.add_argument("--sim-threshold", type=float, default=0.60,
                        help="k-mer Jaccard 相似度去冗余阈值（默认 0.60）")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv= os.path.join(args.outdir, "step4_final_candidates.csv")
    out_fasta  = os.path.join(args.outdir, "step4_final_candidates.fasta")
    funnel_csv = os.path.join(args.outdir, "step3_funnel.csv")
    t0= time.time()

    print("\n" + "="*60)
    print("Step 4: 综合评分 & 多样性感知选取")
    print("="*60)
    print(f"  输入: {args.input}")
    print(f"  目标候选肽数量: {args.target}")
    print(f"  相似度去冗余阈值: {args.sim_threshold}")

    df = pd.read_csv(args.input)
    print(f"  读取肽段数: {len(df):,}")

    if len(df) == 0:
        print("错误：输入文件为空，请检查 step3 是否正常运行。")
        sys.exit(1)

    if "synth_score" not in df.columns:
        print("错误：synth_score 列不存在，请确认 step3 已正常完成。")
        sys.exit(1)

    # 计算生物活性潜力评分
    print("\n[1/3] 计算生物活性潜力评分...")
    df["bio_score"] = df.apply(bio_score, axis=1)

    # 计算综合评分
    print("[2/3] 计算综合评分...")
    df["composite_score"] = df.apply(composite_score, axis=1)

    # 多样性感知选取
    print(f"[3/3] 多样性感知选取（目标 {args.target} 条）...")
    final_df = diversity_select(df, args.target, args.sim_threshold)

    # 添加排名列
    final_df.insert(0, "rank", range(1, len(final_df) + 1))

    # 整理输出列顺序
    priority_cols = [
        "rank", "sequence", "length",
        "molecular_weight", "isoelectric_point", "net_charge_ph7",
        "gravy", "instability_index", "aromaticity",
        "synth_score", "bio_score", "composite_score",
        "source_protein_count", "first_source_protein",]
    existing = [c for c in priority_cols if c in final_df.columns]
    other= [c for c in final_df.columns if c not in existing]
    final_df = final_df[existing + other]

    # 保存
    final_df.to_csv(out_csv, index=False)
    write_fasta(final_df, out_fasta)

    #绘图
    if HAS_PLOT:
        print("\n绘制图表...")
        plot_distributions(final_df, args.outdir)
        plot_funnel(funnel_csv, args.outdir)
    else:
        print("\n跳过绘图（请安装：pip install matplotlib seaborn）")

    #汇总统计
    elapsed = (time.time() - t0) / 60
    print("\n" + "="*60)
    print("最终候选肽统计摘要")
    print("="*60)
    print(f"  最终候选肽数量    : {len(final_df)}")
    print(f"  长度范围          : {final_df['length'].min()}–{final_df['length'].max()} aa")

    if "molecular_weight" in final_df.columns:
        print(f"  分子量范围        : {final_df['molecular_weight'].min():.1f}–"
              f"{final_df['molecular_weight'].max():.1f} Da")
    if "net_charge_ph7" in final_df.columns:
        print(f"  净电荷范围        : {final_df['net_charge_ph7'].min():.1f}–"
              f"{final_df['net_charge_ph7'].max():.1f}")
    if "gravy" in final_df.columns:
        print(f"  GRAVY范围        : {final_df['gravy'].min():.3f}–"
              f"{final_df['gravy'].max():.3f}")
    if "composite_score" in final_df.columns:
        print(f"  综合评分范围      : {final_df['composite_score'].min():.1f}–"
              f"{final_df['composite_score'].max():.1f}")

    print(f"\n  输出 CSV  : {out_csv}")
    print(f"  输出 FASTA: {out_fasta}")
    print(f"  运行时间  : {elapsed:.2f} 分钟")
    print("\n✓ 完成！可将step4_final_candidates.fasta 提交用于分子对接（3.3节）。")


if __name__ == "__main__":
    main()

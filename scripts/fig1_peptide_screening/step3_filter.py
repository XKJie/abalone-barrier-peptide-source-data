#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_filter.py
多阶段严格筛选，记录每步漏斗数量，目标输出 200–500 条。
输出：step3_filtered.csv, step3_funnel.csv
用法：python step3_filter.py -i pipeline_results/step2_with_properties.csv
"""

import argparse
import os
import time
from collections import Counter

import pandas as pd

#──筛选阈值（论文方法章节需如实报告）──────────────────────────
THRESHOLDS = {
    "len_min":5,      # aa，合成成本低、可覆盖结合口袋
    "len_max":         12,     # aa，≤12 aa SPPS成本可控
    "mw_min":          450.0,  # Da，排除极端短肽
    "mw_max":          1500.0, # Da，合成与吸收权衡
    "charge_min":     -1.0,    # pH7.0净电荷
    "charge_max":      4.0,    # 过高正电荷→细胞毒风险
    "gravy_min":      -2.5,    # 过亲水削弱靶蛋白结合
    "gravy_max":       0.3,    # 过疏水→溶解性和毒性风险
    "instab_max":     40.0,    # strict 预设阈值
    "max_pro_ratio":   0.30,   # Pro 比例上限
    "max_same_run":    3,      # 连续相同残基最大数
    "max_hydro_run":   4,      # 连续强疏水残基最大长度
    "synth_min":      60.0,    # 合成可行性评分下限
}

HYDROPHOBIC = frozenset("VILFMWA")
AROMATIC    = frozenset("FWY")

# 毒性规则模式（粗筛，不替代ToxinPred）
TOXIC_PATTERNS = [
    "KKKK","RRRR","KRKR","RRKK",# 极端阳离子
    "LLLL","FFFF","WWWW","IIII",    # 极端疏水
    "LKLL","KLLL","FKFF","WKWK",   # 溶膜 AMP 核心
]

# 过敏原模式（粗筛，不替代 AllerTOP）
ALLERGEN_PATTERNS = [
    "QQPFP","LQPFP","QQQPP",# 麸质/醇溶蛋白
    "RPQQPY","YLQQQ",                # 花生
]


def max_same_run(seq):
    if not seq:
        return 0
    best = cur = 1
    for i in range(1, len(seq)):
        cur = cur + 1 if seq[i] == seq[i-1] else 1
        best = max(best, cur)
    return best


def max_hydro_run(seq):
    best = cur = 0
    for aa in seq:
        if aa in HYDROPHOBIC:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def synth_score(seq):
    """合成可行性评分，满分 100，扣分制。"""
    score = 100.0
    n = len(seq)
    c = Counter(seq)

    score -= 15* c.get("C", 0)           # Cys：二硫键复杂
    score -= 5* c.get("M", 0)           # Met：氧化风险
    for i in range(n-1):
        if seq[i] == "G" and seq[i+1] == "G":
            score -= 5     # GG：聚集倾向
        if seq[i] == "P" and seq[i+1] == "P":
            score -= 10                    # PP：位阻if seq[i] == "D" and seq[i+1] == "P":
            score -= 8                     # DP：酸不稳定
    pro_ratio = c.get("P", 0) / n
    if pro_ratio > 0.20:
        score -= (pro_ratio - 0.20) * 80
    hrun = max_hydro_run(seq)
    if hrun >= 4:
        score -= (hrun - 3) * 10
    if seq[0] in ("Q", "N"):
        score -= 8# N端脱酰胺/环化
    if n >= 11:
        score -= (n - 10) * 3
    return max(0.0, min(100.0, score))


def log_stage(funnel, stage, desc, n):
    funnel.append({"阶段": stage, "筛选条件": desc, "保留数量": n})
    print(f"  [{stage}] {desc:<40s} 保留: {n:>8,}")


def main():
    parser = argparse.ArgumentParser(description="Step3: 多阶段精细筛选")
    parser.add_argument("-i", "--input",  required=True,
                        help="step2_with_properties.csv")
    parser.add_argument("-o", "--outdir", default="pipeline_results")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv    = os.path.join(args.outdir, "step3_filtered.csv")
    funnel_csv = os.path.join(args.outdir, "step3_funnel.csv")
    t0 = time.time()

    print("\n" + "="*60)
    print("Step 3: 多阶段精细筛选")
    print("="*60)

    df = pd.read_csv(args.input)
    funnel = []
    log_stage(funnel, "S0", "读取输入（Stage2输出）", len(df))

    T = THRESHOLDS

    # S1: 长度
    df = df[df["length"].between(T["len_min"], T["len_max"])]
    log_stage(funnel, "S1", f"长度 {T['len_min']}–{T['len_max']} aa", len(df))

    # S2: 分子量
    df = df[df["molecular_weight"].between(T["mw_min"], T["mw_max"])]
    log_stage(funnel, "S2", f"分子量 {T['mw_min']}–{T['mw_max']} Da", len(df))

    # S3: 净电荷
    df = df[df["net_charge_ph7"].between(T["charge_min"], T["charge_max"])]
    log_stage(funnel, "S3", f"净电荷 {T['charge_min']}–{T['charge_max']}", len(df))

    # S4: GRAVY
    df = df[df["gravy"].between(T["gravy_min"], T["gravy_max"])]
    log_stage(funnel, "S4", f"GRAVY {T['gravy_min']}–{T['gravy_max']}", len(df))

    # S5: 不稳定指数
    df = df[df["instability_index"] < T["instab_max"]]
    log_stage(funnel, "S5", f"不稳定指数 < {T['instab_max']}", len(df))

    # S6: 去除含Cys
    df = df[~df["sequence"].str.contains("C")]
    log_stage(funnel, "S6", "不含 Cys（合成复杂）", len(df))

    # S7: Pro 比例
    df = df[df["sequence"].apply(
        lambda s: s.count("P") / len(s) <= T["max_pro_ratio"])]
    log_stage(funnel, "S7", f"Pro 比例 ≤ {T['max_pro_ratio']}", len(df))

    # S8: 连续相同残基
    df = df[df["sequence"].apply(max_same_run) <= T["max_same_run"]]
    log_stage(funnel, "S8", f"连续相同残基 ≤ {T['max_same_run']}", len(df))

    # S9: 连续疏水片段
    df = df[df["sequence"].apply(max_hydro_run) <= T["max_hydro_run"]]
    log_stage(funnel, "S9", f"连续强疏水残基 ≤ {T['max_hydro_run']}", len(df))

    # S10: 毒性模式
    def no_toxic(seq):
        return not any(pat in seq for pat in TOXIC_PATTERNS)
    df = df[df["sequence"].apply(no_toxic)]
    log_stage(funnel, "S10", "无毒性规则模式", len(df))

    # S11: 过敏原模式
    def no_allergen(seq):
        return not any(pat in seq for pat in ALLERGEN_PATTERNS)
    df = df[df["sequence"].apply(no_allergen)]
    log_stage(funnel, "S11", "无过敏原规则模式", len(df))

    # S12: 合成可行性评分
    print("计算合成可行性评分...")
    df["synth_score"] = df["sequence"].apply(synth_score)
    df = df[df["synth_score"] >= T["synth_min"]]
    log_stage(funnel, "S12", f"合成可行性评分 ≥ {T['synth_min']}", len(df))

    # 保存
    df.to_csv(out_csv, index=False)
    pd.DataFrame(funnel).to_csv(funnel_csv, index=False)

    elapsed = (time.time() - t0) / 60
    print(f"\n✓ 完成。")
    print(f"  过滤结果: {out_csv} （{len(df):,} 条）")
    print(f"  漏斗统计: {funnel_csv}")
    print(f"  运行时间: {elapsed:.2f} 分钟")

    if len(df) < 40:
        print("\n⚠ 保留数量不足 40 条，建议适当放宽阈值后重新运行。")
    elif len(df) > 2000:
        print(f"\n提示：保留 {len(df)} 条，step4 会进一步缩小至 40–60 条。")


if __name__ == "__main__":
    main()

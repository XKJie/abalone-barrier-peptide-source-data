#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step4b_dedup.py
对step4 输出的 50 条候选肽进行家族去冗余：
- 检测所有 Jaccard 相似度 > 阈值的序列对
- 每个相似家族只保留 source_protein_count 最高的一条
- 输出去冗余后的最终候选肽表

用法：python step4b_dedup.py -i pipeline_results/step4_final_candidates.csv
"""

import argparse
import os
from collections import defaultdict

import pandas as pd

KMER_K = 2
# 相似度阈值：高于此值认为属于同一家族
# 0.50 = 50% k-mer 重叠，可根据实际情况调整
SIM_THRESHOLD = 0.50


def kmer_set(seq, k=KMER_K):
    return frozenset(seq[i:i+k] for i in range(len(seq) - k + 1))


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def find_families(df, threshold):
    """
    用 Union-Find 将相似度超过阈值的肽段归入同一家族。
    返回：list of list，每个子列表是一个家族的行索引。
    """
    n = len(df)
    sequences = df["sequence"].tolist()
    kmer_sets = [kmer_set(seq) for seq in sequences]

    # Union-Find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # 两两比较（50条，只有50*49/2 = 1225 次，很快）
    for i in range(n):
        for j in range(i + 1, n):
            sim = jaccard(kmer_sets[i], kmer_sets[j])
            if sim > threshold:
                union(i, j)

    # 按根节点分组
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    return list(groups.values())


def main():
    parser = argparse.ArgumentParser(description="Step4b: 候选肽家族去冗余")
    parser.add_argument("-i", "--input",required=True,
                        help="step4_final_candidates.csv")
    parser.add_argument("-o", "--outdir", default="pipeline_results")
    parser.add_argument("--threshold", type=float, default=SIM_THRESHOLD,
                        help=f"Jaccard 相似度阈值（默认 {SIM_THRESHOLD}）")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv= os.path.join(args.outdir, "step4b_final_dedup.csv")
    out_fasta  = os.path.join(args.outdir, "step4b_final_dedup.fasta")
    report_csv = os.path.join(args.outdir, "step4b_family_report.csv")

    print("\n" + "="*60)
    print("Step4b: 候选肽家族去冗余")
    print("="*60)

    df = pd.read_csv(args.input)
    print(f"  输入候选肽数量  : {len(df)}")
    print(f"  k-mer 大小      : {KMER_K}")
    print(f"  相似度阈值      : {args.threshold}")

    #── 检测家族 ──────────────────────────────────────────────
    families = find_families(df, args.threshold)

    singleton_count = sum(1 for f in families if len(f) == 1)
    multi_count     = sum(1 for f in families if len(f) > 1)
    print(f"\n  检测到家族数量  : {len(families)}")
    print(f"  单成员家族      : {singleton_count}（直接保留）")
    print(f"  多成员家族      : {multi_count}（每族保留最佳一条）")

    # ── 打印多成员家族详情 ────────────────────────────────────
    report_rows = []
    kept_indices = []

    print("\n  多成员家族明细：")
    print(f"  {'家族':>4}  {'序列':<14} {'长度':>4} {'来源数':>6} "
          f"{'GRAVY':>7} {'净电荷':>6} {'评分':>7} {'状态'}")
    print("  " + "-"*70)

    family_id = 0
    for group in sorted(families, key=lambda g: -len(g)):
        members = df.iloc[group].copy()

        if len(members) == 1:
            idx = group[0]
            kept_indices.append(idx)
            report_rows.append({
                "family_id": "—",
                "sequence":  members.iloc[0]["sequence"],
                "source_protein_count": members.iloc[0].get("source_protein_count", ""),
                "status":"singleton",
            })
            continue

        family_id += 1

        # 家族内按 source_protein_count 降序，选最高的保留
        members = members.sort_values("source_protein_count", ascending=False)
        best_seq_idx = members.index[0]
        kept_indices.append(best_seq_idx)

        for pos, (idx, row) in enumerate(members.iterrows()):
            status = "KEEP✓" if idx == best_seq_idx else "drop"
            print(f"  {family_id:>4}  {row['sequence']:<14} "
                  f"{row['length']:>4} "
                  f"{int(row.get('source_protein_count', 0)):>6} "
                  f"{row.get('gravy', 0):>7.3f} "
                  f"{row.get('net_charge_ph7', 0):>6.2f} "
                  f"{row.get('composite_score', 0):>7.1f}  "
                  f"{status}")
            report_rows.append({
                "family_id":family_id,
                "sequence":            row["sequence"],
                "source_protein_count": row.get("source_protein_count", ""),
                "gravy":               row.get("gravy", ""),
                "net_charge_ph7":      row.get("net_charge_ph7", ""),
                "composite_score":     row.get("composite_score", ""),
                "status":              status,
            })

    # ── 构建去冗余结果 ────────────────────────────────────────
    final_df = df.loc[sorted(set(kept_indices))].copy()
    final_df = final_df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    if "rank" in final_df.columns:
        final_df = final_df.drop(columns=["rank"])
    final_df.insert(0, "rank", range(1, len(final_df) + 1))

    # ── 保存 CSV ──────────────────────────────────────────────
    final_df.to_csv(out_csv, index=False)

    # ── 保存 FASTA ────────────────────────────────────────────
    with open(out_fasta, "w") as f:
        for i, row in final_df.iterrows():
            rank= row["rank"]
            seq    = row["sequence"]
            mw     = row.get("molecular_weight", "")
            charge = row.get("net_charge_ph7", "")
            gravy  = row.get("gravy", "")
            src    = row.get("source_protein_count", "")
            header = (f">candidate_{rank:03d} | seq={seq} | "
                      f"MW={mw} | charge={charge} | "
                      f"GRAVY={gravy} | source_n={src}")
            f.write(header + "\n")
            f.write(seq + "\n")

    # ── 保存家族报告 ──────────────────────────────────────────
    pd.DataFrame(report_rows).to_csv(report_csv, index=False)

    # ── 汇总 ──────────────────────────────────────────────────
    print("\n" + "="*60)
    print("去冗余结果汇总")
    print("="*60)
    print(f"  输入候选肽      : {len(df)} 条")
    print(f"  去冗余后保留    : {len(final_df)} 条")
    print(f"  筛除冗余肽      : {len(df) - len(final_df)} 条")
    print(f"\n  输出 CSV        : {out_csv}")
    print(f"  输出 FASTA      : {out_fasta}")
    print(f"  家族报告        : {report_csv}")
    print("\n✓ 完成！可将 step4b_final_dedup.fasta 提交用于分子对接（3.3节）。")


if __name__ == "__main__":
    main()

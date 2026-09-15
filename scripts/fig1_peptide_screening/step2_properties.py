#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2_properties.py
读取 step1_unique_peptides.csv，计算全套理化性质。
输出：step2_with_properties.csv
用法：python step2_properties.py -i pipeline_results/step1_unique_peptides.csv
"""

import argparse
import os
import sys
import time
from math import log10

import pandas as pd
from tqdm import tqdm

# 检查BioPython
try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    print("请安装 biopython: pip install biopython")
    sys.exit(1)

# pKa 值 (Henderson-Hasselbalch 计算净电荷用)
PKA = {
    "N_term": 8.0, "C_term": 3.1,
    "K": 10.5, "R": 12.5, "H": 6.0,
    "D": 3.9,  "E": 4.1,
    "C": 8.3,  "Y": 10.1,
}


def net_charge(seq, ph=7.0):
    """Henderson-Hasselbalch 方程计算净电荷。"""
    q = 0.0
    # N端正电
    q += 1.0 / (1.0 + 10** (ph - PKA["N_term"]))
    # C端负电
    q -= 1.0 / (1.0 + 10 ** (PKA["C_term"] - ph))
    #侧链
    pos_res = {"K": 0, "R": 0, "H": 0}
    neg_res = {"D": 0, "E": 0, "C": 0, "Y": 0}
    for aa in seq:
        if aa in pos_res:
            pos_res[aa] += 1
        if aa in neg_res:
            neg_res[aa] += 1
    for aa, cnt in pos_res.items():
        q += cnt / (1.0 + 10 ** (ph - PKA[aa]))
    for aa, cnt in neg_res.items():
        q -= cnt / (1.0 + 10 ** (PKA[aa] - ph))
    return round(q, 3)


def calc_properties(seq):
    """计算单条肽段的理化性质，返回 dict。"""
    try:
        pa = ProteinAnalysis(seq)
        mw         = round(pa.molecular_weight(), 2)
        pi         = round(pa.isoelectric_point(), 2)
        gravy      = round(pa.gravy(), 3)
        instab     = round(pa.instability_index(), 2)
        aromatic   = round(pa.aromaticity(), 3)
    except Exception:
        mw = pi = gravy = instab = aromatic = None

    charge = net_charge(seq,7.0)

    return {
        "molecular_weight":  mw,
        "isoelectric_point": pi,
        "net_charge_ph7":    charge,
        "gravy":             gravy,
        "instability_index": instab,
        "aromaticity":aromatic,
    }


def main():
    parser = argparse.ArgumentParser(description="Step2: 计算理化性质")
    parser.add_argument("-i", "--input",  required=True,
                        help="step1_unique_peptides.csv")
    parser.add_argument("-o", "--outdir", default="pipeline_results")
    parser.add_argument("--batch", type=int, default=10000,
                        help="每批处理数量（控制内存）")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "step2_with_properties.csv")
    t0 = time.time()

    print("\n" + "="*60)
    print("Step 2: 计算理化性质")
    print("="*60)
    print(f"  输入: {args.input}")

    df = pd.read_csv(args.input)
    print(f"  读取肽段数: {len(df):,}")

    tqdm.pandas(desc="  计算中")
    props = df["sequence"].progress_apply(calc_properties)
    props_df = pd.DataFrame(list(props))
    df = pd.concat([df, props_df], axis=1)

    df.to_csv(out_csv, index=False)
    elapsed = (time.time() - t0) / 60
    print(f"\n✓ 完成。输出: {out_csv}")
    print(f"  运行时间: {elapsed:.2f} 分钟")
    print(f"  包含理化性质的肽段数: {len(df):,}")


if __name__ == "__main__":
    main()

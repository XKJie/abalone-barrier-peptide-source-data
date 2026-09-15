#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step1_digest.py
从 protein.faa 进行序贯模拟酶切并去重。
输出：step1_unique_peptides.csv
用法：python step1_digest.py -i protein.faa
"""

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool, cpu_count

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# 酶切规则 (切割位点, 不在其前切割的残基)
# Pepsin pH2.0: 在 F/L/W/Y C端切割，不在P前切割
# Trypsin: 在 K/R C端切割，不在P前切割
# Chymotrypsin高特异性: 在 F/Y/W C端切割，不在P前切割
ENZYME_RULES = {
    "pepsin_ph2":(frozenset("FLWY"), frozenset("P")),
    "trypsin":           (frozenset("KR"),   frozenset("P")),
    "chymotrypsin_high": (frozenset("FYW"),  frozenset("P")),
}
ENZYME_ORDER = ["pepsin_ph2", "trypsin", "chymotrypsin_high"]


def single_digest(seq, sites, not_before, missed=0):
    n = len(seq)
    cut_pos = [0]
    for i in range(n - 1):
        if seq[i] in sites and seq[i + 1] not in not_before:
            cut_pos.append(i + 1)
    cut_pos.append(n)
    frags = [seq[cut_pos[i]:cut_pos[i+1]] for i in range(len(cut_pos)-1)
             if cut_pos[i+1] > cut_pos[i]]
    if missed == 0:
        return frags
    result = []
    nf = len(frags)
    for i in range(nf):
        for j in range(i, min(i + missed + 1, nf)):
            result.append("".join(frags[i:j+1]))
    return result


def sequential_digest(seq, missed=0):
    current = [seq]
    for enzyme in ENZYME_ORDER:
        sites, not_before = ENZYME_RULES[enzyme]
        nxt = []
        for pep in current:
            nxt.extend(single_digest(pep, sites, not_before, missed))
        current = nxt
    return current


def worker(args):
    rid, seq, missed, len_min, len_max = args
    results = []
    for pep in sequential_digest(seq, missed):
        if len_min <= len(pep) <= len_max:
            if all(aa in STANDARD_AA for aa in pep):
                results.append((pep, rid))
    return results


def parse_fasta(filepath):
    records = []
    ns_count = 0
    ns_records = 0
    cur_id = None
    cur_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seq = "".join(cur_seq).upper()
                    ns = sum(1 for aa in seq if aa not in STANDARD_AA)
                    if ns > 0:
                        ns_records += 1
                        ns_count += ns
                    records.append((cur_id, seq))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        seq = "".join(cur_seq).upper()
        ns = sum(1 for aa in seq if aa not in STANDARD_AA)
        if ns > 0:
            ns_records += 1
            ns_count += ns
        records.append((cur_id, seq))
    return records, ns_records, ns_count


def main():
    parser = argparse.ArgumentParser(description="Step1: 模拟酶切 & 去重")
    parser.add_argument("-i", "--input",required=True, help="protein.faa")
    parser.add_argument("-o", "--outdir",  default="pipeline_results")
    parser.add_argument("--missed",  type=int, default=0)
    parser.add_argument("--len-min", type=int, default=3)
    parser.add_argument("--len-max", type=int, default=20)
    parser.add_argument("--workers", type=int, default=min(4, cpu_count()))
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "step1_unique_peptides.csv")
    t0 = time.time()

    print("\n" + "="*60)
    print("Step 1: 模拟序贯酶切 &肽段去重")
    print("="*60)
    print(f"  输入文件: {args.input}")
    print(f"  酶切顺序 : Pepsin(pH2) → Trypsin → Chymotrypsin(高)")
    print(f"  漏切数   : {args.missed}")
    print(f"  长度范围 : {args.len_min}–{args.len_max} aa")
    print(f"  进程数   : {args.workers}")

    print("\n[1/3] 读取 FASTA ...")
    records, ns_recs, ns_chars = parse_fasta(args.input)
    print(f"  蛋白记录总数           : {len(records):>12,}")
    print(f"  含非标准字符的记录数   : {ns_recs:>12,}")
    print(f"  非标准字符总数         : {ns_chars:>12,}")

    print("\n[2/3] 并行酶切中...")
    tasks = [(rid, seq, args.missed, args.len_min, args.len_max)
             for rid, seq in records]

    # pep -> [count, first_source]
    pep_info = {}
    write_count = 0

    with Pool(args.workers) as pool:
        for batch in pool.imap_unordered(worker, tasks, chunksize=200):
            for pep, rid in batch:
                write_count += 1
                if pep not in pep_info:
                    pep_info[pep] = [1, rid]
                else:
                    pep_info[pep][0] += 1

    print(f"  酶切写入次数（含重复） : {write_count:>12,}")
    print(f"  去重后唯一肽段数       : {len(pep_info):>12,}")

    print("\n[3/3] 写入 CSV ...")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence", "length", "source_protein_count", "first_source_protein"])
        for pep, (cnt, src) in pep_info.items():
            writer.writerow([pep, len(pep), cnt, src])

    elapsed = (time.time() - t0) / 60
    print(f"\n✓ 完成。输出: {out_csv}")
    print(f"  运行时间: {elapsed:.2f} 分钟")
    print(f"  唯一肽段: {len(pep_info):,} 条")


if __name__ == "__main__":
    main()

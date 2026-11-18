#!/usr/bin/env python3
"""
00_smoke_test.py
================
REPRO-D 10-minute smoke test for the Neoantigen & Immunogenomics workspace.

Functions:
  1. Validates workspace directory structure (data/raw, data/processed, repro_results, logs).
  2. Deterministically constructs a lightweight synthetic benchmark dataset (seed=20261115):
     - 24 somatic tumor neoantigen candidates (SNV / In-frame Indel / Frameshift)
       with RNA VAF/TPM expression, WT vs. MT MHC-I IC50 (nM), Percentile Rank,
       and Differential Agretopicity Index (DAI = WT_IC50 / MT_IC50).
     - 6 patient / TIL samples with TCR CDR3beta clonotype counts and frequencies
       to verify clonality (1 - normalized Shannon entropy) and top-clone expansion.
  3. Applies standard clinical prioritization thresholds:
     - MT_IC50 <= 500 nM, Percentile Rank <= 2.0%, RNA_TPM >= 1.0, Tumor_RNA_VAF >= 0.05.
  4. Saves cleaned tables to data/processed/ and quantitative summary to repro_results/.
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys
from typing import Dict, List

import numpy as np


SEED = 20261115


def get_project_root() -> pathlib.Path:
    """定位工作区根目录 (reproduction/neoantigen_immunogenomics)，避免依赖当前终端 cwd。"""
    return pathlib.Path(__file__).resolve().parent.parent


def compute_shannon_clonality(counts: np.ndarray) -> float:
    """
    计算免疫组库标准化香农克隆性 (Normalized Shannon Clonality, 范围 0~1)：
      Clonality = 1 - H / ln(N)
    其中 H = -sum(p_i * ln(p_i)) 为香农熵，N 为检测到的克隆总数。
    值越接近 1 表示少数优势 T 细胞克隆高度寡克隆扩增。
    """
    total = float(counts.sum())
    if total <= 0 or len(counts) <= 1:
        return 0.0
    probs = counts[counts > 0] / total
    entropy = -float(np.sum(probs * np.log(probs)))
    max_entropy = math.log(len(counts))
    return max(0.0, min(1.0, 1.0 - (entropy / max_entropy)))


def generate_neoantigen_table(rng: np.random.Generator) -> List[Dict[str, object]]:
    """
    生成 24 条覆盖 missense / in-frame indel / frameshift 的确定性模拟新抗原候选肽，
    按高优先级 (0-7)、边缘灰区 (8-15) 与弱结合/低表达 (16-23) 三梯队构造特征并施加临床过滤阈值。
    """
    genes = [
        "TP53", "KRAS", "PIK3CA", "BRAF", "EGFR", "NRAS",
        "IDH1", "CTNNB1", "ARID1A", "ERBB2", "APC", "SMAD4",
        "PTEN", "CDKN2A", "RB1", "FGFR3", "BRAF", "PIK3R1",
        "SF3B1", "KMT2D", "FAT1", "NF1", "FBXW7", "STK11",
    ]
    var_types = ["missense"] * 16 + ["inframe_indel"] * 4 + ["frameshift"] * 4
    hla_alleles = ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*07:02", "HLA-B*44:02"]
    aa_alphabet = list("ACDEFGHIKLMNPQRSTVWY")

    rows: List[Dict[str, object]] = []
    for i in range(24):
        # 随机生成 9-mer 突变肽段，并仅替换第 5 位模拟野生型对应肽段
        seq = "".join(rng.choice(aa_alphabet, size=9))
        wt_seq = seq[:4] + rng.choice(aa_alphabet) + seq[5:]
        if i < 8:
            # 高优先级强结合子：MT IC50 低、DAI (WT/MT) 高、RNA 表达与 VAF 充足
            mt_ic50 = float(np.round(rng.uniform(12.0, 180.0), 2))
            wt_ic50 = float(np.round(mt_ic50 * rng.uniform(4.0, 25.0), 2))
            rank_pct = float(np.round(rng.uniform(0.05, 1.20), 2))
            rna_tpm = float(np.round(rng.uniform(8.0, 95.0), 2))
            rna_vaf = float(np.round(rng.uniform(0.15, 0.55), 3))
        elif i < 16:
            # 边缘灰区候选：用于测试多维度联合过滤边界（IC50/Rank/TPM/VAF）
            mt_ic50 = float(np.round(rng.uniform(190.0, 650.0), 2))
            wt_ic50 = float(np.round(mt_ic50 * rng.uniform(0.8, 3.0), 2))
            rank_pct = float(np.round(rng.uniform(1.10, 4.50), 2))
            rna_tpm = float(np.round(rng.uniform(0.5, 25.0), 2))
            rna_vaf = float(np.round(rng.uniform(0.03, 0.25), 3))
        else:
            # 阴性对照梯队：弱结合且几乎无表达
            mt_ic50 = float(np.round(rng.uniform(700.0, 4800.0), 2))
            wt_ic50 = float(np.round(rng.uniform(400.0, 3500.0), 2))
            rank_pct = float(np.round(rng.uniform(4.80, 18.0), 2))
            rna_tpm = float(np.round(rng.uniform(0.0, 3.0), 2))
            rna_vaf = float(np.round(rng.uniform(0.0, 0.08), 3))

        # Differential Agretopicity Index (DAI) = WT_IC50 / MT_IC50，>1 表示突变增强 MHC 结合
        dai = float(np.round(wt_ic50 / max(mt_ic50, 1e-6), 3))
        pass_filter = (
            mt_ic50 <= 500.0
            and rank_pct <= 2.0
            and rna_tpm >= 1.0
            and rna_vaf >= 0.05
        )
        rows.append(
            {
                "candidate_id": f"NEO_{i + 1:03d}",
                "gene": genes[i],
                "variant_type": var_types[i],
                "hla_allele": hla_alleles[i % len(hla_alleles)],
                "mt_peptide": seq,
                "wt_peptide": wt_seq,
                "mt_ic50_nM": mt_ic50,
                "wt_ic50_nM": wt_ic50,
                "percentile_rank": rank_pct,
                "dai_ratio": dai,
                "rna_tpm": rna_tpm,
                "rna_vaf": rna_vaf,
                "tier_pass": "PASS" if pass_filter else "FILTERED",
            }
        )
    return rows


def generate_tcr_repertoire(rng: np.random.Generator) -> List[Dict[str, object]]:
    """
    基于 Zipf 幂律分布模拟 6 例样本（3 例治疗前外周血 PBMC_Pre vs 3 例治疗后肿瘤浸润 TIL_Post）：
    TIL_Post 采用更高倾斜度 zipf_alpha，模拟抗原驱动的寡克隆优势扩增。
    """
    samples = [
        ("PBMC_Pre_1", "PBMC_Pre", 1.1),
        ("PBMC_Pre_2", "PBMC_Pre", 1.15),
        ("PBMC_Pre_3", "PBMC_Pre", 1.12),
        ("TIL_Post_1", "TIL_Post", 1.85),
        ("TIL_Post_2", "TIL_Post", 1.95),
        ("TIL_Post_3", "TIL_Post", 1.78),
    ]
    rows: List[Dict[str, object]] = []
    for sample_id, cohort, zipf_alpha in samples:
        n_clones = 80
        ranks = np.arange(1, n_clones + 1)
        weights = 1.0 / np.power(ranks, zipf_alpha)
        weights /= weights.sum()
        # 多项式抽样构造每例样本 2,500 条 TCR 测序读数
        counts = rng.multinomial(2500, weights)
        counts = np.sort(counts)[::-1]
        clonality = compute_shannon_clonality(counts)
        top1_freq = float(counts[0] / counts.sum())
        top5_freq = float(counts[:5].sum() / counts.sum())
        expanded_clones = int(np.sum((counts / counts.sum()) >= 0.01))

        rows.append(
            {
                "sample_id": sample_id,
                "cohort": cohort,
                "total_tcr_reads": int(counts.sum()),
                "unique_clonotypes": int(np.sum(counts > 0)),
                "clonality": round(clonality, 4),
                "top1_clone_freq": round(top1_freq, 4),
                "top5_clone_freq": round(top5_freq, 4),
                "expanded_clones_ge_1pct": expanded_clones,
            }
        )
    return rows


def write_tsv(path: pathlib.Path, rows: List[Dict[str, object]]) -> None:
    """将字典列表按制表符分隔 (TSV, UTF-8) 写出至目标路径。"""
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """冒烟测试主入口：构建模拟数据、施加筛选与单调性断言并写出摘要指标表。"""
    root = get_project_root()
    raw_dir = root / "data" / "raw"
    proc_dir = root / "data" / "processed"
    res_dir = root / "repro_results"
    log_dir = root / "logs"

    for d in (raw_dir, proc_dir, res_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(">>> [1/3] Constructing deterministic neoantigen & TCR smoke datasets...")
    rng = np.random.default_rng(SEED)
    neo_rows = generate_neoantigen_table(rng)
    tcr_rows = generate_tcr_repertoire(rng)

    neo_tsv = proc_dir / "smoke_neoantigen_candidates.tsv"
    tcr_tsv = proc_dir / "smoke_tcr_clonality.tsv"
    write_tsv(neo_tsv, neo_rows)
    write_tsv(tcr_tsv, tcr_rows)

    print(">>> [2/3] Evaluating prioritization filters & repertoire expansion metrics...")
    n_total_neo = len(neo_rows)
    n_pass_neo = sum(1 for r in neo_rows if r["tier_pass"] == "PASS")
    pass_dais = [float(r["dai_ratio"]) for r in neo_rows if r["tier_pass"] == "PASS"]
    median_pass_dai = float(np.median(pass_dais)) if pass_dais else 0.0

    pbmc_clonality = [float(r["clonality"]) for r in tcr_rows if r["cohort"] == "PBMC_Pre"]
    til_clonality = [float(r["clonality"]) for r in tcr_rows if r["cohort"] == "TIL_Post"]

    mean_pbmc = float(np.mean(pbmc_clonality))
    mean_til = float(np.mean(til_clonality))

    print(
        f"    Neoantigen candidates: {n_pass_neo}/{n_total_neo} passed "
        f"(MT_IC50<=500nM, Rank<=2%, TPM>=1, VAF>=0.05) | Median DAI={median_pass_dai:.2f}"
    )
    print(
        f"    TCR Clonality: PBMC_Pre mean={mean_pbmc:.4f} vs TIL_Post mean={mean_til:.4f}"
    )

    assert n_pass_neo >= 5, "Expected >= 5 high-priority neoantigens in smoke test"
    assert mean_til > mean_pbmc, "Expected TIL_Post clonality > PBMC_Pre clonality"

    print(">>> [3/3] Exporting Python smoke verification summary...")
    summary_rows = [
        {"Metric": "Total_Neoantigen_Candidates", "Value": str(n_total_neo)},
        {"Metric": "Prioritized_Neoantigens_PASS", "Value": str(n_pass_neo)},
        {"Metric": "Median_DAI_PASS", "Value": f"{median_pass_dai:.3f}"},
        {"Metric": "TCR_Samples_Evaluated", "Value": str(len(tcr_rows))},
        {"Metric": "PBMC_Pre_Mean_Clonality", "Value": f"{mean_pbmc:.4f}"},
        {"Metric": "TIL_Post_Mean_Clonality", "Value": f"{mean_til:.4f}"},
        {"Metric": "Smoke_Status", "Value": "PASS"},
    ]
    write_tsv(res_dir / "00_smoke_verification_summary.tsv", summary_rows)
    print(f"    Saved: {neo_tsv}")
    print(f"    Saved: {tcr_tsv}")
    print(f"    Saved: {res_dir / '00_smoke_verification_summary.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

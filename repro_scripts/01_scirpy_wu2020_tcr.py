#!/usr/bin/env python3
"""
01_scirpy_wu2020_tcr.py
=======================
Track A full reproduction of scverse/scirpy `tutorial_3k_tcr` (Wu et al. 2020
3,000 tumor/blood/normal-adjacent T cells with paired scRNA-seq + V(D)J TCR).

Workflow:
  1. Load `ir.datasets.wu2020_3k()` MuData (`gex` + `airr` modalities)
  2. Preprocess scRNA-seq (`gex`): gene/cell filtering, normalization, HVG, PCA,
     neighbor graph, cluster annotation mapping, and UMAP visualization
  3. TCR chain indexing & QC (`index_chains`, `chain_qc`), filtering multichain &
     orphan chains
  4. Clonotype definition (nucleotide identity + amino-acid `tcrdist` clusters)
     and clonotype network layout
  5. Clonal expansion, normalized Shannon diversity, repertoire overlap,
     clonotype modularity, and differential clonotype imbalance (`CD8_Teff` vs `CD8_Trm`)
  6. Public antigen-specificity annotation via `ir.datasets.iggytop()`
  7. Export publication-ready figures, summary TSV, and container environment versions
"""

from __future__ import annotations

import csv
import importlib.metadata
import os
import pathlib
import sys
import warnings

# Ensure reproducible cache directories inside data/raw/ before importing scanpy/pooch
PROJ_ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW_CACHE = PROJ_ROOT / "data" / "raw" / "scirpy_cache"
PROC_DIR = PROJ_ROOT / "data" / "processed"
RES_DIR = PROJ_ROOT / "repro_results" / "scirpy_wu2020"
LOG_DIR = PROJ_ROOT / "logs"

for d in (RAW_CACHE, PROC_DIR, RES_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("XDG_CACHE_HOME", str(RAW_CACHE))
os.environ.setdefault("POOCH_BASE_DIR", str(RAW_CACHE))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt
from cycler import cycler
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir

warnings.simplefilter(action="ignore", category=FutureWarning)

SEED = 20261115
np.random.seed(SEED)
sc.settings.verbosity = 2
sc.settings.datasetdir = str(RAW_CACHE)
sc.set_figure_params(figsize=(4.5, 4.2), dpi=150)


def save_fig(fig: plt.Figure, name: str) -> pathlib.Path:
    """将 Matplotlib Figure 导出为高分辨率 PNG 并立即释放内存。"""
    out = RES_DIR / name
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved figure: {out}")
    return out


def main() -> int:
    """Track A 主流程：加载 wu2020_3k 多模态数据、质控、克隆聚类与统计绘图。"""
    print("===================================================================")
    print("  Track A: scverse/scirpy wu2020_3k (scRNA-seq + TCR) Reproduction")
    print(f"  scirpy={ir.__version__}, scanpy={sc.__version__}, muon={mu.__version__}")
    print("===================================================================")

    # 1. 加载 MuData（含 gex 转录组与 airr V(D)J 免疫组库两个模态）
    print("\n>>> [1/7] Loading ir.datasets.wu2020_3k()...")
    mdata = ir.datasets.wu2020_3k()
    raw_gex_cells = int(mdata["gex"].n_obs)
    raw_gex_genes = int(mdata["gex"].n_vars)
    raw_airr_cells = int(mdata["airr"].n_obs)
    print(f"    Raw MuData: gex=({raw_gex_cells} cells x {raw_gex_genes} genes), airr={raw_airr_cells} cells")

    # 2. GEX 转录组标准预处理与 T 细胞亚群名称标准化映射
    print("\n>>> [2/7] Preprocessing GEX modality & mapping T-cell clusters...")
    sc.pp.filter_genes(mdata["gex"], min_cells=10)       # 过滤低检出基因（<10 个细胞表达）
    sc.pp.filter_cells(mdata["gex"], min_genes=100)     # 过滤空液滴/低质量细胞（<100 个基因）
    # [兼容修复] scanpy 1.13 移除了旧版 normalize_per_cell()，统一迁移至 normalize_total()
    sc.pp.normalize_total(mdata["gex"])
    sc.pp.log1p(mdata["gex"])
    sc.pp.highly_variable_genes(mdata["gex"], flavor="cell_ranger", n_top_genes=5000)
    sc.tl.pca(mdata["gex"], random_state=SEED)
    sc.pp.neighbors(mdata["gex"], random_state=SEED)

    # 优先复用原作者发表级 UMAP 坐标（X_umap_orig）以保证亚群拓扑完全对齐；缺失时按固定种子重算
    if "X_umap_orig" in mdata["gex"].obsm:
        mdata["gex"].obsm["X_umap"] = mdata["gex"].obsm["X_umap_orig"]
    else:
        sc.tl.umap(mdata["gex"], random_state=SEED)

    # 将细粒度聚类簇合并为生物学可读亚群（例如 8.3a/b/c-Trm 合并为 CD8_Trm 组织驻留记忆 T 细胞）
    mapping = {
        "nan": "other",
        "3.1-MT": "other",
        "4.1-Trm": "CD4_Trm",
        "4.2-RPL32": "CD4_RPL32",
        "4.3-TCF7": "CD4_TCF7",
        "4.4-FOS": "CD4_FOSS",
        "4.5-IL6ST": "CD4_IL6ST",
        "4.6a-Treg": "CD4_Treg",
        "4.6b-Treg": "CD4_Treg",
        "8.1-Teff": "CD8_Teff",
        "8.2-Tem": "CD8_Tem",
        "8.3a-Trm": "CD8_Trm",
        "8.3b-Trm": "CD8_Trm",
        "8.3c-Trm": "CD8_Trm",
        "8.4-Chrom": "other",
        "8.5-Mitosis": "other",
        "8.6-KLRB1": "other",
    }
    mdata["gex"].obs["cluster"] = mdata["gex"].obs["cluster_orig"].astype(str).map(mapping).fillna("other")
    mdata.update()
    mdata.pull_obs()

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.3), gridspec_kw={"wspace": 0.55})
    for col, ax in zip(["gex:sample", "gex:source", "gex:cluster"], axs):
        mu.pl.embedding(mdata, basis="gex:umap", color=col, ax=ax, show=False)
    save_fig(fig, "01_scirpy_gex_umap_overview.png")

    # 3. TCR V(D)J 链索引与双链配对质控 (Chain QC)
    print("\n>>> [3/7] Running TCR chain indexing & chain QC...")
    ir.pp.index_chains(mdata)
    ir.tl.chain_qc(mdata)

    # 统计含额外链或双胞体风险细胞（extra VJ / extra VDJ / two full chains / multichain）
    multichain_mask = mdata.obs["airr:chain_pairing"].isin(
        ["extra VJ", "extra VDJ", "two full chains", "multichain"]
    )
    extra_or_multi_frac = float(np.sum(multichain_mask) / mdata["airr"].n_obs)
    n_multichain = int(np.sum(mdata.obs["airr:chain_pairing"] == "multichain"))
    n_orphan = int(np.sum(mdata.obs["airr:chain_pairing"].isin(["orphan VDJ", "orphan VJ"])))
    print(f"    Extra/multichain fraction: {extra_or_multi_frac:.3f} | multichain={n_multichain}, orphan={n_orphan}")

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"wspace": 0.4})
    ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col="gex:source", ax=axs[0])
    ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:source", ax=axs[1])
    axs[0].set_title("Receptor Subtype by Tissue Source")
    axs[1].set_title("TCR Chain Pairing QC (Pre-filter)")
    save_fig(fig, "02_scirpy_chain_qc.png")

    # 剔除多链 (multichain=17, 疑似双胞体) 与缺失配对链 (orphan=1138, 无法完整组装 TCR-T)，保留 1,845 个完整配对细胞
    mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
    mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))
    qc_cells = int(mdata.n_obs)
    print(f"    Cells retained after multichain/orphan QC: {qc_cells}")

    # 4. 定义克隆型（两层粒度：核苷酸全同 clone_id vs 氨基酸 tcrdist 相似簇）与克隆网络
    print("\n>>> [4/7] Defining clonotypes (nt identity & aa tcrdist) and network...")
    ir.pp.ir_dist(mdata)
    # 精细克隆型：基于主链 (primary_only) 核苷酸 CDR3 完全一致，追踪同一祖细胞克隆来源
    ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
    ir.tl.clonotype_network(mdata, min_cells=2)

    # 粗粒度相似簇：基于氨基酸 tcrdist 进化距离 (cutoff=15)，捕获趋同抗原识别的相似 TCR
    ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=15)
    ir.tl.define_clonotype_clusters(
        mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any"
    )
    ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")

    n_nt_clonotypes = int(mdata.obs["airr:clone_id"].nunique(dropna=True))
    n_aa_clusters = int(mdata.obs["airr:cc_aa_tcrdist"].nunique(dropna=True))
    print(f"    Unique NT clonotypes: {n_nt_clonotypes} | AA tcrdist clusters: {n_aa_clusters}")

    ax_net = ir.pl.clonotype_network(
        mdata, color="gex:source", base_size=20, label_fontsize=9, panel_size=(7, 6)
    )
    save_fig(ax_net.figure, "03_scirpy_clonotype_network.png")

    # 5. 克隆扩增分级与归一化香农熵多样性计算
    print("\n>>> [5/7] Computing clonal expansion & Shannon diversity...")
    ir.tl.clonal_expansion(mdata)

    fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.8), gridspec_kw={"wspace": 0.45})
    ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:cluster", ax=axs[0])
    ir.pl.alpha_diversity(
        mdata, metric="normalized_shannon_entropy", groupby="gex:cluster", ax=axs[1]
    )
    axs[0].set_title("Clonal Expansion by T-Cell Cluster")
    axs[1].set_title("Normalized Shannon Entropy by Cluster")
    save_fig(fig, "04_scirpy_clonal_expansion_and_diversity.png")

    # 6. 样本间组库重叠、克隆模块化置换检验及 CD8_Teff vs CD8_Trm 差异丰度检验
    print("\n>>> [6/7] Computing repertoire overlap, modularity & cluster imbalance...")
    overlap_df, _, _ = ir.tl.repertoire_overlap(mdata, "gex:sample", inplace=False)
    cg_overlap = ir.pl.repertoire_overlap(
        mdata,
        "gex:sample",
        heatmap_cats=["gex:patient", "gex:source"],
        yticklabels=True,
        xticklabels=True,
    )
    cg_overlap.fig.savefig(
        RES_DIR / "05_scirpy_repertoire_overlap.png", dpi=180, bbox_inches="tight"
    )
    plt.close(cg_overlap.fig)
    print(f"    Saved figure: {RES_DIR / '05_scirpy_repertoire_overlap.png'}")

    # 基于随机子图背景评估哪些 AA tcrdist 克隆簇在转录组状态空间中显著聚集 (FDR 控制)
    ir.tl.clonotype_modularity(mdata, target_col="airr:cc_aa_tcrdist")
    ax_mod = ir.pl.clonotype_modularity(mdata, base_size=20)
    save_fig(ax_mod.figure, "06_scirpy_clonotype_modularity.png")

    # 跨样本统计检验：比较效应 CD8_Teff 与组织驻留记忆 CD8_Trm 间的克隆占比倾斜
    _, imbalance_stat = ir.tl.clonotype_imbalance(
        mdata,
        replicate_col="gex:sample",
        groupby="gex:cluster",
        case_label="CD8_Teff",
        control_label="CD8_Trm",
        inplace=False,
    )
    imbalance_stat.to_csv(RES_DIR / "01_scirpy_cd8_teff_vs_trm_imbalance.tsv", sep="\t", index=False)

    # 7. 克隆趋同性检测与公共抗原表位数据库 (iggytop) 注释（含网络受限优雅降级保护）
    print("\n>>> [7/7] Computing clonotype convergence & querying antigen specificity (iggytop)...")
    # 同一 AA tcrdist 簇内包含多个不同 NT clone_id 的细胞，提示抗原驱动的趋同重组
    ir.tl.clonotype_convergence(mdata, key_coarse="cc_aa_tcrdist", key_fine="clone_id")
    n_convergent = int((mdata.obs["airr:is_convergent"] == "convergent").sum())
    print(f"    Convergent TCR cells (same AA tcrdist, distinct NT clone_id): {n_convergent}")

    iggytop_status = "OK"
    n_annotated_cells: object = 0
    try:
        iggytop = ir.datasets.iggytop()
        ir.pp.ir_dist(mdata, iggytop, metric="identity", sequence="aa")
        ir.tl.ir_query(
            mdata,
            iggytop,
            metric="identity",
            sequence="aa",
            receptor_arms="any",
            dual_ir="any",
        )
        ir.tl.ir_query_annotate(
            mdata,
            iggytop,
            metric="identity",
            sequence="aa",
            include_ref_cols=["antigen_species"],
            strategy="most-frequent",
        )
        n_annotated_cells = int(mdata.obs["airr:antigen_species"].notna().sum())
        print(f"    Cells matched to iggytop antigen_species: {n_annotated_cells}/{qc_cells}")
    except Exception as exc:
        # 容器内若无法访问 GitHub Release Assets CDN (Connection refused)，自动降级展示本地趋同性 UMAP
        iggytop_status = f"SKIPPED ({type(exc).__name__})"
        n_annotated_cells = "NA (iggytop download unreachable)"
        print(f"    [WARN] Skipping iggytop remote query: {exc}")

    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.5), gridspec_kw={"wspace": 0.45})
    mu.pl.embedding(mdata, "gex:umap", color="airr:clonal_expansion", ax=axs[0], show=False)
    axs[0].set_title("Clonal Expansion on GEX UMAP")
    if iggytop_status == "OK":
        mu.pl.embedding(mdata, "gex:umap", color="airr:antigen_species", ax=axs[1], show=False)
        axs[1].set_title("Iggytop Antigen Specificity Match")
    else:
        mu.pl.embedding(mdata, "gex:umap", color="airr:is_convergent", ax=axs[1], show=False)
        axs[1].set_title("TCR Convergent Clonotypes on UMAP")
    save_fig(fig, "07_scirpy_expansion_and_antigen_umap.png")

    # Save processed MuData & quantitative summary
    out_h5mu = PROC_DIR / "wu2020_3k_scirpy_processed.h5mu"
    mdata.write(out_h5mu)
    print(f"    Saved processed MuData: {out_h5mu}")

    summary_rows = [
        ("Raw_GEX_Cells", raw_gex_cells),
        ("Raw_GEX_Genes", raw_gex_genes),
        ("Raw_AIRR_Cells", raw_airr_cells),
        ("Extra_or_Multichain_Fraction", round(extra_or_multi_frac, 4)),
        ("Filtered_Multichain_Cells", n_multichain),
        ("Filtered_Orphan_Cells", n_orphan),
        ("PostQC_Retained_Cells", qc_cells),
        ("Unique_NT_Clonotypes", n_nt_clonotypes),
        ("AA_TCRdist_Clusters", n_aa_clusters),
        ("Convergent_TCR_Cells", n_convergent),
        ("Repertoire_Overlap_Samples", int(overlap_df.shape[0])),
        ("Iggytop_Query_Status", iggytop_status),
        ("Iggytop_Annotated_Cells", n_annotated_cells),
        ("Random_Seed", SEED),
    ]
    summary_path = RES_DIR / "01_scirpy_wu2020_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Metric", "Value"])
        w.writerows(summary_rows)
    print(f"    Saved summary TSV: {summary_path}")

    # Export container package versions
    pkgs = ["scirpy", "scanpy", "anndata", "mudata", "muon", "leidenalg", "igraph", "numpy", "pandas", "scipy", "matplotlib", "parasail", "scikit-bio", "pooch"]
    ver_path = RES_DIR / "scirpy_container_versions.tsv"
    with ver_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Package", "Version"])
        w.writerow(["python", sys.version.split()[0]])
        for p in pkgs:
            try:
                w.writerow([p, importlib.metadata.version(p)])
            except importlib.metadata.PackageNotFoundError:
                w.writerow([p, "Not Installed"])
    print(f"    Saved container versions: {ver_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

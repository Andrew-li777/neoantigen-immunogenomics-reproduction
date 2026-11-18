# 清洗与分析中间产物说明 (`data/processed/`)

> **管理规范**：本目录保存流程生成的标准化中间对象与全量候选新抗原文件。为避免 Git 仓库膨胀，体积较大的 `.h5mu`、`pvacseq_hcc1395/` 与 `pvacfuse_hcc1395/`（含数十 MB 全量表位表与自动生成的 `pVACview` R Shiny 脚手架）均受 [`.gitignore`](../../.gitignore) 屏蔽；精简后的定量摘要表（`*_summary.tsv`）、Top 候选表与发表级图表已统一导出至 [`repro_results/`](../../repro_results/) 供直接查看与版本控制。

---

## 一、冒烟自检中间表（Smoke Test）

由 [00_smoke_test.py](../../repro_scripts/00_smoke_test.py) 基于确定性种子 `seed=20261115` 生成，用于验证 Python → TSV → R 端联调链路：

| 文件名 | 规模 | 核心字段说明 |
|---|---|---|
| `smoke_neoantigen_candidates.tsv` | `24` 条候选肽 | `candidate_id`, `gene`, `variant_type` (`missense`/`inframe_indel`/`frameshift`), `hla_allele`, `mt_peptide`, `wt_peptide`, `mt_ic50_nM`, `wt_ic50_nM`, `percentile_rank`, `dai_ratio`, `rna_tpm`, `rna_vaf`, `tier_pass` (`PASS`/`FILTERED`) |
| `smoke_tcr_clonality.tsv` | `6` 例模拟样本 | `sample_id`, `cohort` (`PBMC_Pre`/`TIL_Post`), `total_tcr_reads`, `unique_clonotypes`, `clonality`, `top1_clone_freq`, `top5_clone_freq`, `expanded_clones_ge_1pct` |

---

## 二、Track A 处理后多模态对象（`scirpy` `wu2020_3k`）

| 文件名 | 体积 | 内部新增模态/观测列 |
|---|---:|---|
| `wu2020_3k_scirpy_processed.h5mu` | `27.75 MB` | 经 [01_scirpy_wu2020_tcr.py](../../repro_scripts/01_scirpy_wu2020_tcr.py) 质控过滤后保留 **`1,845`** 个完整双链配对 T 细胞，新增 `gex:cluster`、`airr:receptor_type/subtype`、`airr:chain_pairing`、`airr:clone_id`（`1,526` 类）、`airr:cc_aa_tcrdist`（`1,490` 簇）、`airr:clonal_expansion`、`airr:clonotype_modularity`、`airr:clonotype_modularity_fdr` 与 `airr:is_convergent`（`174` 个收敛细胞） |

---

## 三、Track B-1 产出：`pvacseq_hcc1395/`（SNV / Indel 新抗原）

由 `neoantigen-pvactools:7.1.5` 容器执行 `pvacseq run` 生成（[02_run_pvacseq_hcc1395.ps1](../../repro_scripts/02_run_pvacseq_hcc1395.ps1)），包含三个子目录：

| 子目录 | 核心产出文件 | 实测规模与说明 |
|---|---|---|
| `combined/` | • `HCC1395_TUMOR_DNA.Combined.all_epitopes.aggregated.tsv`<br>• `HCC1395_TUMOR_DNA.Combined.filtered.tsv`<br>• `HCC1395_TUMOR_DNA.Combined.all_epitopes.tsv`<br>• `app.R` / `server.R` / `ui.R` / `www/` | • `93,411` B：变异级汇总（**`357`** 个编码突变，含 **`48`** 个 `Pass`、`7`** 个 `Subclonal`）<br>• `30,774` B：联合过滤表位（**`52`** 条表位，对应 `50` 个唯一变异；因 9 个变异同时通过 Class I & II 过滤，合并后按 Top Score 优选为 52 行）<br>• `48,226,890` B：未过滤全量肽段预测表<br>• `pVACtools` 自动生成的 `pVACview` 交互式 Shiny 应用代码（本地可直接 `shiny::runApp()` 浏览） |
| `MHC_Class_I/` | • `HCC1395_TUMOR_DNA.MHC_I.filtered.tsv`<br>• `HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv`<br>• `HCC1395_TUMOR_DNA.proximal_variants.tsv`<br>• `HCC1395_TUMOR_DNA.tsv_1-100` ~ `_901-955` | • `25,542` B：**`44`** 条 Class I（8/9-mer）过滤表位（覆盖 `42` 个唯一变异）<br>• `38,918,523` B：Class I 全量表位<br>• `2,052` B：**`20`** 条顺式相邻定相校正记录<br>• 955 个待预测变异分片中间文件 |
| `MHC_Class_II/` | • `HCC1395_TUMOR_DNA.MHC_II.filtered.tsv`<br>• `HCC1395_TUMOR_DNA.MHC_II.all_epitopes.tsv` | • `10,999` B：**`17`** 条 Class II（15-mer, `DRB1*04:05`）过滤表位<br>• `8,306,390` B：Class II 全量表位 |

---

## 四、Track B-2 产出：`pvacfuse_hcc1395/`（基因融合新抗原）

由 `neoantigen-pvactools:7.1.5` 容器执行 `pvacfuse run`（联合 `agfusion_results/` 与 `star-fusion.fusion_predictions.tsv`）生成（[03_run_pvacfuse_hcc1395.ps1](../../repro_scripts/03_run_pvacfuse_hcc1395.ps1)）：

| 子目录 | 核心产出文件 | 实测规模与说明 |
|---|---|---|
| `combined/` | • `HCC1395_TUMOR_RNA.Combined.all_epitopes.aggregated.tsv`<br>• `HCC1395_TUMOR_RNA.Combined.filtered.tsv`<br>• `app.R` / `server.R` / `ui.R` | • **`31`** 条聚合融合转录本异构体（涉及 `5` 个融合基因；其中 **`16`** 条异构体 / `3` 个融合基因达 `Tier Pass`：`TMEM104_RAB37`、`EIF3K_CYP39A1`、`HNRNPUL2_AHNAK`）<br>• **`14`** 条联合过滤融合新抗原表位 |
| `MHC_Class_I/` | • `HCC1395_TUMOR_RNA.MHC_I.filtered.tsv` | **`14`** 条 Class I（8/9-mer）过滤融合表位（如 `YRIYLTIL` `5.15 nM`、`AEMLGDLSA` `16.76 nM`、`FWHQFQRM` `24.93 nM`） |
| `MHC_Class_II/` | • `HCC1395_TUMOR_RNA.MHC_II.filtered.tsv` | **`6`** 条 Class II（15-mer, `DRB1*04:05`）过滤融合表位 |

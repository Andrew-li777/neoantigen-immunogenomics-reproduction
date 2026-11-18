# ==============================================================================
# 02_pvacseq_hcc1395_analysis.R
# Track B (griffithlab/pVACtools v7.1.5) HCC1395 下游统计汇总与临床优先级绘图
#
# 输入:
#   - data/processed/pvacseq_hcc1395/combined/HCC1395_TUMOR_DNA.Combined.all_epitopes.aggregated.tsv
#   - data/processed/pvacseq_hcc1395/combined/HCC1395_TUMOR_DNA.Combined.filtered.tsv
#   - data/processed/pvacseq_hcc1395/MHC_Class_I/HCC1395_TUMOR_DNA.MHC_I.filtered.tsv
#   - data/processed/pvacseq_hcc1395/MHC_Class_II/HCC1395_TUMOR_DNA.MHC_II.filtered.tsv
# 输出:
#   - repro_results/pvacseq_hcc1395/02_pvacseq_hcc1395_summary.tsv
#   - repro_results/pvacseq_hcc1395/02_pvacseq_top_candidates.tsv
#   - repro_results/pvacseq_hcc1395/01_pvacseq_tier_and_allele_distribution.png
#   - repro_results/pvacseq_hcc1395/02_pvacseq_ic50_vs_expression_and_dai.png
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# 定位项目根目录（保证从任意工作目录运行 Rscript 均能找到 data/processed 与 repro_results/）
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  getwd()
}

proj_root <- get_project_root()
pvac_dir <- file.path(proj_root, "data", "processed", "pvacseq_hcc1395")
res_dir <- file.path(proj_root, "repro_results", "pvacseq_hcc1395")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# 4 份 pvacseq 核心产物路径与口径说明：
#   - aggregated.tsv (357 行): 变异级汇总表，记录各变异 Tier 分级与最佳候选肽
#   - Combined.filtered.tsv (52 行 / 50 个唯一变异): 联合 Class I & II 通过过滤的表位
#     （注：44 MHC-I + 17 MHC-II 共涉及 50 个唯一变异，其中 9 个变异同时通过 I/II 类过滤，
#       经 Top Score Filter 按 (variant, transcript) 优选后输出 52 条联合表位）
#   - proximal_variants.tsv (20 条): --phased-proximal-variants-vcf 检出的顺式定相近端变异记录
agg_path <- file.path(pvac_dir, "combined", "HCC1395_TUMOR_DNA.Combined.all_epitopes.aggregated.tsv")
filt_path <- file.path(pvac_dir, "combined", "HCC1395_TUMOR_DNA.Combined.filtered.tsv")
mhc1_path <- file.path(pvac_dir, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.filtered.tsv")
mhc2_path <- file.path(pvac_dir, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.filtered.tsv")
prox_path <- file.path(pvac_dir, "MHC_Class_I", "HCC1395_TUMOR_DNA.proximal_variants.tsv")

# 统一把 "NA" / "None" / "" 识别为缺失值，防止数值列被读成字符型
agg_df <- fread(agg_path, na.strings = c("NA", "None", ""))
filt_df <- fread(filt_path, na.strings = c("NA", "None", ""))
mhc1_df <- fread(mhc1_path, na.strings = c("NA", "None", ""))
mhc2_df <- fread(mhc2_path, na.strings = c("NA", "None", ""))
prox_n <- if (file.exists(prox_path)) max(0L, length(readLines(prox_path)) - 1L) else 0L

cat(sprintf(">>> Aggregated variants: %d | Combined filtered epitopes: %d (MHC-I: %d, MHC-II: %d)\n",
            nrow(agg_df), nrow(filt_df), nrow(mhc1_df), nrow(mhc2_df)))

# 强制转换数值列并计算 Differential Agretopicity Index (DAI = IC50_WT / IC50_MT)
# DAI > 1 表示突变肽比对应野生型肽更强结合患者 MHC（降低中心免疫耐受风险）
num_cols_agg <- c("IC50 MT", "IC50 WT", "%ile MT", "RNA Expr", "RNA VAF", "Allele Expr", "DNA VAF", "Num Passing Peptides")
for (col in num_cols_agg) {
  if (col %in% names(agg_df)) agg_df[, (col) := suppressWarnings(as.numeric(get(col)))]
}
agg_df[, DAI := fifelse(!is.na(`IC50 WT`) & !is.na(`IC50 MT`) & `IC50 MT` > 0, `IC50 WT` / `IC50 MT`, NA_real_)]

# 1. 构建变异级与表位级定量汇总表（记录 Tier Pass/Subclonal 及 GOI 命中数）
tier_counts <- agg_df[, .N, by = Tier][order(-N)]
n_pass_vars <- agg_df[Tier == "Pass", .N]
n_subclonal <- agg_df[Tier == "Subclonal", .N]
n_goi_filt <- filt_df[`Gene of Interest` == TRUE | `Gene of Interest` == "True", .N]
med_ic50_pass <- agg_df[Tier == "Pass", median(`IC50 MT`, na.rm = TRUE)]
med_dai_pass <- agg_df[Tier == "Pass", median(DAI, na.rm = TRUE)]

summary_dt <- data.table(
  Metric = c(
    "Input_Sample",
    "Total_Aggregated_Variants",
    "Tier_Pass_Variants",
    "Tier_Subclonal_Variants",
    "Combined_Filtered_Epitopes",
    "MHC_Class_I_Filtered_Epitopes",
    "MHC_Class_II_Filtered_Epitopes",
    "Filtered_Epitopes_In_Genes_Of_Interest",
    "Phased_Proximal_Variant_Records",
    "Median_Best_MT_IC50_Pass_nM",
    "Median_DAI_Pass"
  ),
  Value = c(
    "HCC1395_TUMOR_DNA",
    as.character(nrow(agg_df)),
    as.character(n_pass_vars),
    as.character(n_subclonal),
    as.character(nrow(filt_df)),
    as.character(nrow(mhc1_df)),
    as.character(nrow(mhc2_df)),
    as.character(n_goi_filt),
    as.character(prox_n),
    sprintf("%.2f", med_ic50_pass),
    sprintf("%.3f", med_dai_pass)
  )
)
fwrite(summary_dt, file.path(res_dir, "02_pvacseq_hcc1395_summary.tsv"), sep = "\t")

# 2. 按最佳突变体亲和力 (IC50 MT) 升序、等位基因特异性表达量降序导出 Top 30 Tier-Pass 候选靶点
top_dt <- agg_df[Tier == "Pass"][order(`IC50 MT`, -`Allele Expr`), .(
  ID, Gene, AA_Change = `AA Change`, Best_Peptide = `Best Peptide`,
  Allele, IC50_MT_nM = `IC50 MT`, IC50_WT_nM = `IC50 WT`,
  Percentile_MT = `%ile MT`, DAI = round(DAI, 3),
  RNA_Expr = `RNA Expr`, RNA_VAF = `RNA VAF`, Allele_Expr = `Allele Expr`,
  DNA_VAF = `DNA VAF`, Num_Passing_Peptides = `Num Passing Peptides`, Tier
)]
fwrite(head(top_dt, 30), file.path(res_dir, "02_pvacseq_top_candidates.tsv"), sep = "\t")

# 3. 图 1：Panel A 统计临床分级 (Tier) 数量；Panel B 统计过滤表位在各患者 HLA 等位基因与肽段长度 (8/9/15aa) 分布
tier_counts[, Tier := factor(Tier, levels = rev(Tier))]
p1a <- ggplot(tier_counts, aes(x = Tier, y = N, fill = Tier == "Pass")) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.15, size = 3.5) +
  coord_flip(ylim = c(0, max(tier_counts$N) * 1.15)) +
  scale_fill_manual(values = c("TRUE" = "#d95f02", "FALSE" = "#4c78a8"), guide = "none") +
  labs(
    title = "A. Aggregated Variant Tiers (HCC1395)",
    subtitle = sprintf("Total: %d variants | Tier Pass: %d | Subclonal: %d",
                       nrow(agg_df), n_pass_vars, n_subclonal),
    x = "pVACseq Clinical Tier",
    y = "Variant Count"
  ) +
  theme_bw(base_size = 11)

allele_len <- filt_df[, .N, by = .(`HLA Allele`, `Peptide Length` = factor(`Peptide Length`))]
p1b <- ggplot(allele_len, aes(x = `HLA Allele`, y = N, fill = `Peptide Length`)) +
  geom_col(position = "stack", width = 0.65) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "B. Filtered Neoepitopes by HLA Allele & Length",
    subtitle = sprintf("Combined filtered: %d epitopes (MHC-I: %d, MHC-II: %d)",
                       nrow(filt_df), nrow(mhc1_df), nrow(mhc2_df)),
    x = "Patient HLA Allele",
    y = "Filtered Epitope Count",
    fill = "Peptide Length (aa)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig1 <- p1a + p1b + plot_layout(widths = c(1, 1.15))
ggsave(file.path(res_dir, "01_pvacseq_tier_and_allele_distribution.png"),
       fig1, width = 12.2, height = 4.8, dpi = 170)

# 4. 图 2：Panel A 展示聚合变异的 MHC 结合力 vs 等位基因特异性表达量；Panel B 展示过滤表位 IC50 vs 差异免疫原性 (DAI)
#    [视觉修复] 此前强结合靶点 FARP1/ADPRHL1 位于 log10(IC50) 最左侧易贴边截断，故使用 ggrepel + 左侧 14% 留白防裁剪
plot_agg <- agg_df[!is.na(`IC50 MT`) & `IC50 MT` > 0 & !is.na(`Allele Expr`)]
plot_agg[, Is_Pass := fifelse(Tier == "Pass", "Pass", "Other Tiers")]
top_labels <- head(plot_agg[Tier == "Pass"][order(`IC50 MT`)], 8)

label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = paste0(Gene, " (", `AA Change`, ")")),
    size = 2.9, color = "#b30000", box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0.1, max.overlaps = 20, seed = 20261115,
    show.legend = FALSE
  )
} else {
  geom_text(
    data = top_labels,
    aes(label = paste0(Gene, " (", `AA Change`, ")")),
    hjust = -0.08, vjust = 0.3, size = 2.9, color = "#b30000",
    check_overlap = TRUE, show.legend = FALSE
  )
}

p2a <- ggplot(plot_agg, aes(x = log10(`IC50 MT`), y = log10(`Allele Expr` + 1), color = Is_Pass)) +
  geom_vline(xintercept = log10(500), linetype = "dashed", color = "#555555") +
  geom_point(aes(alpha = Is_Pass, size = `DNA VAF`)) +
  label_layer +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.06))) +
  scale_color_manual(values = c("Pass" = "#d95f02", "Other Tiers" = "#8da0cb")) +
  scale_alpha_manual(values = c("Pass" = 0.9, "Other Tiers" = 0.45), guide = "none") +
  scale_size_continuous(range = c(1.5, 4.2)) +
  labs(
    title = "A. MHC Binding vs Allele-Specific Expression",
    subtitle = "Top Tier-Pass clonal neoantigens labeled",
    x = "log10(Best Mutant IC50 [nM])",
    y = "log10(Allele-Specific RNA Expr + 1)",
    color = "Clinical Tier",
    size = "DNA VAF"
  ) +
  theme_bw(base_size = 11)

filt_plot <- copy(filt_df)
filt_plot[, MT_IC50 := suppressWarnings(as.numeric(`Best MT IC50 Score`))]
filt_plot[, Fold_Change := suppressWarnings(as.numeric(`Corresponding Fold Change`))]
filt_plot[, RNA_Expr := suppressWarnings(as.numeric(`Gene Expression`))]
filt_plot <- filt_plot[!is.na(MT_IC50) & MT_IC50 > 0 & !is.na(Fold_Change) & Fold_Change > 0]

p2b <- ggplot(filt_plot, aes(x = log10(MT_IC50), y = log2(Fold_Change),
                             color = `HLA Allele`, size = log10(RNA_Expr + 1))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "#666666") +
  geom_vline(xintercept = log10(500), linetype = "dashed", color = "#555555") +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "B. Filtered Epitopes: Binding vs Agretopicity (DAI)",
    subtitle = "Points above y=0 bind stronger as mutant than wild-type",
    x = "log10(Best MT IC50 [nM])",
    y = "log2(DAI = WT_IC50 / MT_IC50)",
    color = "HLA Allele",
    size = "log10(RNA Expr + 1)"
  ) +
  theme_bw(base_size = 11)

fig2 <- p2a + p2b + plot_layout(widths = c(1.05, 1.1))
ggsave(file.path(res_dir, "02_pvacseq_ic50_vs_expression_and_dai.png"),
       fig2, width = 12.6, height = 5.0, dpi = 170)

cat(">>> Exported Track B pVACseq figures & summary tables to:", res_dir, "\n")

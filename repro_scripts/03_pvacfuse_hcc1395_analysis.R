# ==============================================================================
# 03_pvacfuse_hcc1395_analysis.R
# Track B-2 (griffithlab/pVACtools v7.1.5 - pvacfuse) HCC1395 融合新抗原下游分析
#
# 功能说明:
#   1. 读取 pvacfuse 输出的聚合转录本表 (31 行)、联合过滤表位 (14 行) 与 I/II 类过滤表
#   2. 统计各临床分级 (Pass / PoorBinder / LowExpr)、融合基因数与读数/表达量支持度
#   3. 导出汇总 TSV 与 Top 融合新抗原候选表，并绘制融合事件分级与 IC50 vs 读数支持度双联图
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# 自适应定位工作区根目录（兼容命令行与 RStudio）
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
fuse_dir <- file.path(proj_root, "data", "processed", "pvacfuse_hcc1395")
res_dir <- file.path(proj_root, "repro_results", "pvacfuse_hcc1395")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

agg_path <- file.path(fuse_dir, "combined", "HCC1395_TUMOR_RNA.Combined.all_epitopes.aggregated.tsv")
filt_path <- file.path(fuse_dir, "combined", "HCC1395_TUMOR_RNA.Combined.filtered.tsv")
mhc1_path <- file.path(fuse_dir, "MHC_Class_I", "HCC1395_TUMOR_RNA.MHC_I.filtered.tsv")
mhc2_path <- file.path(fuse_dir, "MHC_Class_II", "HCC1395_TUMOR_RNA.MHC_II.filtered.tsv")

# 统一把 "NA" / "None" / "" 识别为缺失值，避免数值列被强制读成字符型
agg_df <- fread(agg_path, na.strings = c("NA", "None", ""))
filt_df <- fread(filt_path, na.strings = c("NA", "None", ""))
mhc1_df <- fread(mhc1_path, na.strings = c("NA", "None", ""))
mhc2_df <- fread(mhc2_path, na.strings = c("NA", "None", ""))

num_cols <- c("IC50 MT", "%ile MT", "Expr", "Read Support", "Num Passing Peptides")
for (col in num_cols) {
  if (col %in% names(agg_df)) agg_df[, (col) := suppressWarnings(as.numeric(get(col)))]
}

cat(sprintf(">>> Fusion aggregated isoforms: %d | Combined filtered epitopes: %d (MHC-I: %d, MHC-II: %d)\n",
            nrow(agg_df), nrow(filt_df), nrow(mhc1_df), nrow(mhc2_df)))

# 1. 导出融合新抗原汇总指标表
n_pass_iso <- agg_df[Tier == "Pass", .N]
n_pass_genes <- uniqueN(agg_df[Tier == "Pass", Gene])
med_ic50_pass <- agg_df[Tier == "Pass", median(`IC50 MT`, na.rm = TRUE)]

summary_dt <- data.table(
  Metric = c(
    "Input_Sample",
    "Total_Aggregated_Fusion_Isoforms",
    "Unique_Aggregated_Fusion_Genes",
    "Tier_Pass_Fusion_Isoforms",
    "Tier_Pass_Unique_Fusion_Genes",
    "Combined_Filtered_Fusion_Epitopes",
    "MHC_Class_I_Filtered_Epitopes",
    "MHC_Class_II_Filtered_Epitopes",
    "Median_Best_MT_IC50_Pass_nM"
  ),
  Value = c(
    "HCC1395_TUMOR_RNA",
    as.character(nrow(agg_df)),
    as.character(uniqueN(agg_df$Gene)),
    as.character(n_pass_iso),
    as.character(n_pass_genes),
    as.character(nrow(filt_df)),
    as.character(nrow(mhc1_df)),
    as.character(nrow(mhc2_df)),
    sprintf("%.2f", med_ic50_pass)
  )
)
fwrite(summary_dt, file.path(res_dir, "03_pvacfuse_hcc1395_summary.tsv"), sep = "\t")

# 2. 按最佳突变体 IC50 与 Read Support 排序，导出融合新抗原候选明细
top_dt <- agg_df[order(`IC50 MT`, -`Read Support`), .(
  ID, Gene, Best_Transcript = `Best Transcript`, Best_Peptide = `Best Peptide`,
  Allele, IC50_MT_nM = `IC50 MT`, Percentile_MT = `%ile MT`,
  Expr_FFPM = Expr, Read_Support = `Read Support`,
  Num_Passing_Peptides = `Num Passing Peptides`, Tier
)]
fwrite(top_dt, file.path(res_dir, "03_pvacfuse_top_candidates.tsv"), sep = "\t")

# 3. 双联图：Panel A 为临床 Tier 分层计数；Panel B 为结合力 vs 融合表达/读数支持散点图
tier_counts <- agg_df[, .N, by = Tier][order(-N)]
tier_counts[, Tier := factor(Tier, levels = rev(Tier))]

p1 <- ggplot(tier_counts, aes(x = Tier, y = N, fill = Tier == "Pass")) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.8) +
  coord_flip(ylim = c(0, max(tier_counts$N) * 1.2)) +
  scale_fill_manual(values = c("TRUE" = "#d95f02", "FALSE" = "#4c78a8"), guide = "none") +
  labs(
    title = "A. Fusion Isoform Clinical Tiers (HCC1395)",
    subtitle = sprintf("Total isoforms: %d (%d genes) | Pass: %d (%d genes)",
                       nrow(agg_df), uniqueN(agg_df$Gene), n_pass_iso, n_pass_genes),
    x = "pVACfuse Clinical Tier",
    y = "Aggregated Fusion Isoform Count"
  ) +
  theme_bw(base_size = 11)

# 整理过滤后融合表位坐标数据，对缺失 Read Support 的异构体做安全标记避免 log 报错
plot_df <- copy(agg_df[!is.na(`IC50 MT`) & `IC50 MT` > 0])
plot_df[, Read_Plot := fifelse(is.na(`Read Support`) | `Read Support` <= 0, 1, `Read Support`)]
plot_df[, Has_StarFusion := fifelse(!is.na(`Read Support`), "STAR-Fusion Quantified", "No Read Support")]

# 每个融合基因取一条最佳异构体标注标签，保持图面清晰不重复
gene_labels <- plot_df[Tier == "Pass", .SD[which.min(`IC50 MT`)], by = Gene]

label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(
    data = gene_labels,
    aes(label = paste0(Gene, "\n(", `Best Peptide`, ", ", sprintf("%.1f nM", `IC50 MT`), ")")),
    size = 2.9, color = "#b30000", box.padding = 0.75, point.padding = 0.5,
    nudge_x = -0.12, nudge_y = 0.12, min.segment.length = 0.05, seed = 20261115,
    show.legend = FALSE
  )
} else {
  geom_text(
    data = gene_labels,
    aes(label = paste0(Gene, " (", `Best Peptide`, ")")),
    vjust = -0.8, size = 3.0, color = "#b30000", check_overlap = TRUE, show.legend = FALSE
  )
}

p2 <- ggplot(plot_df, aes(x = log10(`IC50 MT`), y = log10(Read_Plot + 1), color = Tier, shape = Has_StarFusion)) +
  geom_vline(xintercept = log10(500), linetype = "dashed", color = "#555555") +
  geom_point(size = 3.2, alpha = 0.85) +
  label_layer +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.16))) +
  scale_color_manual(values = c("Pass" = "#d95f02", "PoorBinder" = "#7570b3", "LowExpr" = "#4c78a8")) +
  labs(
    title = "B. Fusion Neoepitopes: MHC Binding vs Read Support",
    subtitle = sprintf("Combined filtered: %d epitopes (MHC-I: %d, MHC-II: %d)",
                       nrow(filt_df), nrow(mhc1_df), nrow(mhc2_df)),
    x = "log10(Best Mutant IC50 [nM])",
    y = "log10(Junction+Spanning Read Support + 1)",
    color = "Clinical Tier",
    shape = "STAR-Fusion Link"
  ) +
  theme_bw(base_size = 11)

fig <- p1 + p2 + plot_layout(widths = c(1, 1.35))
ggsave(file.path(res_dir, "03_pvacfuse_fusion_neoepitopes.png"), fig, width = 12.6, height = 5.0, dpi = 170)
cat(">>> Exported Track B-2 pVACfuse summary & figure to:", res_dir, "\n")

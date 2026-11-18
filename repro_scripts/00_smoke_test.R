# ==============================================================================
# 00_smoke_test.R
# REPRO-D 冒烟测试 R 端可视化与环境校验脚本
#
# 功能:
#   1. 读取 00_smoke_test.py 输出的新抗原优先级表与 TCR 克隆性汇总表
#   2. 生成双联质控图:
#      - Panel A: 突变肽段 MT_IC50 vs DAI 免疫原性散点筛选图 (标注 Tier PASS 阈值线)
#      - Panel B: PBMC_Pre vs TIL_Post TCR 克隆性 (Clonality) 与扩增克隆数对比图
#   3. 导出 repro_results/00_smoke_neoantigen_tcr_qc.png 与 logs/sessionInfo.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) {
  suppressPackageStartupMessages(library(patchwork))
}

# 动态解析脚本所在位置并定位工作区根目录（兼容从仓库根目录或 scripts 目录执行）
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
proc_dir <- file.path(proj_root, "data", "processed")
res_dir <- file.path(proj_root, "repro_results")
log_dir <- file.path(proj_root, "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 1. 读取 Python 冒烟脚本生成的 TSV 输入并校验存在性
neo_path <- file.path(proc_dir, "smoke_neoantigen_candidates.tsv")
tcr_path <- file.path(proc_dir, "smoke_tcr_clonality.tsv")

if (!file.exists(neo_path) || !file.exists(tcr_path)) {
  stop("Smoke test TSV inputs not found. Run 00_smoke_test.py first.")
}

neo_df <- fread(neo_path)
tcr_df <- fread(tcr_path)

# 固定因子水平顺序，保证图例与坐标轴分组顺序稳定
neo_df[, tier_pass := factor(tier_pass, levels = c("PASS", "FILTERED"))]
tcr_df[, cohort := factor(cohort, levels = c("PBMC_Pre", "TIL_Post"))]

# 2. Panel A：IC50 (log10) vs DAI (log2) 筛选象限图
#    竖虚线标注 500 nM 经典强/中结合门槛；横点线标注 DAI=1.0（突变/野生型亲和力相等线）
p1 <- ggplot(neo_df, aes(x = log10(mt_ic50_nM), y = log2(dai_ratio + 0.01), color = tier_pass, size = rna_tpm)) +
  geom_vline(xintercept = log10(500), linetype = "dashed", color = "#555555") +
  geom_hline(yintercept = log2(1.0 + 0.01), linetype = "dotted", color = "#888888") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("PASS" = "#d95f02", "FILTERED" = "#7570b3")) +
  labs(
    title = "A. Neoantigen Candidate Prioritization",
    subtitle = sprintf("PASS: %d/%d (MT IC50 <= 500 nM, Rank <= 2%%, TPM >= 1)",
                       sum(neo_df$tier_pass == "PASS"), nrow(neo_df)),
    x = "log10(Mutant MHC-I IC50 [nM])",
    y = "log2(DAI = WT_IC50 / MT_IC50)",
    color = "Filter Tier",
    size = "RNA TPM"
  ) +
  theme_bw(base_size = 11)

# 3. Panel B：PBMC_Pre vs TIL_Post 克隆性对比箱线+散点图（验证寡克隆扩增趋势）
p2 <- ggplot(tcr_df, aes(x = cohort, y = clonality, fill = cohort)) +
  geom_boxplot(width = 0.45, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, color = "#222222") +
  scale_fill_manual(values = c("PBMC_Pre" = "#2b8cbe", "TIL_Post" = "#e34a33")) +
  labs(
    title = "B. TCR Repertoire Clonality (TIL vs PBMC)",
    subtitle = "Higher clonality indicates antigen-driven T-cell expansion",
    x = "Sample Cohort",
    y = "Normalized Shannon Clonality (0-1)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# 4. 拼图导出 PNG（若未安装 patchwork 则优雅回退至单图输出）并归档 R sessionInfo()
out_png <- file.path(res_dir, "00_smoke_neoantigen_tcr_qc.png")
if (has_patchwork) {
  combined <- p1 + p2 + plot_layout(widths = c(1.25, 1))
  ggsave(out_png, combined, width = 11.5, height = 4.8, dpi = 160)
} else {
  ggsave(out_png, p1, width = 6.8, height = 4.8, dpi = 160)
}

cat("Saved QC figure:", out_png, "\n")

sink(file.path(log_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

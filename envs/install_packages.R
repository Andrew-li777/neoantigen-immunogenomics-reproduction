# ==============================================================================
# 新抗原与肿瘤免疫复现项目：R 环境最小依赖检查与安装脚本
#
# 用法:
#   Rscript envs/install_packages.R --dry-run            # 仅检查缺失包，不执行安装
#   Rscript envs/install_packages.R                      # 安装核心 CRAN 绘图与表格包
#   Rscript envs/install_packages.R --profile ressler    # 检查/补充 ressler2024_neobcc 扩展包
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

profile <- "core"
if ("--profile" %in% args) {
  idx <- which(args == "--profile")
  if (idx < length(args)) profile <- args[idx + 1]
} else if (any(grepl("^--profile=", args))) {
  profile <- sub("^--profile=", "", args[grep("^--profile=", args)][1])
}

# 核心 CRAN 包：表格读写 (data.table/dplyr)、发表级绘图 (ggplot2/patchwork/ggrepel/pheatmap) 与单细胞生态 (Seurat)
cran_core <- c(
  "BiocManager", "ggplot2", "data.table", "dplyr",
  "patchwork", "ggrepel", "pheatmap", "RColorBrewer", "Seurat"
)

# 核心 Bioconductor 基础容器包
bioc_core <- c(
  "SummarizedExperiment"
)

profile_extras <- switch(
  profile,
  "ressler" = list(
    cran = c("SoupX", "yaml", "rmarkdown", "knitr", "ggrepel"),
    bioc = c("scRepertoire", "fgsea", "AUCell")
  ),
  "neotcr" = list(
    cran = c("glmnet", "pROC", "survival", "survminer"),
    bioc = character(0)
  ),
  list(cran = character(0), bioc = character(0))
)

target_cran <- unique(c(cran_core, profile_extras$cran))
target_bioc <- unique(c(bioc_core, profile_extras$bioc))

installed <- rownames(installed.packages())
missing_cran <- setdiff(target_cran, installed)
missing_bioc <- setdiff(target_bioc, installed)

cat("=== 新抗原与肿瘤免疫 R 环境依赖检查 (Profile:", profile, ") ===\n")
cat("R 版本:", R.version.string, "\n")
cat("CRAN 目标:", length(target_cran), "个 | 已装:", length(target_cran) - length(missing_cran),
    "| 缺失:", length(missing_cran), "\n")
if (length(missing_cran) > 0) cat("  缺失 CRAN:", paste(missing_cran, collapse = ", "), "\n")

cat("Bioc 目标:", length(target_bioc), "个 | 已装:", length(target_bioc) - length(missing_bioc),
    "| 缺失:", length(missing_bioc), "\n")
if (length(missing_bioc) > 0) cat("  缺失 Bioc:", paste(missing_bioc, collapse = ", "), "\n")

if (dry_run) {
  cat("\n[--dry-run] 仅完成检查，未执行任何包安装。\n")
  quit(status = 0)
}

if (length(missing_cran) > 0) {
  cat("\n>>> 安装缺失 CRAN 包...\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (length(missing_bioc) > 0) {
  cat("\n>>> 安装缺失 Bioconductor 包...\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

cat("\n=== R 依赖检查/安装流程结束 ===\n")

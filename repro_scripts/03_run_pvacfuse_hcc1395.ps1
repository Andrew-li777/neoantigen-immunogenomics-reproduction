<#
.SYNOPSIS
    Track B-2 可复现调度脚本：在 neoantigen-pvactools:7.1.5 Docker 容器内运行 pvacfuse
    对 HCC1395 基因融合数据（AGFusion + STAR-Fusion）进行 MHC-I/II 融合新抗原预测，
    并调用 R 下游脚本输出摘要指标、Top 候选表与汇总可视化图。
#>

$ProjRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$OutDir = Join-Path $ProjRoot "data\processed\pvacfuse_hcc1395"
$ResDir = Join-Path $ProjRoot "repro_results\pvacfuse_hcc1395"
$LogFile = Join-Path $ProjRoot "logs\03_pvacfuse_hcc1395.log"
New-Item -ItemType Directory -Force -Path $OutDir, $ResDir | Out-Null

# [1/2] 启动轻量离线容器执行 pvacfuse run：
#   - 输入1: AGFusion 输出目录（含 20 个基因融合事件的 cDNA/CDS/蛋白序列注释）
#   - 样本名: HCC1395_TUMOR_RNA（作为输出文件名前缀）
#   - HLA 等位基因: 与 HCC1395 pvacseq 保持一致（OptiType HLA-I + DRB1*04:05）
#   - 预测引擎: 自带离线权重的 MHCnuggetsI (8/9-mer) 与 MHCnuggetsII (15-mer)
#   - --starfusion-file: 关联 Junction/Spanning reads 支持度与 FFPM 表达量，用于 Tier Pass 判定
Write-Host ">>> [1/2] Running pvacfuse inside neoantigen-pvactools:7.1.5..."
docker run --rm -v "${ProjRoot}:/workspace" -w /workspace neoantigen-pvactools:7.1.5 `
  pvacfuse run `
  /workspace/data/raw/HCC1395/HCC1395_inputs/agfusion_results `
  HCC1395_TUMOR_RNA `
  "HLA-A*29:02,HLA-B*45:01,HLA-B*82:02,HLA-C*06:02,DRB1*04:05" `
  MHCnuggetsI MHCnuggetsII `
  /workspace/data/processed/pvacfuse_hcc1395 `
  --starfusion-file /workspace/data/raw/HCC1395/HCC1395_inputs/star-fusion.fusion_predictions.tsv `
  -e1 8,9 -e2 15 `
  -t 4 2>&1 | Tee-Object -FilePath $LogFile

if ($LASTEXITCODE -ne 0) { throw "pvacfuse Docker run failed!" }

# [2/2] 调用本机 R 脚本解析 pvacfuse 聚合与过滤表，输出图表与摘要表
Write-Host ">>> [2/2] Running R downstream fusion neoantigen analysis..."
Rscript (Join-Path $PSScriptRoot "03_pvacfuse_hcc1395_analysis.R")
if ($LASTEXITCODE -ne 0) { throw "03_pvacfuse_hcc1395_analysis.R failed!" }

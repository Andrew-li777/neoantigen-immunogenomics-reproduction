<#
.SYNOPSIS
    Reproducible runner for Track B: griffithlab/pVACtools v7.1.5 on HCC1395_inputs
    Runs pvacseq inside neoantigen-pvactools:7.1.5 Docker container, exports
    container versions, and runs R downstream summary & plotting.
#>

$ProjRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$OutDir = Join-Path $ProjRoot "data\processed\pvacseq_hcc1395"
$ResDir = Join-Path $ProjRoot "repro_results\pvacseq_hcc1395"
$LogFile = Join-Path $ProjRoot "logs\02_pvacseq_hcc1395.log"
New-Item -ItemType Directory -Force -Path $OutDir, $ResDir | Out-Null

# [1/3] 启动 Track B 隔离容器执行 pvacseq run（体细胞 SNV/Indel 新抗原预测）：
#   - 输入1: VEP 注释且含 RNA 表达/VAF (GX/TX/RAD/RAF) 的 VCF (1,470 条变异)
#   - 样本名: HCC1395_TUMOR_DNA（需与 VCF Header 肿瘤列严格一致）
#   - HLA 列表: OptiType HLA-I (A*29:02, B*45:01, B*82:02, C*06:02) + HLA-II DRB1*04:05
#   - 预测方法: 自带离线权重的 MHCnuggetsI (-e1 8,9-mer) 与 MHCnuggetsII (-e2 15-mer)
#   - --normal-sample-name: 配对正常对照列 HCC1395_NORMAL_DNA，过滤胚系假阳性
#   - --phased-proximal-variants-vcf: 启用顺式定相近端变异氨基酸序列校正 (输出 20 条记录)
#   - -t 4: 容器内并行 4 线程加速预测
Write-Host ">>> [1/3] Running pvacseq inside neoantigen-pvactools:7.1.5..."
docker run --rm -v "${ProjRoot}:/workspace" -w /workspace neoantigen-pvactools:7.1.5 `
  pvacseq run `
  /workspace/data/raw/HCC1395/HCC1395_inputs/annotated.expression.vcf.gz `
  HCC1395_TUMOR_DNA `
  "HLA-A*29:02,HLA-B*45:01,HLA-B*82:02,HLA-C*06:02,DRB1*04:05" `
  MHCnuggetsI MHCnuggetsII `
  /workspace/data/processed/pvacseq_hcc1395 `
  -e1 8,9 -e2 15 `
  --normal-sample-name HCC1395_NORMAL_DNA `
  --phased-proximal-variants-vcf /workspace/data/raw/HCC1395/HCC1395_inputs/phased.vcf.gz `
  -t 4 2>&1 | Tee-Object -FilePath $LogFile

if ($LASTEXITCODE -ne 0) { throw "pvacseq Docker run failed!" }

# [2/3] 在同一容器镜像内导出精确 Python 包版本清单，保障环境可审计性
Write-Host ">>> [2/3] Exporting pVACtools container package versions..."
docker run --rm -v "${ProjRoot}:/workspace" -w /workspace neoantigen-pvactools:7.1.5 python -c "
import csv, importlib.metadata, sys
pkgs = ['pvactools', 'mhcnuggets', 'mhcflurry', 'tensorflow', 'vaxrank', 'varcode', 'isovar', 'vcfpy', 'pysam', 'numpy', 'pandas', 'pyarrow', 'polars']
with open('/workspace/repro_results/pvacseq_hcc1395/pvactools_container_versions.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['Package', 'Version'])
    w.writerow(['python', sys.version.split()[0]])
    for p in pkgs:
        w.writerow([p, importlib.metadata.version(p)])
"

# [3/3] 调用宿主机 R 脚本完成下游分层统计、DAI 计算与双联可视化绘图
Write-Host ">>> [3/3] Running R downstream prioritization summary & plots..."
Rscript (Join-Path $PSScriptRoot "02_pvacseq_hcc1395_analysis.R")
if ($LASTEXITCODE -ne 0) { throw "02_pvacseq_hcc1395_analysis.R failed!" }

<#
.SYNOPSIS
    Pre-flight check for Docker Desktop, WSL2, Nextflow, Python, and R readiness
    under Windows when local user WSL2 distro (vm) is unavailable.
#>

Write-Host "==================================================================="
Write-Host "  Neoantigen & Immunogenomics Pre-Flight Environment Check"
Write-Host "==================================================================="

# 1. Docker Engine check
Write-Host "`n[1/4] Checking Docker Desktop..."
try {
    $dockerVer = docker version --format "{{.Client.Version}} (Server: {{.Server.Version}}, OS: {{.Server.Os}}/{{.Server.Arch}})" 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  [OK] Docker: $dockerVer" -ForegroundColor Green
    } else {
        Write-Host "  [WARN] Docker daemon not reachable: $dockerVer" -ForegroundColor Yellow
    }
} catch {
    Write-Host "  [ERR] Docker CLI not found." -ForegroundColor Red
}

Write-Host "`n  Local cached Docker images:"
$imgs = docker images --format "  - {{.Repository}}:{{.Tag}} ({{.Size}})" 2>$null
if ($imgs) { $imgs } else { Write-Host "  (none)" }

# 2. WSL2 status
Write-Host "`n[2/4] Checking WSL2 distributions..."
$wslOut = wsl -l -v 2>$null
if ($wslOut) {
    $wslOut | ForEach-Object { Write-Host "  $_" }
} else {
    Write-Host "  [INFO] No WSL2 output returned."
}

# 3. Java & Nextflow native check
Write-Host "`n[3/4] Checking native Java / Nextflow on Windows PATH..."
$hasJava = Get-Command java -ErrorAction SilentlyContinue
$hasNf = Get-Command nextflow -ErrorAction SilentlyContinue
Write-Host ("  java     : " + ($(if ($hasJava) { $hasJava.Source } else { "Not installed on Windows PATH" })))
Write-Host ("  nextflow : " + ($(if ($hasNf) { $hasNf.Source } else { "Not installed on Windows PATH" })))

# 4. 结合当前环境（无可用用户 WSL2 vm）输出已落地的双容器推荐策略
Write-Host "`n[4/4] Recommended Execution Strategy (No active user WSL2 VM):"
Write-Host "  - Track A (scirpy / TCR scRNA-seq):"
Write-Host "      Use isolated Linux container neoantigen-scirpy:py313 (envs/Dockerfile.scirpy)"
Write-Host "      or an isolated .venv with envs/requirements.txt."
Write-Host "  - Track B (Neoantigen Prediction: pvacseq + pvacfuse):"
Write-Host "      Use lightweight offline container neoantigen-pvactools:7.1.5 (envs/Dockerfile.pvactools)"
Write-Host "      via repro_scripts/02_run_pvacseq_hcc1395.ps1 & 03_run_pvacfuse_hcc1395.ps1."
Write-Host "  - Optional Nextflow DSL2 (nf-core/epitopeprediction):"
Write-Host "      Install a clean Ubuntu WSL2 distro (wsl --install -d Ubuntu) OR run Nextflow"
Write-Host "      inside a Linux Docker container with matching host/container workdir mount."
Write-Host "==================================================================="

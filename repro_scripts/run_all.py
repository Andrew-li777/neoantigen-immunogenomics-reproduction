#!/usr/bin/env python3
"""
run_all.py
==========
Unified automation orchestrator & output integrity checker for
`reproduction/neoantigen_immunogenomics`.

Supported modes:
  - smoke      : Run synthetic REPRO-D smoke pipeline (00_smoke_test.py + 00_smoke_test.R)
  - track-a    : Run Track A inside neoantigen-scirpy:py313 Docker container
  - track-b    : Run Track B-1 (pvacseq) & B-2 (pvacfuse) inside neoantigen-pvactools:7.1.5 + R
  - downstream : Re-run only host R downstream analyses (02_pvacseq + 03_pvacfuse) on existing TSVs
  - all        : Execute smoke + track-a + track-b end-to-end and verify all deliverables

Usage:
  python reproduction/neoantigen_immunogenomics/repro_scripts/run_all.py --mode smoke
  python reproduction/neoantigen_immunogenomics/repro_scripts/run_all.py --mode downstream
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import pathlib
import subprocess
import sys
import time
from typing import Dict, List


def get_project_root() -> pathlib.Path:
    """定位项目根目录 (reproduction/neoantigen_immunogenomics)。"""
    return pathlib.Path(__file__).resolve().parent.parent


def run_step(cmd: List[str], cwd: pathlib.Path, log_lines: List[str]) -> float:
    """执行子进程步骤，实时回显标准输出并收集耗时与退出状态至运行日志。"""
    t0 = time.perf_counter()
    print(f">>> Running: {' '.join(cmd)}")
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    elapsed = round(time.perf_counter() - t0, 2)
    log_lines.append(f"=== CMD: {' '.join(cmd)} (exit={proc.returncode}, {elapsed}s) ===")
    if proc.stdout:
        print(proc.stdout.rstrip())
        log_lines.append(proc.stdout.rstrip())
    if proc.stderr:
        print(proc.stderr.rstrip(), file=sys.stderr)
        log_lines.append("[STDERR]\n" + proc.stderr.rstrip())
    if proc.returncode != 0:
        raise RuntimeError(f"Step failed ({proc.returncode}): {' '.join(cmd)}")
    print(f"<<< Done in {elapsed}s\n")
    return elapsed


def collect_environment_versions(root: pathlib.Path) -> pathlib.Path:
    """采集宿主机 Python、R（含 CRAN/Bioc 核心包）与 Docker Desktop 精确版本至 TSV。"""
    py_pkgs = [
        "numpy", "scipy", "pandas", "matplotlib", "seaborn",
        "biopython", "anndata", "scanpy", "mudata", "scirpy", "pytest",
    ]
    rows: List[Dict[str, str]] = [
        {
            "Component": "Python",
            "Version": sys.version.split()[0],
            "Source": "Runtime",
        }
    ]
    for pkg in py_pkgs:
        try:
            ver = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            ver = "Not Installed"
        rows.append({"Component": f"py:{pkg}", "Version": ver, "Source": "PyPI/venv"})

    # 采集宿主机 R 版本及实际使用的下游统计绘图包（含 ggrepel 0.9.8）
    r_expr = (
        "pkgs <- c('ggplot2','data.table','dplyr','patchwork','ggrepel','Seurat','BiocManager','pheatmap');"
        "cat('R\\t', paste0(R.version$major,'.',R.version$minor), '\\tBase R\\n', sep='');"
        "for (p in pkgs) {"
        "  v <- tryCatch(as.character(packageVersion(p)), error=function(e) 'Not Installed');"
        "  cat('R:', p, '\\t', v, '\\tCRAN/Bioc\\n', sep='')"
        "}"
    )
    r_proc = subprocess.run(
        ["Rscript", "-e", r_expr],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if r_proc.returncode == 0:
        for line in r_proc.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 3:
                rows.append({"Component": parts[0], "Version": parts[1], "Source": parts[2]})

    # 采集宿主机 Docker Desktop 客户端版本号
    d_proc = subprocess.run(
        ["docker", "version", "--format", "{{.Client.Version}}"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    docker_ver = d_proc.stdout.strip() if d_proc.returncode == 0 else "Unavailable"
    rows.append({"Component": "Docker", "Version": docker_ver, "Source": "Docker Desktop"})

    out_tsv = root / "repro_results" / "environment_versions.tsv"
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Component", "Version", "Source"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return out_tsv


def verify_expected_outputs(root: pathlib.Path, mode: str) -> None:
    """对已执行阶段的关键 TSV 摘要表与 PNG 图表执行存在性及非空断言校验。"""
    res = root / "repro_results"
    required: List[pathlib.Path] = []
    if mode in ("smoke", "all"):
        required += [
            res / "00_smoke_verification_summary.tsv",
            res / "00_smoke_neoantigen_tcr_qc.png",
        ]
    if mode in ("track-a", "all"):
        required += [
            res / "scirpy_wu2020" / "01_scirpy_wu2020_summary.tsv",
            res / "scirpy_wu2020" / "01_scirpy_cd8_teff_vs_trm_imbalance.tsv",
            res / "scirpy_wu2020" / "07_scirpy_expansion_and_antigen_umap.png",
        ]
    if mode in ("track-b", "downstream", "all"):
        required += [
            res / "pvacseq_hcc1395" / "02_pvacseq_hcc1395_summary.tsv",
            res / "pvacseq_hcc1395" / "02_pvacseq_top_candidates.tsv",
            res / "pvacseq_hcc1395" / "01_pvacseq_tier_and_allele_distribution.png",
            res / "pvacseq_hcc1395" / "02_pvacseq_ic50_vs_expression_and_dai.png",
            res / "pvacfuse_hcc1395" / "03_pvacfuse_hcc1395_summary.tsv",
            res / "pvacfuse_hcc1395" / "03_pvacfuse_top_candidates.tsv",
            res / "pvacfuse_hcc1395" / "03_pvacfuse_fusion_neoepitopes.png",
        ]
    for path in required:
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"Expected output missing or empty: {path}")
    print(f">>> Verified {len(required)} output artifact(s) for mode={mode} [OK]")


def main() -> int:
    """解析运行模式并按序调度冒烟测试、Docker 双轨流程或下游 R 分析。"""
    parser = argparse.ArgumentParser(description="Neoantigen & Immunogenomics Reproduction Runner")
    parser.add_argument(
        "--mode",
        choices=["smoke", "track-a", "track-b", "downstream", "all"],
        default="smoke",
        help="Execution mode (default: smoke)",
    )
    args = parser.parse_args()

    root = get_project_root()
    log_dir = root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_lines: List[str] = []

    print("===================================================================")
    print("  Neoantigen & Immunogenomics Reproducible Workflow Orchestrator")
    print(f"  Workspace: {root}")
    print(f"  Mode     : {args.mode.upper()}")
    print("===================================================================\n")

    if args.mode in ("smoke", "all"):
        run_step([sys.executable, "repro_scripts/00_smoke_test.py"], root, log_lines)
        run_step(["Rscript", "repro_scripts/00_smoke_test.R"], root, log_lines)

    if args.mode in ("track-a", "all"):
        run_step(
            [
                "docker", "run", "--rm",
                "-v", f"{root}:/workspace", "-w", "/workspace",
                "neoantigen-scirpy:py313",
                "python", "/workspace/repro_scripts/01_scirpy_wu2020_tcr.py",
            ],
            root,
            log_lines,
        )

    if args.mode in ("track-b", "all"):
        run_step(
            ["powershell", "-ExecutionPolicy", "Bypass", "-File", "repro_scripts/02_run_pvacseq_hcc1395.ps1"],
            root,
            log_lines,
        )
        run_step(
            ["powershell", "-ExecutionPolicy", "Bypass", "-File", "repro_scripts/03_run_pvacfuse_hcc1395.ps1"],
            root,
            log_lines,
        )

    if args.mode == "downstream":
        run_step(["Rscript", "repro_scripts/02_pvacseq_hcc1395_analysis.R"], root, log_lines)
        run_step(["Rscript", "repro_scripts/03_pvacfuse_hcc1395_analysis.R"], root, log_lines)

    verify_expected_outputs(root, args.mode)
    env_tsv = collect_environment_versions(root)
    print(f"Environment snapshot exported: {env_tsv}")

    log_name = "smoke_run.log" if args.mode == "smoke" else f"orchestrator_{args.mode}.log"
    log_file = log_dir / log_name
    log_file.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"Run log written: {log_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

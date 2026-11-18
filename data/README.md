# 新抗原与肿瘤免疫复现项目数据字典 (`data/`)

> 本目录遵循 `REPRO-D` 工程规范：原始测序 FASTQ、比对 BAM、原始 VCF/MAF、参考基因组与外部下载测试包严格只读存放于 `data/raw/`（受 `.gitignore` 屏蔽，不提交至 Git）；清洗后的新抗原候选表、样本表（`samplesheet.csv`）及 `AnnData`/`MuData`/`RDS` 中间对象存放于 `data/processed/`。

---

## 目录结构约定

```tree
data/
├── raw/            # 原始测序 / VCF / 融合 / V(D)J 对象只读区（详单见 data/raw/README.md）
├── processed/      # 清洗后的新抗原优先级表、TCR 克隆表、.h5mu 中间对象（详单见 data/processed/README.md）
└── README.md       # 本数据总台账（登记来源、样本规模、HLA 分型与体积）
```

> 分级细则文档：[**`data/raw/README.md`**](./raw/README.md) ｜ [**`data/processed/README.md`**](./processed/README.md)

---

## 数据登记表

| 数据集名称 | 来源平台 / Accession / 仓库入口 | 样本规模与组学模态 | 关键文件格式 | 磁盘占用 | 存放路径 |
|---|---|---|---|---|---|
| `smoke_synthetic_neo_tcr` | 本地 `00_smoke_test.py` 确定性生成（`seed=20261115`） | 24 条体细胞候选新抗原肽 + 6 例样本 TCR 克隆型表 | `.tsv` | < 10 KB | `data/processed/smoke_*.tsv` |
| `epitopeprediction_test` *(未拉取 / 架构对照)* | `nf-core/epitopeprediction` `-profile test` | GRCh38 chr4/chr19 子集（4 个变异；本工作区因无可用 WSL2 `vm` 改由 `pVACtools` 容器直接验证） | `.vcf` / `.fa.gz` / `.tar.gz` | ~15 MB *(未下载)* | *未拉取* |
| `scirpy_wu2020_3k` *(已实跑)* | `scverse/scirpy` (`ir.datasets.wu2020_3k()`) | 3,000 个泛癌种（肺/肠/子宫内膜/肾）肿瘤/癌旁/外周血 T 细胞 scRNA-seq + V(D)J TCR（质控后保留 1,845 个完整配对细胞） | `.h5mu` (`MuData`) | 原始 `16.48 MB` / 处理后 `27.75 MB` | `data/raw/scirpy_cache/` & `data/processed/wu2020_3k_scirpy_processed.h5mu` |
| `pVACtools_HCC1395` *(已实跑 Seq + Fuse)* | `griffithlab/pVACtools` (`HCC1395_inputs.zip`, `pvactools v7.1.5`) | HCC1395 乳腺癌配对 WES+RNA-Seq（1,470 变异 + 20 融合；HLA：A\*29:02, B\*45:01, B\*82:02, C\*06:02, DRB1\*04:05；`pvacseq` 357 变异/52 过滤表位；`pvacfuse` 31 融合异构体/14 过滤表位） | `.vcf.gz` / `.fa.gz` / `.tsv` | 原始 `24.22 MiB`（zip）/ 产物见 `processed/README.md` | `data/raw/HCC1395/HCC1395_inputs/` & `data/processed/pvacseq_hcc1395/` & `data/processed/pvacfuse_hcc1395/` |

# 原始只读输入数据清单 (`data/raw/`)

> **管理规范**：本目录严格保持**只读**（`REPRO-D` 规范），所有原始体细胞 VCF、定相参考、参考蛋白组、融合预测结果及缓存数据集均受工作区 [`.gitignore`](../../.gitignore) 屏蔽，不提交大文件至 Git；仅本文档随仓记录全部原始文件来源、规模与字段契约。

---

## 一、`HCC1395/` —— 乳腺癌细胞系多组学配对基准包（Track B）

- **来源**：[griffithlab/pVACtools](https://github.com/griffithlab/pVACtools) 官方基准测试包 `HCC1395_inputs.zip`（下载包体积 `25,399,820` bytes ≈ `24.22 MiB`）
- **生物学背景**：经典三阴性乳腺癌细胞系 `HCC1395`（肿瘤）及其配对正常淋巴母细胞系 `HCC1395BL`（Normal）的配对 WES + RNA-Seq 多组学注释输入
- **存放路径**：`data/raw/HCC1395/HCC1395_inputs/`

| 文件 / 目录 | 实测字节大小 | 用途与关键数据特征 | 被哪个脚本使用 |
|---|---:|---|---|
| `annotated.expression.vcf.gz` (+ `.tbi`) | `2,050,827` B (`+21,977` B) | 经 Ensembl VEP（`CSQ`）注释并整合 RNA-Seq 基因/转录本表达量（`GX`/`TX`）与等位基因读数支持（`RAD`/`RAF`）的体细胞 VCF；共 **`1,470`** 条变异，含配对样本列 `HCC1395_NORMAL_DNA` 与 `HCC1395_TUMOR_DNA` | [02_run_pvacseq_hcc1395.ps1](../../repro_scripts/02_run_pvacseq_hcc1395.ps1) |
| `phased.vcf.gz` (+ `.tbi`) | `8,024,980` B (`+45,374` B) | 同染色单体顺式定相 VCF（包含体细胞突变与邻近胚系杂合 SNP），供 `pvacseq --phased-proximal-variants-vcf` 修正长肽翻译序列 | [02_run_pvacseq_hcc1395.ps1](../../repro_scripts/02_run_pvacseq_hcc1395.ps1) |
| `optitype_normal_result.tsv` | `111` B | `OptiType` 对正常样本实测的 HLA-I 四位分型：`A*29:02`、`B*45:01`、`B*82:02`、`C*06:02`（共 1,748 支持 reads；另有 Class II 常用参考位点 `DRB1*04:05`） | Track B `pvacseq` & `pvacfuse` HLA 参数输入 |
| `star-fusion.fusion_predictions.tsv` | `46,172` B | `STAR-Fusion` 预测融合基因汇总表（`1` 行表头 + **`20`** 条融合事件），记录 `JunctionReadCount`、`SpanningFragCount`、断点坐标与标准化表达量 `FFPM`（例如 `EIF3K--CYP39A1`、`TMEM104--RAB37`） | [03_run_pvacfuse_hcc1395.ps1](../../repro_scripts/03_run_pvacfuse_hcc1395.ps1) |
| `agfusion_results/` | 20 个子目录 | `AGFusion` 对上述 20 个融合事件生成的剪切后融合 cDNA（`*_cdna.fa`）、CDS（`*_cds.fa`）、蛋白（`*_protein.fa`）、外显子/域注释表及结构示意图 | [03_run_pvacfuse_hcc1395.ps1](../../repro_scripts/03_run_pvacfuse_hcc1395.ps1) |
| `Homo_sapiens.GRCh38.pep.all.fa.gz` | `14,376,754` B | Ensembl GRCh38 人类全蛋白组参考 FASTA，可用于参考蛋白组相似性比对过滤（`--peptide-fasta`） | 扩展参考库 |
| `.gz` | `28` B | 上游 `HCC1395_inputs.zip` 解压自带空占位文件（可忽略，未参与计算） | 无 |

---

## 二、`scirpy_cache/` —— 单细胞转录组与 V(D)J TCR 基准数据（Track A）

- **来源**：`scverse/scirpy` 官方内置接口 `ir.datasets.wu2020_3k()`（通过 `pooch` 自动拉取并重定向缓存至此目录，源自 Wu et al. *Nature* 2020 泛癌种（肺、结直肠、子宫内膜、肾癌肿瘤/癌旁/外周血）T 细胞数据集子集）
- **存放路径**：`data/raw/scirpy_cache/`

| 文件 | 体积 | 数据模态与关键维度 | 被哪个脚本使用 |
|---|---:|---|---|
| `scirpy/0.26.1/wu2020_3k.h5mu` | `16.48 MB` | 标准 `MuData` 多模态对象：<br>• `gex`（转录组）：`3,000` 个 T 细胞 × `30,727` 个基因，含 `X_umap_orig`、`sample`、`source`（tumor / normal_adjacent / blood）、`patient`、`cluster_orig`<br>• `airr`（V(D)J TCR）：`3,000` 个细胞符合 AIRR 标准的 α/β 链 `awkward` 结构矩阵 | [01_scirpy_wu2020_tcr.py](../../repro_scripts/01_scirpy_wu2020_tcr.py) |
| `matplotlib/fontlist-v3.11.0.json` | 缓存索引 | 容器内首次启动 Matplotlib 渲染生成的字体检索缓存（重定向 `XDG_CACHE_HOME` 避免污染宿主机） | Track A 容器运行时 |

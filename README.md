# Sepsis scRNA-seq Analysis Scripts

> **单细胞 RNA 测序分析脚本整理版**：基于上传的五份 R 脚本，分别整理了数据质控/注释、CellChat 细胞通讯、Monocle 2 拟时序、scTenifoldKnk 虚拟敲除与 DoRothEA/VIPER 转录因子活性分析。
>
> **English:** Curated R scripts for a sepsis-related scRNA-seq workflow. Scripts are modular and require user-provided input data and dataset-specific metadata.

## 📁 Project structure / 项目结构

```text
.
├── README.md
├── .gitignore
├── config/
│   ├── cluster_annotations.csv          # 原始人工注释映射（需自行审核）
│   └── sample_metadata.example.csv      # 样本分组示例
├── data/
│   └── README.md                        # 输入文件说明；数据不会提交到 Git
├── scripts/
│   ├── main.R                           # QC / Harmony / clustering / annotation
│   ├── cell.R                           # CellChat cell-cell communication
│   ├── monocle_main.R                   # Monocle 2 (DDRTree) pseudotime
│   ├── scTenifoldKnk.R                  # In-silico single-gene knockout
│   └── tf.R                             # DoRothEA / VIPER TF activity
├── tests/
│   └── validate_repository.py            # 辅助静态检查，不替代 R/数据测试
├── docs/
│   ├── main.md
│   ├── cell.md
│   ├── monocle_main.md
│   ├── scTenifoldKnk.md
│   └── tf.md
└── results/                             # 运行结果（已 gitignore）
```

## 🔬 Module overview / 脚本说明

| Script | Function | Required input | Main output |
| --- | --- | --- | --- |
| [`main.R`](scripts/main.R) | 10X 读取、质控、标准化、Harmony、聚类与人工注释 | `data/output_folder/<sample>/` | `data/af.rds`, `results/main/` |
| [`cell.R`](scripts/cell.R) | CellChat 配体–受体与信号通路分析 | `data/af_Activated_Macrophages.rds`，`Celltype` 列 | `results/cellchat/` |
| [`monocle_main.R`](scripts/monocle_main.R) | Monocle **2** / DDRTree 拟时序与基因动态 | 同上，需原始 `RNA` counts | `results/monocle2/` |
| [`scTenifoldKnk.R`](scripts/scTenifoldKnk.R) | scTenifoldKnk 单基因虚拟敲除与显著扰动基因 | 同上，含目标细胞类型 | `results/scTenifoldKnk/` |
| [`tf.R`](scripts/tf.R) | DoRothEA ABC regulon + VIPER 转录因子活性 | `data/af1_Figure7_FINAL_reannotated_noNA.rds`，`Celltype_revised` 列 | `results/dorothea/` |

每个模块的参数和注意事项见 [`docs/`](docs/main.md)。

## 🧰 Software and dependencies / 环境依赖

建议使用 **R 4.x + Seurat v5**；需要根据各包与 R / Bioconductor 版本的兼容性安装。以下是安装入口示例（不是锁定版本的环境）：

```r
install.packages(c(
  "Seurat", "SeuratObject", "ggplot2", "dplyr", "patchwork",
  "harmony", "igraph", "ggraph", "ggrepel", "pheatmap",
  "remotes", "clustree"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("monocle", "Biobase", "dorothea", "ComplexHeatmap"))

# CellChat v1/v2 API used in scripts/cell.R:
remotes::install_github("jinworks/CellChat")

# scTenifoldKnk (confirm installed version's API):
install.packages("scTenifoldKnk")
```

以上不保证所有操作系统和包版本均可直接安装；安装失败时请查看对应官方文档。CellChat 的新一代 Spatial CellChat v3 **不是**这个脚本的直接替代品。

## ▶️ How to run / 运行方法

1. 在仓库根目录工作，准备 `data/output_folder/<sample>/` 中的 Cell Ranger 10X 文件（矩阵、基因、barcode），或对应 `.rds` 文件。数据文件不在 Git 仓库中。
2. 审核 [`config/cluster_annotations.csv`](config/cluster_annotations.csv) 的原始手工标签：它仅对应上传代码的 cluster 0–24，不保证适用于新数据。
3. 如需 **HC/Sepsis** 组别比例图，将 `config/sample_metadata.example.csv` 复制为 `config/sample_metadata.csv`，填写真实**样本文件夹名**与对应 `condition`（否则该图自动跳过）。
4. 依需求执行脚本，例如：

```bash
Rscript scripts/main.R
Rscript scripts/cell.R
Rscript scripts/monocle_main.R
Rscript scripts/scTenifoldKnk.R
Rscript scripts/tf.R
```

**重要：这五份脚本不是无缝串行管线。** `main.R` 保存 `data/af.rds`；原始其他三个模块依赖额外准备的 `af_Activated_Macrophages.rds`，`tf.R` 依赖另一个重新注释的 RDS。若希望构建完整串行流程，须先定义这些对象是如何筛选和生成的，不能直接把 `af.rds` 当成它们。

## ⚙️ Dataset-specific settings / 需要自行检查的参数

- **质控与聚类**（`main.R`）：保留原始范围 `nCount_RNA = 1000–40000`，`nFeature_RNA = 200–5000`，`percent.MT ≤ 25`，`n_pcs = 9`，`resolution = 1.0`。阈值应按项目重新评估。
- **样本标签**（`main.R`）：不再假设“文件夹排序后前两组是 HC”，改用 `sample_metadata.csv` 显式映射。
- **细胞身份**（`main.R`）：人工注释映射仅为原始代码的 0–24 簇；若聚类结果不同，脚本将明确停止并提示修改。
- **物种**：当前 CellChat DB 与 DoRothEA regulon 为 **human**，线粒体基因前缀使用 `^MT-`。
- **Monocle**：原始脚本实际使用 `monocle` / DDRTree（Monocle 2），并非 `monocle3`；可设置 `root_state` 指定具有生物学依据的轨迹起点。
- **虚拟敲除**：原始注释写“五个”，但真正启用的是 `KIF1B`, `CTSO`, `BLOC1S1` **三个**。脚本按实际三基因运行，可自行加入 `SORT1`, `M6PR`；显著性规则保持 `|Z| ≥ 2` 且 adjusted `P < 0.01`。
- **计算资源**：VIPER 需要把表达矩阵转成稠密矩阵；已加内存阈值保护。scTenifoldKnk 运算时间和内存消耗较高。

## 🛠️ Changes from the original files / 整理内容

- 统一 UTF-8 注释、代码分区、路径命名、`results/` 输出与 `sessionInfo()` 记录。
- 删除机器专用的 `G:\...` 绝对路径、重复包加载与环境清空命令，不自动安装依赖。
- 修正原始可确认的问题：`percent.MT`/`percent.mt` 不一致、未定义的 `de`、UMAP 的对象/细胞身份混用、Monocle 输出覆盖与拟时序基因排序、DoRothEA 数据加载方式等。
- 为 scTenifoldKnk 新旧版本的 QC 参数名提供兼容分支（具体见 [`docs/scTenifoldKnk.md`](docs/scTenifoldKnk.md)）。
- 尽量保留原有分析目的、主要阈值、候选基因与图表；部分探索性或未执行注释段落没有照搬。

## ✅ Reproducibility and limitations / 复现声明

仓库只包含脚本、文档和示例配置，**没有提供原始 10X 数据或 RDS 对象**，因此不能凭仓库直接重复出原项目的图或统计结果。已检查文本结构、括号配对、配置数据及部分边界条件；仍需在具备实际数据和相应 R 包的环境中进行完整语法解析和运行验证。人工细胞类型、分组和伪时间方向都需要独立的生物学审查。

无需额外 Python 库的辅助检查（**不是 R 解析器**）：

```bash
python3 tests/validate_repository.py
```

请不要未经授权提交临床数据、样本标识或包含敏感信息的分析结果。`.gitignore` 已默认排除 `data/` 与 `results/` 的数据内容；新建文件上传前仍应检查 `git status`。

## 📚 Official references / 官方参考

- [SeuratObject: LayerData / Layers](https://satijalab.github.io/seurat-object/reference/Layers.html)
- [SeuratObject: JoinLayers](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
- [CellChat (v1/v2)](https://github.com/jinworks/CellChat)
- [Monocle 2 (Bioconductor)](https://bioconductor.org/packages/monocle)
- [scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk)
- [DoRothEA (Bioconductor)](https://bioconductor.org/packages/dorothea)

### GitHub upload / 上传建议

新建空仓库后，在本目录执行：

```bash
git init
git add README.md .gitignore scripts docs tests config data/README.md results/.gitkeep
git commit -m "Organize sepsis scRNA-seq R analysis scripts"
git branch -M main
git remote add origin https://github.com/<your-username>/<your-repository>.git
git push -u origin main
```

请将 `<your-username>` 和 `<your-repository>` 换为你的 GitHub 账户/仓库名，并先确认没有不应公开的数据。

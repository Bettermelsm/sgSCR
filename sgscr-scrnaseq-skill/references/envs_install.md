# 隔离环境可复现安装说明

配套脚本：`scripts/common/setup_envs.sh`（远程 Linux 服务器上运行，mamba 优先，conda 回退）。

## 四个环境一览

| 环境 | 用途 | 主要依赖 | 建立方式 |
|---|---|---|---|
| `sgscr_main` | S0–S6, S10, S10b, S11（Seurat 主环境） | R ≥ 4.2, Seurat v5, SeuratObject, SingleR, clustree, yaml, jsonlite, ggplot2, hedgehog, testthat | conda-forge 建 R 基础 env + CRAN/Bioc 补装 |
| `sgscr_cellchat` | S7 细胞通讯 | CellChat（BiocManager；conda-forge 无 r-cellchat，2026-10 实测 solve 失败） | 基础 env + BiocManager |
| `sgscr_trajectory` | S8 拟时序 | slingshot（+ 可选 monocle3） | conda |
| `sgscr_scenic` | S9 调控子 | pyscenic CLI（grn/ctx/aucell 三步） | pip/conda (python) |

环境名以用户 Config `envs.*` 为准；脚本默认上述名字，可用
`bash setup_envs.sh --cellchat <name> --trajectory <name> --scenic <name>` 覆盖。

## 快速开始

```bash
# 服务器上，一次性创建全部三个隔离环境（main 缺失时一并创建）
bash scripts/common/setup_envs.sh

# 只补装主环境 R 包
bash scripts/common/setup_envs.sh --main-pkgs-only
```

安装完成后 S0 环境盘点会校验四个 env 与关键包，缺项写进 Quality_Gate FAIL detail。

## 可选：S2b 环境 RNA 校正后端（默认不需要）

S2b **默认关闭**，不开就不需要任何额外包。开启（Config `ambient.correction.enabled: true`）需在
`sgscr_main` 中安装以下任一后端：

- **decontX**（推荐，来自 celda 包）：
  ```R
  # BiocManager::install("celda")   # 提供 decontX
  ```
- 脚本行为约定：开启但两个后端都缺失时 **显式失败并报告缺失包名**（quit status 1），
  绝不静默跳过、绝不伪造校正结果；校正产物写 `rds/s2b_obj_corrected.rds`（assay `decontX`），不覆盖 S2 产物。

## 版本锚点

不锁死小版本；以 Manifest `package_versions` 记录每次运行实际版本。
已知坑：SeuratObject 5.2.0 上 v3 Assay 的 JoinLayers 无方法（见 anti_patterns.md）。

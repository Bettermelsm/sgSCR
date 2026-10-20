---
name: sgscr-scrnaseq-skill
description: >-
  sgSCR scRNA-seq 标准化分析流程 — 10X 数据导入与批次整合、QC、归一化/HVG、PCA（定量选 PC 数）、
  多分辨率聚类 + clustree + ARI 辅助分辨率选择、SingleR 两层注释、Compartment 亚聚类，
  以及 CellChat 细胞通讯 / slingshot 拟时序 / pySCENIC 转录因子调控三个隔离环境模块。
  当用户提出 scRNA-seq / 单细胞测序 / Seurat / 细胞聚类注释 / CellChat / 细胞通讯 /
  拟时序 / Monocle3 / slingshot / SCENIC / regulon / 10X 三联文件相关请求时触发。
license: MIT
compatibility: "Linux（远程服务器）；R 主环境 Seurat v5；conda/mamba 隔离环境 ×3"
metadata:
  version: "0.2.0-beta"
  status: "BETA"
  owner: "sea"
  created: "2026-10-02"
  last_updated: "2026-10-20"
  last_review: ""
  risk_level: "medium"
  visibility: "public"
  triggers: "scRNA-seq, 单细胞, Seurat, 聚类注释, SingleR, CellChat, 细胞通讯, 拟时序, Monocle3, slingshot, SCENIC, pySCENIC, regulon, 10X, UMAP"
  platforms: "claude-code, kiro"
  locales: "zh-CN"
  dependencies: "yaml, jsonlite; 远程: Seurat, SingleR, clustree, CellChat, slingshot, pySCENIC"
  requires_user_action: "true"
  estimated_tokens: ""
  risk_notes: "人工确认点：QC 阈值 / 聚类分辨率 / cluster→主群映射表"
---
# sgSCR — scRNA-seq 标准化分析流程 Skill

## 适用场景（触发条件）

任一关键词或情境出现即触发：
- 用户提供/讨论 10X scRNA-seq 数据（barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz 三联文件）需要导入、整合、聚类、注释
- 需要从单细胞对象出 UMAP、marker 基因表、细胞类型注释（SingleR）或主群收敛合并
- 用户询问细胞通讯（CellChat / macrophage→fibroblast / TGFb / CCL / MIF 通路）
- 用户询问拟时序（trajectory / Monocle3 / slingshot / pseudotime）
- 用户询问转录因子调控（SCENIC / pySCENIC / regulon）
- 需要对既有 sgSCR 分析做断点续跑（`--from sN`）或检查 manifest 溯源记录

不适用：bulk RNA-seq（走 sgFLW）、空间转录组（暂未覆盖）、_SMART-seq2 全长数据（未验证）。

## 核心结论（先看这个）

> **流程 = S0–S11 十二个独立 Step，每个 Step 输入/输出落盘为文件， Quality_Gate 三态（PASS/FAIL/NEEDS_HUMAN_REVIEW）**。
> 三个新模块（CellChat/拟时序/SCENIC）各跑在独立 conda 环境，与 Seurat 主环境靠 rds/csv 文件交接，互不共享内存状态。
> 所有参数集中在一份 YAML Config；每次运行的种子/包版本/参数快照/耗时记入 JSON Manifest。

## 平台基线

- **计算端**：远程 Linux 服务器（数据留在服务器本地，本机不搬运大数据）
- **连接**：优先 ssh 直连（Config `remote.host_alias`），不可用走 frp 隧道（Config `remote.frp_fallback`）
- **环境**：R 主环境（S0–S6, S10, S11）+ 三个隔离 conda env：`sgscr_cellchat`(S7)、`sgscr_trajectory`(S8)、`sgscr_scenic`(S9)
- **配置**：`config/config.example.yaml` → 用户填成 `config.yaml`；`<TBD, 用户确认>` 占位符未填完不允许运行

## 标准工作流（60 秒上手）

```bash
# 1. 同步脚本与配置到远程
bash scripts/common/remote_exec.sh sync --config config/config.yaml

# 2. S0 环境与数据盘点（自动，无人工确认）
bash scripts/common/run_step.sh s0 --config config/config.yaml

# 3. 依次执行 S1–S4（S2/S4 有 NEEDS_HUMAN_REVIEW 暂停点，复查图后 --from 续跑）
bash scripts/common/run_step.sh s2 --config config/config.yaml   # → 复查 QC 小提琴图
bash scripts/common/run_step.sh --from s3 --config config/config.yaml

# 4. S5 两层注释 → 人工填写 cluster_mapping_confirmed.csv → 续跑 S6–S11
bash scripts/common/run_step.sh --from s6 --config config/config.yaml
```

## Step 速查表

| Step | 内容 | 环境 | 人工确认 |
|---|---|---|---|
| S0 | 环境盘点 + 三联文件校验 | main | 无 |
| S1 | Read10X 批量导入 + FindIntegrationAnchors/IntegrateData | main | 无 |
| S2 | QC（percent.mt/percent.HB 过滤） | main | ✅ QC 小提琴图 |
| S3 | 归一化/HVG/PCA + 定量选 PC 数（累计方差≥阈值 且 单PC<阈值 + 拐点） | main | 无 |
| S4 | 多分辨率聚类 + clustree + ARI 辅助选分辨率 | main | ✅ clustree 图 |
| S5 | SingleR 自动注释 + FindAllMarkers → 人工合并主群映射表 | main | ✅ 映射表 |
| S6 | Compartment 亚聚类（myeloid/fibroblast） | main | 无 |
| S6b | 亚群注释与保留过滤（标签 CSV + keep 语义 yes/sens/no，main/strict 双模式） | main | 无 |
| S7 | CellChat 细胞通讯（NC vs UC 比较） | sgscr_cellchat | 无 |
| S8 | 拟时序（默认 slingshot） | sgscr_trajectory | 无 |
| S9 | SCENIC（默认 pySCENIC，grn/ctx/aucell 三步） | sgscr_scenic | 无 |
| S10 | 目标基因追踪（SLC26A3/GUCA2A） | main | 无 |
| S11 | 出版级出图 + Manifest 最终汇总 | main | 无 |

## 决策矩阵

| 情境 | 动作 |
|---|---|
| ssh 直连不通 | 自动走 `remote.frp_fallback` 隧道；两者皆不通则终止并分别报告失败详情 |
| Species 与 SingleR 参考集物种不一致 | S5 终止，先改 Config（历史上曾因人数据误用小鼠参考集踩坑） |
| 某步 FAIL | 看该步 Manifest 条目的 `quality_gate.detail` 与原始报错，修完 `--from sN` 续跑 |
| NEEDS_HUMAN_REVIEW | 复查指定产物 → 人工确认 → `--from sN` 续跑，不是失败 |
| 前序产物 rds 缺失 | 终止并报告缺失路径，不要空状态硬跑 |

## 踩坑清单（详见 anti_patterns.md）

- **不要**在脚本内硬编码阈值/路径——一切以 Config 为准
- **不要**覆盖 `obj$SingleR_label_raw`——人工合并只新增 `major_cluster` 列
- **不要**把 CellChat/Monocle3/pySCENIC 装进 R 主环境——依赖会互相冲突
- **不要**在 Manifest 中改写历史 Step 条目——只追加
- 群内过碎（37 类）≠ 失败——用 marker 驱动的人工合并收敛到 8–10 主群

## 验证清单

- [ ] S0 Quality_Gate：ssh/frp 任一可达 + 全部样本三联文件齐备 + 三个隔离环境存在
- [ ] 每个 Step 的 Manifest 条目 `quality_gate.status` 非空
- [ ] 选定 PC 数与分辨率均已连同依据写入 Manifest
- [ ] `figures/` 产物命名匹配 `Fig<X>_Panel<Y>_<描述>.<格式>`
- [ ] 最终 `manifest.json` 字段齐备（参照 `references/cross_skill_convention.md` 字段表）

## 参考资料

- `references/sgSCR_source.md` — 来源仓库（Bettermelsm/sgSCR）与 PART1–PART8 改造对照
- `references/cross_skill_convention.md` — 目录布局/命名/Manifest 字段/日志格式约定
- `references/singleR_reference_matrix.md` — SingleR 参考集选择矩阵

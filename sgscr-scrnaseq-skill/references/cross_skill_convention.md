# Cross-Skill Convention — sgscr-scrnaseq-skill v1.1

一页跨 skill 最小约定。约定优于框架：本文件是文档，不是共享代码库。其他 skill 可引用本文件作为权威版本。

## 1. 目录布局（一次运行的远程工作目录）

```
{remote_work_dir}/
├── rds/            # 各 Step 中间对象 sN_*.rds（含可选 s2b_obj_corrected.rds，不覆盖 s2 产物）
├── figures/        # 全部图表 Fig<X>_Panel<Y>_<描述>.<格式> + 同名 .source.csv
├── markers/        # FindAllMarkers 结果
├── annotation/     # SingleR 预测 + cluster_mapping 模板/确认表
├── cellchat/       # S7 产物
├── trajectory/     # S8 产物
├── scenic/         # S9 产物
├── target_genes/   # S10 产物
├── stats/          # S10b 产物：sample_level_summary.csv / *_per_sample_detail.csv / design_note.txt
├── diagnostics/    # S10 诊断：ambient_rna.csv / sample_representation.csv
├── module_score/   # S10 模块评分：by_cell_group.csv / by_sample.csv / by_cluster_sample.csv / <set>_genes_used.csv
├── metrics/        # ari_by_resolution.csv 等指标表
└── manifest/       # manifest.json + config 快照 + 各步报告
```

各 Step 只写入自己对应的子目录。

## 2. 命名规则

| 对象 | 规则 | 示例 |
|---|---|---|
| rds | `s<N>_<短描述>.rds` | `s3_obj_pca.rds` |
| 图表 | `Fig<X>_Panel<Y>_<描述>.<格式>` | `Fig3_PanelB_UMAP_major_cluster.pdf` |
| 图源数据 | 与图同名 `.source.csv` | `Fig3_PanelB_UMAP_major_cluster.source.csv` |
| 样本 ID | 原文不转换 | 见 samples.manifest |
| 分组 | 与 Config `group` 字段原文一致 | 见 design.groups |

## 3. Manifest 字段名（登记表）

顶层：`run_id`, `spec_feature`, `created_at`, `random_seed`, `species`, `singleR_reference`, `config_snapshot_path`, `input_checksums`, `excluded_samples`, `steps`, `figures_panel_map`, `warnings`, `convention_version`

steps[] 条目：`step_id`, `started_at`, `finished_at`, `duration_sec`, `env_name`, `package_versions`, `quality_gate{status, detail}`, `human_confirmation{confirmed_by, confirmed_at, file}`, `resolved_params`, `outputs`

figures_panel_map[] 条目：`fig_id`, `panel_id`, `file`, `science_question`, `source_csv`, `font_family`

warnings：`ambient_rna[]`（环境 RNA 告警文本）、`sample_representation[]`（单样本主导告警文本）

- `quality_gate.status` 取值固定：`PASS` | `FAIL` | `NEEDS_HUMAN_REVIEW`
- **新增字段须先在本表登记**，禁止私自扩字段

## 4. 日志格式

```
[<ISO8601时间戳>] [<LEVEL>] [<step_id>] <message>
```

- LEVEL ∈ `INFO` | `WARN` | `ERROR`
- 统一经 `common/log_utils.R` 的 `log_info()/log_warn()/log_error()` 输出，脚本内禁止手写 `cat()`/`print()` 当日志
- 示例：`[2026-10-02T09:00:05+08:00] [INFO] [s3_dimreduce] selected pc.use=22`

## 5. 出图与源数据契约（v1.1 新增）

- 所有出版级图必须经 `common/save_plot.R` 的 `sgs_save_plot()`/`sgs_save_heatmap()` 导出；
  字体（`figure.font_family`）与基准字号（`figure.base_size`）只来自 Config
- `source_df` 必填：含样本标识列（sample_id/sample/orig.ident 之一）与细胞数列（n_cells/n_cell/cell_count 之一），
  缺任一即错误——不存在"只出图不出源数据"的调用点
- 每张图同时落盘同名 `.source.csv`；`source_csv` 与 `font_family` 记入 `figures_panel_map`

## 6. 样本水平统计契约（v1.1 新增，S10b）

- 组均值 = 样本均值之均值；细胞数不进权重
- 完全分离 = min(一组) > max(另一组)；margin 取成立方向的 min−max；未分离 margin ≤0 照样输出
- `separation_robust` = separation 且 margin ≥ `stats.min_separation_margin`；margin_ratio = margin/Δ（Δ=0 → NA）
- 不提供细胞级 p 值；Manifest `resolved_params.min_attainable_two_sided_p` = 2/choose(n0+n1, n0)
- `stats/sample_level_summary.csv` 列契约见 `scripts/s10b_sample_stats.R` 表头与 design.md 需求 23


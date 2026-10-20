# Cross-Skill Convention — sgscr-scrnaseq-skill v1.0

一页跨 skill 最小约定。约定优于框架：本文件是文档，不是共享代码库。其他 skill 可引用本文件作为权威版本。

## 1. 目录布局（一次运行的远程工作目录）

```
{remote_work_dir}/
├── rds/            # 各 Step 中间对象 sN_*.rds
├── figures/        # 全部图表 Fig<X>_Panel<Y>_<描述>.<格式>
├── markers/        # FindAllMarkers 结果
├── annotation/     # SingleR 预测 + cluster_mapping 模板/确认表
├── cellchat/       # S7 产物
├── trajectory/     # S8 产物
├── scenic/         # S9 产物
├── target_genes/   # S10 产物
├── metrics/        # ari_by_resolution.csv 等指标表
└── manifest/       # manifest.json + config 快照 + 各步报告
```

各 Step 只写入自己对应的子目录。

## 2. 命名规则

| 对象 | 规则 | 示例 |
|---|---|---|
| rds | `s<N>_<短描述>.rds` | `s3_obj_pca.rds` |
| 图表 | `Fig<X>_Panel<Y>_<描述>.<格式>` | `Fig3_PanelB_UMAP_major_cluster.pdf` |
| 样本 ID | 原文不转换 | `NC1`, `UC4` |
| 分组 | 与 Config `group` 字段原文一致 | `NC`, `UC` |

## 3. Manifest 字段名（登记表）

顶层：`run_id`, `spec_feature`, `created_at`, `random_seed`, `species`, `singleR_reference`, `config_snapshot_path`, `input_checksums`, `excluded_samples`, `steps`, `figures_panel_map`, `convention_version`

steps[] 条目：`step_id`, `started_at`, `finished_at`, `duration_sec`, `env_name`, `package_versions`, `quality_gate{status, detail}`, `human_confirmation{confirmed_by, confirmed_at, file}`, `resolved_params`, `outputs`

- `quality_gate.status` 取值固定：`PASS` | `FAIL` | `NEEDS_HUMAN_REVIEW`
- **新增字段须先在本表登记**，禁止私自扩字段

## 4. 日志格式

```
[<ISO8601时间戳>] [<LEVEL>] [<step_id>] <message>
```

- LEVEL ∈ `INFO` | `WARN` | `ERROR`
- 统一经 `common/log_utils.R` 的 `log_info()/log_warn()/log_error()` 输出，脚本内禁止手写 `cat()`/`print()` 当日志
- 示例：`[2026-10-02T09:00:05+08:00] [INFO] [s3_dimreduce] selected pc.use=22`

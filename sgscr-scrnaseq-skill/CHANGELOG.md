# Changelog

遵循语义化版本 V{Major}.{Minor}.{Patch}（-Alpha/Beta/RC 后缀），推送前与用户确认版本号。

## [0.3.0-beta] - 2026-10-28

### Added（阶段 12–17：样本水平统计、环境 RNA、代表性告警、模块评分、图源数据契约、文档通用化）
- **S10b 样本水平分离统计**（`scripts/s10b_sample_stats.R` + `scripts/common/sample_stats_utils.R`）：
  组均值只由样本均值构成（M_g = mean of m_s，细胞数不进权重）；完全分离判据 = min(一组) > max(另一组)，
  margin 取成立方向 min−max（未分离 ≤0 照样输出）；margin_ratio = margin/Δ（Δ=0 → NA）；
  margin < `stats.min_separation_margin` → `separation_robust=FALSE`；LOO 数据集（ALL + LOO_<sample_id>）；
  低细胞数（< `stats.min_cells_per_sample`）样本明细标 `low_n_flag`、汇总排除、`n_samples_used` 按组；
  不设细胞级 p 值，Manifest 记录最小可达双侧 p = 2/choose(n0+n1, n0)；
  输出列契约（metric/cell_group/dataset/mean_group0/1/delta/min/max/separation/margin/margin_ratio/robust/n_used/value_<sample_id>）
- **S10 环境 RNA 诊断**（需求 24）：`ambient.lineage_markers` 谱系→本不应表达基因，按 (谱系,基因,样本)
  输出 `diagnostics/ambient_rna.csv` 检出比例；超 `ambient.warn_threshold` 告警入 Manifest `warnings.ambient_rna`
- **S2b 环境 RNA 校正**（`scripts/s2b_ambient_correct.R`，需求 24.4）：**默认关闭**照常跳过；
  开启而后端（celda/decontX）缺失必须**显式失败并报告缺失包名**，不静默、不伪造校正；
  开启时产出 `rds/s2b_obj_corrected.rds`（assay "decontX"，不覆盖 s2 产物）。
  **注：校正后端（celda/decontX）未在开发环境安装，校正路径未经实际运行验证（遵守"不安装新包"约束）；
  已验证的是默认关闭路径的逐元素相等（`tests/unit/test_s2b_skip_matrix_equal.R`）与开启-缺后端的显式失败分支**
- **S10 样本代表性告警**（需求 25）：(组,细胞群) 最大单样本占比 > `diagnostics.single_sample_dominance_threshold`
  （默认 0.70）→ 告警 + `single_sample_driven=TRUE` + 主导样本 ID；经 `diagnostics/sample_representation.csv`
  传递到 S10b 汇总表与 S11 图注
- **S10 模块特征评分**（需求 29）：`module_scores` 配置集合名→基因列表；AddModuleScore 缺失基因
  只标记不静默丢弃（`<set>_genes_used.csv`）；输出 by_cell_group/by_sample/by_cluster_sample；
  只写 `<set>_score1` 数值列，静态检查禁止离散标签回写
- **出图源数据契约**（需求 18.5–18.7）：`sgs_save_plot`/`sgs_save_heatmap` 统一套用
  `figure.font_family`/`figure.base_size`（字体来自 Config 非脚本字面量）；强制 `source_df`
  （含样本标识列与细胞数列，缺一即错误）；同名 `.source.csv` 同目录落盘；源数据路径+字体记入
  Manifest `figures_panel_map`；s2/s3/s4/s6/s7/s8/s10 全部调用点改造，无"只出图不出源数据"遗留
- S11 质量门控纳入 s10b；`figures_panel_map` 记录 source_csv 与 font_family；图注注入 `single_sample_driven` 注记
- 测试：属性测试 ×5（23.1 组均值/23.3-23.9 分离判据与实测量级/25.1-25.2 代表性/29.3 模块缺失标记/16 出图契约）
  + 单元测试 ×2（23.2/23.5/23.6 LOO 与最小可达 p/24.1-24.4 与 18.5-18.7），全部 Rscript 直跑通过
- 文档通用化：SKILL.md 移除案例细节改以配置契约描述；案例迁入 `examples/case-UC-10X-2025.md`
  （7 样本纳入决策、57,230 细胞、决议参数、科学结论摘要）

### Fixed
- s10 环境 RNA 块引用不存在变量 `clusters` → `cluster`（mock 全链暴露）
- 测试侧三处判据错误（分离方向取错分支、检出比例计数错、hedgehog 双参数生成器 API 误用）

### Fixed（P12-d，Kiro 0053 PR per-file review）
- `margin_ratio = margin/|Δ|`（原 margin/Δ 在 case 低于 control 时给出误导性负比例）
- 未分离时 margin 取两个方向 gap 中更接近 0 的（语义 = "距分离还差多少"），不再依赖分支顺序
- `sgs_separation_stats` 的 `min_separation_margin` **取消默认值**——缺省 0 会让任何已分离指标恒判 robust
- `sgs_sample_level_columns` 列契约与实际表头对齐（base + value_<sid> + low_n_<sid> + single_sample_driven）
- **本底对照群双口径并列输出**（24.3 科学判定项）：auto `nonexpressing_background`（检出≤background_max_detection，
  局限：会把环境 RNA 污染重的群排除出本底→本底差值被低估→目标谱系特异性偏强，偏差方向有利于预期结论）
  + 谱系 `lineage_background`（`ambient.background_exclude_lineages` 之外的细胞），口径选择记入 Manifest
  `resolved_params.background_exclude_lineages`；论文报告以谱系口径为主（最终口径待用户/Kiro 确认）
- s10b `n_cells` 改为行构造时随细胞集合一并计算（mock 重跑暴露：本底口径伪群不在
  `table(clusters, samples)` 中，事后合并回算得 0 → 全行被 low_n_flag 排除 → 汇总全 NA）

### Note（P12-c 验证路径披露，Kiro 0052 要求）
- s10 主导告警端到端验证经**调低阈值**（single_sample_dominance_threshold 0.70→0.25）的 mock 场景完成，
  未遇天然超阈数据；默认阈值 0.70 下该路径逻辑相同（同一 `sgs_single_sample_dominance` 判据），仅触发条件不同。

### Known issues
- s11_finalize 对"复用上游产物的续跑 run"按 11 步完整计数会 FAIL（修复方案待定，同 0.2.0）

## [0.2.0-beta] - 2026-10-20

### Added
- **S6b 亚群注释与保留过滤**（`scripts/s6b_subcluster_filter.R` + `scripts/common/subcluster_filter_utils.R`）：
  标签 CSV（compartment/subcluster/keep/proposed_label）与 s6 对象一一对应校验（缺失/多余/重复即 FAIL）；
  keep 语义 yes/sens/no，filter_mode main|strict；亚群列 sprintf 精确匹配（禁止 grep 首列）；
  Manifest 记录 CSV MD5、各亚群/各样本保留数
- s7/s8/s9 优先消费 S6b 清洗产物：s7 在通讯计算前剔除被过滤细胞；s8 优先读 s6b + 支持
  `trajectory.start_cluster` 指定起点簇（如 myeloid=classical_monocyte）；SCENIC 导出改用清洗后子集
- CellChat liftCellChat 组集合对齐（各条件 cell group 并集，修复 mergeCellChat "非整合陈列"）
- s8 multi-compartment：target_compartment 支持列表，primary/backup 角色标注，FAIL 隔离不互相阻断
- Seurat v5 layer 清理工具 `sgs_drop_scale_layers`（CreateAssay5Object 重建 + Key 保留），
  `sgs_assert_layers_clean` 落盘前断言，`sgs_check_marker_coverage` marker 覆盖门检
- 属性测试 `tests/property/test_subcluster_filter.R`；集成测试 `mock_cellchat_lift.R`、`mock_s8_multicomp.R`
- 真实数据（UC 结肠炎 7 样本 57,230 细胞）全链验证：S6b→s7–s10 全 PASS；strict 敏感性全 PASS

### Fixed
- FindAllMarkers 静默零输出（S1 逐样本 ScaleData 残留 scale.data.<sample> layers）
- v3 Assay 在 SeuratObject 5.2.0 上 JoinLayers 无方法；CreateAssay5Object 重建后 Key 丢失
- s8 成功判定误用 detail 字符串匹配 → 显式 comp_ok 标志；UMAP 列名硬编码 → 动态取用

### Known issues
- s11_finalize 对"复用上游产物的续跑 run"按 11 步完整计数会 FAIL（copy-run 的 Manifest 无 s0–s5 step 条目），
  修复方案待定（占位 step 条目 / 接受 provenance / 动态步数）

## [0.1.0-Alpha] - 2026-10-02

### Added（阶段 0：Skill 骨架与 Config/Manifest 基础设施）
- 目录骨架（SKILL.md / config/ / scripts/common/ / references/ / examples/ / tests/）
- SKILL.md 能力声明与触发条件
- references/sgSCR_source.md（来源引用 + PART1–PART8 改造对照表）
- references/cross_skill_convention.md v1.0（目录布局/命名/Manifest 字段表/日志格式）
- references/singleR_reference_matrix.md（占位，阶段 4 补全）
- config/config.example.yaml（全字段分组，TBD 占位符原样保留）
- scripts/common/config_loader.R（一次性列全缺失字段；TBD 视为未填）
- scripts/common/manifest_utils.R（只追加 Manifest，JSON）
- scripts/common/log_utils.R（统一日志格式）
- tests/property/test_config_loader.R（P2/P3）、test_manifest_utils.R（P13/P14）— 待 R 环境执行

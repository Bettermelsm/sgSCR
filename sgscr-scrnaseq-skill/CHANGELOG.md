# Changelog

遵循语义化版本 V{Major}.{Minor}.{Patch}（-Alpha/Beta/RC 后缀），推送前与用户确认版本号。

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

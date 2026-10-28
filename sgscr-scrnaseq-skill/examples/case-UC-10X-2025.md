# 案例：人溃疡性结肠炎结肠黏膜 10X（2025）

> 本文件是该技能首个真实案例的**决议记录与科学结论摘要**，案例细节从 SKILL.md 迁入此处（SKILL.md 保持项目无关）。
> 涉及样本 ID、分组名与服务器路径仅作案例记录，不构成技能默认值。

## 数据与纳入决策

- 组织：人溃疡性结肠炎（UC）结肠黏膜，10X 三联文件
- 样本：7 样本 = 3 对照（NC1/NC3/NC4）+ 4 UC（UC1/UC2/UC3/UC4）
- **7 样本纳入决策**：UC4 因 QC 门控不合格被排除（excluded 列表 + 排除原因写 Manifest），最终 7 样本（3+4）进入整合分析
- QC 后细胞数：**57,230**；QC 阈值与各步决议参数以 Manifest `resolved_params` 为准

## 各步实际决议参数（摘要）

| 决议点 | 实际取值 | 依据 |
|---|---|---|
| QC 阈值 | min_features 200 / max_features 6000 / max_percent_mt 20（UC 组织线粒体偏高，放宽到 20） | QC 小提琴图人工复查 |
| PC 数 | 定量选 PC（累计方差 + 单 PC 贡献 + 拐点三判据） | S3 Manifest |
| 聚类分辨率 | ARI 辅助 + clustree 人工确认 | S4 Manifest |
| 主群合并 | marker 驱动人工合并为 8–10 主群（含 epi_absorptive/epi_secretory/myeloid/fibroblast/endothelial 等） | cluster_mapping_confirmed.csv |
| 目标基因 | 离子转运轴 + 杯状/分泌谱系基因 + 炎症小体模块评分 | 项目科学问题 |

## 样本水平统计（s10b）关键结论

- **炎症小体模块评分**：样本水平**完全分离**——3 个 NC 样本（−0.214/−0.249/−0.299）全部低于 4 个 UC 样本（−0.037/−0.092/−0.176/−0.193），即 min(NC) > max(UC)，分离稳健
- **NOS2**：检出细胞占比仅 ~2.4%，如实报告"接近零表达"，不做过度解读
- **离子转运轴本底警示**：CA2（本底群 margin ≈ 0.0006）与 BEST4（吸收上皮 margin ≈ 0.0015）虽形式上"完全分离"，但 margin 低于 `stats.min_separation_margin`（0.01），一律标记 `separation_robust=FALSE`——量级属噪声
- 小样本设计（3 vs 4）：不设细胞级 p 值；最小可达双侧 p = 2/choose(7,3) ≈ 0.057，写入报告说明

## 小鼠验证加测清单（按修正后优先级）

1. **首选**：GUCA2A（abs，效应量最大）+ GUCA2B（abs/sec 两隔室均分离）
2. **次选**：AQP8 (abs) > CA2 (abs) > CA2 (sec) > AQP8 (sec) > BEST4 (abs)
3. **不进清单**：GUCY2C、CFTR、SLC26A3（分离不成立或不稳健）

> 清单以 `human_STable_iontransport_separation.csv`（含 ctrl/case/bg min-max 审计列与判定定义表）为准；bg_diff 列为谱系本底差值。

## 经验教训（本项目贡献给技能的通用规则）

- 分离判定必须先算**谱系本底群**的分离再做"上皮特异"结论——曾因本底分离被误判为 TRUE 导致 GUCA2A 结论翻转（P13-d 修正）
- margin 量级 1e-3 级别的"完全分离"在基因表达本底下属于噪声，必须配 `separation_robust` 双重标记
- 组均值只由样本均值构成：细胞数不进权重，避免样本内细胞数失衡扭曲组间比较
- 出图强制源数据（.source.csv）：所有发表图可从源数据完整复算

# SingleR 参考集选择矩阵

> 占位文档：阶段 4（两层注释）执行时补全。物种校验逻辑见 `common/species_check`（P5）。

## 已知参考集（待补充验证）

| 参考集 | 物种 | 组织覆盖 | 本案例适配度 | 备注 |
|---|---|---|---|---|
| HumanPrimaryCellAtlasData | human | 广泛 | 候选 | DesignR 包，需下载 |
| BlueprintEncodeData | human | 血液/基质 | 候选 | |
| MouseRNAseqData | **mouse** | — | ❌ 禁用于本案例 | 历史缺陷：人数据误用小鼠参考集，已加 P5 校验拦截 |

## 校验规则（P5）

Config `species.organism` 必须与所选参考集所属物种一致，不一致 S5 直接终止并报告双方物种名。
参考集→物种映射表在本文件补全后，由校验函数读取本文件（或 Config 中显式 `reference_species` 字段）。

# Anti-Patterns — sgscr-scrnaseq-skill

已踩过/明确禁止的做法。发现新反模式请追加并注明日期。

1. **硬编码阈值/路径进 R 脚本**（2026-10-02）——一切参数以 Config 为准；config_loader 会把 TBD 占位符当缺失拦截，不要为绕过校验把默认值写回脚本。
2. **覆盖 `obj$SingleR_label_raw`**（2026-10-02）——人工合并只新增 `major_cluster` 列；历史缺陷"人数据误用小鼠参考集"已由 P5 物种校验拦截，不要删校验图省事。
3. **把 CellChat/Monocle3/pySCENIC 装进 R 主环境**（2026-10-02）——Matrix/Seurat/Bioc/numpy 依赖互相冲突，一律走隔离 conda env。
4. **改写 Manifest 历史 Step 条目**（2026-10-02）——只追加；发现前一步写错了，追加更正条目而不是原地改。
5. **日志用手写 cat()/print()**（2026-10-02）——统一 log_utils.R，否则日志无法按约定正则解析（P17）。
6. **人工确认点被跳过直接续跑**（2026-10-02）——QC 图/clustree/映射表三处 NEEDS_HUMAN_REVIEW 必须等用户确认。

# sgSCR 来源引用与改造对照

**来源仓库**：https://github.com/Bettermelsm/sgSCR
**核心脚本**：`SG_scRNA_seq_FLOW_V2.1.R`（单体 R 脚本，PART1–PART8，Seurat v5 工作流）
**首个验证案例**：远程服务器上人 UC 肠黏膜 10X 数据（NC1/NC3/NC4 vs UC1/UC2/UC3/UC4，7 样本，QC 后 57,230 细胞；UC4 保留，裁决依据见 Manifest sample_decisions）

## PART1–PART8 → S0–S11 改造对照表

| 原脚本 | 原功能 | Skill 归属 | 改造内容 |
|---|---|---|---|
| PART1 | 环境加载/路径硬编码 | S0 + `common/config_loader.R` | 路径与参数全部入 Config；新增三联文件校验、隔离环境盘点 |
| PART2 | 批量 Read10X 导入 | S1 | 新增 manifest 溯源（输入校验和）、Quality_Gate |
| PART3 | QC 过滤（percent.mt/percent.HB） | S2 | 阈值入 Config；新增 QC 小提琴图人工复查暂停点 |
| PART4 | 归一化/HVG/PCA；**个性化改造：定量选 PC 数** | S3 | 规则固化为「累计方差≥阈值 且 单PC贡献<阈值 + 拐点检测」，阈值入 Config，决议值写入 Manifest（P8） |
| PART5 | 多分辨率聚类 + clustree；**个性化改造：ARI 辅助选分辨率** | S4 | ARI 规则确定性化，决议值写入 Manifest（P9）；新增 clustree 人工复查暂停点 |
| PART6 | SingleR 注释 + FindAllMarkers | S5 | 新增物种-参考集一致性校验（P5，修复人数据误用小鼠参考集缺陷）；新增两层注释（人工合并映射表 → major_cluster，P10） |
| PART7 | 结果保存 rds | S1–S11 | 固化到各 Step 落盘约定（`rds/sN_*.rds`） |
| PART8 | 出图 | S11 | 统一 `save_plot()`，`Fig<X>_Panel<Y>_<描述>` 命名，格式/dpi 入 Config（P16） |
| （无） | Compartment 亚聚类 | S6 | 新增 |
| （无） | 细胞通讯 | S7 | 新增，CellChat，隔离环境 `sgscr_cellchat` |
| （无） | 拟时序 | S8 | 新增，默认 slingshot（Monocle3 为 Config 可选项），隔离环境 `sgscr_trajectory` |
| （无） | 转录因子调控 | S9 | 新增，默认 pySCENIC（R 版 SCENIC 为可选项），隔离环境 `sgscr_scenic` |
| （无） | 目标基因追踪 | S10 | 新增（SLC26A3/GUCA2A 定位与 UC vs NC 差异） |
| （无） | 溯源 | 全程 | 新增 JSON Manifest（只追加），随机种子贯穿全程（P15） |

## 保留承诺（需求 20）

原脚本的全部既有能力（整合、QC、注释、聚类、clustree、rds 保存）在改造中保留为回归验收项，对应 `tests/unit/` 集成测试。

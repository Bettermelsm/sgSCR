# sgscr-scrnaseq-skill

sgSCR scRNA-seq 标准化分析流程（Agent Skill）。能力声明与触发条件见 [SKILL.md](SKILL.md)。

面向初级用户的快速说明：
1. 本 skill 驱动 SGlcl01 上的 R/conda 环境做单细胞分析，本机只发命令、取图表
2. 使用前把 `config/config.example.yaml` 复制为 `config/config.yaml`，填完所有 `<TBD, 用户确认>` 占位符
3. 按 SKILL.md「标准工作流」逐步执行；每步结果都落盘，中断可 `--from sN` 续跑
4. 流程文档（需求/设计/任务清单）：`D:\01SG_FILES\20250904UC2025\.kiro\specs\sgscr-scrnaseq-skill\`

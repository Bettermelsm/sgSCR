# config_loader.R — Config_File 加载与必需字段校验
# Validates: Requirements 2.3, 2.4, 6.3, 6.4 (P2, P3)
#
# 设计要点：
# - 一次性读入 YAML，解析为参数对象，全程无交互式输入（P2）
# - 缺必需字段时一次性列出【全部】缺失字段名并终止，而非只报第一个（P3）
# - 需要人工决定的占位符 <TBD, 用户确认> 视为缺失（防止占位值流入生产运行）
# - 脚本内不出现任何硬编码阈值/路径

suppressPackageStartupMessages(library(yaml))

`%||%` <- function(a, b) if (is.null(a)) b else a   # 全链通用空值兜底

# 必需字段表（点路径表示嵌套层级）。新增字段须同步登记 references/cross_skill_convention.md。
.sgs_REQUIRED_FIELDS <- c(
  "run.run_id", "run.random_seed",
  "species.organism", "species.singleR_reference", "species.reference_species",
  "samples.manifest",
  "remote.host_alias", "remote.ssh_target", "remote.frp_fallback.host",
  "remote.frp_fallback.port", "remote.data_dir", "remote.remote_work_dir",
  "remote.triplet_files",
  "envs.main", "envs.cellchat", "envs.trajectory", "envs.scenic",
  "qc.min_features", "qc.max_features", "qc.min_cells_per_gene",
  "qc.max_percent_mt", "qc.max_percent_hb", "qc.min_cells_pass_gate",
  "integration.dims", "integration.reduction", "integration.nfeatures_hvg",
  "pca.npcs_compute", "pca.cum_variance_threshold_pct",
  "pca.single_pc_contrib_threshold_pct", "pca.elbow_diff_threshold",
  "pca.fallback_npc",
  "clustering.candidate_resolutions", "clustering.ari_drop_threshold",
  "compartments.targets", "compartments.subcluster_resolution",
  "cellchat.pathways_of_interest", "cellchat.min_cells_per_group",
  "trajectory.tool", "trajectory.target_compartment",
  "scenic.backend", "scenic.species_db_path", "scenic.target_subset",
  "target_genes.genes", "target_genes.science_questions",
  "figure.dpi", "figure.width_in", "figure.height_in",
  "figure.formats", "figure.panel_map",
  "paths.figures_dir", "paths.rds_dir", "paths.manifest_path"
)

# 沿点路径逐层取值；任一层不存在返回 NULL
.sgs_get_at_path <- function(x, path) {
  parts <- strsplit(path, ".", fixed = TRUE)[[1]]
  cur <- x
  for (p in parts) {
    if (is.null(cur) || !p %in% names(cur)) return(NULL)
    cur <- cur[[p]]
  }
  cur
}

# 判定占位符/空值：TBD 文案、NULL、长度 0、空串均视为"未填"
.sgs_is_unresolved <- function(val) {
  if (is.null(val)) return(TRUE)
  if (length(val) == 0) return(TRUE)
  if (is.character(val)) {
    joined <- paste(val, collapse = " ")
    if (all(nzchar(trimws(val))) == FALSE && !is.data.frame(val)) return(TRUE)
    if (grepl("TBD", joined, fixed = TRUE)) return(TRUE)
  }
  FALSE
}

# 返回缺失（或仍为占位符）的必需字段点路径向量；合法时返回 character(0)
sgs_missing_fields <- function(cfg, required = .sgs_REQUIRED_FIELDS) {
  missing <- character(0)
  for (f in required) {
    val <- .sgs_get_at_path(cfg, f)
    if (.sgs_is_unresolved(val)) missing <- c(missing, f)
  }
  missing
}

# 主入口：读文件 → 校验 → 返回参数对象。缺字段直接停，一次性列全。
sgs_load_config <- function(config_path, required = .sgs_REQUIRED_FIELDS) {
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path), call. = FALSE)
  }
  cfg <- yaml::read_yaml(config_path)
  missing <- sgs_missing_fields(cfg, required)
  if (length(missing) > 0) {
    msg <- c(
      sprintf("Config incomplete: %d required field(s) missing or unresolved (TBD):", length(missing)),
      paste0("  - ", missing)
    )
    stop(paste(msg, collapse = "\n"), call. = FALSE)
  }
  cfg
}

# s10b_sample_stats.R — 样本级统计与敏感性分析（阶段12 固化，需求 23，主环境）
# Validates: Requirements 23.1-23.9, 24.3（本底对照群）, 25.2（代表性标记传递）
#
# 计算顺序（design.md）：每样本均值 m_s → 组均值 M_g → delta；细胞数不进权重
# 指标来源：target genes（s5 表达）+ module scores（s10 产物 module_score/by_cluster_sample.csv）
# 输出：stats/sample_level_summary.csv + stats/sensitivity_leave_one_out.csv（列契约见 design.md）
#       明细 stats/per_sample_detail.csv（低细胞数样本 low_n_flag 标记）
#       本设计最小可达双侧 p 值记入 Manifest 与 stats/design_note.txt（23.5）
# 不输出任何细胞级 p 值（23.5）

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s10b_sample_stats.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R",
            "sample_stats_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s10b_sample_stats"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

in_rds <- file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds")
run_root <- normalizePath(file.path(dirname(cfg$paths$manifest_path), ".."))  # = {remote_work_dir}/run
mod_file <- file.path(run_root, "module_score", "by_cluster_sample.csv")
repr_file <- file.path(run_root, "diagnostics", "sample_representation.csv")
out_dir <- file.path(run_root, "stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pre <- tryCatch(sgs_resume_plan("s10b", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

problems <- c()
if (!file.exists(in_rds)) problems <- c(problems, "s5_obj_annotated.rds")
if (!file.exists(mod_file)) problems <- c(problems, "module_score/by_cluster_sample.csv（s10 先行）")
if (length(problems) > 0) sgs_preflight(STEP, problems, character(0))
stale_out <- file.path(out_dir, "sample_level_summary.csv")
if (file.exists(stale_out)) sgs_check_freshness(stale_out, in_rds)

# ---- 参数（23.8/23.6，均来自 Config） ----
min_margin <- if (!is.null(cfg$stats$min_separation_margin)) cfg$stats$min_separation_margin else 0.01
min_cells  <- if (!is.null(cfg$stats$min_cells_per_sample)) cfg$stats$min_cells_per_sample else 10
bg_max_det <- if (!is.null(cfg$ambient$background_max_detection)) cfg$ambient$background_max_detection else 0.05

grp_cfg <- if (!is.null(cfg$design$groups)) as.character(unlist(cfg$design$groups)) else c("ctrl", "case")
g0 <- grp_cfg[1]; g1 <- if (length(grp_cfg) > 1) grp_cfg[2] else grp_cfg[1]

obj <- readRDS(in_rds)
obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
meta <- obj@meta.data
samples <- as.character(meta$sample_id)
groups  <- as.character(meta$group)
clusters <- as.character(meta$major_cluster)

# ---- 指标矩阵：target genes + module scores（每细胞一值） ----
tgt <- as.character(cfg$target_genes$genes)
tgt_kept <- intersect(tgt, rownames(obj))
if (length(setdiff(tgt, tgt_kept))) {
  log_warn(STEP, sprintf("目标基因缺失（只标记不阻断）: %s", paste(setdiff(tgt, tgt_kept), collapse = ",")))
}
expr_m <- if (length(tgt_kept)) as.matrix(FetchData(obj, vars = tgt_kept, layer = "data")) else NULL

mod <- read.csv(mod_file, check.names = FALSE)
# by_cluster_sample.csv 契约：key 列 cell_group/sample_id/group/n_cells + 每模块评分一列（列名=<set>_score）
mod_sets <- setdiff(names(mod), c("cell_group", "sample_id", "group", "n_cells"))
if (!length(mod_sets)) stop("module_score/by_cluster_sample.csv 无模块评分列", call. = FALSE)

# ---- 组装每指标 × 每细胞群 × 每样本 均值长表 ----
# detail 行：(metric, cell_group, sample_id, group, mean, n_cells, low_n_flag)
make_detail <- function(vals, sid, grp, clu, metric, cell_group) {
  data.frame(metric = metric, cell_group = cell_group, sample_id = sid, group = grp,
             value = as.numeric(vals), stringsAsFactors = FALSE)
}
det_list <- list()
# 目标基因 × 各 major_cluster
if (!is.null(expr_m)) {
  for (gcol in colnames(expr_m)) {
    for (cl in sort(unique(clusters))) {
      i <- clusters == cl
      det_list[[length(det_list) + 1]] <- data.frame(
        metric = gcol, cell_group = cl,
        sample_id = samples[i], group = groups[i], value = expr_m[i, gcol],
        stringsAsFactors = FALSE)
    }
    # 24.3 本底对照群：该基因按主群检出判定为不表达的细胞群合并
    det_pct <- tapply(expr_m[, gcol] > 0, clusters, mean)
    bg_cl <- sgs_background_clusters(det_pct, max_detection = bg_max_det)
    if (length(bg_cl)) {
      i <- clusters %in% bg_cl
      det_list[[length(det_list) + 1]] <- data.frame(
        metric = gcol, cell_group = "nonexpressing_background",
        sample_id = samples[i], group = groups[i], value = expr_m[i, gcol],
        stringsAsFactors = FALSE)
    }
  }
}
# 模块评分 × 细胞群（来自 s10 交叉聚合表的逐样本均值；列名即 <set>_score）
if (nrow(mod)) {
  det_list[[length(det_list) + 1]] <- do.call(rbind, lapply(mod_sets, function(ms) {
    data.frame(metric = ms, cell_group = mod$cell_group,
               sample_id = mod$sample_id, group = mod$group, value = mod[[ms]],
               stringsAsFactors = FALSE)
  }))
}
detail <- do.call(rbind, det_list)
if (is.null(detail) || !nrow(detail)) stop("无可用指标（target genes 与 module scores 均为空）", call. = FALSE)
detail$n_cells <- 1L  # 占位，下面按 (cell_group, sample) 重算真实细胞数

# ---- 低细胞数标记（23.6）：按 (cell_group, sample) 细胞数 ----
cell_n <- as.data.frame(table(cell_group = clusters, sample_id = samples), responseName = "n")
detail <- merge(detail, cell_n, by = c("cell_group", "sample_id"), all.x = TRUE)
detail$n_cells[is.na(detail$n_cells)] <- 0L
detail$low_n_flag <- detail$n_cells < min_cells

write.csv(detail[, c("metric", "cell_group", "sample_id", "group", "n_cells", "value", "low_n_flag")],
          file.path(out_dir, "per_sample_detail.csv"), row.names = FALSE)

# ---- 单样本代表性标记传递（25.2：从 s10 诊断表读到汇总行） ----
driven_key <- character(0)
if (file.exists(repr_file)) {
  repr <- read.csv(repr_file, check.names = FALSE)
  driven_key <- sprintf("%s|%s", repr$group, repr$cell_group)[repr$single_sample_driven %in% TRUE]
}
driven_for_row <- function(grp_v, clu) {
  any(driven_key == sprintf("%s|%s", grp_v, clu))
}

# ---- 汇总 + LOO（23.1/23.2/23.3/23.7/23.8/23.9） ----
sid_universe <- sort(unique(detail$sample_id))   # 全样本列集：ALL 与 LOO 行保持同列契约
run_one <- function(d, dataset_label) {
  m_s_all <- tapply(d$value, d$sample_id, mean)
  grp_of <- stats::setNames(d$group[match(names(m_s_all), d$sample_id)], names(m_s_all))
  keep <- !d$low_n_flag
  m_s <- tapply(d$value[keep], d$sample_id[keep], mean)
  grp_of <- grp_of[names(m_s)]
  st <- sgs_separation_stats(m_s, grp_of, g0, g1, min_separation_margin = min_margin)
  out <- data.frame(
    metric = d$metric[1], cell_group = d$cell_group[1], dataset = dataset_label,
    mean_group0 = st$mean_group0, mean_group1 = st$mean_group1, delta = st$delta,
    min_group0 = st$min_group0, max_group0 = st$max_group0,
    min_group1 = st$min_group1, max_group1 = st$max_group1,
    separation = st$separation, separation_margin = st$separation_margin,
    margin_ratio = st$margin_ratio, separation_robust = st$separation_robust,
    n_samples_used_group0 = st$n_samples_used_group0,
    n_samples_used_group1 = st$n_samples_used_group1, stringsAsFactors = FALSE)
  # 每样本逐值列（value_<sid>，全样本统一列集）；低细胞数样本列置 NA 并附 low_n_<sid>=TRUE 标记；
  # LOO 行中被剔除的样本置 NA（low_n_<sid> 保持其真实标记）
  for (s in sid_universe) {
    i_row <- which(d$sample_id == s)[1]
    low_n_s <- if (is.na(i_row)) FALSE else isTRUE(d$low_n_flag[i_row])
    out[[paste0("value_", s)]] <- if (is.na(i_row) || low_n_s) NA_real_ else m_s_all[[s]]
    out[[paste0("low_n_", s)]] <- low_n_s
  }
  out$single_sample_driven <- driven_for_row(d$group[1], d$cell_group[1])
  out
}

key_cols <- c("metric", "cell_group")
summary_rows <- list(); loo_rows <- list()
for (k in split(seq_len(nrow(detail)), detail[key_cols])) {
  d <- detail[k, ]
  if (!nrow(d)) next
  summary_rows[[length(summary_rows) + 1]] <- run_one(d, "ALL")
  for (loo_lab in names(sgs_loo_samples(unique(d$sample_id)))) {
    keep_s <- sgs_loo_samples(unique(d$sample_id))[[loo_lab]]
    d2 <- d[d$sample_id %in% keep_s, ]
    if (nrow(d2)) loo_rows[[length(loo_rows) + 1]] <- run_one(d2, loo_lab)
  }
}
summary_df <- do.call(rbind, summary_rows)
loo_df <- do.call(rbind, loo_rows)

write.csv(summary_df, file.path(out_dir, "sample_level_summary.csv"), row.names = FALSE)
write.csv(loo_df, file.path(out_dir, "sensitivity_leave_one_out.csv"), row.names = FALSE)

# ---- 23.5：最小可达双侧 p 值（按各组最大样本数计） ----
n0 <- length(unique(samples[groups == g0])); n1 <- length(unique(samples[groups == g1]))
p_min <- sgs_min_attainable_pvalue(n0, n1)
note <- sprintf(paste("小样本设计声明（需求23.5）：不输出细胞级 p 值（细胞非独立生物学重复）。",
  "本设计每组样本数 n0=%d, n1=%d，精确秩检验下最小可达双侧 p = %.4g。",
  "stats.min_separation_margin=%g, stats.min_cells_per_sample=%d。"), n0, n1, p_min, min_margin, min_cells)
writeLines(note, file.path(out_dir, "design_note.txt"))
log_info(STEP, note)

low_n_list <- sort(unique(detail$sample_id[detail$low_n_flag]))
gate <- sgs_qc_gate(auto_pass = TRUE,
  detail = sprintf("metrics=%d cell_groups=%d; summary %d rows, LOO %d rows; low_n samples: %s; min attainable two-sided p=%.4g",
                   length(unique(detail$metric)), length(unique(detail$cell_group)),
                   nrow(summary_df), nrow(loo_df),
                   if (length(low_n_list)) paste(low_n_list, collapse = ",") else "(none)", p_min))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(
    min_cells_per_sample = min_cells,
    min_separation_margin = min_margin,
    small_sample_design = list(max_samples_per_group = max(n0, n1),
                               omit_cell_level_pvalue = TRUE),
    min_attainable_two_sided_p = p_min,
    n_samples_used = list(group0 = n0, group1 = n1),
    low_n_samples = as.list(low_n_list),
    background_max_detection = bg_max_det),
  outputs = as.list(c(file.path(out_dir, "sample_level_summary.csv"),
                      file.path(out_dir, "sensitivity_leave_one_out.csv"),
                      file.path(out_dir, "per_sample_detail.csv"),
                      file.path(out_dir, "design_note.txt"))),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
log_info(STEP, "S10b complete")

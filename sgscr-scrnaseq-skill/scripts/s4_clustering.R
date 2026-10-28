# s4_clustering.R — 多分辨率聚类 + clustree + ARI 辅助分辨率选择（PART5 固化）
# Validates: Requirements 10.1-10.3, 8.4, 20.3
# 输入：s3_obj_pca.rds（使用 Config resolutions 的 pc.use 由 S3 resolutions 或重算）
# 输出：rds/s4_obj_cluster.rds、figures/P4.1_Resolution_clustree.*、metrics/ari_by_resolution.csv
# Quality_Gate：自动=ARI 曲线产出且选定分辨率有效；人工=复查 clustree 图

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s4_clustering.R --config <config.yaml> [--from s4]", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resolution_selection.R", "resume_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s4_clustering"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# ---- 断点续跑前置检查（P6 真实路径） ----
pre <- tryCatch(sgs_resume_plan("s4", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

in_rds <- file.path(cfg$paths$rds_dir, "s3_obj_pca.rds")
obj <- readRDS(in_rds)

# S3 已决议的 pc.use（读 Manifest resolutions；缺则按 S3 同规则重算，双保险）
pc_use <- sgs_manifest_get_resolved(cfg$paths$manifest_path, "pc_used")
if (is.null(pc_use)) {
  stdev <- obj@reductions$pca@stdev
  vp <- stdev^2 / sum(stdev^2) * 100
  pc_use <- sgs_select_pcs(vp, cfg$pca$cum_variance_threshold_pct,
                           cfg$pca$single_pc_contrib_threshold_pct,
                           cfg$pca$elbow_diff_threshold,
                           cfg$pca$fallback_npc)$pc_use
}
obj <- Seurat::RunUMAP(obj, dims = 1:pc_use, verbose = FALSE)

# ---- SNN graph + 多分辨率聚类（需求 10.1） ----
obj <- Seurat::FindNeighbors(obj, dims = 1:pc_use, reduction = "pca", verbose = FALSE)
resolutions <- as.numeric(cfg$clustering$candidate_resolutions)
for (r in resolutions) {
  obj <- Seurat::FindClusters(obj, resolution = r, verbose = FALSE)
}
clustree_plt <- clustree::clustree(obj, prefix = "RNA_snn_res.")
# 阶段15：源数据 = 各分辨率 × 簇 × 样本细胞数
res_cols <- sprintf("RNA_snn_res.%g", resolutions)
clu_src <- do.call(rbind, lapply(res_cols, function(rc) {
  tb <- table(obj@meta.data[[rc]], obj$sample_id)
  df <- as.data.frame(tb, responseName = "n_cells")
  names(df)[1] <- "cluster"
  data.frame(resolution = rc, cluster = as.character(df$cluster),
             sample_id = as.character(df$Var2), n_cells = df$n_cells,
             stringsAsFactors = FALSE)
}))
figres <- sgs_save_plot(cfg, "P4.1", "PanelA", "Resolution_clustree", clustree_plt, clu_src, STEP)
figs <- unlist(figres$files)

# ---- ARI 相邻分辨率曲线（需求 10.2，纯 R ARI） ----
ari_seq <- vapply(2:length(resolutions), function(i) {
  c1 <- obj@meta.data[[sprintf("RNA_snn_res.%g", resolutions[i - 1])]]
  c2 <- obj@meta.data[[sprintf("RNA_snn_res.%g", resolutions[i])]]
  sgs_ari(c1, c2)
}, numeric(1))
metrics_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "metrics")
dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)
ari_df <- data.frame(
  resolution = resolutions,
  n_clusters = vapply(sprintf("RNA_snn_res.%g", resolutions),
                       function(cn) length(unique(obj@meta.data[[cn]])), numeric(1)),
  ari_to_next = c(ari_seq, NA))
write.csv(ari_df, file.path(metrics_dir, "ari_by_resolution.csv"), row.names = FALSE)

# ---- 分辨率选择（P9） ----
sel <- sgs_select_resolution(resolutions, ari_seq, cfg$clustering$ari_drop_threshold)
sel_col <- sprintf("RNA_snn_res.%g", sel$resolution)
obj$selected_resolution <- obj@meta.data[[sel_col]]
obj$seurat_clusters <- obj@meta.data[[sel_col]]
log_info(STEP, sprintf("resolution selected: %g (%s)", sel$resolution, sel$basis))

# ---- Quality_Gate（P7 通用函数） ----
gate <- sgs_qc_gate(
  auto_pass = file.exists(file.path(metrics_dir, "ari_by_resolution.csv")) &&
              !is.na(sel$resolution) && length(unique(obj$selected_resolution)) > 1,
  detail = sprintf("selected res=%g (%s); pc.use=%s", sel$resolution, sel$basis, pc_use),
  human_review_file = "figures/P4.1_Resolution_clustree (confirm resolution, then rerun --from s5)")
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
out_rds <- file.path(cfg$paths$rds_dir, "s4_obj_cluster.rds")
saveRDS(obj, out_rds)
sgs_manifest_set_resolved(cfg$paths$manifest_path, "resolution", sel$resolution)
sgs_manifest_set_resolved(cfg$paths$manifest_path, "resolution_basis", sel$basis)
sgs_manifest_append_step(
  cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat")),
                          clustree = as.character(packageVersion("clustree"))),
  quality_gate = gate,
  resolved_params = list(resolution = sel$resolution, basis = sel$basis,
                         pc_used = pc_use,
                         ari_by_resolution = as.list(ari_seq)),
  outputs = c(list(out_rds, file.path(metrics_dir, "ari_by_resolution.csv")), as.list(figs)),
  started_at = started_at, finished_at = finished_at)
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S4 complete (awaiting human review of clustree)")

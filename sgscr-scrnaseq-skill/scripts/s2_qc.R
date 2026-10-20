# s2_qc.R — QC 过滤（percent.mt / percent.HB 等，阈值全部来自 Config）
# Validates: Requirements 6.2, 8.4, 20.2
# 输入：s1_obj_int.rds；输出：rds/s2_obj_qc.rds + figures/P2.1_QC_violin.*
# Quality_Gate：自动=过滤后细胞数 > qc.min_cells_pass_gate；
#               人工=复查 QC 小提琴图 → 状态 NEEDS_HUMAN_REVIEW（确认后 --from s3 续跑）

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s2_qc.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s2_qc"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

in_rds <- file.path(cfg$paths$rds_dir, "s1_obj_int.rds")
if (!file.exists(in_rds)) {
  log_error(STEP, paste("missing prerequisite:", in_rds)); quit(status = 1)
}
obj <- readRDS(in_rds)
n_before <- ncol(obj)

# ---- QC 指标与过滤（阈值全部来自 Config，需求 6.2） ----
# percent.mt/HB 必须算在 RNA assay 的 counts 上（integrated assay 无 counts 层）
# P12-a A3：线粒体前缀随物种（human "^MT-" / mouse "^mt-"），config qc.mt_pattern 可覆盖
mt_pattern <- if (!is.null(cfg$qc$mt_pattern)) cfg$qc$mt_pattern else {
  switch(cfg$species$organism, human = "^MT-", mouse = "^mt-",
         stop(sprintf("qc.mt_pattern 未配置且物种 %s 无默认线粒体前缀", cfg$species$organism), call. = FALSE))
}
Seurat::DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
obj[["percent.HB"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^HB[^(P)]")

obj <- subset(obj,
  subset = nFeature_RNA > cfg$qc$min_features &
           nFeature_RNA < cfg$qc$max_features &
           percent.mt < cfg$qc$max_percent_mt &
           percent.HB < cfg$qc$max_percent_hb)
n_after <- ncol(obj)
log_info(STEP, sprintf("QC: %d -> %d cells", n_before, n_after))

# ---- QC 小提琴图 ----
plt <- Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA",
                                         "percent.mt", "percent.HB"),
                       group.by = "sample_id", ncol = 2) &
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
figs <- sgs_save_plot(cfg, "P2.1", "PanelA", "QC_violin", plt, STEP)

# ---- Quality_Gate：自动判定 + 人工复查标记 ----
cells_ok <- n_after > cfg$qc$min_cells_pass_gate
auto_status <- if (cells_ok) "PASS" else "FAIL"
detail <- sprintf("cells %d -> %d (gate>%d); thresholds mt<%s hb<%s feats(%s,%s)",
                  n_before, n_after, cfg$qc$min_cells_pass_gate,
                  cfg$qc$max_percent_mt, cfg$qc$max_percent_hb,
                  cfg$qc$min_features, cfg$qc$max_features)
status <- if (auto_status == "PASS") "NEEDS_HUMAN_REVIEW" else "FAIL"
detail <- paste0(detail, "; human review required: check P2.1_QC_violin, then rerun with --from s3")
log_info(STEP, paste("Quality_Gate:", status, "-", detail))

finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
out_rds <- file.path(cfg$paths$rds_dir, "s2_obj_qc.rds")
saveRDS(obj, out_rds)
sgs_manifest_append_step(
  cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = list(status = status, detail = detail),
  resolved_params = list(n_before = n_before, n_after = n_after,
                         qc_thresholds = list(min_features = cfg$qc$min_features,
                                              max_features = cfg$qc$max_features,
                                              max_percent_mt = cfg$qc$max_percent_mt,
                                              max_percent_hb = cfg$qc$max_percent_hb)),
  outputs = c(list(out_rds), as.list(figs)),
  started_at = started_at, finished_at = finished_at)
if (status == "FAIL") quit(status = 1)   # NEEDS_HUMAN_REVIEW 正常退出，等人工确认
log_info(STEP, "S2 complete (awaiting human review of QC violin)")

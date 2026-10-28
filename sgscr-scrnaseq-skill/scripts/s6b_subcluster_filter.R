# s6b_subcluster_filter.R — S6b 亚群注释与保留过滤（P11D-c 派工项 1，需求 12 扩展）
# Validates: 0018 派工：读 compartments.subcluster_labels CSV，按 filter_mode 过滤 s6 子对象
#
# 输入：rds/s6_<comp>_subcluster.rds（S6 产物）+ 标签 CSV（compartment,subcluster,keep,proposed_label）
# 校验：compartment/subcluster 与 s6 对象一一对应（缺失或多余即 FAIL）；
#       亚群列用 sprintf("%s_snn_res.%g", assay, res) 精确匹配（禁止 grep 首列——父 S4 扫描列在前的实测坑）
# 输出：rds/s6b_<comp>_subcluster.rds（含 subcluster_label 元数据列）
# Manifest：标签 CSV 的 MD5、filter_mode、各亚群保留数、各样本保留数
# Quality_Gate：校验通过且各 comp 过滤后细胞数 > 0 且与保留亚群细胞数之和一致

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(digest))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s6b_subcluster_filter.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R",
            "subcluster_filter_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s6b_subcluster_filter"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

csv_path <- cfg$compartments$subcluster_labels
mode <- if (is.null(cfg$compartments$filter_mode)) "main" else cfg$compartments$filter_mode
targets <- as.character(cfg$compartments$targets)
res_sub <- cfg$compartments$subcluster_resolution
res_col <- sprintf("%s_snn_res.%g", "RNA", res_sub)   # 精确构造，禁止 grep 首列（P11C 实测坑）

# ---- 前置校验（先于一切计算） ----
problems <- c()
if (is.null(csv_path) || !file.exists(csv_path)) problems <- c(problems, sprintf("标签 CSV: %s", csv_path))
for (comp in targets) {
  f <- file.path(cfg$paths$rds_dir, sprintf("s6_%s_subcluster.rds", comp))
  if (!file.exists(f)) problems <- c(problems, basename(f))
}
if (!mode %in% c("main", "strict")) problems <- c(problems, sprintf("filter_mode 非法: %s", mode))
if (length(problems) > 0) sgs_preflight(STEP, problems, character(0))

csv_df <- read.csv(csv_path, stringsAsFactors = FALSE)
csv_md5 <- digest::digest(file = csv_path, algo = "md5")

# ---- 逐 comp：读对象 → 校验 → 过滤 → 落盘 ----
objs <- list(); ids_by_comp <- list()
for (comp in targets) {
  o <- readRDS(file.path(cfg$paths$rds_dir, sprintf("s6_%s_subcluster.rds", comp)))
  if (!res_col %in% colnames(o@meta.data)) {
    log_error(STEP, sprintf("%s: 对象无精确亚群列 %s（保留列: %s）", comp, res_col,
                            paste(grep("snn_res", colnames(o@meta.data), value = TRUE), collapse = ",")))
    quit(status = 1)
  }
  objs[[comp]] <- o
  ids_by_comp[[comp]] <- as.character(sort(as.integer(unique(as.character(o@meta.data[[res_col]])))))
}

val <- sgs_subcluster_validate(csv_df, targets, ids_by_comp)
if (!val$ok) {
  log_error(STEP, paste("标签 CSV 校验 FAIL:", paste(val$problems, collapse = "; ")))
  gate <- sgs_qc_gate(auto_pass = FALSE, detail = paste(val$problems, collapse = "; "))
  sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
    package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
    quality_gate = gate, resolved_params = list(labels_csv = csv_path, filter_mode = mode),
    outputs = list(), started_at = started_at,
    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  quit(status = 1)
}

out_files <- list(); details <- character(0)
per_sub <- list(); per_sample <- list()
for (comp in targets) {
  keep_ids <- sgs_subcluster_keep_ids(val$keep_map[[comp]], mode)
  if (length(keep_ids) == 0) {
    log_error(STEP, sprintf("%s: mode=%s 下保留亚群为空", comp, mode)); quit(status = 1)
  }
  o <- objs[[comp]]
  n_before <- ncol(o)
  n_keep_expected <- sum(as.character(o@meta.data[[res_col]]) %in% keep_ids)
  filt <- sgs_filter_subobject(o, res_col, keep_ids, val$labels_map[[comp]])
  chk <- sgs_filter_count_check(filt, res_col, keep_ids, n_keep_expected)
  if (!chk$ok) {
    log_error(STEP, sprintf("%s: 过滤后细胞数 %d != 保留亚群细胞数之和 %d", comp, chk$n_actual, chk$n_expected))
    quit(status = 1)
  }
  out_rds <- file.path(cfg$paths$rds_dir, sprintf("s6b_%s_subcluster.rds", comp))
  saveRDS(filt, out_rds)
  out_files[[comp]] <- out_rds
  # Manifest 明细：各亚群保留数（含被剔除的，保留数记 0）
  all_ids <- ids_by_comp[[comp]]
  sub_tab <- table(factor(as.character(o@meta.data[[res_col]]), levels = all_ids))
  kept_flag <- all_ids %in% keep_ids
  per_sub[[comp]] <- stats::setNames(as.integer(ifelse(kept_flag, sub_tab, 0)), all_ids)
  per_sample[[comp]] <- as.list(table(factor(as.character(filt$orig.ident))))
  details <- c(details, sprintf("%s[%s]: %d/%d cells kept (%d/%d subclusters: %s)",
    comp, mode, ncol(filt), n_before, length(keep_ids), length(all_ids),
    paste(keep_ids, collapse = ",")))
  log_info(STEP, details[length(details)])
}

gate <- sgs_qc_gate(auto_pass = length(out_files) == length(targets),
                    detail = paste(details, collapse = "; "))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))
sgs_manifest_set_resolved(cfg$paths$manifest_path, "subcluster_filter_labels_md5", csv_md5)
sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(labels_csv = csv_path, labels_md5 = csv_md5, filter_mode = mode,
                         res_col = res_col,
                         kept_per_subcluster = per_sub, kept_per_sample = per_sample),
  outputs = as.list(unlist(out_files)),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S6b complete")

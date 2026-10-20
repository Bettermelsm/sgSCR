# s5_annotation.R — 两层注释（SingleR 自动 + FindAllMarkers + 人工合并映射应用）
# Validates: Requirements 11.1-11.4, 5.3-5.4 (P5), P10
#
# Phase A（默认）: 物种校验(P5) → SingleR 自动注释 → FindAllMarkers →
#   生成 annotation/cluster_mapping_template.csv → Gate=NEEDS_HUMAN_REVIEW（等待人工填写）
# Phase B（--mapping <confirmed.csv>）: 校验覆盖全部簇 → 应用映射写 obj$major_cluster
#   （保留 obj$SingleR_label_raw 不覆盖）→ rds/s5_obj_annotated.rds → Gate=PASS
#
# 用法:
#   Rscript s5_annotation.R --config <yaml>                    # Phase A
#   Rscript s5_annotation.R --config <yaml> --mapping <csv>    # Phase B（人工确认后）

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleR))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
map_idx <- grep("^--mapping$", args)
map_path <- if (length(map_idx) > 0) args[map_idx + 1] else NA_character_
if (is.na(config_path)) stop("Usage: Rscript s5_annotation.R --config <yaml> [--mapping <csv>]", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "annotation_utils.R", "assay_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s5_annotation"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

pre <- tryCatch(sgs_resume_plan("s5", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

obj <- readRDS(file.path(cfg$paths$rds_dir, "s4_obj_cluster.rds"))
obj <- SeuratObject::JoinLayers(obj, assay = "RNA")   # S1 按样本拆层 → 合并后才能取单层
cluster_col <- "selected_resolution"
cluster_ids <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
phaseA_rds <- file.path(cfg$paths$rds_dir, "s5_phaseA_annot.rds")
anno_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "annotation")
dir.create(anno_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Phase B: 应用人工确认的映射（P10） ----
if (!is.na(map_path)) {
  phaseA_rds <- file.path(cfg$paths$rds_dir, "s5_phaseA_annot.rds")
  if (!file.exists(phaseA_rds)) {
    log_error(STEP, sprintf("未找到 Phase A 产物 %s — 请先不带 --mapping 运行 S5", phaseA_rds))
    quit(status = 1)
  }
  obj <- readRDS(phaseA_rds)
  map_df <- read.csv(map_path, stringsAsFactors = FALSE, colClasses = "character")
  mapping <- tryCatch(sgs_validate_confirmed_mapping(map_df, cluster_ids),
                      error = function(e) e)
  if (inherits(mapping, "error")) { log_error(STEP, conditionMessage(mapping)); quit(status = 1) }
  raw_before <- obj$SingleR_label_raw
  obj$major_cluster <- sgs_apply_cluster_mapping(as.character(obj@meta.data[[cluster_col]]), mapping)
  stopifnot(identical(obj$SingleR_label_raw, raw_before))   # P10: 原始注释不可变
  obj <- sgs_drop_scale_layers(obj, assay = "RNA")   # P11B 复用 Phase A 对象时兜底清理
  sgs_assert_layers_clean(obj, assay = "RNA")
  out_rds <- file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds")
  saveRDS(obj, out_rds)
  sgs_manifest_set_resolved(cfg$paths$manifest_path, "cluster_to_major_cluster", as.list(mapping))
  gate <- sgs_qc_gate(auto_pass = TRUE,
                      detail = sprintf("mapping applied: %d clusters -> %d major clusters",
                                       length(mapping), length(unique(obj$major_cluster))))
  sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
    package_versions = list(Seurat = as.character(packageVersion("Seurat")),
                            SingleR = as.character(packageVersion("SingleR"))),
    quality_gate = gate, resolved_params = list(mapping = as.list(mapping)),
    outputs = as.list(out_rds), started_at = started_at,
    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))
  log_info(STEP, "S5 complete")
  quit(status = 0)
}

# ---- Phase A: SingleR 自动注释（P5 物种校验在参考集构建处强制） ----
ref_name <- cfg$species$singleR_reference
ref_species <- cfg$species$reference_species
sgs_check_species(cfg$species$organism, ref_species)   # P5: 不一致直接终止

if (identical(ref_name, "mock")) {
  # 烟雾测试用 mock 参考：按簇聚合表达构造 (genes × types) 矩阵，直接以矩阵作参考
  ids <- as.character(obj@meta.data[[cluster_col]])
  mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  ref_mat <- sapply(cluster_ids, function(cl) rowMeans(mat[, ids == cl, drop = FALSE]))
  ref <- ref_mat
  ref_labels <- paste0("MockType_", colnames(ref_mat))
  ref_lab_vec <- ref_labels
} else {
  if (!requireNamespace("celldex", quietly = TRUE)) {
    stop("celldex 未安装，无法加载参考集 ", ref_name, call. = FALSE)
  }
  ref <- get(ref_name, envir = asNamespace("celldex"))()   # 如 HumanPrimaryCellAtlasData
  ref_labels <- ref$label.main
  ref_lab_vec <- NULL
}

pred <- SingleR(test = GetAssayData(obj, assay = "RNA", layer = "data"),
                ref = ref, labels = ref_labels)
obj$SingleR_label_raw <- as.character(pred$labels)
# 残留 layer 清理提前到 Phase A 落盘之前（P11A）：phaseA_rds 是 Phase B 与下游 s6/s7 的
# 共同输入，落盘时就必须是干净的 counts/data，而不是依赖下游各自清理
obj <- sgs_drop_scale_layers(obj, assay = "RNA")
sgs_assert_layers_clean(obj, assay = "RNA")   # 落盘前断言
write.csv(data.frame(cell_id = colnames(obj), SingleR_label_raw = obj$SingleR_label_raw,
                     pruned_score = pred$pruned.scores[, 1] %||% NA),
          file.path(anno_dir, "singleR_pred.csv"), row.names = FALSE)
saveRDS(obj, phaseA_rds)   # Phase B 从此对象应用映射

# ---- FindAllMarkers（需求 11.2） ----
Idents(obj) <- cluster_col
markers <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
markers_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "markers")
dir.create(markers_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(markers, file.path(markers_dir, "Markers_all_clusters.csv"), row.names = FALSE)
# Quality_Gate 新检查项：marker 覆盖簇数必须等于实际簇数，缺失即 FAIL 并点名
cov <- sgs_check_marker_coverage(markers, cluster_ids)

# ---- 人工合并模板（需求 11.3；人工填写 confirmed_major_cluster 列） ----
top_n <- 10
top_markers <- vapply(split(markers, markers$cluster), function(m) {
  paste(head(m$gene[order(m$avg_log2FC, decreasing = TRUE)], top_n), collapse = ";")
}, character(1))
maj <- vapply(cluster_ids, function(cl) {
  sgs_majority_label(obj$SingleR_label_raw[obj@meta.data[[cluster_col]] == cl])
}, character(1))
template <- data.frame(
  cluster_id = cluster_ids,
  n_cells = as.integer(table(obj@meta.data[[cluster_col]])[cluster_ids]),
  top_markers = unname(top_markers[cluster_ids]),
  singleR_majority_label = unname(maj[cluster_ids]),
  suggested_major_cluster = "",
  confirmed_major_cluster = "", stringsAsFactors = FALSE)
template_path <- file.path(anno_dir, "cluster_mapping_template.csv")
write.csv(template, template_path, row.names = FALSE)

# ---- Quality_Gate: 人工确认点（设计明确要求停下等待，不自行填写） ----
gate <- sgs_qc_gate(
  auto_pass = nrow(markers) > 0 && length(cluster_ids) > 0 && cov$ok,
  detail = sprintf("%d clusters; marker 覆盖 %d/%d%s; reference=%s (%s)",
                   length(cluster_ids), cov$n_covered, cov$n_expected,
                   if (!cov$ok) sprintf("（缺: %s）", paste(cov$missing, collapse = ",")) else "",
                   ref_name, ref_species),
  human_review_file = template_path)
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat")),
                          SingleR = as.character(packageVersion("SingleR"))),
  quality_gate = gate,
  resolved_params = list(reference = ref_name, reference_species = ref_species),
  outputs = list(template_path, file.path(anno_dir, "singleR_pred.csv"),
                 file.path(markers_dir, "Markers_all_clusters.csv"), phaseA_rds),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S5 Phase A complete — fill confirmed_major_cluster, then rerun with --mapping")

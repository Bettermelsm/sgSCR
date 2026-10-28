# s1_import_integrate.R — 批量 Read10X 导入 + FindIntegrationAnchors/IntegrateData
# Validates: Requirements 3.4, 20.1
# 输入：S0 通过；Config samples.manifest（跳过 excluded）指向的 10X 三联文件
# 输出：{paths.rds_dir}/s1_obj_int.rds；Manifest 追加
# Quality_Gate（自动）：细胞总数 > qc.min_cells_pass_gate 且 整合后 assay 为 "integrated"

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s1_import_integrate.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R",
            "triplet_check.R", "save_plot.R")) source(file.path(HERE, "common", f))

STEP <- "s1_import_integrate"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

dir.create(cfg$paths$rds_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 前置：S0 产物存在（缺则终止，需求 7.4 的最简形式） ----
if (!file.exists(cfg$paths$manifest_path)) {
  log_error(STEP, paste("S0 manifest missing, run S0 first:", cfg$paths$manifest_path))
  quit(status = 1)
}

# ---- 批量导入（excluded 自动跳过） ----
excluded_ids <- vapply(cfg$samples$excluded, function(e) e$sample_id, character(1))
obj_list <- list()
for (s in cfg$samples$manifest) {
  sid <- s$sample_id
  if (sid %in% excluded_ids) {
    log_info(STEP, sprintf("skip excluded sample: %s", sid)); next
  }
  d <- file.path(cfg$remote$data_dir, sid)
  # features.tsv 可能只有 1 列（纯基因符号，如 PISA 生成的矩阵），此时须 gene.column=1
  feat <- list.files(d, pattern = "^features\\.tsv(\\.gz)?$")
  n_cols <- if (length(feat) == 1) {
    con <- gzfile(file.path(d, feat[1]), "rt")
    n <- max(vapply(head(readLines(con, n = 20), 1),
                    function(l) length(strsplit(l, "\t", fixed = TRUE)[[1]]), integer(1)))
    close(con); n
  } else 2L
  counts <- Seurat::Read10X(data.dir = d, gene.column = min(2L, n_cols))
  # 注：不用 min.features —— Seurat 5.3.0 该参数触发 LogMap "Rownames must be
  # supplied" bug（2026-10-02 实测）；特征过滤统一放 S2 的 subset 完成
  obj_list[[sid]] <- Seurat::CreateSeuratObject(
    counts, project = sid,
    min.cells = cfg$qc$min_cells_per_gene)
  obj_list[[sid]]$sample_id <- sid
  obj_list[[sid]]$group <- s$group
  log_info(STEP, sprintf("%s: %d cells imported", sid, ncol(obj_list[[sid]])))
}

total_cells <- sum(vapply(obj_list, ncol, numeric(1)))
log_info(STEP, sprintf("total cells before integration: %d", total_cells))

# ---- 整合前置：每样本 NormalizeData / FindVariableFeatures / ScaleData ----
# Seurat v5 FindIntegrationAnchors 要求每个对象先完成这三步（2026-10-02 实测补齐）
for (sid in names(obj_list)) {
  obj_list[[sid]] <- Seurat::NormalizeData(obj_list[[sid]], verbose = FALSE)
  obj_list[[sid]] <- Seurat::FindVariableFeatures(
    obj_list[[sid]], selection.method = "vst",
    nfeatures = cfg$integration$nfeatures_hvg, verbose = FALSE)
  obj_list[[sid]] <- Seurat::ScaleData(obj_list[[sid]], verbose = FALSE)
}

# ---- 批次整合（需求 20.1 保留能力） ----
anchors <- Seurat::FindIntegrationAnchors(
  object.list = obj_list,
  dims = 1:cfg$integration$dims,
  reduction = cfg$integration$reduction)
obj <- Seurat::IntegrateData(
  anchorset = anchors,
  dims = 1:cfg$integration$dims)
rm(obj_list, anchors); invisible(gc())

Seurat::DefaultAssay(obj) <- "integrated"
obj[["RNA"]] <- obj[["RNA"]]

# ---- Quality_Gate ----
assay_ok <- identical(Seurat::DefaultAssay(obj), "integrated")
cells_ok <- total_cells > cfg$qc$min_cells_pass_gate
status <- if (assay_ok && cells_ok) "PASS" else "FAIL"
detail <- sprintf("cells=%d (gate>%d); assay=%s", total_cells,
                  cfg$qc$min_cells_pass_gate, Seurat::DefaultAssay(obj))
log_info(STEP, paste("Quality_Gate:", status, "-", detail))

# ---- 落盘 ----
finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
out_rds <- file.path(cfg$paths$rds_dir, "s1_obj_int.rds")
saveRDS(obj, out_rds)
sgs_manifest_append_step(
  cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = list(status = status, detail = detail),
  resolved_params = list(total_cells = total_cells,
                         integration = list(dims = cfg$integration$dims,
                                            reduction = cfg$integration$reduction)),
  outputs = list(out_rds),
  started_at = started_at, finished_at = finished_at)
if (status == "FAIL") quit(status = 1)
log_info(STEP, "S1 complete")

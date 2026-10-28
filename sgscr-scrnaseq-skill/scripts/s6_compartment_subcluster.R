# s6_compartment_subcluster.R — Compartment 亚聚类（需求 12）
# Validates: Requirements 12.1-12.3
#
# 对 Config compartments.targets 中每个 Compartment：
#   从 s5_obj_annotated.rds 按 major_cluster 子集化 → 独立 PCA → SNN →
#   subcluster_resolution 聚类 → rds/s6_<compartment>_subcluster.rds
#   + figures/P6_<compartment>_UMAP.*
# Quality_Gate（自动）：各目标 Compartment 均产出子对象且细胞数 > 0，否则 FAIL

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s6_compartment_subcluster.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R",
            "subcluster_filter_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s6_compartment_subcluster"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

pre <- tryCatch(sgs_resume_plan("s6", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

# 产物新鲜度（P11 增强）：已有 s6 亚聚类产物不得早于 s5（如 s5 换映射表重跑后 s6 必须跟跑）
s5_rds <- file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds")
stale_s6 <- list.files(cfg$paths$rds_dir, pattern = "^s6_.*_subcluster\\.rds$", full.names = TRUE)
for (f in stale_s6) sgs_check_freshness(f, s5_rds)

obj <- readRDS(s5_rds)
obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
targets <- as.character(cfg$compartments$targets)
res_sub <- cfg$compartments$subcluster_resolution
# P11A：子聚类 PC 上限不再硬编码 10（anti_patterns 第 1 条：阈值一律进 Config）
sub_npcs <- if (is.null(cfg$compartments$subcluster_npcs)) 10L else as.integer(cfg$compartments$subcluster_npcs)
if (is.null(targets) || length(targets) == 0) {
  log_error(STEP, "compartments.targets 为空 — Config 必须以可配置参数指定亚聚类 Compartment 列表（需求 12.3）")
  quit(status = 1)
}

out_files <- list()
gate_problems <- character(0)
for (comp in targets) {
  n_cells <- sum(as.character(obj$major_cluster) == comp)
  if (n_cells == 0) {
    gate_problems <- c(gate_problems, sprintf("%s: 0 cells (major_cluster 无此标签)", comp))
    log_error(STEP, gate_problems[length(gate_problems)])
    next
  }
  sub <- subset(obj, subset = major_cluster == comp)
  log_info(STEP, sprintf("%s: %d cells, subclustering at res=%s", comp, n_cells, res_sub))
  sub <- Seurat::NormalizeData(sub, verbose = FALSE)
  sub <- Seurat::FindVariableFeatures(sub, verbose = FALSE)
  sub <- Seurat::ScaleData(sub, verbose = FALSE)
  npcs_use <- min(sub_npcs, n_cells - 1)
  log_info(STEP, sprintf("%s: npc=%d (config 上限 %d, n_cells-1=%d)", comp, npcs_use, sub_npcs, n_cells - 1))
  sub <- Seurat::RunPCA(sub, npcs = npcs_use, verbose = FALSE)
  sub <- Seurat::FindNeighbors(sub, dims = 1:npcs_use, verbose = FALSE)
  sub <- Seurat::FindClusters(sub, resolution = res_sub, verbose = FALSE)
  sub <- Seurat::RunUMAP(sub, dims = 1:npcs_use, verbose = FALSE)

  # P12-a D5：删除继承自父对象的 *_snn_res.* 扫描列，只保留本分辨率的亚群列
  # （防消费方 grep 首列误取——P11C 实测坑）
  sub <- sgs_drop_inherited_snn_cols(sub, sprintf("RNA_snn_res.%g", res_sub))

  # P6_<compartment>_UMAP（需求产物命名约定）
  umap_df <- as.data.frame(Embeddings(sub, "umap"))
  umap_df$subcluster <- as.character(sub@meta.data[[sprintf("RNA_snn_res.%g", res_sub)]])
  plt <- ggplot(umap_df, aes(umap_1, umap_2, color = subcluster)) +
    geom_point(size = 0.6) + theme_classic() +
    ggtitle(sprintf("%s subclusters (n=%d)", comp, n_cells))
  # 阶段15：源数据 = 逐细胞 UMAP + 亚群 + 样本 + 各(亚群,样本)细胞数
  umap_src <- data.frame(
    subcluster = umap_df$subcluster,
    sample_id = as.character(sub$sample_id),
    umap_1 = umap_df$umap_1, umap_2 = umap_df$umap_2,
    stringsAsFactors = FALSE)
  umap_src$n_cells <- as.integer(table(paste(umap_src$subcluster, umap_src$sample_id))[
    paste(umap_src$subcluster, umap_src$sample_id)])
  figres <- sgs_save_plot(cfg, "P6", "PanelA", paste0(comp, "_subcluster_UMAP"), plt, umap_src, STEP)
  figs <- unlist(figres$files)

  out_rds <- file.path(cfg$paths$rds_dir, sprintf("s6_%s_subcluster.rds", comp))
  saveRDS(sub, out_rds)
  out_files[[comp]] <- c(list(out_rds), as.list(figs))
  sgs_manifest_set_resolved(cfg$paths$manifest_path,
                            sprintf("subcluster_cells_%s", comp), n_cells)
  sgs_manifest_set_resolved(cfg$paths$manifest_path,
                            sprintf("subcluster_resolution_%s", comp),
                            length(unique(sub@meta.data[[sprintf("RNA_snn_res.%g", res_sub)]])))
}

# ---- Quality_Gate（自动）：每个目标 Compartment 均产出子对象且细胞数 > 0 ----
gate <- sgs_qc_gate(
  auto_pass = length(out_files) == length(targets),
  detail = if (length(gate_problems)) paste(gate_problems, collapse = "; ")
           else sprintf("%d compartments subclustered: %s", length(out_files),
                        paste(names(out_files), collapse = ", ")))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(targets = as.list(targets), subcluster_resolution = res_sub,
                         subcluster_npcs = sub_npcs),
  outputs = as.list(unlist(out_files)),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S6 complete")

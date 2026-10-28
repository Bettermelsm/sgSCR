# s3_dimreduce.R — 归一化/HVG/Scale/PCA + PART4 定量选 PC 数（P8）
# Validates: Requirements 9.1-9.3, 20.3
# 输入：s2_obj_qc.rds（人工已确认 QC）；输出：rds/s3_obj_pca.rds + Elbow/PCA 图
# Manifest 记录选定 PC 数与选择依据（需求 9.3 / P13 resolutions）

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
from_arg <- args[match("--from", args) + 1]   # 预留：--from s3 时校验 S2 已确认

if (is.na(config_path)) stop("Usage: Rscript s3_dimreduce.R --config <config.yaml> [--from s3]", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R",
            "save_plot.R", "pc_selection.R")) source(file.path(HERE, "common", f))

STEP <- "s3_dimreduce"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

in_rds <- file.path(cfg$paths$rds_dir, "s2_obj_qc.rds")
if (!file.exists(in_rds)) {
  log_error(STEP, paste("missing prerequisite:", in_rds)); quit(status = 1)
}
obj <- readRDS(in_rds)

# ---- 归一化 / HVG / Scale / PCA（需求 20.3 保留能力，参数来自 Config） ----
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                    nfeatures = cfg$integration$nfeatures_hvg)
obj <- Seurat::ScaleData(obj, vars.to.regress = "percent.mt")
obj <- Seurat::RunPCA(obj, npcs = cfg$pca$npcs_compute, approx = FALSE)

# ---- PART4 定量选 PC 数（P8） ----
stdev <- obj@reductions$pca@stdev
var_pct <- stdev^2 / sum(stdev^2) * 100
sel <- sgs_select_pcs(var_pct,
                      cum_th = cfg$pca$cum_variance_threshold_pct,
                      single_th = cfg$pca$single_pc_contrib_threshold_pct,
                      elbow_diff_threshold = cfg$pca$elbow_diff_threshold,
                      fallback_npc = cfg$pca$fallback_npc)
# 人工覆盖通道（阶段 11 Reviewer 指示）：pca.manual_pc_override 为正整数时以人工值替换
# 自动算法结果；自动算法建议值保留在 Manifest 的 auto 字段供对比追溯（basis 如实标 manual_override）
sel <- sgs_apply_pc_override(sel, cfg$pca$manual_pc_override, cfg$pca$npcs_compute)
log_info(STEP, sprintf("PC selection: pc.use=%d (basis=%s, co1=%s, co2=%d%s)",
                       sel$pc_use, sel$basis, sel$co1, sel$co2,
                       if (!is.null(sel$auto)) sprintf(
                         "; 自动算法建议 pc.use=%d (basis=%s) 未采用", sel$auto$pc_use,
                         sel$auto$basis) else ""))

# ---- 图 ----
elbow_df <- data.frame(PC = seq_along(var_pct), var_exp = var_pct,
                       cum_exp = cumsum(var_pct))
plt_elbow <- ggplot2::ggplot(elbow_df, ggplot2::aes(PC, var_exp)) +
  ggplot2::geom_col(fill = "grey70") +
  ggplot2::geom_line(ggplot2::aes(y = cum_exp / max(cum_exp) * max(var_exp)),
                     color = "steelblue") +
  ggplot2::geom_vline(xintercept = sel$pc_use, linetype = "dashed", color = "red") +
  ggplot2::annotate("text", x = sel$pc_use + 1, y = max(elbow_df$var_exp) * 0.9,
                    label = sprintf("pc.use=%d (%s)", sel$pc_use, sel$basis),
                    hjust = 0, color = "red") +
  ggplot2::labs(title = "Elbow plot with selected PC number",
                x = "PC", y = "% variance explained (bars) / cumulative (line)")
# 阶段15：源数据（PC 方差表；非样本级图附全量细胞数与全样本标识，逐字面满足契约）
elbow_src <- elbow_df
elbow_src$sample_id <- "(all)"
elbow_src$n_cells <- ncol(obj)
figres <- sgs_save_plot(cfg, "P3.2", "PanelA", "ElbowPlot_annotated", plt_elbow, elbow_src, STEP)
figs <- unlist(figres$files)

plt_pca <- Seurat::DimPlot(obj, reduction = "pca", group.by = "sample_id")
pca_emb <- Seurat::Embeddings(obj, reduction = "pca")[, 1:2, drop = FALSE]
pca_src <- data.frame(sample_id = as.character(obj$sample_id),
                      PC1 = pca_emb[, 1], PC2 = pca_emb[, 2], stringsAsFactors = FALSE)
pca_src$n_cells <- as.integer(table(pca_src$sample_id)[pca_src$sample_id])
figres2 <- sgs_save_plot(cfg, "P3.1", "PanelA", "PCA_by_sample", plt_pca, pca_src, STEP)
figs <- c(figs, unlist(figres2$files))

# ---- Quality_Gate ----
pcs_ok <- !is.na(sel$pc_use) && sel$pc_use >= 2 && sel$pc_use <= cfg$pca$npcs_compute
status <- if (pcs_ok) "PASS" else "FAIL"
detail <- sprintf("pc.use=%d within [2,%d]; rule=%s; co1=%s co2=%d",
                  sel$pc_use, cfg$pca$npcs_compute, sel$basis, sel$co1, sel$co2)
log_info(STEP, paste("Quality_Gate:", status, "-", detail))

finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
out_rds <- file.path(cfg$paths$rds_dir, "s3_obj_pca.rds")
saveRDS(obj, out_rds)

# P13：决议值回写 Manifest（resolutions + step resolved_params 双写）
sgs_manifest_set_resolved(cfg$paths$manifest_path, "pc_used", sel$pc_use)
sgs_manifest_set_resolved(cfg$paths$manifest_path, "pc_selection_basis", sel$basis)
sgs_manifest_append_step(
  cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = list(status = status, detail = detail),
  resolved_params = list(pc_used = sel$pc_use,
                         pc_selection_basis = sel$basis,
                         auto_pc_suggestion = if (!is.null(sel$auto)) sel$auto else NULL,
                         manual_override = !is.null(sel$auto),
                         selection_rule = sprintf("co1=cumvar>=%s&single<%s; co2=elbow; pc.use=min",
                                                  cfg$pca$cum_variance_threshold_pct,
                                                  cfg$pca$single_pc_contrib_threshold_pct)),
  outputs = c(list(out_rds), as.list(figs)),
  started_at = started_at, finished_at = finished_at)
if (status == "FAIL") quit(status = 1)
log_info(STEP, "S3 complete")

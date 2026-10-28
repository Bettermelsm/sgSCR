# export_scenic_input.R — S9 输入导出（主环境执行，供 pySCENIC 消费）
# 从 s5_obj_annotated.rds 取 scenic.target_subset 对应子集，
# 输出表达矩阵 CSV（genes × cells）与细胞元数据 CSV。

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
outdir <- args[grep("^--outdir$", args) + 1]
if (is.na(config_path) || is.na(outdir)) {
  stop("Usage: Rscript export_scenic_input.R --config <yaml> --outdir <dir>", call. = FALSE)
}

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "preflight_utils.R")) source(file.path(HERE, f))

suppressPackageStartupMessages(library(Seurat))
cfg <- sgs_load_config(config_path)
obj <- readRDS(file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds"))
obj <- SeuratObject::JoinLayers(obj, assay = "RNA")

target <- cfg$scenic$target_subset
# P11D-c：target_subset 命中 compartments.targets 且存在 S6b 清洗产物时，直接用清洗后子对象
s6b_f <- if (!is.null(target) && target %in% as.character(cfg$compartments$targets)) {
  file.path(cfg$paths$rds_dir, sprintf("s6b_%s_subcluster.rds", target))
} else NULL
if (!is.null(s6b_f) && file.exists(s6b_f)) {
  obj <- readRDS(s6b_f)
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  log_info("s9_export", sprintf("subset %s via S6b: %d cells (清洗后)", target, ncol(obj)))
} else if (!is.null(target) && target %in% unique(as.character(obj$major_cluster))) {
  obj <- subset(obj, subset = major_cluster == target)
  log_info("s9_export", sprintf("subset %s: %d cells", target, ncol(obj)))
} else {
  log_info("s9_export", sprintf("target_subset '%s' 不在 major_cluster，用全量 %d cells",
                                target, ncol(obj)))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# P12-a D7：freshness——已有 SCENIC 链产物不得早于实际输入（s5 或 S6b 清洗产物），
# 否则 s9 复用旧 adjacencies/regulons 会静默用陈旧输入
s9_ref <- if (exists("s6b_f") && !is.null(s6b_f) && file.exists(s6b_f)) s6b_f else {
  p <- file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds"); if (file.exists(p)) p else NULL
}
if (!is.null(s9_ref)) {
  for (f in c("expr.csv", "adjacencies.tsv", "regulons.csv", "auc_matrix.csv", "regulon_activity.csv")) {
    fp <- file.path(outdir, f)
    if (file.exists(fp)) sgs_check_freshness(fp, s9_ref)
  }
}
mat <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
write.csv(data.frame(gene = rownames(mat), mat, check.names = FALSE),
          file.path(outdir, "expr.csv"), row.names = FALSE)
write.csv(as.data.frame(obj@meta.data[, c("major_cluster", "group"), drop = FALSE]),
          file.path(outdir, "cell_meta.csv"))
log_info("s9_export", sprintf("exported %d genes x %d cells", nrow(mat), ncol(mat)))

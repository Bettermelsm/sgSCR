# s10_target_genes.R — 目标基因追踪（需求 16，主环境）
# Validates: Requirements 16.1-16.3 (P11, P12)
#
# 输入：rds/s5_obj_annotated.rds + config target_genes（基因清单 + science_questions）
# P12 过滤：目标基因 ∩ 对象基因（纯函数 sgs_target_gene_filter）；交集为空 → 告警不终止
# 输出：target_genes/<genes>_expr.csv（各细胞表达 + 分组）、target_genes/target_gene_stats.csv
#       （各 major_cluster 均值 + UC/NC 差异）、figures/P10_target_gene_dotplot.*、
#       manifest 追加基因→科学问题对应表

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s10_target_genes.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R",
            "target_gene_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s10_target_genes"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

in_rds <- file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds")
targets <- as.character(cfg$target_genes$genes)
science_q <- cfg$target_genes$science_questions

# 断点续跑统一入口（P6）：与其他 Step 相同的 sgs_resume_plan 检查
pre <- tryCatch(sgs_resume_plan("s10", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

# ---- P11 前置校验 ----
problems <- c()
if (!file.exists(in_rds)) problems <- c(problems, "s5_obj_annotated.rds")
if (length(targets) == 0) problems <- c(problems, "target_genes.genes")
if (length(problems) > 0) sgs_preflight(STEP, problems, character(0))
# 产物新鲜度（P11 增强）：已有 s10 产物不得早于上游 s5 产物（如 s5 换映射表重跑后 s10 必须跟跑）
stale_out <- file.path(dirname(cfg$paths$manifest_path), "..", "target_genes",
                       "target_gene_stats.csv")
if (file.exists(stale_out)) sgs_check_freshness(stale_out, in_rds)

obj <- readRDS(in_rds)
obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
kept <- sgs_target_gene_filter(targets, rownames(obj))
missing <- setdiff(targets, kept)

if (length(kept) == 0) {
  # P12/需求16.2：交集空 → 告警 + 空结果表，非终止
  log_warn(STEP, sprintf("目标基因与对象基因交集为空（清单: %s）——输出空结果表",
                         paste(targets, collapse = ", ")))
  out_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "target_genes")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  empty_stats <- data.frame(gene = character(0), major_cluster = character(0),
                            mean_expr = numeric(0), mean_UC = numeric(0),
                            mean_NC = numeric(0), diff_UC_minus_NC = numeric(0),
                            p_value = numeric(0), stringsAsFactors = FALSE)
  write.csv(empty_stats, file.path(out_dir, "target_gene_stats.csv"), row.names = FALSE)
  gate <- sgs_qc_gate(auto_pass = TRUE,
                      detail = sprintf("WARN: 目标基因交集为空（%d 个均不在对象中）", length(targets)))
  sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
    package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
    quality_gate = gate,
    resolved_params = list(targets = targets, kept = character(0), warning = "基因交集为空"),
    outputs = list(file.path(out_dir, "target_gene_stats.csv")),
    started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  log_info(STEP, "S10 complete (empty intersection)")
  quit(status = 0)
}

# ---- 各 major_cluster 定位 + UC/NC 差异 ----
expr <- FetchData(obj, vars = kept, slot = "data")
meta <- obj@meta.data
group <- as.character(meta$group)
cluster <- as.character(meta$major_cluster)

stats <- do.call(rbind, lapply(kept, function(g) {
  v <- expr[[g]]
  uc <- v[group == "UC"]; nc <- v[group == "NC"]
  p <- tryCatch(wilcox.test(uc, nc)$p.value, error = function(e) NA_real_)
  data.frame(gene = g,
             major_cluster = paste(sort(unique(cluster[v > 0])), collapse = ";"),
             mean_expr = mean(v),
             mean_UC = if (length(uc)) mean(uc) else NA_real_,
             mean_NC = if (length(nc)) mean(nc) else NA_real_,
             diff_UC_minus_NC = if (length(uc) && length(nc)) mean(uc) - mean(nc) else NA_real_,
             p_value = p, stringsAsFactors = FALSE)
}))

out_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "target_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(stats, file.path(out_dir, "target_gene_stats.csv"), row.names = FALSE)

# 各细胞表达表（基因→科学问题随 manifest 记录，文件按基因命名）
expr_df <- cbind(data.frame(cell = colnames(obj), group = group, major_cluster = cluster),
                 expr[, kept, drop = FALSE])
expr_file <- file.path(out_dir, sprintf("%s_expr.csv", paste(kept, collapse = "_")))
write.csv(expr_df, expr_file, row.names = FALSE)

# ---- Dotplot（目标基因 × major_cluster） ----
plt <- DotPlot(obj, features = kept, group.by = "major_cluster") +
  coord_flip() + theme_classic() + RotatedAxis() +
  ggtitle(sprintf("Target genes (kept %d/%d)", length(kept), length(targets)))
figs <- sgs_save_plot(cfg, "P10", "PanelA", "target_gene_dotplot", plt, STEP)

gate <- sgs_qc_gate(
  auto_pass = TRUE,
  detail = sprintf("kept %d/%d 目标基因 (%s)%s；stats+%s 落盘",
                   length(kept), length(targets), paste(kept, collapse = ","),
                   if (length(missing)) sprintf("，缺失: %s", paste(missing, collapse = ",")) else "",
                   basename(expr_file)))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

# 基因→科学问题对应表记入 manifest（需求16.3）
sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(
    targets = targets, kept = kept,
    science_questions = as.list(stats::setNames(
      science_q[intersect(kept, names(science_q))],
      intersect(kept, names(science_q)))),
    missing_genes = as.list(missing)),
  outputs = as.list(c(file.path(out_dir, "target_gene_stats.csv"), expr_file, figs)),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
log_info(STEP, "S10 complete")

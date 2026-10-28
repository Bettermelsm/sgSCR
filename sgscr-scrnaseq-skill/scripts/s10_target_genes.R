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
            "target_gene_utils.R", "sample_stats_utils.R")) {
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
                            mean_expr = numeric(0), mean_case = numeric(0),
                            mean_ctrl = numeric(0), diff_case_minus_ctrl = numeric(0),
                            p_value = numeric(0), stringsAsFactors = FALSE)
  write.csv(empty_stats, file.path(out_dir, "target_gene_stats.csv"), row.names = FALSE)
  # 空交集分支也产出 s10b 依赖的占位与诊断空表（阶段13/16：诊断不因目标基因缺失而消失）
  diag_dir0 <- file.path(dirname(cfg$paths$manifest_path), "diagnostics")
  ms_dir0 <- file.path(dirname(cfg$paths$manifest_path), "module_score")
  dir.create(diag_dir0, recursive = TRUE, showWarnings = FALSE)
  dir.create(ms_dir0, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(lineage = character(0), gene = character(0), sample_id = character(0),
                       n_cells = integer(0), pct_detected = numeric(0)),
            file.path(diag_dir0, "ambient_rna.csv"), row.names = FALSE)
  write.csv(data.frame(group = character(0), cell_group = character(0),
                       dominant_sample = character(0), max_fraction = numeric(0),
                       n_cells = integer(0), single_sample_driven = logical(0)),
            file.path(diag_dir0, "sample_representation.csv"), row.names = FALSE)
  write.csv(data.frame(cell_group = character(0), sample_id = character(0), group = character(0),
                       n_cells = integer(0)),
            file.path(ms_dir0, "by_cluster_sample.csv"), row.names = FALSE)
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

# ---- 各 major_cluster 定位 + 组间差异（P12-a A1：组别从 config design 读取，代码零项目字面量） ----
grp_cfg <- if (!is.null(cfg$design$groups)) as.character(unlist(cfg$design$groups)) else c("ctrl", "case")
ctrl_grp <- grp_cfg[1]
case_grp <- if (length(grp_cfg) > 1) grp_cfg[2] else grp_cfg[1]
expr <- FetchData(obj, vars = kept, slot = "data")
meta <- obj@meta.data
group <- as.character(meta$group)
cluster <- as.character(meta$major_cluster)
samples <- as.character(meta$sample_id)

stats <- do.call(rbind, lapply(kept, function(g) {
  v <- expr[[g]]
  case_v <- v[group == case_grp]; ctrl_v <- v[group == ctrl_grp]
  p <- tryCatch(wilcox.test(case_v, ctrl_v)$p.value, error = function(e) NA_real_)
  data.frame(gene = g,
             major_cluster = paste(sort(unique(cluster[v > 0])), collapse = ";"),
             mean_expr = mean(v),
             mean_case = if (length(case_v)) mean(case_v) else NA_real_,
             mean_ctrl = if (length(ctrl_v)) mean(ctrl_v) else NA_real_,
             diff_case_minus_ctrl = if (length(case_v) && length(ctrl_v)) mean(case_v) - mean(ctrl_v) else NA_real_,
             p_value = p, stringsAsFactors = FALSE)
}))

# P12-a D1：per-major 拆分（需求 16.2）——基因 × major_cluster 的分组均值与差值
per_major <- do.call(rbind, lapply(kept, function(g) {
  v <- expr[[g]]
  do.call(rbind, lapply(sort(unique(cluster)), function(cl) {
    i <- cluster == cl
    case_v <- v[i & group == case_grp]; ctrl_v <- v[i & group == ctrl_grp]
    data.frame(gene = g, major_cluster = cl,
               mean_case = if (length(case_v)) mean(case_v) else NA_real_,
               mean_ctrl = if (length(ctrl_v)) mean(ctrl_v) else NA_real_,
               diff_case_minus_ctrl = if (length(case_v) && length(ctrl_v)) mean(case_v) - mean(ctrl_v) else NA_real_,
               n_case = length(case_v), n_ctrl = length(ctrl_v), stringsAsFactors = FALSE)
  }))
}))

out_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "target_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(stats, file.path(out_dir, "target_gene_stats.csv"), row.names = FALSE)
write.csv(per_major, file.path(out_dir, "target_gene_stats_per_major.csv"), row.names = FALSE)

# ---- 阶段13 固化（需求 24.1/24.2）：环境 RNA 诊断 ----
# ambient.lineage_markers: 谱系 -> 该谱系本不应表达的基因；按 (谱系,基因,样本) 输出检出比例
run_root10 <- normalizePath(file.path(dirname(cfg$paths$manifest_path), ".."))  # = {remote_work_dir}/run
diag_dir <- file.path(run_root10, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
ambient_warn_threshold <- if (!is.null(cfg$ambient$warn_threshold)) cfg$ambient$warn_threshold else 0.1
ambient_warnings <- list()
if (!is.null(cfg$ambient$lineage_markers)) {
  lms <- cfg$ambient$lineage_markers
  amb_rows <- list()
  for (lin in names(lms)) {
    i_lin <- cluster == lin
    if (!any(i_lin)) {
      log_warn(STEP, sprintf("ambient.lineage_markers 谱系 %s 在对象中无细胞，跳过", lin)); next
    }
    for (gn in as.character(unlist(lms[[lin]]))) {
      if (!(gn %in% rownames(obj))) next
      v <- FetchData(obj, gn, layer = "data")[[1]]
      for (s in sort(unique(samples[i_lin]))) {
        i <- i_lin & samples == s
        pct <- mean(v[i] > 0)
        amb_rows[[length(amb_rows) + 1]] <- data.frame(
          lineage = lin, gene = gn, sample_id = s,
          n_cells = sum(i), pct_detected = pct, stringsAsFactors = FALSE)
        if (pct > ambient_warn_threshold) {
          w <- sprintf("环境RNA告警: 谱系 %s 中 %s 检出 %.1f%%（样本 %s）> 阈值 %.1f%%",
                       lin, gn, 100 * pct, s, 100 * ambient_warn_threshold)
          log_warn(STEP, w)
          ambient_warnings[[length(ambient_warnings) + 1]] <- w
        }
      }
    }
  }
  if (length(amb_rows)) {
    write.csv(do.call(rbind, amb_rows), file.path(diag_dir, "ambient_rna.csv"), row.names = FALSE)
  }
} else {
  log_info(STEP, "ambient.lineage_markers 未配置，环境RNA诊断输出空（不告警）")
}

# ---- 阶段14 固化（需求 25.1/25.2）：样本代表性告警 ----
dom_thr <- if (!is.null(cfg$diagnostics$single_sample_dominance_threshold))
  cfg$diagnostics$single_sample_dominance_threshold else 0.70
repr_rows <- list(); repr_warnings <- list()
for (cl in sort(unique(cluster))) {
  for (gp in sort(unique(group))) {
    i <- cluster == cl & group == gp
    if (!any(i)) next
    cnt <- table(samples[i])
    dom <- sgs_single_sample_dominance(cnt, threshold = dom_thr)
    repr_rows[[length(repr_rows) + 1]] <- data.frame(
      group = gp, cell_group = cl, dominant_sample = dom$dominant_sample,
      max_fraction = dom$max_fraction, n_cells = sum(i),
      single_sample_driven = dom$driven, stringsAsFactors = FALSE)
    if (dom$driven) {
      w <- sprintf("样本代表性告警: (%s, %s) %.1f%% 细胞来自单一样本 %s",
                   gp, cl, 100 * dom$max_fraction, dom$dominant_sample)
      log_warn(STEP, w)
      repr_warnings[[length(repr_warnings) + 1]] <- w
    }
  }
}
repr_df <- do.call(rbind, repr_rows)
write.csv(repr_df, file.path(diag_dir, "sample_representation.csv"), row.names = FALSE)

# ---- 阶段16 固化（需求 29.1-29.5）：模块特征评分（可配置步骤） ----
# module_scores: 集合名 -> 基因列表；只写 <set>_score 数值型 metadata，不新增离散标签（29.4/29.5）
ms_dir <- file.path(run_root10, "module_score")
dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)
ms_meta_notes <- list(); cbg_rows <- list(); bs_rows <- list()
if (!is.null(cfg$module_scores) && length(cfg$module_scores)) {
  genes_all10 <- rownames(obj)
  for (set_nm in names(cfg$module_scores)) {
    gset <- as.character(unlist(cfg$module_scores[[set_nm]]))
    usage <- sgs_module_gene_usage(gset, genes_all10)   # 29.3 缺失只标记不静默丢弃
    used <- usage$used; missing <- usage$missing
    write.csv(data.frame(gene = gset,
                         missing = gset %in% missing,
                         used_in_calculation = gset %in% used),
              file.path(ms_dir, sprintf("%s_genes_used.csv", set_nm)), row.names = FALSE)
    if (length(missing)) {
      log_warn(STEP, sprintf("模块 %s 缺失基因（已标记）: %s", set_nm, paste(missing, collapse = ",")))
    }
    if (length(used) < 2) {
      log_warn(STEP, sprintf("模块 %s 可用基因 %d 个（<2），跳过评分", set_nm, length(used)))
      next
    }
    # 小对象自适应：AddModuleScore 按"feature 所在表达 bin"抽 ctrl 个对照基因，
    # bin 内基因数（≈总基因数/nbin）< ctrl 时 sample.int 越界，故按基因池大小收敛 nbin/ctrl
    n_genes_pool <- length(rownames(obj))
    if (n_genes_pool < 12) {
      log_warn(STEP, sprintf("模块 %s：基因池仅 %d 个（<12），跳过评分", set_nm, n_genes_pool))
      next
    }
    ms_nbin <- max(2L, min(24L, as.integer(floor(n_genes_pool / 6))))
    ms_ctrl <- max(1L, min(100L, as.integer(floor(n_genes_pool / ms_nbin)) - 1L))
    obj <- Seurat::AddModuleScore(obj, features = list(used), name = paste0(set_nm, "_score"),
                                  nbin = ms_nbin, ctrl = ms_ctrl)
    sc <- obj@meta.data[[paste0(set_nm, "_score1")]]
    stopifnot(is.numeric(sc))   # 29.5 静态防线：评分必须是数值列
    cbg_rows[[set_nm]] <- data.frame(set = set_nm, cell_group = cluster, score = sc,
                                     stringsAsFactors = FALSE)
    bs_rows[[set_nm]] <- data.frame(set = set_nm, cell_group = cluster, sample_id = samples,
                                    group = group, score = sc, stringsAsFactors = FALSE)
    ms_meta_notes[[set_nm]] <- list(genes_used = used, genes_missing = as.list(missing))
  }
  if (length(cbg_rows)) {
    all_cbg <- do.call(rbind, cbg_rows)
    agg_cbg <- aggregate(score ~ set + cell_group, data = all_cbg, FUN = mean)
    names(agg_cbg)[3] <- "mean_score"
    agg_cbg$n_cells <- as.integer(table(all_cbg$cell_group)[agg_cbg$cell_group])
    write.csv(agg_cbg, file.path(ms_dir, "by_cell_group.csv"), row.names = FALSE)
    all_bs <- do.call(rbind, bs_rows)
    agg_bs <- aggregate(score ~ set + sample_id + group, data = all_bs, FUN = mean)
    names(agg_bs)[4] <- "mean_score"
    agg_bs$n_cells <- as.integer(table(all_bs$sample_id)[agg_bs$sample_id])
    write.csv(agg_bs, file.path(ms_dir, "by_sample.csv"), row.names = FALSE)
    # 交叉聚合：每模块一列（<set>_score），供 s10b 直接消费
    wide_cs <- NULL
    for (st in names(bs_rows)) {
      d <- bs_rows[[st]]   # 细胞级：set/sample_id/group/score
      agg <- aggregate(score ~ cell_group + sample_id + group, data = d, FUN = mean)
      names(agg)[4] <- paste0(st, "_score")
      wide_cs <- if (is.null(wide_cs)) agg else merge(wide_cs, agg,
        by = c("cell_group", "sample_id", "group"), all = TRUE)
    }
    score_cols <- grep("_score$", names(wide_cs), value = TRUE)
    wide_cs$n_cells <- as.integer(table(paste(cluster, samples))[
      paste(wide_cs$cell_group, wide_cs$sample_id)])
    write.csv(wide_cs[, c("cell_group", "sample_id", "group", score_cols, "n_cells")],
              file.path(ms_dir, "by_cluster_sample.csv"), row.names = FALSE)
  }
  # 29.5 单元测试伴随检查：对象 metadata 中本步新增列全部数值型
  new_cols <- grep("_score1$", colnames(obj@meta.data), value = TRUE)
  if (length(new_cols)) {
    stopifnot(all(vapply(new_cols, function(cc) is.numeric(obj@meta.data[[cc]]), logical(1))))
  }
} else {
  log_info(STEP, "module_scores 未配置，模块评分步骤空转（by_cluster_sample.csv 写空表占位）")
  write.csv(data.frame(cell_group = character(0), sample_id = character(0), group = character(0),
                       n_cells = integer(0)),
            file.path(ms_dir, "by_cluster_sample.csv"), row.names = FALSE)
}

# 各细胞表达表（基因→科学问题随 manifest 记录，文件按基因命名）
expr_df <- cbind(data.frame(cell = colnames(obj), group = group, major_cluster = cluster),
                 expr[, kept, drop = FALSE])
expr_file <- file.path(out_dir, sprintf("%s_expr.csv", paste(kept, collapse = "_")))
write.csv(expr_df, expr_file, row.names = FALSE)

# ---- Dotplot（目标基因 × major_cluster；阶段15：图与源数据同出） ----
plt <- DotPlot(obj, features = kept, group.by = "major_cluster") +
  coord_flip() + theme_classic() + RotatedAxis() +
  ggtitle(sprintf("Target genes (kept %d/%d)", length(kept), length(targets)))
dot_src <- do.call(rbind, lapply(kept, function(g) {
  data.frame(gene = g, major_cluster = cluster, sample_id = samples, group = group,
             expr = expr[[g]], stringsAsFactors = FALSE)
}))
dot_src$n_cells <- as.integer(table(paste(cluster, samples))[
  paste(dot_src$major_cluster, dot_src$sample_id)])
figres <- sgs_save_plot(cfg, "P10", "PanelA", "target_gene_dotplot", plt, dot_src, STEP)
figs <- unlist(figres$files)

gate <- sgs_qc_gate(
  auto_pass = TRUE,
  detail = sprintf("kept %d/%d 目标基因 (%s)%s；stats+%s 落盘；ambient告警 %d 项；代表性告警 %d 项；模块 %s",
                   length(kept), length(targets), paste(kept, collapse = ","),
                   if (length(missing)) sprintf("，缺失: %s", paste(missing, collapse = ",")) else "",
                   basename(expr_file),
                   length(ambient_warnings), length(repr_warnings),
                   if (length(ms_meta_notes)) paste(names(ms_meta_notes), collapse = ",") else "(none)"))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

# 基因→科学问题对应表记入 manifest（需求16.3）；24.2/25.1 告警与 29 记录同条目入 Manifest
sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(
    targets = targets, kept = kept,
    science_questions = as.list(stats::setNames(
      science_q[intersect(kept, names(science_q))],
      intersect(kept, names(science_q)))),
    missing_genes = as.list(missing),
    warnings = list(ambient_rna = as.list(ambient_warnings),
                    sample_representation = as.list(repr_warnings)),
    module_scores = ms_meta_notes),
  outputs = as.list(c(file.path(out_dir, "target_gene_stats.csv"),
                      file.path(out_dir, "target_gene_stats_per_major.csv"), expr_file, figs,
                      if (file.exists(file.path(diag_dir, "ambient_rna.csv")))
                        file.path(diag_dir, "ambient_rna.csv") else NULL,
                      file.path(diag_dir, "sample_representation.csv"),
                      if (length(ms_meta_notes))
                        file.path(ms_dir, "by_cluster_sample.csv") else NULL)),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
log_info(STEP, "S10 complete")

# s8_trajectory.R — 拟时序模块（需求 14，隔离环境 sgscr_trajectory，slingshot）
# Validates: Requirements 14.1-14.2 (P11)；P11A 补充派发：multi-compartment
#
# 输入：rds/s6_<comp>_subcluster.rds（S6 产物），comp 遍历 trajectory.target_compartment
#   （支持列表，向后兼容单值；真实运行 ["myeloid","fibroblast"]，myeloid=primary，
#   fibroblast=backup，标注写入 Manifest）
# 前置校验（P11）：文件存在、slingshot 可加载、对象含聚类列与 UMAP
# 执行语义：各 compartment 独立执行；某一 comp FAIL 记 Gate=FAIL 但不阻断其余 comp；
#   全部 comp 均 FAIL 时进程退出码 1
# 输出：trajectory/pseudotime_<comp>.rds、figures/P8_*

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s8_trajectory.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s8_trajectory"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# ---- 断点续跑统一入口（P6） ----
pre <- tryCatch(sgs_resume_plan("s8", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }

# ---- target_compartment：列表或单值（向后兼容） ----
tc <- cfg$trajectory$target_compartment
if (is.null(tc)) { log_error(STEP, "trajectory.target_compartment 未配置"); quit(status = 1) }
comps <- as.character(if (length(tc) > 1 || is.list(tc)) unlist(tc) else tc)
comps <- unique(comps)

# P11A：primary=第一个，其余为 backup（论文主线/候补存档），写入 Manifest
role_of <- setNames(rep("backup", length(comps)), comps); role_of[[1]] <- "primary"

if (!requireNamespace("slingshot", quietly = TRUE)) {
  log_error(STEP, "slingshot包不可加载"); quit(status = 1)
}
cluster_col <- sprintf("RNA_snn_res.%g", cfg$compartments$subcluster_resolution)
out_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 各 compartment 独立执行，FAIL 隔离（P11A 派发） ----
gate_detail <- character(0)
out_all <- list()
n_ok <- 0L
comp_ok <- setNames(logical(length(comps)), comps)   # 显式成功标志（P11A 修正：不得以 detail 文本判成败）
for (comp in comps) {
  comp_started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")   # P12-a D3：每 comp 条目独立 started_at
  # P11D-c：优先读 S6b 清洗产物（存在时）；freshness 检查随之以实际输入为准
  in_rds <- file.path(cfg$paths$rds_dir, sprintf("s6b_%s_subcluster.rds", comp))
  input_src <- "s6b"
  if (!file.exists(in_rds)) {
    in_rds <- file.path(cfg$paths$rds_dir, sprintf("s6_%s_subcluster.rds", comp))
    input_src <- "s6"
  }
  g <- tryCatch({
    if (!file.exists(in_rds)) stop(sprintf("%s 不存在（需先跑 S6/S6b）", basename(in_rds)), call. = FALSE)
    # 产物新鲜度（P11 增强）：已有 s8 产物不得早于上游 s6 产物
    stale_out <- file.path(out_dir, sprintf("pseudotime_%s.rds", comp))
    if (file.exists(stale_out)) sgs_check_freshness(stale_out, in_rds)

    obj <- readRDS(in_rds)
    sgs_preflight(STEP, c(cluster_col, "UMAP"),
                  c(cluster_col, if (!is.null(Embeddings(obj, "umap"))) "UMAP"))

    # P11D-c：myeloid 指定起点簇（trajectory.start_cluster: <comp>: <簇id>；未配置则不指定）
    start_clus <- NULL
    if (!is.null(cfg$trajectory$start_cluster[[comp]])) {
      start_clus <- as.character(cfg$trajectory$start_cluster[[comp]])
      if (!start_clus %in% as.character(obj@meta.data[[cluster_col]])) {
        stop(sprintf("start_cluster %s 不在 %s 的亚群集合中", start_clus, comp), call. = FALSE)
      }
      log_info(STEP, sprintf("%s: slingshot 起点簇 = %s", comp, start_clus))
    }
    sds <- slingshot(Embeddings(obj, "umap"), clusterLabels = obj@meta.data[[cluster_col]],
                     start.clus = start_clus)
    pt <- slingPseudotime(sds)
    if (is.null(pt) || all(is.na(pt))) stop("slingshot 产出空拟时序", call. = FALSE)
    obj$pseudotime <- rowMeans(pt, na.rm = TRUE)

    out_rds <- file.path(out_dir, sprintf("pseudotime_%s.rds", comp))
    saveRDS(list(sds = sds, pseudotime_per_lineage = pt,
                 cell_pseudotime = obj$pseudotime, cluster_col = cluster_col), out_rds)

    df <- as.data.frame(Embeddings(obj, "umap"))
    df$pseudotime <- obj$pseudotime
    emb_cols <- colnames(df)[1:2]   # 列名随环境为 UMAP_1/umap_1，动态取用
    plt <- ggplot(df, aes(x = .data[[emb_cols[1]]], y = .data[[emb_cols[2]]], color = pseudotime)) +
      geom_point(size = 0.7) + scale_color_viridis_c() + theme_classic() +
      ggtitle(sprintf("%s pseudotime (n=%d, %s)", comp, ncol(obj), role_of[[comp]]))
    figs <- sgs_save_plot(cfg, "P8", "PanelA", paste0(comp, "_trajectory"), plt, STEP)

    gate <- sgs_qc_gate(
      auto_pass = sum(!is.na(obj$pseudotime)) > 0 &&
                  length(unique(obj@meta.data[[cluster_col]])) >= 2,
      detail = sprintf("%s[%s]: %d cells (%.0f%% with pseudotime), %d lineages, range [%.2f, %.2f]",
                       comp, role_of[[comp]], ncol(obj), 100 * mean(!is.na(obj$pseudotime)),
                       ncol(pt), min(obj$pseudotime, na.rm = TRUE), max(obj$pseudotime, na.rm = TRUE)))
    sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$trajectory,
      package_versions = list(slingshot = as.character(packageVersion("slingshot"))),
      quality_gate = gate,
      resolved_params = list(target_compartment = comp, role = role_of[[comp]],
                             cluster_col = cluster_col, input = input_src,
                             start_cluster = if (!is.null(start_clus)) start_clus else NULL),
      outputs = as.list(c(out_rds, figs)),
      started_at = comp_started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))
    out_all[[comp]] <- c(list(out_rds), as.list(figs))
    comp_ok[[comp]] <- identical(gate$status, "PASS")   # gate FAIL 也计入该 comp 未通过
    gate$detail
  }, error = function(e) {
    # FAIL 隔离：记录该 comp 失败原因，继续下一个 comp
    g <- sgs_qc_gate(auto_pass = FALSE, detail = sprintf("%s: %s", comp, conditionMessage(e)))
    sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$trajectory,
      package_versions = list(slingshot = as.character(packageVersion("slingshot"))),
      quality_gate = g,
      resolved_params = list(target_compartment = comp, role = role_of[[comp]]),
      outputs = list(),
      started_at = comp_started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    log_error(STEP, paste("Quality_Gate:", g$status, "-", g$detail))
    g$detail
  })
  gate_detail <- c(gate_detail, g)
  if (comp_ok[[comp]]) n_ok <- n_ok + 1L   # 显式成功标志，不靠 detail 字符串判别
}

# 总 Gate 语义（P11A）：backup FAIL 不阻断——primary（第一个 comp）PASS 且至少一个
# comp PASS 即整体 PASS；primary FAIL（或全部 FAIL）才 FAIL 退出阻断下游
gate <- sgs_qc_gate(auto_pass = comp_ok[[comps[1]]] && n_ok > 0,
                    detail = paste(gate_detail, collapse = "; "))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))
sgs_manifest_set_resolved(cfg$paths$manifest_path, "trajectory_compartment_roles", as.list(role_of))
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S8 complete")

# s7_cellchat.R — CellChat 细胞通讯模块（需求 13，隔离环境 sgscr_cellchat）
# Validates: Requirements 13.1-13.2 (P11)
#
# 前置校验（P11）：输入对象需含 major_cluster 与 group 元数据列；CellChat 可加载；
# 各组细胞数 >= cellchat.min_cells_per_group。任一缺失 → 终止并报告恰好那组缺失项。
# 实跑：createCellChat → computeCommunProb → computeCommunProbPathway → aggregateNet
#   → 按 Config group 拆分 → mergeCellChat 差异比较 → 通路气泡图
# 产物：cellchat/cellchat_<group>.rds、cellchat/comm_prob_<group>.csv、差异图

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s7_cellchat.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s7_cellchat"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# ---- P11 前置校验（先于一切计算） ----
problems <- c()
if (!file.exists(file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds"))) {
  problems <- c(problems, "s5_obj_annotated.rds")
}
if (!requireNamespace("CellChat", quietly = TRUE)) {
  problems <- c(problems, "CellChat包(隔离环境sgscr_cellchat未就绪)")
}
if (length(problems) > 0) {
  sgs_preflight(STEP, problems, character(0))   # 恰好报告缺失集并终止（P11）
}

obj <- readRDS(file.path(cfg$paths$rds_dir, "s5_obj_annotated.rds"))
meta_cols <- c("major_cluster", "group")
sgs_preflight(STEP, meta_cols, meta_cols[meta_cols %in% colnames(obj@meta.data)])

groups <- unique(as.character(obj$group))
min_cells <- cfg$cellchat$min_cells_per_group
small <- groups[vapply(groups, function(g) sum(obj$group == g) < min_cells, logical(1))]
if (length(small) > 0) {
  log_error(STEP, sprintf("组细胞数不足 (< %s): %s", min_cells, paste(small, collapse = ", ")))
  quit(status = 1)
}

# ---- S6b 清洗（P11D-c）：存在 s6b_<comp>_subcluster.rds 时，剔除其中被过滤掉的细胞 ----
# s6b 是 s6 子对象按标签 CSV 过滤后的产物；其与 s6 的细胞差集 = 本 compartment 清洗剔除细胞。
# 通讯计算在清洗后的全对象上进行（major_cluster 标签不变）。
cleaned_n <- 0L
for (comp in as.character(cfg$compartments$targets)) {
  s6_f <- file.path(cfg$paths$rds_dir, sprintf("s6_%s_subcluster.rds", comp))
  s6b_f <- file.path(cfg$paths$rds_dir, sprintf("s6b_%s_subcluster.rds", comp))
  if (file.exists(s6_f) && file.exists(s6b_f)) {
    s6_cells <- colnames(readRDS(s6_f))
    s6b_cells <- colnames(readRDS(s6b_f))
    drop_cells <- setdiff(s6_cells, s6b_cells)
    drop_cells <- intersect(drop_cells, colnames(obj))
    if (length(drop_cells) > 0) {
      obj <- subset(obj, cells = setdiff(colnames(obj), drop_cells))
      cleaned_n <- cleaned_n + length(drop_cells)
    }
  }
}
if (cleaned_n > 0) {
  log_info(STEP, sprintf("S6b 清洗：剔除 %d 个细胞（myeloid/fibroblast 亚群过滤），剩余 %d",
                         cleaned_n, ncol(obj)))
}

# 组集合对齐（2026-10-07 real-run s7 FAIL 固化，P11A 按 Reviewer 派发改用官方 liftCellChat）：
# mergeCellChat/netVisual_heatmap 要求各条件 cell group 集合一致，否则对比矩阵维度不匹配
# （"obj2 - obj1: 非整合陈列"）。两级处理：
# (1) createCellChat 之前剔除 exclude_groups（如 low_quality）——这部分细胞不参与任何计算；
# (2) merge 之前用 liftCellChat 把各条件对象提升到分组并集（官方比较流程做法），
#     使某条件缺簇（真实数据：UC 中 dendritic 与 epi_other 为 0 细胞）后仍可安全比较。
excluded <- if (is.null(cfg$cellchat$exclude_groups)) character(0) else as.character(cfg$cellchat$exclude_groups)
if (length(excluded) > 0) {
  drop_cells <- as.character(obj$major_cluster) %in% excluded
  log_info(STEP, sprintf("exclude_groups 剔除 %d 个 cluster 共 %d 细胞: %s",
                         length(excluded), sum(drop_cells), paste(excluded, collapse = ", ")))
  obj <- subset(obj, subset = major_cluster %in% setdiff(unique(as.character(obj$major_cluster)), excluded))
}
if (length(unique(as.character(obj$major_cluster))) < 2) {
  log_error(STEP, "剔除 exclude_groups 后仅剩不足 2 个 cluster，无法做条件间比较")
  quit(status = 1)
}
sgs_preflight(STEP, "cellchat.pathways_of_interest",
              if (!is.null(cfg$cellchat$pathways_of_interest)) "cellchat.pathways_of_interest")

suppressPackageStartupMessages(library(CellChat))

cc_dir <- file.path(dirname(cfg$paths$manifest_path), "..", "cellchat")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)
cc_objs <- list()
for (g in groups) {
  sub <- subset(obj, subset = group == g)
  data.input <- GetAssayData(sub, assay = "RNA", layer = "data")
  meta <- data.frame(labels = as.character(sub$major_cluster),
                     group = as.character(sub$group), row.names = colnames(sub))
  cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cc <- CellChat::setIdent(object = cc, ident.use = "labels")
  db_obj <- get(sprintf("CellChatDB.%s", cfg$species$organism),
                envir = asNamespace("CellChat"))   # 数据库物种随 Config（P5 一致性承接）
  cc@DB <- db_obj
  cc <- subsetData(cc)
  if (sum(cc@data.signaling) == 0) {
    log_error(STEP, sprintf("group %s: subsetData 后信令基因为空 — 输入基因名与 %s CellChatDB 无交集（冒烟数据请用 make_mock_10x.R 生成真实基因名）",
                            g, cfg$species$organism))
    next
  }
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # 无 presto 时用标准 Wilcoxon
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc_objs[[g]] <- cc
  saveRDS(cc, file.path(cc_dir, sprintf("cellchat_%s.rds", g)))
  # 通讯概率矩阵导出（长表 CSV，base R，avoid reshape2 依赖）
  prob <- cc@net$prob[, , 1]   # source × target
  df <- data.frame(source = rep(rownames(prob), ncol(prob)),
                   target = rep(colnames(prob), each = nrow(prob)),
                   prob = as.vector(prob))
  write.csv(df, file.path(cc_dir, sprintf("comm_prob_%s.csv", g)), row.names = FALSE)
  log_info(STEP, sprintf("group %s: %d interactions computed", g, nrow(df)))
}

# ---- 组间差异（>=2 组时） ----
figs <- c()
if (length(cc_objs) >= 2) {
  # merge 前 liftCellChat 对齐到分组并集（各条件 filterCommunication 剔除小簇后集合可能不一致）
  label_sets <- lapply(cc_objs, function(cc) sort(as.character(levels(cc@idents))))
  all_labels <- sort(unique(unlist(label_sets)))
  lifted <- length(unique(label_sets)) != 1L
  if (lifted) {
    log_info(STEP, sprintf("各条件 cell group 集合不一致（%s），liftCellChat 对齐到并集 %d 组",
                           paste(vapply(label_sets, length, integer(1)), collapse = "/"),
                           length(all_labels)))
    cc_objs <- lapply(cc_objs, function(cc) CellChat::liftCellChat(cc, group.new = all_labels))
  }
  merged <- mergeCellChat(cc_objs, add.names = groups)
  saveRDS(merged, file.path(cc_dir, "cellchat_merged.rds"))
  # netVisual_heatmap 返回 ComplexHeatmap 对象（非 ggplot），用其原生设备落盘
  hm <- CellChat::netVisual_heatmap(merged, comparison = c(1, 2))
  figs <- sgs_save_heatmap(cfg, "P7", "PanelA", "CellChat_diff_interaction", hm, STEP)
  pathways <- cfg$cellchat$pathways_of_interest
  log_info(STEP, sprintf("关注通路: %s（差异分析图见 P7 产物）",
                         paste(pathways, collapse = ", ")))
}

gate <- sgs_qc_gate(auto_pass = length(cc_objs) == length(groups),
                    detail = sprintf("CellChat computed for %d/%d groups",
                                     length(cc_objs), length(groups)))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))
sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$cellchat,
  package_versions = list(CellChat = as.character(packageVersion("CellChat"))),
  quality_gate = gate,
  resolved_params = list(groups = as.list(groups), min_cells = min_cells,
                         exclude_groups = as.list(excluded),
                         cellchat_lifted = lifted,
                         label_union = as.list(if (length(cc_objs) >= 2) all_labels else
                           as.character(levels(cc_objs[[1]]@idents))),
                         s6b_cells_removed = cleaned_n,
                         pathways = as.list(cfg$cellchat$pathways_of_interest)),
  outputs = as.list(c(file.path(cc_dir, sprintf("cellchat_%s.rds", groups)), figs)),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (gate$status == "FAIL") quit(status = 1)
log_info(STEP, "S7 complete")

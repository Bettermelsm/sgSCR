# subcluster_filter_utils.R — S6b 亚群注释过滤（P11D-c 派工项 1）
# Validates: 0018 派工：亚群标签 CSV 校验、keep 语义、过滤纯函数
#
# keep 语义（0018）：yes=保留；sens=主分析保留、敏感性分析剔除；no=剔除。
# filter_mode: main → keep ∈ {yes, sens}；strict → keep == "yes"。
# 校验为纯函数，属性测试直接消费；任何缺失/多余亚群 → FAIL（问题逐条列出）。

# 校验标签 CSV 与 s6 对象亚群一一对应。csv_df 需含 compartment, subcluster, keep 列。
# 返回 list(ok, problems, labels_map(comp→subcluster→label), keep_map(comp→subcluster→keep))
sgs_subcluster_validate <- function(csv_df, targets, ids_by_comp) {
  problems <- character(0)
  need_cols <- c("compartment", "subcluster", "keep")
  miss_cols <- setdiff(need_cols, colnames(csv_df))
  if (length(miss_cols)) {
    return(list(ok = FALSE,
                problems = sprintf("标签 CSV 缺列: %s", paste(miss_cols, collapse = ",")),
                labels_map = NULL, keep_map = NULL))
  }
  labels_map <- list(); keep_map <- list()
  for (comp in targets) {
    if (!comp %in% csv_df$compartment) {
      problems <- c(problems, sprintf("标签 CSV 无 compartment=%s 的任何行", comp)); next
    }
    rows <- csv_df[csv_df$compartment == comp, ]
    ids <- as.character(ids_by_comp[[comp]])
    csv_ids <- as.character(rows$subcluster)
    miss <- setdiff(ids, csv_ids)     # 对象有、CSV 无 → 漏
    extra <- setdiff(csv_ids, ids)    # CSV 有、对象无 → 多
    if (length(miss)) problems <- c(problems, sprintf("%s: CSV 缺亚群 %s", comp, paste(miss, collapse = ",")))
    if (length(extra)) problems <- c(problems, sprintf("%s: CSV 多出亚群 %s", comp, paste(extra, collapse = ",")))
    dup <- csv_ids[duplicated(csv_ids)]
    if (length(dup)) problems <- c(problems, sprintf("%s: CSV 亚群重复 %s", comp, paste(unique(dup), collapse = ",")))
    bad_keep <- setdiff(unique(rows$keep), c("yes", "sens", "no"))
    if (length(bad_keep)) problems <- c(problems, sprintf("%s: 非法 keep 值 %s（只允许 yes/sens/no）", comp, paste(bad_keep, collapse = ",")))
    m <- rows$subcluster; names(m) <- if ("proposed_label" %in% colnames(rows)) rows$proposed_label else rows$subcluster
    labels_map[[comp]] <- stats::setNames(as.character(rows$proposed_label), as.character(rows$subcluster))
    keep_map[[comp]] <- stats::setNames(as.character(rows$keep), as.character(rows$subcluster))
  }
  list(ok = length(problems) == 0, problems = problems,
       labels_map = labels_map, keep_map = keep_map)
}

# keep 语义 → 保留亚群 ID 集（纯函数）
sgs_subcluster_keep_ids <- function(keep_map_comp, mode) {
  stopifnot(mode %in% c("main", "strict"))
  allowed <- if (mode == "main") c("yes", "sens") else "yes"
  ids <- names(keep_map_comp)[keep_map_comp %in% allowed]
  as.character(sort(as.integer(ids)))
}

# 过滤纯函数：按保留亚群子集化对象，写 subcluster_label 元数据列。
# res_col 必须精确匹配（调用方用 sprintf("%s_snn_res.%g", assay, res) 构造，禁止 grep 首列）。
sgs_filter_subobject <- function(obj, res_col, keep_ids, labels_map_comp) {
  cl <- as.character(obj@meta.data[[res_col]])
  keep_cells <- cl %in% keep_ids
  if (!any(keep_cells)) stop("过滤后 0 细胞（keep_ids 与对象无交集）", call. = FALSE)
  out <- subset(obj, cells = colnames(obj)[keep_cells])
  cl2 <- as.character(out@meta.data[[res_col]])
  out$subcluster_label <- unname(labels_map_comp[cl2])
  out
}

# 落盘前自检：过滤后细胞数 == 保留亚群细胞数之和（纯断言，供 S6b 门控与属性测试）
sgs_filter_count_check <- function(obj, res_col, keep_ids, n_expected) {
  cl <- as.character(obj@meta.data[[res_col]])
  n_actual <- sum(cl %in% keep_ids)
  list(ok = identical(as.integer(n_actual), as.integer(n_expected)),
       n_actual = n_actual, n_expected = as.integer(n_expected))
}

# P12-a D5：删除子对象继承自父对象的 *_snn_res.* 扫描列（保留指定亚群列）。
# 背景：s6 子对象带着父对象 S4 分辨率扫描列（如 RNA_snn_res.0.1..1.1），
# 消费方 grep 首列会误取（P11C 实测）；根治 = 落盘前删除 + 取列一律 sprintf 精确匹配。
sgs_drop_inherited_snn_cols <- function(obj, keep_col) {
  snn <- grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE)
  drop <- setdiff(snn, keep_col)
  if (length(drop)) obj@meta.data[drop] <- NULL
  obj
}

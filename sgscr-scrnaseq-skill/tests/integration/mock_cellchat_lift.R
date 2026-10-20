# 集成测试（mock）：复现"两条件分组集合不同 → non-conformable"并验证 liftCellChat 修复
# Validates: P11A 派发 1（2026-10-07 real-run s7 FAIL 固化）
#
# 用 mock 数据构造 cond1={A,B,C} / cond2={A,B}（C 在 cond2 不存在，
# 模拟真实 UC 中 dendritic/epi_other 为 0 细胞），分别验证：
#   (1) 不对齐直接 merge+heatmap 复现原错误（贴原文）；
#   (2) liftCellChat 提升到分组并集后 merge+heatmap 成功。
# Run: Rscript tests/integration/mock_cellchat_lift.R   （需 CellChat 环境，如 sgscr_cellchat）

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) package <- args[1] else package <- NA
if (!is.na(package) && package == "mock") {
  # 供无 CellChat 环境的 CI 跳过
  cat("SKIP: CellChat not requested\n"); quit(status = 0)
}
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(Matrix))

set.seed(42)
db <- CellChatDB.human
db_interaction <- if (!is.null(db$interaction)) db$interaction else db@interaction   # 2.x 为 list
pairs <- head(db_interaction, 60)   # 注意：ligand/receptor 列是复合物名，基因须用 symbol 列
pair_lig_genes <- lapply(pairs$ligand.symbol, function(s) unlist(strsplit(s, ", ")))
pair_rec_genes <- lapply(pairs$receptor.symbol, function(s) unlist(strsplit(s, ", ")))
genes <- sort(unique(c(unlist(pair_lig_genes), unlist(pair_rec_genes))))

mk_cc <- function(labels_set, prefix) {
  cells_per <- 60
  cell_names <- unlist(lapply(labels_set, function(l) sprintf("%s_%s_%d", prefix, l, 1:cells_per)))
  lab_vec <- rep(labels_set, each = cells_per)
  # 每个 cluster 高表达一组真实 DB 配对的配体+受体，保证 Wilcoxon/triMean 有可检出的 L-R 信号
  mat <- sapply(seq_along(lab_vec), function(i) {
    m <- match(lab_vec[i], labels_set)
    sel <- ((m - 1) * 10 + 1):(m * 10)   # 每组独占 20 个 DB 配对（互不重叠，需 length(labels_set)*10 <= 60）
    base <- stats::rnorm(length(genes), mean = 0.1, sd = 0.5)
    base[genes %in% unlist(pair_lig_genes[sel])] <- 4
    base[genes %in% unlist(pair_rec_genes[sel])] <- 4
    stats::rnbinom(length(genes), size = 5, mu = exp(base))
  })
  dimnames(mat) <- list(genes, cell_names)
  # CellChat 要求 log-normalized 数据（官方流程 NormalizeData+log1p 后传入）
  mat <- log1p(t(t(mat) / colSums(mat)) * 10000)
  meta <- data.frame(labels = factor(lab_vec, levels = labels_set), row.names = cell_names)
  cc <- createCellChat(object = as(mat, "dgCMatrix"), meta = meta, group.by = "labels")
  cc <- CellChat::setIdent(object = cc, ident.use = "labels")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # 无 presto 时用标准 Wilcoxon
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 30)   # >60 不会剔组；组差异由下面构造
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# cond1: A,B,C；cond2: A,B —— C 在 cond2 不存在（对齐真实 UC 中 dendritic/epi_other 为 0 细胞）
cond1 <- mk_cc(c("A", "B", "C"), "c1")
cond2 <- mk_cc(c("A", "B"), "c2")

cat("label sets: cond1=", paste(levels(cond1@idents), collapse = ","),
    " cond2=", paste(levels(cond2@idents), collapse = ","), "\n", sep = "")

## (1) 复现：不对齐直接 merge + heatmap
merged_raw <- tryCatch({
  mergeCellChat(list(cond1, cond2), add.names = c("NC", "UC"))
}, error = function(e) { cat("mergeCellChat 原始报错: ", conditionMessage(e), "\n"); NULL })
hm_err <- NULL
if (!is.null(merged_raw)) {
  hm_err <- tryCatch({
    CellChat::netVisual_heatmap(merged_raw, comparison = c(1, 2))
    "NO_ERROR(heatmap succeeded)"
  }, error = function(e) conditionMessage(e))
  cat("netVisual_heatmap 原始报错: ", hm_err, "\n")
}

## (2) 修复：liftCellChat 到并集
all_labels <- sort(union(levels(cond1@idents), levels(cond2@idents)))
lift1 <- liftCellChat(cond1, group.new = all_labels)
lift2 <- liftCellChat(cond2, group.new = all_labels)
cat("lifted sets: ", paste(levels(lift1@idents), collapse = ","),
    " | ", paste(levels(lift2@idents), collapse = ","), "\n", sep = "")
merged <- mergeCellChat(list(lift1, lift2), add.names = c("NC", "UC"))
hm_ok <- tryCatch({
  CellChat::netVisual_heatmap(merged, comparison = c(1, 2))
  "HEATMAP_OK"
}, error = function(e) paste("STILL_FAILS:", conditionMessage(e)))
cat("lift 后 netVisual_heatmap: ", hm_ok, "\n")

ok <- identical(hm_ok, "HEATMAP_OK")
cat(if (ok) "MOCK_LIFT_TEST_PASS\n" else "MOCK_LIFT_TEST_FAIL\n")
quit(status = if (ok) 0 else 1)

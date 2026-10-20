# Property tests for marker_panel_utils（P12-a A4）与 sgs_drop_inherited_snn_cols（P12-a D5）
#   sgs_load_marker_panel：真实面板 human_colon.yaml 可加载、66 基因、跨谱系零重复；
#     缺 config 键 / 面板不存在 / 面板结构非法 → 明确报错
#   sgs_drop_inherited_snn_cols：任意继承 SNN 列组合下，只保留 keep_col；其余列原样保留
# Run: Rscript tests/property/test_marker_panel_and_snn.R


suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(yaml))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
SKILL <- file.path(HERE, "..", "..")
source(file.path(SKILL, "scripts", "common", "marker_panel_utils.R"))
source(file.path(SKILL, "scripts", "common", "subcluster_filter_utils.R"))

test_that("human_colon panel loads: 66 genes, zero cross-lineage duplicates", {
  p <- sgs_load_marker_panel(list(annotation = list(marker_panel = "human_colon")),
                             skill_dir = SKILL)
  expect_identical(p$species, "human")
  expect_identical(p$tissue, "colon")
  expect_gte(p$n_genes, 60L)                    # P11D 实测 66；防面板被悄悄删减
  expect_equal(length(p$genes), length(unique(p$genes)))
  all_g <- unlist(p$lineages)
  expect_false(any(duplicated(all_g)))
  # 每个基因都在全集中（lineages 与 genes 一致）
  expect_setequal(all_g, p$genes)
})

test_that("loader errors are explicit on missing key / missing file / bad structure", {
  expect_error(sgs_load_marker_panel(list(annotation = NULL), skill_dir = SKILL),
               "marker_panel")
  expect_error(sgs_load_marker_panel(list(annotation = list(marker_panel = "no_such_panel")),
                                     skill_dir = SKILL),
               "不存在")
})

test_that("bad-structure panel fails loudly", {
  root <- file.path(tempdir(), sprintf("sgscr_panel_test_%d", as.integer(Sys.time())))
  dir.create(file.path(root, "references", "marker_panels"), recursive = TRUE)
  writeLines("species: human\ntissue: colon\nlineages: {}",
             file.path(root, "references", "marker_panels", "empty_panel.yaml"))
  writeLines("species: human\ntissue: colon\nlineages:\n  a: [X, Y]\n  b: [Y, Z]",
             file.path(root, "references", "marker_panels", "dup_panel.yaml"))
  expect_error(sgs_load_marker_panel(list(annotation = list(marker_panel = "empty_panel")),
                                     skill_dir = root), "lineages")
  expect_error(sgs_load_marker_panel(list(annotation = list(marker_panel = "dup_panel")),
                                     skill_dir = root), "重复")
})

test_that("D5: only keep_col survives; arbitrary inherited SNN columns are dropped", {
  pool <- sprintf("RNA_snn_res.%s", c("0.1", "0.3", "0.5", "0.7", "0.9", "1.1"))
  n <- 30
  results <- logical(2^length(pool))
  for (i in seq_along(results)) {
    mask <- as.logical(intToBits(i - 1L)[seq_along(pool)])
    inherited <- pool[mask]
    m <- Matrix::rsparsematrix(20, n, 0.2)
    o <- CreateSeuratObject(counts = m)
    for (cl in inherited) o@meta.data[[cl]] <- as.character(sample(0:2, n, replace = TRUE))
    keep_col <- "RNA_snn_res.1"
    o@meta.data[[keep_col]] <- as.character(rep(0:2, each = 10))
    o$batch <- "b1"
    out <- sgs_drop_inherited_snn_cols(o, keep_col)
    snn_left <- grep("_snn_res\\.", colnames(out@meta.data), value = TRUE)
    results[i] <- setequal(snn_left, keep_col) &&
                  all(inherited %in% colnames(o@meta.data)) &&          # 原对象不被就地修改
                  identical(as.character(out$batch), rep("b1", n)) &&
                  identical(as.character(out@meta.data[[keep_col]]), as.character(rep(0:2, each = 10)))
  }
  expect_true(all(results),
              info = sprintf("失败组合索引: %s",
                             paste(which(!results), collapse = ",")))
})

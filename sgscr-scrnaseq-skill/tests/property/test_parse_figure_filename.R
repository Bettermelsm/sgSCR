# Property tests for sgs_parse_figure_filename（阶段 11：Reviewer 登记隐患的补测）
#   覆盖：Fig/P 两种前缀、点分编号（P2.1）、无前缀纯数字（回退归入 P 系列）、
#         含特殊字符描述、不可解析文件名返回 ok=FALSE。
# Run: Rscript tests/property/test_parse_figure_filename.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "save_plot.R"))

test_that("parse: Fig/P prefix, dotted and plain ids all resolve", {
  property <- forall(
    gen.bind(function(pre) gen.bind(function(num) gen.map(
      function(panel) list(pre = pre, num = num, panel = panel),
      gen.element(c("PanelA", "PanelB", "PanelC1"))),
      gen.element(c("8", "2.1", "10"))),
      gen.element(c("Fig", "P"))),
    function(x) {
      f <- sprintf("%s%s_%s_umap_major.pdf", x$pre, x$num, x$panel)
      fp <- sgs_parse_figure_filename(f)
      expect_true(fp$ok)
      expect_identical(fp$fig_id, paste0(x$pre, x$num))
      expect_identical(fp$panel_id, x$panel)
      expect_identical(fp$desc, "umap_major")
    }
  )
})

test_that("parse: bare-number fallback is classified as P-series (documented default)", {
  fp <- sgs_parse_figure_filename("3_PanelB_target_gene_dotplot.csv")
  expect_true(fp$ok)
  expect_identical(fp$fig_id, "P3")   # 无前缀 → 补 "P"（约定仅 Fig/P 两种合法前缀）
  fp2 <- sgs_parse_figure_filename("P10_PanelA_x.tiff")
  expect_identical(fp2$fig_id, "P10")
})

test_that("parse: unparseable names return ok=FALSE and keep basename as desc", {
  for (f in c("heatmap_final.pdf", "P8_PanelA.pdf", "Fig3_umap.png",
              "P8_PanelA_.pdf")) {
    fp <- sgs_parse_figure_filename(f)
    expect_false(fp$ok, info = f)
    expect_identical(fp$desc, f, info = f)
    expect_true(is.na(fp$fig_id), info = f)
  }
})

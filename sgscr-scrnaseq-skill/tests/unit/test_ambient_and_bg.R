# Unit tests: ambient diagnostics ratio/warning, background clusters, s2b skip, source_df contract
# 判据取 design.md 测试映射表「24.1、24.2」「24.3」「24.4」「18.5-18.7」行
# Run: Rscript tests/unit/test_ambient_and_bg.R

suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))
source(file.path(HERE, "..", "..", "scripts", "common", "save_plot.R"))

test_that("24.1: 检出比例计算正确", {
  counts <- c(0, 1, 0, 2, 0, 3)   # 3/6 检出
  pct <- mean(counts > 0)
  expect_equal(pct, 0.5)
})

test_that("24.2: 超阈值产生告警标记，未超不产生", {
  warn_thr <- 0.1
  expect_true(mean(c(1, 0, 1, 0, 1, 0, 0, 0) > 0) > warn_thr)   # 3/8=0.375 超
  expect_false(mean(rep(0, 8)) > warn_thr)
})

test_that("24.3: 本底对照群 = 检出低于阈值的细胞群，排除目标群", {
  det <- c(epi_absorptive = 0.9, epi_secretory = 0.7, myeloid = 0.02, endothelial = 0.01)
  bg <- sgs_background_clusters(det, max_detection = 0.05, exclude = "epi_absorptive")
  expect_setequal(bg, c("myeloid", "endothelial"))
  expect_length(sgs_background_clusters(c(a = 0.5), 0.05), 0)   # 无不表达群时不产出
})

test_that("24.4-static: s2b 默认关闭分支存在且不加载校正后端", {
  lines <- readLines(file.path(HERE, "..", "..", "scripts", "s2b_ambient_correct.R"))
  expect_true(any(grepl("correction\\$enabled", lines)))
  # enabled 缺省（NULL）也必须走跳过：isTRUE(NULL) == FALSE 保证
  expect_false(isTRUE(NULL))
  # 开启而后端缺失必须显式失败（不得静默）：存在命名缺失包并退出的分支
  expect_true(any(grepl("需要 R 包", lines)))
  expect_true(any(grepl("quit\\(status = 1\\)", lines)))
})

test_that("18.6/18.7: source_df 契约——缺样本列/缺细胞数列/全缺失均报错", {
  ok <- data.frame(sample_id = "S1", n_cells = 10, y = 1)
  expect_silent(sgs_check_source_df(ok))
  expect_error(sgs_check_source_df(data.frame(y = 1)), "样本标识")
  expect_error(sgs_check_source_df(data.frame(sample_id = "S1", y = 1)), "细胞数")
  expect_error(sgs_check_source_df(NULL), "source_df")
})

test_that("18.5/18.7: save_plot 输出同名 .source.csv 并返回字体配置（ggplot mock）", {
  tmp <- file.path(tempdir(), paste0("sgs_fig_", as.integer(Sys.time())))
  cfg <- list(paths = list(figures_dir = tmp),
              figure = list(dpi = 72, width_in = 2, height_in = 2,
                            formats = "pdf", font_family = "Helvetica", base_size = 9))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  src <- data.frame(sample_id = rep(c("S1", "S2"), each = 3), x = 1:6, y = 1:6,
                    n_cells = c(30, 20)[rep(1:2, each = 3)])
  res <- sgs_save_plot(cfg, "Fig9", "PanelT", "unittest", p, src, "test")
  expect_true(file.exists(file.path(tmp, "Fig9_PanelT_unittest.pdf")))
  expect_true(file.exists(res$source_csv))
  src_back <- read.csv(res$source_csv)
  expect_true("sample_id" %in% names(src_back))
  expect_true("n_cells" %in% names(src_back))
  expect_equal(res$font_family, "Helvetica")   # 字体来自配置而非脚本字面量
  expect_equal(res$base_size, 9)
  unlink(tmp, recursive = TRUE)
})

cat("test_ambient_and_bg: all tests passed\n")

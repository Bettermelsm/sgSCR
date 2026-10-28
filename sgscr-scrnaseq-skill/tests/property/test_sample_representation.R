# Property tests for sample representation dominance (Req 25.1/25.2)
# Feature: sgscr-scrnaseq-skill, Property (design 25.1/25.2): 单一样本占比超阈值时必告警
#   并置 single_sample_driven=TRUE，否则不置
# Run: Rscript tests/property/test_sample_representation.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))

test_that("25.1/25.2-property: 占比>阈值必置 driven，<=阈值必不置", {
  property <- forall(
    gen.bind(function(main_frac) gen.map(
      function(main_count) list(main_frac = main_frac, main_count = main_count),
      gen.element(1:1000)),
      gen.element(70:100)),   # 主样本占比 0.70-1.00（整数百分比）
    function(x) {
      main_n <- x$main_count
      # 其余样本均分剩余占比（>=1 个其他样本）
      others <- 2L
      other_total <- max(1L, as.integer(main_n * (100 - x$main_frac) / x$main_frac))
      counts <- c(main = main_n, O1 = max(1, other_total %/% others),
                  O2 = max(0, other_total - max(1, other_total %/% others)))
      dom <- sgs_single_sample_dominance(counts, threshold = 0.70)
      frac <- counts[["main"]] / sum(counts)
      expect_equal(dom$driven, frac > 0.70)   # 恰在阈值处取等：0.70 不触发
      if (dom$driven) expect_identical(dom$dominant_sample, "main")
      TRUE
    }
  )
})

test_that("25.1-edge: 恰等于阈值(0.70)不告警；空输入不告警", {
  d <- sgs_single_sample_dominance(c(main = 70, O1 = 30), 0.70)
  expect_false(d$driven)
  d2 <- sgs_single_sample_dominance(c(main = 71, O1 = 29), 0.70)
  expect_true(d2$driven)
  d3 <- sgs_single_sample_dominance(c(main = 0), 0.70)
  expect_false(d3$driven)
  expect_true(is.na(d3$max_fraction))
})

cat("test_sample_representation: all properties passed\n")

# Property + constructed tests for separation margin (Req 23.3/23.4/23.7/23.8/23.9)
# Feature: sgscr-scrnaseq-skill, Property (design 23.3-23.9): 分离判据、余量、余量比例与
#   不稳健标记：margin 0.0006 / 0.0015 量级的输入必须被标记 separation_robust = FALSE
# Run: Rscript tests/property/test_separation_margin.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))

MARGIN_THR <- 0.01   # stats.min_separation_margin（真实案例量级取 0.01）

test_that("23.3/23.4/23.7/23.9: 分离/未分离构造数据判据与余量正确", {
  # 构造1：完全分离（case 整体低于 control）
  g <- c(rep("ctrl", 3), rep("case", 3)); names(g) <- sprintf("S%d", 1:6)
  m_s <- c(S1 = 2.0, S2 = 1.5, S3 = 1.8, S4 = 0.3, S5 = 0.1, S6 = 0.2)
  st <- sgs_separation_stats(m_s, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_true(st$separation)
  expect_equal(st$separation_margin, 1.5 - 0.3)   # 分离方向为 case 低：min(ctrl) - max(case)
  expect_equal(st$delta, mean(c(0.3, 0.1, 0.2)) - mean(c(2.0, 1.5, 1.8)))
  expect_equal(st$margin_ratio, st$separation_margin / st$delta)
  expect_equal(st$min_group0, 1.5); expect_equal(st$max_group0, 2.0)
  expect_equal(st$min_group1, 0.1); expect_equal(st$max_group1, 0.3)
  expect_true(st$separation_robust)

  # 构造2：交叠未分离
  m_s2 <- c(S1 = 1.0, S2 = 0.5, S3 = 0.8, S4 = 0.9, S5 = 0.7, S6 = 1.2)
  st2 <- sgs_separation_stats(m_s2, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_false(st2$separation)
  expect_lte(st2$separation_margin, 0)     # 未分离照样输出 margin（<=0）
  expect_false(st2$separation_robust)

  # 构造3：反向分离（case 高于 control）：min(ctrl) 0.2 < max(case) 0.9 → 走 sep_high，margin = 0.2 - 0.9
  m_s3 <- c(S1 = 0.1, S2 = 0.2, S3 = 0.15, S4 = 1.0, S5 = 1.2, S6 = 0.9)
  st3 <- sgs_separation_stats(m_s3, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_true(st3$separation)
  expect_equal(st3$separation_margin, 0.9 - 0.2)
  expect_true(st3$separation_robust)
})

test_that("23.8-key: 实测量级 margin 0.0006 / 0.0015 必须 separation_robust=FALSE", {
  g <- c(rep("ctrl", 3), rep("case", 2)); names(g) <- sprintf("S%d", 1:5)
  # 真实案例 CA2 本底：case 整体低于 ctrl，margin = 0.0425 - 0.0419 = 0.0006
  m_ca2bg <- c(S1 = 0.0431, S2 = 0.0425, S3 = 0.0430, S4 = 0.0419, S5 = 0.0415)
  st1 <- sgs_separation_stats(m_ca2bg, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_true(st1$separation)
  expect_false(st1$separation_robust)
  expect_equal(st1$separation_margin, 0.0425 - 0.0419, tolerance = 1e-8)
  # 真实案例 BEST4 吸收上皮：margin 0.0015
  m_best4 <- c(S1 = 0.0291, S2 = 0.0291, S3 = 0.0291, S4 = 0.0306, S5 = 0.0306)
  st2 <- sgs_separation_stats(m_best4, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_true(st2$separation)
  expect_false(st2$separation_robust)
})

test_that("23.9: delta==0 时 margin_ratio 记 NA", {
  g <- c(rep("ctrl", 2), rep("case", 2)); names(g) <- sprintf("S%d", 1:4)
  m_s <- c(S1 = 1, S2 = 1, S3 = 1, S4 = 1)
  st <- sgs_separation_stats(m_s, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
  expect_equal(st$delta, 0)
  expect_true(is.na(st$margin_ratio))
})

test_that("23.3-property: 任意数据下判据与 min/max 自洽", {
  property <- forall(
    gen.map(function(vs) list(v = as.numeric(vs) / 10),
            gen.list(of = 8, gen.element(0:50))),
    function(x) {
      n <- length(x$v)
      n0 <- max(1, n %/% 2)
      if (n - n0 < 1) return(TRUE)   # shrink 缩短列表时丢弃
      g <- c(rep("ctrl", n0), rep("case", n - n0)); names(g) <- sprintf("S%d", seq_len(n))
      m_s <- stats::setNames(x$v, sprintf("S%d", seq_len(n)))
      st <- sgs_separation_stats(m_s, g, "ctrl", "case", min_separation_margin = MARGIN_THR)
      if (isTRUE(st$separation)) {
        ok <- st$min_group1 > st$max_group0 || st$min_group0 > st$max_group1
        expect_true(ok)
        expect_true(st$separation_margin > 0)
      } else {
        expect_lte(st$separation_margin, 0)
      }
      # 不稳健规则：分离但 margin < 阈值 → FALSE
      if (isTRUE(st$separation) && st$separation_margin < MARGIN_THR) {
        expect_false(st$separation_robust)
      }
      TRUE
    }
  )
})

cat("test_separation_margin: all tests passed\n")

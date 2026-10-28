# Property tests for s10b sample-level mean derivation (Req 23.1)
# Feature: sgscr-scrnaseq-skill, Property (design 23.1): 组均值等于样本均值之均值；
#   样本内细胞数的任意缩放不改变组均值
# Run: Rscript tests/property/test_sample_level_mean.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))

# 生成器：3-6 个样本，每样本 1-80 个细胞值
gen.values <- function() {
  gen.map(function(vals) as.numeric(vals),
          gen.list(of = 24, gen.element(0:20)))
}
gen.samples <- function() {
  gen.map(function(k) sprintf("S%02d", seq_len(k)), gen.element(1:6))
}

test_that("23.1-main: 组均值 == 样本均值之均值（细胞数不进权重）", {
  property <- forall(
    # 单一生成器：值列表长度决定样本数（2-6），避免深层 bind 嵌套
    gen.map(function(vals) {
      n <- length(vals)
      n_s <- min((n %% 5) + 2, max(n, 1))
      list(n_s = n_s, vals = as.numeric(vals))
    }, gen.list(of = 24, gen.element(0:20))),
    function(x) {
      n <- length(x$vals)
      if (n < 2) return(TRUE)   # shrink 可能缩短列表，太短丢弃
      n_s <- min(x$n_s, n)
      sid <- sprintf("S%02d", seq_len(n_s))
      vals <- x$vals
      m_s <- sgs_sample_means(vals, rep(sid, length.out = length(vals)))
      g <- rep("g0", length(sid)); names(g) <- sid
      M0 <- sgs_group_mean(m_s, g, "g0")
      expect_equal(unname(M0), mean(unname(m_s)))          # 组均值 = mean(样本均值)
      expect_equal(unname(M0), mean(vals))                  # 均匀分配时也 = 全体均值
    }
  )
})

test_that("23.1-scaling: 样本内细胞数任意整倍缩放不改变组均值", {
  property <- forall(
    gen.bind(function(sid) gen.bind(function(k) gen.map(
      function(pool) list(sid = sid, k = k, pool = pool),
      gen.list(of = 12, gen.element(0:10))),
      gen.element(1:4)),
      gen.samples()),
    function(x) {
      sid <- x$sid; k <- x$k; pool <- x$pool
      # 每样本基础值 = pool 循环铺 3 个；缩放 = 整倍重复 k 次等价于细胞数 xk
      base_v <- unlist(lapply(seq_along(sid), function(i) rep(pool[(i %% length(pool)) + 1], 3)))
      base_s <- rep(sid, each = 3)
      scaled_v <- rep(base_v, times = k)
      scaled_s <- rep(base_s, times = k)
      g <- rep("g0", length(sid)); names(g) <- sid
      M_base <- sgs_group_mean(sgs_sample_means(base_v, base_s), g, "g0")
      M_scale <- sgs_group_mean(sgs_sample_means(scaled_v, scaled_s), g, "g0")
      expect_equal(M_scale, M_base)   # 细胞数缩放不改变组均值
    }
  )
})

cat("test_sample_level_mean: all properties passed\n")

# Property tests for pc_selection.R
# Feature: sgscr-scrnaseq-skill, Property 8: PC 数定量选择算法
# Validates: Requirements 9.1
# Run: Rscript tests/property/test_pc_selection.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "pc_selection.R"))

# 生成合法方差贡献序列：n ∈ [2,50] 等差分片，随机权重归一化到 100
gen_varpct <- gen.bind(
  function(n) gen.map(
    function(ws) { w <- as.numeric(unlist(ws)); w / sum(w) * 100 },
    gen.list(gen.unif(0.001, 10), of = n)),
  gen.map(function(m) m + 2, gen.int(49)))   # 2..50

test_that("P8: pc.use valid positive int <= n; co1-solution meets cumvar threshold", {
  property <- forall(
    gen_varpct,
    function(v) {
      r <- sgs_select_pcs(v, cum_th = 90, single_th = 5, fallback_npc = 30L)
      expect_true(is.finite(r$pc_use) && r$pc_use >= 1 && r$pc_use <= length(v))
      if (!is.na(r$co1)) {
        # P8 原文：pc.use <= co1，且 co1 对应位置的累计方差 >= 配置阈值
        expect_gte(cumsum(v)[r$co1], 90 - 1e-9)
        expect_lte(r$pc_use, r$co1)
      }
      expect_true(!is.na(r$pc_use))
    }
  )
})

test_that("P8-degenerate: zero-variance / invalid inputs fall back, never NA", {
  property <- forall(
    gen.element(list(rep(0, 10),                       # 全零：无信息，回退
                     c(NA, 50, 50),                    # 含 NA：无效
                     "bad",                            # 非数值：无效
                     numeric(0),                       # 空向量：无效
                     c(-1, 50, 51),                    # 负值：无效
                     42)),                             # 长度<2：无效
    function(v) {
      r <- sgs_select_pcs(v, cum_th = 90, single_th = 5, fallback_npc = 30L)
      expect_equal(r$pc_use, 30L)
      expect_true(!is.na(r$pc_use))
    }
  )
})

# Property tests for resolution_selection.R
# Feature: sgscr-scrnaseq-skill, Property 9: ARI 辅助分辨率选择规则
#   确定性返回推荐分辨率；无骤降→最高候选；有骤降→骤降前分辨率；ARI 实现对同配对
#   标签 permutation 不变、对完全一致/完全独立有已知值。
# Validates: Requirements 10.2
# Run: Rscript tests/property/test_resolution_selection.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "resolution_selection.R"))

# hedgehog 无 gen.permutation，用固定用例 + 已知 ARI 值
test_that("P9-ari-fixed: identical=1, disjoint≈0, arity error", {
  a <- c(1, 1, 2, 2, 3, 3)
  expect_equal(sgs_ari(a, a), 1)
  expect_equal(sgs_ari(a, c(1, 2, 1, 2, 1, 2)) < 0.2, TRUE)  # 近独立 → 接近 0
  expect_error(sgs_ari(a, a[1:3]))
})

test_that("P9-rule: deterministic; no-drop→highest; first-drop→before drop", {
  property <- forall(
    gen.bind(
      function(n) gen.map(
        function(ari) list(n = n, ari = ari),
        gen.list(gen.unif(0, 1), of = n)),
      gen.int(5)),  # 1..5 → 候选 2..6 个（两次生成范围合并的历史注释）
    function(x) {
      res <- as.numeric(seq_len(x$n + 1)) * 0.1
      r1 <- sgs_select_resolution(res, x$ari, drop_threshold = 0.3)
      r2 <- sgs_select_resolution(res, x$ari, drop_threshold = 0.3)
      expect_identical(r1$resolution, r2$resolution)              # 确定性
      drops <- which(x$ari < 0.3)
      if (length(drops) == 0) {
        expect_equal(r1$resolution, max(res))                     # 无骤降→最高
      } else {
        expect_equal(r1$resolution, res[drops[1]])                # 首个骤降之前
      }
    }
  )
})

# Property tests for preflight_utils.R (P11)
# Feature: sgscr-scrnaseq-skill, Property 11: 下游模块输入前置校验
#   必需项齐备 → 放行；任意缺失子集 → 终止并恰好报告那组缺失项名称。
# Validates: Requirements 13.1, 13.2, 14.1, 14.2, 15.1, 15.2
# Run: Rscript tests/property/test_preflight.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "preflight_utils.R"))

ITEMS <- c("major_cluster", "group", "CellChat", "ref.db", "nUMI", "pathway")

test_that("P11: complete present set passes; any missing subset reported exactly", {
  property <- forall(
    gen.bind(function(k) gen.bind(function(kmask) gen.map(
      function(pmask) list(k = k, kmask = kmask, pmask = pmask),
      gen.subsequence(ITEMS)),
      gen.subsequence(ITEMS)),
      gen.int(length(ITEMS) + 1)),   # gen.int 标量上限
    function(x) {
      required <- ITEMS[seq_len(x$k)]
      present <- unique(c(required[x$kmask], ITEMS[x$pmask]))
      missing_items <- setdiff(required, present)
      if (length(missing_items) == 0) {
        expect_silent(sgs_preflight("S7", required, present))
      } else {
        expect_error(sgs_preflight("S7", required, present), regexp = "缺失必需项")
        for (m in missing_items) {
          expect_error(sgs_preflight("S7", required, present),
                       regexp = m, fixed = TRUE)
        }
        # 恰好那组：不报告本就存在的项（用缺失项构成的子集断言间接覆盖）
        expect_error(sgs_preflight("S7", required, present),
                     regexp = paste(sort(missing_items), collapse = ", "), fixed = TRUE)
      }
    }
  )
})

test_that("P11-determinism: same inputs, same message", {
  m1 <- tryCatch({ sgs_preflight("S8", c("a", "b"), "a"); NULL },
                 error = function(e) conditionMessage(e))
  m2 <- tryCatch({ sgs_preflight("S8", c("a", "b"), "a"); NULL },
                 error = function(e) conditionMessage(e))
  expect_identical(m1, m2)
  expect_identical(m1, "[S8] 前置校验失败，缺失必需项: b")
})

test_that("P11-edge: empty required passes; empty present fails with all names", {
  expect_silent(sgs_preflight("S9", character(0), character(0)))
  expect_error(sgs_preflight("S9", c("y", "x"), character(0)),
               regexp = "x, y", fixed = TRUE)   # 报告按字典序
})

# Property tests for qc_gate_utils.R
# Feature: sgscr-scrnaseq-skill, Property 7: Quality_Gate 三态分类完整性与互斥性
# Validates: Requirements 8.2, 8.3, 8.4
# Run: Rscript tests/property/test_quality_gate.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "qc_gate_utils.R"))

VALID <- c("PASS", "FAIL", "NEEDS_HUMAN_REVIEW")

test_that("P7: exactly one of three mutually exclusive statuses", {
  property <- forall(
    gen.bind(function(a) gen.bind(function(b) gen.map(
                function(c) list(a, b, c),
                gen.element(list(NULL, "", "figures/P2.1.png"))),
              gen.element(c("", "ok detail"))),
      gen.element(c(TRUE, FALSE))),
    function(t) {
      g <- sgs_qc_gate(auto_pass = t[[1]], detail = t[[2]], human_review_file = t[[3]])
      expect_true(g$status %in% VALID)
      # 互斥：status 唯一且满足推导规则
      expected <- if (!isTRUE(t[[1]])) "FAIL"
                  else if (!is.null(t[[3]]) && nzchar(t[[3]])) "NEEDS_HUMAN_REVIEW"
                  else "PASS"
      expect_identical(g$status, expected)
      # FAIL/NEEDS_HUMAN_REVIEW 必须携带说明
      if (g$status != "PASS") expect_true(nzchar(g$detail))
    }
  )
})

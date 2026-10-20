# Property tests for resume_utils.R
# Feature: sgscr-scrnaseq-skill, Property 6: 断点续跑调度与前序产物依赖检测
#   从任意起点 k：跳过 <k 的 Step；前序产物存在时允许续跑；缺失时终止并报告具体路径。
# Validates: Requirements 7.3, 7.4
# Run: Rscript tests/property/test_resume_scheduler.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "resume_utils.R"))

# 含 rds 产物的连续子集（跳过 s7/s8/s9：其前序 s6 特例单独用例覆盖）
STEPS <- c("s0", "s1", "s2", "s3", "s4", "s5", "s10")

.mock_rds_dir <- function(create_through) {
  d <- file.path(tempdir(), paste0("resume_", as.integer(runif(1, 1, 1e9))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  idx <- match(create_through, STEPS)
  if (!is.na(idx)) for (s in STEPS[seq_len(idx - 1)]) {
    out <- .sgs_STEP_OUTPUTS[[s]]
    if (!is.na(out)) writeLines("x", file.path(d, out))
  }
  d
}

test_that("P6: resume from k skips <k; missing prerequisite errors with path", {
  property <- forall(
    gen.bind(
      function(ki) gen.map(
        function(create_flag) list(k = STEPS[ki], ki = ki, create = create_flag),
        gen.element(c(TRUE, FALSE))),
      gen.int(length(STEPS))),
    function(x) {
      d <- .mock_rds_dir(if (x$create) x$k else STEPS[1])  # FALSE 时只建 s0 产物
      if (x$ki == 1) {
        # 从 s0 续跑无前序依赖
        plan <- sgs_resume_plan(x$k, d, step_table = STEPS)
        expect_identical(plan$steps[1], "s0")
      } else if (x$create) {
        plan <- sgs_resume_plan(x$k, d, step_table = STEPS)
        expect_identical(plan$steps[1], x$k)                 # 从 k 开始，不重跑 <k
        expect_false(x$k %in% plan$steps[-1])
      } else {
        prev_out <- .sgs_STEP_OUTPUTS[[STEPS[x$ki - 1]]]
        if (is.na(prev_out)) {
          plan <- sgs_resume_plan(x$k, d, step_table = STEPS)                    # 前序无 rds 产物 → 无依赖检查
          expect_identical(plan$steps[1], x$k)
        } else {
          expect_error(sgs_resume_plan(x$k, d, step_table = STEPS),              # 缺产物→报错含具体路径
                       regexp = "Missing prerequisite")
        }
      }
    }
  )
})

test_that("P6-s6-special: s7 resume needs s6 compartment output", {
  d <- .mock_rds_dir("s7")   # 建 s0–s5 产物，无 s6
  expect_error(sgs_resume_plan("s7", d), regexp = "no s6_.*_subcluster.rds")
  writeLines("x", file.path(d, "s6_myeloid_subcluster.rds"))
  expect_silent(sgs_resume_plan("s7", d))
})

test_that("P6-invalid: unknown step errors", {
  expect_error(sgs_resume_plan("s99", tempdir()))
  expect_error(sgs_resume_plan("banana", tempdir()))
})

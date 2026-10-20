# Property tests for manifest_utils.R
# Requires: hedgehog, jsonlite
# Run on an R host (remote main env or any R with the two packages):
#   Rscript tests/property/test_manifest_utils.R
#
# Feature: sgscr-scrnaseq-skill, Property 13: Manifest 决议值回写一致性
#   写入 resolutions 后读出的值与实际使用的值完全相等；字段名与 convention 登记表一致。
#   Validates: Requirements 5.5, 9.3, 10.4, 11.5, 17.2, 17.5
#
# Feature: sgscr-scrnaseq-skill, Property 14: Manifest 追加的单调性
#   steps 条目数 == 已完成 Step 数；顺序与完成顺序一致；历史条目不被修改或删除。
#   Validates: Requirements 17.3

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(jsonlite))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "manifest_utils.R"))

gen_value <- gen.choice(
  gen.element(c("Epithelial", "Myeloid", "Fibroblast")),
  gen.int(100),
  gen.unif(0, 100))

.test_manifest <- function() {
  f <- tempfile(fileext = ".json")
  sgs_manifest_init(f, run_id = "test_run", random_seed = 42,
                    species = "human", singleR_reference = "HPCA",
                    config_snapshot_path = "cfg.yaml")
  f
}

test_that("P13: resolved value round-trips exactly", {
  kv_gen <- gen.bind(
    function(k) gen.map(function(v) list(key = k, value = v), gen_value),
    gen.element(c("pc_used", "resolution", "reference")))
  property <- forall(
    kv_gen,
    function(kv) {
      f <- .test_manifest()
      sgs_manifest_set_resolved(f, kv$key, kv$value)
      m <- jsonlite::fromJSON(f, simplifyVector = FALSE)
      expect_equal(m$resolutions[[kv$key]], kv$value)   # P13语义为值相等；JSON往返整型会变双精度
    }
  )
})

test_that("P14: steps append-only, ordered, immutable", {
  steps_gen <- gen.bind(
    function(n) gen.list(gen.element(c("s0_env_check", "s1_import", "s2_qc", "s3_dimreduce")), of = n),
    gen.int(10))
  property <- forall(
    steps_gen,
    function(step_ids) {
      f <- .test_manifest()
      before <- character(0)
      ok <- TRUE
      for (i in seq_along(step_ids)) {
        sid <- step_ids[[i]]
        sgs_manifest_append_step(
          f, step_id = sid, env_name = "main",
          quality_gate = list(status = "PASS", detail = "t"),
          started_at = "2026-10-02T09:00:00+0800",
          finished_at = "2026-10-02T09:00:01+0800")
        m <- jsonlite::fromJSON(f, simplifyVector = FALSE)
        n <- length(m$steps)
        if (n != i) ok <- FALSE
        if (m$steps[[n]]$step_id != sid) ok <- FALSE
        # 历史条目未被改写
        if (i > 1 && m$steps[[i - 1]]$step_id != before[i - 1]) ok <- FALSE
        before[i] <- sid
      }
      expect_true(ok)
    }
  )
})

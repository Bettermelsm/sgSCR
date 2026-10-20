# Property tests for env_resolver.R
# Feature: sgscr-scrnaseq-skill, Property 4: 模块名到隔离环境的解析与缺失检测
#   对任意已注册模块标识，解析返回 Config 配置的环境名；未注册模块 / 环境不存在时
#   终止并报告缺失环境名，不静默回退主环境。
#   Validates: Requirements 4.4, 4.5
# Run: Rscript tests/property/test_env_resolver.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "env_resolver.R"))

STEPS <- names(.sgs_STEP_ENV)   # s0..s11
ENVS_CFG <- list(main = "sgscr_main", cellchat = "sgscr_cellchat",
                 trajectory = "sgscr_trajectory", scenic = "sgscr_scenic")

test_that("P4-registered: every registered step resolves to its configured env", {
  property <- forall(
    gen.element(STEPS),
    function(st) {
      expect_equal(sgs_resolve_env(st, ENVS_CFG), ENVS_CFG[[.sgs_STEP_ENV[[st]]]])
    }
  )
})

test_that("P4-unregistered: unknown step errors, never falls back to main", {
  property <- forall(
    gen.element(c("s12", "s99", "sx", "s7x_bogus", "step3")),
    function(bad) {
      expect_error(sgs_resolve_env(bad, ENVS_CFG), regmatches = FALSE)
    }
  )
})

test_that("P4-missing-env: nonexistent env errors with env name, no silent fallback", {
  property <- forall(
    gen.element(c("s7", "s8", "s9")),
    function(st) {
      env <- sgs_resolve_env(st, ENVS_CFG)
      expect_error(sgs_require_env(env, function(e) FALSE),
                   regexp = env, fixed = TRUE)   # 报错信息含缺失环境名
      expect_silent(sgs_require_env(env, function(e) TRUE))
    }
  )
})

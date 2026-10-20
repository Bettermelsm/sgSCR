# Property tests for P15: 随机种子可复现性
#   P15: 同一 random_seed 下，轻量 mock Step 的随机依赖输出完全一致；
#        不同 seed 产生不同输出（在足够大的样本下必然不同）。
#   按 design.md：P15 用轻量 mock 逻辑验证，不触发真实 Seurat/SCENIC 大计算。
# Run: Rscript tests/property/test_seed_reproducibility.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

# 轻量 mock Step：消费 cfg$run$random_seed 的随机依赖计算（模拟 set.seed + 采样管线）
mock_step_run <- function(seed, n = 500L) {
  set.seed(seed)
  vals <- sample(1e6, n, replace = TRUE)
  clusters <- sample(c("A", "B", "C"), n, replace = TRUE)
  list(sum = sum(vals), first5 = vals[1:5], tab = as.list(table(clusters)))
}

test_that("P15: same seed -> bit-identical mock step output (arbitrary seeds)", {
  property <- forall(
    gen.bind(function(seed) gen.map(
      function(n) list(seed = seed, n = n),
      gen.element(100L:2000L)),
      gen.int(100000L)),
    function(x) {
      a <- mock_step_run(x$seed, x$n)
      b <- mock_step_run(x$seed, x$n)
      expect_identical(a, b)
    }
  )
})

test_that("P15: different seeds -> different mock step output", {
  property <- forall(
    gen.element(1L:1000L),
    function(seed) {
      a <- mock_step_run(seed)
      b <- mock_step_run(seed + 7919L)   # 质数步进，保证样本差异
      expect_false(identical(a$first5, b$first5))
    }
  )
})

test_that("P15-static: every step script sets the configured seed as its RNG source", {
  # testthat 的求值环境不一定继承顶层 HERE，这里独立重算
  here <- dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
  scripts <- list.files(file.path(here, "..", "..", "scripts"),
                        pattern = "^s[0-9]+_.*\\.R$", full.names = TRUE)
  expect_true(length(scripts) >= 8)
  for (f in scripts) {
    txt <- paste(readLines(f), collapse = "\n")
    expect_match(txt, "set\\.seed\\(cfg\\$run\\$random_seed\\)",
                 fixed = FALSE,
                 info = sprintf("%s 未使用 cfg$run$random_seed 作为种子源", basename(f)))
  }
})

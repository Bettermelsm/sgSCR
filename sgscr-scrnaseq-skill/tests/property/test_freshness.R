# Property tests for freshness check (P11 增强, 阶段 11)
#   对任意 (上游 mtime, 下游 mtime) 偏移组合：下游不早于上游 → 放行；
#   下游早于任一上游 → 终止并点名恰好那组过期上游；文件缺失 → 终止并点名。
# Run: Rscript tests/property/test_freshness.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "preflight_utils.R"))

new_tmp <- function() { f <- tempfile(fileext = ".csv"); writeLines("x", f); f }

test_that("freshness: downstream newer than all upstreams passes", {
  property <- forall(
    gen.bind(function(up_off) gen.map(
      function(n_up) list(up_off = up_off, n_up = n_up),
      gen.element(1L:3L)),
      gen.element(1L:3600L)),   # 下游比每个上游新 up_off 秒
    function(x) {
      ups <- replicate(x$n_up, new_tmp())
      dn <- new_tmp()
      base <- as.numeric(Sys.time()) - 7200
      for (u in ups) Sys.setFileTime(u, base)
      Sys.setFileTime(dn, base + x$up_off + 1)
      expect_silent(sgs_check_freshness(dn, ups))
      for (u in c(ups, dn)) unlink(u)
    }
  )
})

test_that("freshness: stale upstream is named exactly (any position among several)", {
  property <- forall(
    gen.bind(function(n_up) gen.map(
      function(stale_i) list(stale_i = stale_i, n_up = n_up),
      gen.element(seq_len(n_up))),   # 第 stale_i 个上游比下游新（stale_i ≤ n_up）
      gen.element(1L:3L)),
    function(x) {
      ups <- replicate(x$n_up, new_tmp())
      dn <- new_tmp()
      base <- as.numeric(Sys.time()) - 7200
      Sys.setFileTime(dn, base)
      for (i in seq_along(ups)) {
        Sys.setFileTime(ups[i], if (i == x$stale_i) base + 100 else base - 100)
      }
      expect_error(sgs_check_freshness(dn, ups, tol_sec = 0),
                   class = "error")
      err <- tryCatch(sgs_check_freshness(dn, ups, tol_sec = 0),
                      error = function(e) conditionMessage(e))
      expect_match(err, basename(ups[x$stale_i]), fixed = TRUE)     # 点名过期上游
      for (u in c(ups, dn)) unlink(u)
    }
  )
})

test_that("freshness: missing file is named; fresh chain OK", {
  good <- new_tmp(); dn <- new_tmp()
  Sys.setFileTime(good, 1000); Sys.setFileTime(dn, 2000)
  expect_silent(sgs_check_freshness(dn, good))
  expect_error(sgs_check_freshness(dn, file.path(tempdir(), "不存在_x.csv")), "不存在")
  unlink(c(good, dn))
})

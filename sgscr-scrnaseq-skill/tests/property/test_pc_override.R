# Property tests for sgs_apply_pc_override（阶段 11：人工 PC 覆盖分支）
#   对任意 (合法方差序列, 合法 override)：pc.use == override 且 basis == "manual_override"，
#   自动算法建议保留在 auto 字段且 applied=FALSE；override 缺失时行为与纯自动一致；
#   非法 override（非整数/越界）必须报错，不得静默采用。
# Run: Rscript tests/property/test_pc_override.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "pc_selection.R"))

test_that("override: valid override replaces auto result, auto kept for trace", {
  property <- forall(
    gen.bind(function(pc) gen.map(
      function(ov) list(pc = pc, ov = ov),
      gen.element(2L:50L)),
      gen.element(3L:50L)),
    function(x) {
      var_pct <- rep(100 / (x$pc + 5), x$pc + 5)   # 平坦序列：自动算法必有确定解
      auto <- sgs_select_pcs(var_pct, cum_th = 90, single_th = 5,
                             elbow_diff_threshold = 0.1, fallback_npc = 30L)
      sel <- sgs_apply_pc_override(auto, x$ov, npcs_compute = 50L)
      expect_identical(as.integer(sel$pc_use), as.integer(x$ov))
      expect_identical(sel$basis, "manual_override")
      expect_false(is.null(sel$auto))
      expect_identical(as.integer(sel$auto$pc_use), as.integer(auto$pc_use))
      expect_identical(sel$auto$basis, auto$basis)
      expect_false(sel$auto$applied)
    }
  )
})

test_that("override: null/NA override returns pure auto result unchanged", {
  property <- forall(
    gen.element(2L:60L),
    function(n) {
      var_pct <- rep(100 / n, n)
      auto <- sgs_select_pcs(var_pct, cum_th = 90, single_th = 5,
                             elbow_diff_threshold = 0.1, fallback_npc = 30L)
      s1 <- sgs_apply_pc_override(auto, NULL, npcs_compute = 50L)
      s2 <- sgs_apply_pc_override(auto, NA, npcs_compute = 50L)
      expect_identical(s1, auto)
      expect_identical(s2, auto)
      expect_true(is.null(s1$auto))   # 未覆盖时不应出现 auto 对比字段
    }
  )
})

test_that("override: invalid override must error, never silently applied", {
  auto <- list(pc_use = 3L, basis = "elbow", co1 = 10L, co2 = 3L)
  for (bad in list(0, 1, 51, 3.5, "abc", c(10, 20))) {
    expect_error(sgs_apply_pc_override(auto, bad, npcs_compute = 50L),
                 class = "error", info = paste(bad, collapse = ","))
  }
})

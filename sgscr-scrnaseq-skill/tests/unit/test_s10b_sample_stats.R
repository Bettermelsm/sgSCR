# Unit tests for s10b building blocks: LOO set, low-n handling, min attainable p (Req 23.2/23.5/23.6)
# 判据取 design.md 测试映射表「23.2」「23.5」「23.6」行
# Run: Rscript tests/unit/test_s10b_sample_stats.R

suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))

test_that("23.2: LOO 结果集条目数等于样本数，每条恰好剔除一个样本", {
  sids <- c("NC1", "NC3", "UC1", "UC2")
  loo <- sgs_loo_samples(sids)
  expect_length(loo, length(sids))
  expect_setequal(names(loo), sprintf("LOO_%s", sids))
  for (i in seq_along(loo)) {
    dropped <- sub("^LOO_", "", names(loo)[i])
    expect_false(dropped %in% loo[[i]])
    expect_setequal(loo[[i]], setdiff(sids, dropped))   # 其余样本全保留
  }
})

test_that("23.5: 输出契约不含细胞级 p 值列；最小可达双侧 p 值正确", {
  cols <- sgs_sample_level_columns(c("NC1", "UC1"))
  expect_false(any(grepl("p_value|pvalue|p_value_cell", cols, ignore.case = TRUE)))
  # 最小双侧 p = 2 / choose(n0+n1, n0)：n0=n1=3 → 2/20=0.1；n0=3,n1=4 → 2/35
  expect_equal(sgs_min_attainable_pvalue(3, 3), 0.1)
  expect_equal(sgs_min_attainable_pvalue(3, 4), 2 / 35)
  expect_true(is.na(sgs_min_attainable_pvalue(0, 3)))
})

test_that("23.6: 低细胞数样本在明细中被标记、汇总中被排除，n_samples_used 相应减少", {
  # (cell_group, sample) 细胞数：NC3 只有 5 个细胞（< min_cells=10）
  n_cells <- c(NC1 = 50, NC3 = 5, UC1 = 40, UC2 = 30)
  flag <- sgs_low_n_flag(n_cells, min_cells_per_sample = 10)
  expect_true(flag[["NC3"]])
  expect_false(flag[["NC1"]])
  # 汇总：排除 NC3 后该组样本数减少
  g <- c(NC1 = "ctrl", NC3 = "ctrl", UC1 = "case", UC2 = "case")
  m_s <- c(NC1 = 1.0, NC3 = 99.0, UC1 = 0.5, UC2 = 0.4)   # NC3 若入汇总必拉爆组均值
  st <- sgs_separation_stats(m_s, g, "ctrl", "case", min_separation_margin = 0.01)
  st_excl <- sgs_separation_stats(m_s[names(m_s) != "NC3"], g[names(g) != "NC3"],
                                  "ctrl", "case", min_separation_margin = 0.01)
  expect_equal(st$n_samples_used_group0, 2)   # 全量参考：含 NC3
  expect_equal(st_excl$n_samples_used_group0, 1)
  expect_false(isTRUE(all.equal(st_excl$mean_group0, 50)))  # 排除后不被 99 拉爆
  expect_equal(st_excl$mean_group0, 1.0)
})

cat("test_s10b_sample_stats: all tests passed\n")

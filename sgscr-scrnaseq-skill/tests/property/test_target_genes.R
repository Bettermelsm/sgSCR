# Property tests for target_gene_utils.R (P12)
# Feature: sgscr-scrnaseq-skill
#   P12: 目标基因交集过滤——过滤结果恰等于清单∩结果表；交集空返回空表而非报错或全表 (Req 16.2)
# Run: Rscript tests/property/test_target_genes.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "target_gene_utils.R"))

gen.gene <- function() gen.map(function(i) paste0("G", i), gen.int(1e6))

test_that("P12-main: filtered set == intersection, order preserved, duplicates dropped", {
  property <- forall(
    gen.bind(function(pool) gen.bind(function(idx) gen.map(
      function(dup) {
        targets <- unique(pool[as.integer(unlist(idx)) %% length(pool) + 1])   # 目标清单为基因池的随机子集
        list(pool = pool, targets = rep(targets, times = dup))
      },
      gen.element(1:3)),
      gen.list(of = 8, gen.int(1e6))),
      gen.list(of = 10, gen.gene())),
    function(x) {
      x$targets <- as.character(unlist(x$targets)); x$pool <- as.character(unlist(x$pool))
      kept <- sgs_target_gene_filter(x$targets, x$pool)
      expect_setequal(kept, intersect(unique(x$targets), x$pool))
      # 顺序保持 targets 首次出现顺序
      expect_identical(kept, unique(x$targets[x$targets %in% x$pool]))
    }
  )
})

test_that("P12-empty: empty intersection -> empty result, no error, never full table", {
  property <- forall(
    gen.bind(function(a) gen.bind(function(b) gen.map(
      function(i) list(a = a, b = b, pick = i),
      gen.element(c("a_only", "b_only"))),
      gen.list(of = 10, gen.gene())),
      gen.list(of = 10, gen.gene())),
    function(x) {
      disjoint_a <- setdiff(x$a, x$b)
      disjoint_b <- setdiff(x$b, x$a)
      if (length(disjoint_a) == 0 || length(disjoint_b) == 0) expect_success()
      else if (x$pick == "a_only") {
        expect_identical(sgs_target_gene_filter(disjoint_a, x$b), character(0))
      } else {
        expect_identical(sgs_target_gene_filter(x$a, disjoint_b), character(0))
      }
    }
  )
  # 确定性：空交集显式用例
  expect_identical(sgs_target_gene_filter(c("A", "B"), c("C", "D")), character(0))
  # 不报错
  expect_silent(sgs_target_gene_filter(c("A"), c("C")))
})

test_that("P12-edge: NA/blank targets ignored; empty inputs handled", {
  expect_identical(sgs_target_gene_filter(c("A", NA, ""), c("A", "B")), "A")
  expect_identical(sgs_target_gene_filter(character(0), c("A")), character(0))
  expect_identical(sgs_target_gene_filter(c("A"), character(0)), character(0))
})

test_that("P12-determinism: repeated runs identical", {
  r1 <- sgs_target_gene_filter(c("B", "A", "B"), c("A", "B", "C"))
  r2 <- sgs_target_gene_filter(c("B", "A", "B"), c("A", "B", "C"))
  expect_identical(r1, r2)
  expect_identical(r1, c("B", "A"))   # 保持清单顺序
})

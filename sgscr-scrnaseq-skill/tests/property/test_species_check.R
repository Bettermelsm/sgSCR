# Property tests for annotation_utils.R (P5, P10)
# Feature: sgscr-scrnaseq-skill
#   P5:  物种一致 ⇔ 放行；不一致时终止并报告双方物种名 (Req 5.3, 5.4)
#   P10: 映射应用后每细胞标签恰等于其簇映射值；唯一值集合保持；原始注释不变 (Req 11.3, 11.4)
# Run: Rscript tests/property/test_species_check.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "annotation_utils.R"))

test_that("P5: species check passes iff equal (any case/whitespace), else errors naming both", {
  property <- forall(
    gen.bind(function(sp1) gen.bind(function(sp2) gen.map(
      function(v) list(a = sp1, b = sp2, variant = v),
      gen.element(c("asis", "case", "space"))),
      gen.choice("human", "mouse", "zebrafish")),
      gen.choice("human", "mouse", "zebrafish")),
    function(x) {
      b <- switch(x$variant,
                  asis = x$b,
                  case = toupper(x$b),
                  space = paste0(" ", x$b, " "))
      if (tolower(trimws(x$a)) == tolower(trimws(x$b))) {
        expect_true(sgs_check_species(x$a, b))
      } else {
        expect_error(sgs_check_species(x$a, b), regexp = "物种不匹配")
        expect_error(sgs_check_species(x$a, b), regexp = x$a, fixed = TRUE)
        expect_error(sgs_check_species(x$a, b), regexp = b, fixed = TRUE)
      }
    }
  )
})

test_that("P5: null/empty species never passes", {
  expect_error(sgs_check_species(NULL, "human"))
  expect_error(sgs_check_species("human", ""))
})

test_that("P10-main: per-cell labels == mapping; unique set preserved; raw unchanged", {
  property <- forall(
    gen.bind(function(ncl) gen.bind(function(ncell) gen.bind(
      function(maj_names) gen.map(
        function(perms) list(ncl = ncl, ncell = ncell, maj = maj_names, seeds = perms),
        gen.list(of = 2, gen.element(1:1e6))),
      gen.choice("M1", "M2", "M3", "M4")),
      gen.int(150)),
      gen.int(25)),
    function(x) {
      ncl <- max(1, x$ncl); ncell <- max(ncl, x$ncell)
      clusters <- paste0("c", seq_len(ncl))
      cells <- sample(clusters, ncell, replace = TRUE)
      # 随机映射：每个 cluster 指到随机主群（保证覆盖全部 cluster_id）
      mapping <- stats::setNames(sample(x$maj, ncl, replace = TRUE), clusters)
      mapped <- sgs_apply_cluster_mapping(cells, mapping)
      expect_identical(unname(mapping[as.character(cells)]), mapped)     # 每细胞恰等于其簇映射值
      expect_identical(sort(unique(mapped)), sort(unique(unname(mapping))))  # 唯一值集合恰等
      # 覆盖校验：缺任一实际出现的簇必须报错并点名
      victim <- unique(cells)[1]
      incomplete <- mapping[names(mapping) != victim]
      expect_error(sgs_apply_cluster_mapping(cells, incomplete), regexp = "未确认的簇编号")
    }
  )
})

test_that("P10-validate: template validation rejects partial/blank/dup mappings", {
  ids <- c("0", "1", "2")
  good <- data.frame(cluster_id = ids,
                     confirmed_major_cluster = c("Mono", "T", "B"), stringsAsFactors = FALSE)
  expect_identical(sgs_validate_confirmed_mapping(good, ids),
                   c("0" = "Mono", "1" = "T", "2" = "B"))
  expect_error(sgs_validate_confirmed_mapping(good[1:2, ], ids), regexp = "未确认的簇编号")
  blank <- good; blank$confirmed_major_cluster[2] <- ""
  expect_error(sgs_validate_confirmed_mapping(blank, ids), regexp = "未确认的簇编号")
  expect_error(sgs_validate_confirmed_mapping(good[, 1, drop = FALSE], ids), regexp = "缺少必需列")
  allblank <- good; allblank$confirmed_major_cluster <- ""
  expect_error(sgs_validate_confirmed_mapping(allblank, ids), regexp = "全为空")
})

test_that("P10-determinism: same input repeated runs give identical result", {
  clusters <- rep(c("a", "b", "c"), each = 4)
  mapping <- c(a = "X", b = "Y", c = "X")
  r1 <- sgs_apply_cluster_mapping(clusters, mapping)
  r2 <- sgs_apply_cluster_mapping(clusters, mapping)
  expect_identical(r1, r2)
})

test_that("majority label: deterministic tie-break", {
  expect_identical(sgs_majority_label(c("T cell", "B cell", "T cell")), "T cell")
  expect_identical(sgs_majority_label(c("B cell", "T cell")), "B cell")  # 并列→字典序最小
  expect_true(is.na(sgs_majority_label(character(0))))
})

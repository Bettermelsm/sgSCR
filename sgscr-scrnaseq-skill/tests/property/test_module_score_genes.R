# Property tests for module score gene usage marking (Req 29.3) + 29.1/29.5 checks
# Feature: sgscr-scrnaseq-skill, Property (design 29.3): 对任意「基因集合 × 表达矩阵行名」
#   组合，缺失基因必须出现在 missing 标记中而非被静默丢弃；参与计算基因列表等于交集
# Run: Rscript tests/property/test_module_score_genes.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "sample_stats_utils.R"))

gen.gene <- function() gen.map(function(i) paste0("G", i), gen.int(1e6))

test_that("29.3-property: missing 全标记、used == 交集、二者不相交、并集 == 集合", {
  property <- forall(
    gen.bind(function(pool) gen.bind(function(extra) gen.map(
      function(dup) {
        gset <- unique(c(sample(pool, 4), extra))
        gset <- rep(gset, times = dup)   # 集合含重复也应收敛
        list(pool = pool, gset = gset)
      },
      gen.element(1:2)),
      gen.list(of = 3, gen.gene())),
      gen.list(of = 12, gen.gene())),
    function(x) {
      pool <- unique(as.character(x$pool))
      gset <- unique(as.character(x$gset))
      u <- sgs_module_gene_usage(gset, pool)
      expect_setequal(u$missing, setdiff(gset, pool))   # 缺失必须全标记
      expect_setequal(u$used, intersect(gset, pool))    # 参与计算 == 交集
      expect_length(intersect(u$used, u$missing), 0)
      expect_setequal(c(u$used, u$missing), gset)
      TRUE
    }
  )
})

test_that("29.1-static+29.5-static: s10 评分路径只写数值列、无离散标签回写", {
  lines <- readLines(file.path(HERE, "..", "..", "scripts", "s10_target_genes.R"))
  seg <- paste(lines, collapse = "\n")
  # 静态检查：评分写回仅 AddModuleScore 的 *_score1 数值列
  expect_true(grepl("AddModuleScore", seg))
  expect_true(grepl("is.numeric\\(sc\\)", seg))
  # 不得出现按阈值切分评分并写回字符/因子标签的分支（需求 29.4/29.5）
  expect_false(grepl("cut\\(.*score|score.*<-\\s*ifelse|score_label|as.factor\\(.*score",
                     seg, ignore.case = TRUE))
  # 缺失基因走 genes_used.csv 落盘（只标记不静默）
  expect_true(grepl("_genes_used\\.csv", seg))
})

cat("test_module_score_genes: all tests passed\n")

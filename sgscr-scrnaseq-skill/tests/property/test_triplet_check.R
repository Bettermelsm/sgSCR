# Property tests for triplet_check.R
# Feature: sgscr-scrnaseq-skill, Property 1: 三联文件完整性校验
#   对任意样本目录与任意"缺失文件子集"（{barcodes,features,matrix}.tsv.gz 的所有子集，含空/全集），
#   校验函数准确报告恰好那组缺失文件名；子集为空时判定通过。
#   Validates: Requirements 3.4, 3.5
# Run: Rscript tests/property/test_triplet_check.R  (on R host with hedgehog/testthat)

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "triplet_check.R"))

TRIPLET <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

test_that("P1: missing subset reported exactly, empty subset passes", {
  # 用子集的位掩码枚举全部 8 种缺失组合（全集空间小，穷举优于随机）
  ok <- TRUE
  for (mask in 0:7) {
    dir <- file.path(tempdir(), sprintf("triplet_case_%d", mask))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(TRIPLET)) {
      if (!bitwAnd(mask, 2^(i - 1))) writeLines("x", file.path(dir, TRIPLET[i]))
    }
    missing <- sgs_check_triplet(dir, TRIPLET)
    expected <- TRIPLET[as.logical(bitwAnd(mask, c(1, 2, 4)))]
    if (!setequal(missing, expected)) ok <- FALSE
    gate <- sgs_triplet_gate(list(S1 = missing))
    if (!identical(gate$pass, mask == 0)) ok <- FALSE   # 仅空缺失子集通过 gate
  }
  expect_true(ok)
})

test_that("P1-batch: excluded samples skipped, others checked", {
  dir <- file.path(tempdir(), "triplet_batch")
  dir.create(dir, showWarnings = FALSE)
  dir.create(file.path(dir, "NC1"), showWarnings = FALSE)
  dir.create(file.path(dir, "UC4"), showWarnings = FALSE)  # excluded 但目录在
  writeLines("x", file.path(dir, "NC1", "barcodes.tsv.gz"))
  # NC1 缺 2 个文件；UC4 在 excluded 中应被跳过
  res <- sgs_check_all_samples(
    dir,
    list(list(sample_id = "NC1", group = "NC"),
         list(sample_id = "UC4", group = "UC")),
    TRIPLET,
    excluded = list(list(sample_id = "UC4", reason = "test")))
  expect_length(res, 1)                      # UC4 未被检查
  expect_length(res[["NC1"]], 2)             # NC1 报告恰好 2 个缺失
})

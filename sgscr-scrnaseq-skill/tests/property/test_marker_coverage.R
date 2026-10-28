# Property tests for assay_utils（阶段 11：s5 缺陷固化修复）
#   sgs_drop_scale_layers：含分样本残留 layer 时重建为 counts/data；无残留时原样返回；
#   sgs_check_marker_coverage：任意 (markers 覆盖子集, 实际簇集) —— 缺失集合被精确点名，
#   覆盖数/期望数正确；空 marker 表时全部簇判缺失。
# Run: Rscript tests/property/test_marker_coverage.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "assay_utils.R"))

mk_obj <- function(n_s) {
  # 用与 S1 完全相同的成因构造残留：逐样本 Normalize+Scale 后 merge，
  # merge 会把各样本 scale.data 以 <layer>.<样本名> 分层拼接（真实残留机制）
  per <- lapply(seq_len(n_s), function(i) {
    m <- Matrix::rsparsematrix(20, 10, 0.3)
    dimnames(m) <- list(sprintf("G%d", 1:20), sprintf("s%dc%d", i, 1:10))
    o <- CreateSeuratObject(counts = m, project = sprintf("S%d", i))
    o <- NormalizeData(o, verbose = FALSE)
    ScaleData(o, verbose = FALSE)
  })
  if (n_s == 1) return(per[[1]])
  merge(per[[1]], y = per[-1])
}

test_that("drop_scale_layers: split scale layers removed, counts/data intact", {
  property <- forall(
    gen.bind(function(n_cl) gen.map(
      function(n_s) list(n_cl = n_cl, n_s = n_s),
      gen.element(2L:4L)),
      gen.element(2L:5L)),
    function(x) {
      obj <- mk_obj(x$n_s)
      lay <- Layers(obj[["RNA"]])
      expect_true(any(grepl("^scale.data\\.", lay)))   # 造出残留
      ref <- JoinLayers(obj, assay = "RNA")   # 参照：仅合并 counts/data 后取 counts 总和
      counts_before <- sum(GetAssayData(ref[["RNA"]], layer = "counts"))
      obj2 <- sgs_drop_scale_layers(obj, assay = "RNA")
      expect_setequal(Layers(obj2[["RNA"]]), c("counts", "data"))
      expect_true(inherits(obj2[["RNA"]], "Assay5"))              # (a) 输出必须 Assay5（0008）
      expect_error(SeuratObject::JoinLayers(obj2, assay = "RNA"), NA)  # (b) 下游 JoinLayers 不报错（0007 回归）
      obj3 <- sgs_drop_scale_layers(obj2, assay = "RNA")          # (c) 幂等：重复调用不报错
      expect_setequal(Layers(obj3[["RNA"]]), c("counts", "data"))
      expect_identical(sum(GetAssayData(obj3[["RNA"]], layer = "counts")),
                       sum(GetAssayData(obj2[["RNA"]], layer = "counts")))
      expect_identical(sum(GetAssayData(obj2[["RNA"]], layer = "counts")), counts_before)  # (d)
      expect_identical(sum(GetAssayData(obj2[["RNA"]], layer = "data")),
                       sum(GetAssayData(ref[["RNA"]], layer = "data")))   # (d) data 同样保持
      expect_identical(as.integer(dim(obj2)[2]), as.integer(dim(obj)[2]))      # (e) 细胞不丢
      expect_identical(as.integer(dim(obj2)[1]), as.integer(dim(obj)[1]))      # (e) 基因不丢
    }
  )
})

test_that("drop_scale_layers: no residue -> object returned unchanged", {
  m <- Matrix::rsparsematrix(20, 25, 0.3)
  dimnames(m) <- list(sprintf("G%d", 1:20), sprintf("c%d", 1:25))
  obj <- CreateSeuratObject(counts = m)
  obj <- NormalizeData(obj, verbose = FALSE)   # 仅 counts/data，无任何 scale 残留
  lay0 <- Layers(obj[["RNA"]])
  obj2 <- sgs_drop_scale_layers(obj, assay = "RNA")
  expect_setequal(Layers(obj2[["RNA"]]), lay0)
  # 不存在的 assay：原样返回，不报错
  expect_identical(sgs_drop_scale_layers(obj, assay = "不存在的assay"), obj)
})

test_that("marker coverage: missing clusters named exactly, any subset", {
  property <- forall(
    gen.bind(function(n_cl) gen.map(
      function(k) list(n_cl = n_cl, k = k),
      gen.element(0L:n_cl)),
      gen.element(3L:8L)),
    function(x) {
      ids <- as.character(seq_len(x$n_cl) - 1)         # 0..n-1，模拟真实簇号
      missing_ids <- if (x$k == 0) ids else ids[sample.int(x$n_cl, x$k)]
      covered <- setdiff(ids, missing_ids)
      mk <- function(cl) data.frame(gene = "X", cluster = cl, stringsAsFactors = FALSE)
      markers <- if (length(covered))
        do.call(rbind, lapply(covered, mk)) else mk(NA)[0, ]
      cov <- sgs_check_marker_coverage(markers, ids)
      expect_setequal(cov$missing, missing_ids)   # 恰好点名缺失集
      expect_identical(as.integer(cov$n_covered), length(covered))
      expect_identical(as.integer(cov$n_expected), length(ids))
      expect_identical(cov$ok, length(missing_ids) == 0)
    }
  )
})

test_that("marker coverage: empty marker table -> all clusters missing, full row set", {
  m0 <- data.frame(gene = character(0), cluster = character(0))
  cov <- sgs_check_marker_coverage(m0, c("0", "1", "2"))
  expect_false(cov$ok)
  expect_setequal(cov$missing, c("0", "1", "2"))
  cov2 <- sgs_check_marker_coverage(NULL, c("0", "1"))
  expect_false(cov2$ok)
})

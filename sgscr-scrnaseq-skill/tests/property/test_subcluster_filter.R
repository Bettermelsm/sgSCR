# Property tests for subcluster_filter_utils（P11D-c 派工项 1）
#   sgs_subcluster_validate：CSV 漏任一亚群必 FAIL 且精确点名；多余亚群 FAIL；重复/非法 keep FAIL
#   sgs_subcluster_keep_ids：main={yes,sens}、strict={yes}（任意 keep 组合）
#   sgs_filter_subobject + sgs_filter_count_check：随机 keep 组合下，
#     输出细胞数 == 保留亚群细胞数之和；subcluster_label 全部非 NA 且可回溯
# Run: Rscript tests/property/test_subcluster_filter.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "subcluster_filter_utils.R"))

mk_obj <- function(ids, prefix = "c") {
  n <- length(ids) * 20
  m <- Matrix::rsparsematrix(50, n, 0.2)
  m@x <- as.numeric(m@x > 0)
  dimnames(m) <- list(sprintf("G%d", 1:50), sprintf("%s%d_%d", prefix, rep(as.integer(ids), each = 20), 1:n))
  o <- CreateSeuratObject(counts = m)
  o@meta.data[[res_col]] <- as.character(rep(ids, each = 20))
  o
}
res_col <- "RNA_snn_res.1"

test_that("validate: missing/extra/duplicate/illegal-keep subclusters are each named exactly", {
  property <- forall(
    gen.bind(function(n_cl) gen.map(
      function(k) list(n_cl = n_cl, k = k),
      gen.element(1L:(n_cl - 1L))),
      gen.element(2L:5L)),
    function(x) {
      ids <- as.character(seq_len(x$n_cl) - 1)
      missing_ids <- sort(ids[sample.int(x$n_cl, x$k)])   # 与 validate 的 setdiff 顺序（ids 升序）一致
      covered <- setdiff(ids, missing_ids)
      csv <- data.frame(compartment = "myeloid", subcluster = covered,
                        keep = "yes", proposed_label = "L", stringsAsFactors = FALSE)
      v <- sgs_subcluster_validate(csv, "myeloid", list(myeloid = ids))
      expect_false(v$ok)
      expect_setequal(v$problems, sprintf("myeloid: CSV 缺亚群 %s", paste(missing_ids, collapse = ",")))
      # 多余亚群
      csv_extra <- rbind(csv, data.frame(compartment = "myeloid", subcluster = "99",
                                         keep = "yes", proposed_label = "X"))
      v2 <- sgs_subcluster_validate(csv_extra, "myeloid", list(myeloid = ids))
      expect_false(v2$ok)
      expect_true(any(grepl("多出亚群 99", v2$problems)))
      # 非法 keep
      csv_bad <- csv; csv_bad$keep[1] <- "maybe"
      v3 <- sgs_subcluster_validate(csv_bad, "myeloid", list(myeloid = ids))
      expect_false(v3$ok)
      expect_true(any(grepl("非法 keep", v3$problems)))
    }
  )
})

test_that("keep_ids: main keeps yes+sens, strict keeps only yes, any combination", {
  property <- forall(
    gen.element(2L:6L),
    function(n) {
      ids <- as.character(seq_len(n) - 1)
      keeps <- sample(c("yes", "sens", "no"), n, replace = TRUE)
      km <- stats::setNames(keeps, ids)
      expect_setequal(sgs_subcluster_keep_ids(km, "main"), ids[keeps %in% c("yes", "sens")])
      expect_setequal(sgs_subcluster_keep_ids(km, "strict"), ids[keeps == "yes"])
    }
  )
})

test_that("filter: cell count == sum of kept subclusters; labels traceable", {
  property <- forall(
    gen.bind(function(n_cl) gen.map(
      function(sel_mask) list(n_cl = n_cl, sel = sel_mask),
      gen.element(1L:(2^n_cl))),
      gen.element(3L:5L)),
    function(x) {
      ids <- as.character(seq_len(x$n_cl) - 1)
      obj <- mk_obj(ids)
      km <- stats::setNames(rep("no", x$n_cl), ids)
      km[sample.int(x$n_cl, max(1, x$sel %% x$n_cl))] <- "yes"
      keep_ids <- sgs_subcluster_keep_ids(km, "strict")
      labels <- stats::setNames(paste0("L", ids), ids)
      n_expected <- sum(as.character(obj@meta.data[[res_col]]) %in% keep_ids)
      filt <- sgs_filter_subobject(obj, res_col, keep_ids, labels)
      chk <- sgs_filter_count_check(filt, res_col, keep_ids, n_expected)
      expect_true(chk$ok)
      expect_identical(as.integer(ncol(filt)), as.integer(n_expected))   # 细胞数恒等
      expect_false(any(is.na(filt$subcluster_label)))                    # 标签全部落上
      expect_setequal(as.character(unique(filt@meta.data[[res_col]])), keep_ids)
      # 可回溯：每个保留亚群的标签与 CSV 映射一致
      for (id in keep_ids) {
        expect_identical(unique(as.character(
          filt$subcluster_label[filt@meta.data[[res_col]] == id])), labels[[id]])
      }
    }
  )
})

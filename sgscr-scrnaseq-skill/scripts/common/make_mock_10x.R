# make_mock_10x.R — 生成 10X 格式 mock 数据（真实基因名，支持 CellChat 冒烟）
# 用法: Rscript make_mock_10x.R --out <dir> [--samples 6] [--cells 120]
# 基因：CellChatDB.human 配体/受体等真实基因 + 常驻表达基因；
# 每个样本内按 10 个伪簇注入 marker 富集，保证 S4 聚类与 S7 通讯分析有信号。

args <- commandArgs(trailingOnly = TRUE)
geta <- function(f, d) { i <- grep(paste0("^--", f, "$"), args) + 1; if (length(i)) args[i] else d }
out_dir <- normalizePath(geta("out", stop("--out required")), mustWork = FALSE)
n_samples <- as.integer(geta("samples", "6"))
n_cells <- as.integer(geta("cells", "120"))
ids <- strsplit(geta("ids", paste(sprintf("S%d", seq_len(n_samples)), collapse=",")), ",")[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(CellChat))
db <- CellChatDB.human$interaction
db_genes <- unique(c(db$ligand.symbol, db$receptor.symbol))
set.seed(42)
tf_path <- geta("tfs", "")
tf_genes <- if (nzchar(tf_path) && file.exists(tf_path)) {
  head(readLines(tf_path), 60)   # 注入真实 TF，保证 S9 grn 的 gene∩TF 非空
} else character(0)
genes <- c(sample(db_genes, 250 - length(tf_genes)), tf_genes, sprintf("HK%i", 1:50))   # 300 基因
# 可选：用真实共享 motif 的基因替换 HK 占位基因，保证 S9 ctx pruning 有真实可检出的
# motif 富集信号（随机基因对任何 motif 的 NES 都过不了阈值，属 mock 缺陷而非管线缺陷）
mg_path <- geta("motif-genes", "")
mg_genes <- if (nzchar(mg_path) && file.exists(mg_path)) head(readLines(mg_path), 300) else character(0)
if (length(mg_genes) > 0) {
  genes <- setdiff(genes, sprintf("HK%i", 1:50))
  genes <- unique(c(genes, sample(mg_genes)))
}
# 可选：通过 --target-genes 注入 config 的目标基因（基因清单来自配置，不在本脚本内置），保证 S10 正向路径可冒烟
tg_path <- geta("target-genes", "")
if (nzchar(tg_path) && file.exists(tg_path)) {
  genes <- unique(c(genes, readLines(tg_path)))
}

n_clusters <- 10
markers <- matrix(sample(genes), ncol = n_clusters)        # 每伪簇 30 个富集基因
for (s in seq_len(n_samples)) {
  sid <- ids[s]
  d <- file.path(out_dir, sid); dir.create(d, showWarnings = FALSE)
  counts <- matrix(rpois(n_cells * length(genes), lambda = 2),
                   nrow = length(genes), dimnames = list(genes, NULL))
  cl <- sample(seq_len(n_clusters), n_cells, replace = TRUE)
  for (k in seq_len(n_clusters)) {
    idx <- which(cl == k)
    counts[markers[, k], idx] <- counts[markers[, k], idx, drop = FALSE] +
      matrix(rpois(length(markers[, k]) * length(idx), 6), nrow = length(markers[, k]))
  }
  colnames(counts) <- sprintf("%s_CELL_%d", sid, seq_len(n_cells))
  writeMM(as(counts, "dgCMatrix"), file.path(d, "matrix.mtx"))
  writeLines(colnames(counts), file.path(d, "barcodes.tsv"))
  writeLines(sprintf("ID%d\t%s\tGene Expression", seq_along(genes), genes), file.path(d, "features.tsv"))
  # 10X 标准是 .gz 后缀：生成 gzip 副本（triplet check 按 .gz 校验）
  for (f in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
    fin <- file.path(d, f); fout <- file.path(d, paste0(f, ".gz"))
    zin <- gzfile(fout, "wb"); writeLines(readLines(fin), zin); close(zin)
  }
  cat(sid, ":", n_cells, "cells,", nrow(counts), "genes\n")
}

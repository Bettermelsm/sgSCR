# triplet_check.R — 10X 三联文件完整性校验
# Validates: Requirements 3.4, 3.5 (P1)
#
# Feature: sgscr-scrnaseq-skill, Property 1: 对任意样本目录与任意"缺失文件子集"，
# 校验函数准确报告恰好那组缺失的文件名；子集为空（三文件齐备）时判定通过。
#
# 纯文件系统逻辑，不依赖 Seurat，可在测试中用临时目录 mock。

# 对单个样本目录校验三联文件；返回缺失文件名向量（齐备时 character(0)）
sgs_check_triplet <- function(sample_dir, triplet_files) {
  missing <- character(0)
  for (f in triplet_files) {
    if (!file.exists(file.path(sample_dir, f))) missing <- c(missing, f)
  }
  missing
}

# 批量校验：对 manifest 中每个样本（跳过 excluded）检查 <data_dir>/<sample_id>/ 下三联文件
# 返回 list(sample_id = 缺失文件向量)；全部齐备时各元素为 character(0)
sgs_check_all_samples <- function(data_dir, sample_manifest, triplet_files,
                                  excluded = list()) {
  excluded_ids <- vapply(excluded, function(e) e$sample_id, character(1))
  result <- list()
  for (s in sample_manifest) {
    if (s$sample_id %in% excluded_ids) next   # 剔除样本不触发错误（设计 Error Handling 节）
    result[[s$sample_id]] <- sgs_check_triplet(
      file.path(data_dir, s$sample_id), triplet_files)
  }
  result
}

# 汇总判定：任一样本有缺失即 FALSE，并返回人类可读的缺失清单
sgs_triplet_gate <- function(check_result) {
  bad <- Filter(function(x) length(x) > 0, check_result)
  if (length(bad) == 0) {
    list(pass = TRUE, detail = "all sample triplets complete")
  } else {
    lines <- vapply(names(bad), function(id)
      sprintf("%s: missing %s", id, paste(bad[[id]], collapse = ", ")), character(1))
    list(pass = FALSE, detail = paste(lines, collapse = "; "))
  }
}

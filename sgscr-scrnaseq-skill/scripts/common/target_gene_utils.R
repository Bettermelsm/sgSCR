# target_gene_utils.R — S10 目标基因工具（P12 交集过滤，纯函数）
# Validates: Requirements 16.2 (P12)

# P12: 返回目标基因清单与结果表基因集合的交集（保持 targets 原顺序、去重）；
# 交集为空时返回空向量——不报错、不返回全表。
sgs_target_gene_filter <- function(targets, result_genes) {
  targets <- as.character(targets)
  result_genes <- as.character(result_genes)
  targets <- targets[!is.na(targets) & nzchar(targets)]
  unique(targets[targets %in% result_genes])
}

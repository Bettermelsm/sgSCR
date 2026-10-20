# resolution_selection.R — ARI 辅助聚类分辨率选择（PART5 改造固化）
# Validates: Requirements 10.2 (P9)
#
# Feature: sgscr-scrnaseq-skill, Property 9: 对任意候选分辨率列表与相邻分辨率 ARI 序列，
# 规则确定性地返回一个推荐分辨率；无骤降时返回最高候选，有骤降时返回骤降发生前
# 的分辨率；同一输入重复运行结果相同（确定性）。

# ARI（adjusted Rand index）纯 R 实现，不引入 mclust 依赖
sgs_ari <- function(a, b) {
  if (length(a) != length(b) || length(a) == 0) {
    stop("sgs_ari: inputs must be equal-length non-empty vectors", call. = FALSE)
  }
  tab <- table(a, b)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  comb <- function(x) sum(choose(x[x >= 2], 2))
  sum_ij <- comb(as.vector(tab))
  sum_ik <- comb(rowSums(tab))
  sum_jk <- comb(colSums(tab))
  expected <- sum_ik * sum_jk / choose(n, 2)
  max_index <- (sum_ik + sum_jk) / 2
  (sum_ij - expected) / (max_index - expected)
}

# 分辨率选择规则（确定性）：
#   ari_seq[k] = ARI(clusters(res[k]), clusters(res[k+1]))，长度 = length(resolutions)-1
#   - 无骤降（所有 ari_seq >= drop_threshold）→ 返回最高候选
#   - 第一个 ari_seq[k] < drop_threshold → 返回 res[k]（骤降发生前的分辨率）
sgs_select_resolution <- function(resolutions, ari_seq, drop_threshold) {
  n <- length(resolutions)
  if (n < 1 || length(ari_seq) != n - 1 || anyNA(ari_seq) ||
      !is.numeric(drop_threshold)) {
    stop("sgs_select_resolution: invalid inputs (need >=1 resolution, ari_seq length n-1)",
         call. = FALSE)
  }
  if (n == 1) return(list(resolution = resolutions[1], basis = "single candidate",
                          first_drop_idx = NA_integer_))
  drops <- which(ari_seq < drop_threshold)
  if (length(drops) == 0) {
    list(resolution = resolutions[n], basis = "no ARI drop; highest candidate",
         first_drop_idx = NA_integer_)
  } else {
    k <- drops[1]
    list(resolution = resolutions[k], basis = sprintf("ARI drop at res[%d]->[%d]", k, k + 1),
         first_drop_idx = as.integer(k))
  }
}

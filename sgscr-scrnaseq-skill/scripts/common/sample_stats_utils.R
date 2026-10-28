# sample_stats_utils.R — s10b 样本级统计纯函数库（阶段12/13/14 计算核心）
# Validates: Requirements 23.1-23.9, 24.3, 25.1-25.2
# 设计依据 design.md「样本级统计与敏感性分析：s10b 设计」：
#   - 组均值只由样本均值构成（细胞数不进入权重，23.1）
#   - 完全分离 = min(一组样本均值) > max(另一组)（23.3）
#   - 分离余量 margin = 成立方向上的 min - max；未分离时 margin<=0 照样输出（23.7）
#   - margin_ratio = margin / delta，delta==0 记 NA（23.9）
#   - margin < stats.min_separation_margin → separation_robust = FALSE（23.8）
#   - 本底对照群：该基因按主群检出判定为不表达的细胞群合并（24.3）
#   - 单一样本主导告警：某(组,细胞群)最大样本占比 > 阈值（25.1/25.2）

# ---- 23.1：每样本均值 → 组均值（细胞数不进权重） ----
sgs_sample_means <- function(values, sample_vec) {
  stopifnot(length(values) == length(sample_vec))
  ok <- !is.na(values) & !is.na(sample_vec)
  if (!any(ok)) return(stats::setNames(numeric(0), character(0)))
  tapply(values[ok], sample_vec[ok], mean)
}

sgs_group_mean <- function(m_s, grp_of_sample, g) {
  i <- names(m_s)[grp_of_sample[names(m_s)] == g]
  if (!length(i)) return(NA_real_)
  mean(m_s[i])
}

# ---- 23.3/23.4/23.7/23.8/23.9：分离判据 + 余量（对一组每样本均值 m_s 计算） ----
# grp_of_sample: 命名向量 sample_id -> group；g0=对照（design.groups 首项），g1=case
sgs_separation_stats <- function(m_s, grp_of_sample, g0, g1,
                                 min_separation_margin) {
  # min_separation_margin 无默认值（P12-d 必修3）：稳健性判据必须显式给阈值，
  # 缺省 0 会让任何已分离指标恒判 robust
  i0 <- names(m_s)[grp_of_sample[names(m_s)] == g0]
  i1 <- names(m_s)[grp_of_sample[names(m_s)] == g1]
  v0 <- unname(m_s[i0]); v1 <- unname(m_s[i1])
  if (!length(v0) || !length(v1)) {
    return(list(mean_group0 = NA_real_, mean_group1 = NA_real_, delta = NA_real_,
                min_group0 = NA_real_, max_group0 = NA_real_,
                min_group1 = NA_real_, max_group1 = NA_real_,
                separation = NA, separation_margin = NA_real_,
                margin_ratio = NA_real_, separation_robust = FALSE,
                n_samples_used_group0 = length(v0), n_samples_used_group1 = length(v1)))
  }
  M0 <- mean(v0); M1 <- mean(v1); delta <- M1 - M0
  min0 <- min(v0); max0 <- max(v0); min1 <- min(v1); max1 <- max(v1)
  # 完全分离：minA > maxB 或 minB > maxA（A=case, B=control）
  sep_low  <- min1 > max0    # case 整体低于 control
  sep_high <- min0 > max1    # case 整体高于 control
  separation <- sep_low || sep_high
  # 未分离时两项皆 ≤0：取更接近 0 的（P12-d 必修2）——margin 语义是"距分离还差多少"
  margin <- if (sep_low) min1 - max0 else if (sep_high) min0 - max1 else max(min1 - max0, min0 - max1)
  # ratio 用 |delta|（P12-d 必修1）：case 低于 control 时 delta 为负，
  # margin/abs(delta) 保持正值，避免表中出现误导性负比例
  margin_ratio <- if (isTRUE(delta != 0)) margin / abs(delta) else NA_real_
  separation_robust <- isTRUE(separation) && (margin >= min_separation_margin)
  list(mean_group0 = M0, mean_group1 = M1, delta = delta,
       min_group0 = min0, max_group0 = max0, min_group1 = min1, max_group1 = max1,
       separation = separation, separation_margin = margin,
       margin_ratio = margin_ratio, separation_robust = separation_robust,
       n_samples_used_group0 = length(v0), n_samples_used_group1 = length(v1))
}

# ---- 23.2：LOO 数据集构造（返回剔除单个样本后的 sample_id 向量列表） ----
sgs_loo_samples <- function(sample_ids) {
  stats::setNames(lapply(sample_ids, function(s) setdiff(sample_ids, s)),
                  sprintf("LOO_%s", sample_ids))
}

# ---- 23.5：小样本设计下可达到的最小双侧 p 值 ----
# 精确秩（置换）检验：全部 choose(n0+n1, n0) 种标签排列中，最极端的 1 种（及其对称）
# 给出最小双侧 p = 2 / choose(n0+n1, n0)；任一组样本数为 0 时 NA
sgs_min_attainable_pvalue <- function(n0, n1) {
  if (is.na(n0) || is.na(n1) || n0 < 1 || n1 < 1) return(NA_real_)
  2 / choose(n0 + n1, n0)
}

# ---- 23.6：低细胞数样本过滤（明细标记 + 汇总排除） ----
# 返回 list(keep = 逻辑向量(按输入顺序), low_n_flag = ...)
sgs_low_n_flag <- function(n_cells_per_sample, min_cells_per_sample) {
  n_cells_per_sample < min_cells_per_sample
}

# ---- 25.1/25.2：单一样本主导判定 ----
# counts: 命名向量 sample_id -> 细胞数（某组某细胞群内）；返回 list(driven, dominant_sample, max_fraction)
sgs_single_sample_dominance <- function(counts, threshold = 0.70) {
  total <- sum(counts)
  if (is.na(total) || total <= 0) {
    return(list(driven = FALSE, dominant_sample = NA_character_, max_fraction = NA_real_))
  }
  frac <- counts / total
  top <- names(which.max(frac))
  list(driven = isTRUE(unname(frac[top]) > threshold), dominant_sample = top,
       max_fraction = unname(frac[top]))
}

# ---- 29.3：模块基因集合 ∩ 表达矩阵行名（缺失只标记不静默丢弃） ----
sgs_module_gene_usage <- function(gene_set, matrix_rownames) {
  gene_set <- unique(as.character(gene_set))
  used <- intersect(gene_set, as.character(matrix_rownames))
  missing <- setdiff(gene_set, as.character(matrix_rownames))
  list(used = used, missing = missing)
}

# ---- 24.3：本底对照群自动选取 ----
# detection_by_cluster: 命名向量 cluster -> 该基因检出比例；不表达 = 检出 <= max_detection
sgs_background_clusters <- function(detection_by_cluster, max_detection = 0.05,
                                    exclude = character(0)) {
  cand <- names(detection_by_cluster)[
    detection_by_cluster <= max_detection & !(names(detection_by_cluster) %in% exclude)]
  cand
}

# ---- s10b 输出表列契约（design.md「输出表列契约」，顺序即列顺序） ----
# ---- 23.5/23.6/25.2：s10b 输出表列契约（与实际表头一一对应，P12-d 必修4） ----
sgs_sample_level_columns <- function(sample_ids) {
  base <- c("metric", "cell_group", "dataset", "mean_group0", "mean_group1", "delta",
            "min_group0", "max_group0", "min_group1", "max_group1",
            "separation", "separation_margin", "margin_ratio", "separation_robust",
            "n_samples_used_group0", "n_samples_used_group1")
  c(base,
    sprintf("value_%s", sample_ids),      # 每样本逐值列（低细胞数样本为 NA）
    sprintf("low_n_%s", sample_ids),      # 每样本低细胞数标记列
    "single_sample_driven")               # 代表性告警传递列（25.2）
}

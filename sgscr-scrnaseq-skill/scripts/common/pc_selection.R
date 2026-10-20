# pc_selection.R — PART4 定量选 PC 数算法（固化自 SG_scRNA_seq_FLOW_V2.1.R）
# Validates: Requirements 9.1 (P8)
#
# Feature: sgscr-scrnaseq-skill, Property 8: 对任意合法方差贡献百分比序列（非负、和为100、
# 长度>=2），返回的 pc.use 为正整数且不超过序列长度；规则有解时 pc.use 处累计方差 >= 配置
# 阈值；全零方差差异等退化输入不返回 NA（回退到配置的 fallback 值）。
#
# 规则（与原 PART4 一致，阈值入 Config）：
#   co1：累计方差 >= cum_th 且 单PC贡献 < single_th 的最早位置
#   co2：拐点检测——相邻方差贡献差最大的位置 +1
#   pc.use = min(co1, co2)，夹取到 [1, n]；退化输入用 fallback_npc

sgs_select_pcs <- function(var_pct, cum_th, single_th, elbow_diff_threshold = 0,
                           fallback_npc = NA_integer_) {
  n <- length(var_pct)
  # 输入合法性：非负数值向量
  if (!is.numeric(var_pct) || n < 2 || anyNA(var_pct) || any(var_pct < 0) ||
      sum(var_pct) == 0) {   # 全零方差：无信息可依，按 P8 回退
    return(list(pc_use = fallback_npc, basis = "fallback:invalid_input",
                co1 = NA_integer_, co2 = NA_integer_))
  }

  cumvar <- cumsum(var_pct)

  # co1：累计方差达标 且 该 PC 自身贡献已低于阈值（继续加 PC 收益低于阈值线）
  idx <- which(cumvar >= cum_th & var_pct < single_th)
  co1 <- if (length(idx) > 0) as.integer(idx[1]) else NA_integer_

  # co2：拐点 = 方差贡献相邻差分（降幅）最大的位置 +1；全零方差时 which.max 返回 1 → co2=2，
  # 不会 NA（P8 退化输入要求）
  drops <- -diff(var_pct)          # 正值 = 降幅
  co2 <- as.integer(which.max(drops) + 1)

  candidates <- c(co1, co2)
  candidates <- candidates[!is.na(candidates)]
  if (length(candidates) == 0) {
    return(list(pc_use = fallback_npc, basis = "fallback:no_rule_solution",
                co1 = co1, co2 = co2))
  }

  pc_use <- min(candidates)
  pc_use <- max(1L, min(as.integer(pc_use), n))   # 夹取到合法区间
  basis <- if (!is.na(co1) && co1 <= co2) "cumvar+singlePC rule" else "elbow"
  list(pc_use = pc_use, basis = basis, co1 = co1, co2 = co2)
}

# 人工覆盖（阶段 11 Reviewer 指示）：cfg$pca$manual_pc_override 为正整数时跳过自动算法。
# 保留自动算法结果于 auto 字段，留"人工覆盖了什么、算法原本会给什么"的对比痕迹（P13 决议可追溯）。
sgs_apply_pc_override <- function(sel, override, npcs_compute) {
  if (is.null(override) || is.na(override)) return(sel)
  if (length(override) != 1 || !is.numeric(override) || is.na(override) ||
      override != floor(override)) {
    stop(sprintf("pca.manual_pc_override 非法: %s（须为 [2, npcs_compute=%d] 内的整数）",
                 paste(override, collapse = ","), npcs_compute), call. = FALSE)
  }
  manual <- as.integer(override)
  if (manual < 2 || manual > npcs_compute) {
    stop(sprintf("pca.manual_pc_override 非法: %s（须为 [2, npcs_compute=%d] 内的整数）",
                 paste(override, collapse = ","), npcs_compute), call. = FALSE)
  }
  list(pc_use = manual, basis = "manual_override", co1 = sel$co1, co2 = sel$co2,
       auto = list(pc_use = sel$pc_use, basis = sel$basis,
                   co1 = sel$co1, co2 = sel$co2, applied = FALSE))
}

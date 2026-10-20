# annotation_utils.R — 两层注释通用函数
# Validates: Requirements 5.3, 5.4 (P5), 11.3, 11.4 (P10)
#
# Property 5: 物种与参考集一致性校验 — 当且仅当两者一致时放行；不一致时终止
#   并报告「物种不匹配」及双方具体物种名。
# Property 10: 映射应用 — 应用后每细胞 major_cluster 恰等于其 cluster_id 的映射值；
#   major_cluster 唯一值集合恰等于映射表确认列唯一值集合；原始注释列不变。

# ---- P5: 物种一致性校验（放行返回 TRUE，否则 stop 含双方物种名） ----
sgs_check_species <- function(config_species, reference_species) {
  if (missing(config_species) || missing(reference_species) ||
      is.null(config_species) || is.null(reference_species) ||
      !nzchar(config_species) || !nzchar(reference_species)) {
    stop("物种不匹配校验失败: config 或 reference 物种为空", call. = FALSE)
  }
  if (tolower(trimws(config_species)) != tolower(trimws(reference_species))) {
    stop(sprintf("物种不匹配: config=%s, reference=%s",
                 config_species, reference_species), call. = FALSE)
  }
  TRUE
}

# ---- SingleR 标签多数票（确定性：并列时取字典序最小） ----
sgs_majority_label <- function(labels) {
  if (length(labels) == 0) return(NA_character_)
  tab <- table(labels)
  names(tab)[which.max(tab)][1]   # which.max 并列时取首个；table 按字典序排序保证确定
}

# ---- 确认映射表校验：列齐、非空、覆盖全部 cluster_id ----
# 返回以 cluster_id 为名、confirmed_major_cluster 为值的命名向量
sgs_validate_confirmed_mapping <- function(mapping_df, cluster_ids) {
  need <- c("cluster_id", "confirmed_major_cluster")
  if (is.null(mapping_df) || !all(need %in% names(mapping_df))) {
    stop(sprintf("映射表缺少必需列: %s",
                 paste(setdiff(need, names(mapping_df)), collapse = ", ")), call. = FALSE)
  }
  conf <- trimws(as.character(mapping_df$confirmed_major_cluster))
  if (all(is.na(conf) | conf == "")) {
    stop("映射表 confirmed_major_cluster 列全为空（人工确认未完成）", call. = FALSE)
  }
  ids <- as.character(mapping_df$cluster_id)
  missing <- setdiff(as.character(cluster_ids), ids[!is.na(conf) & conf != ""])
  if (length(missing) > 0) {
    stop(sprintf("映射表未确认的簇编号: %s",
                 paste(sort(missing), collapse = ", ")), call. = FALSE)
  }
  keep <- !is.na(conf) & conf != ""
  stats::setNames(conf[keep], ids[keep])
}

# ---- P10: 应用映射（cluster_ids = 每细胞簇标签向量，mapping = 命名向量） ----
sgs_apply_cluster_mapping <- function(cluster_ids, mapping) {
  if (is.null(names(mapping)) || anyDuplicated(names(mapping))) {
    stop("sgs_apply_cluster_mapping: mapping 必须是以 cluster_id 为名的无重复命名向量",
         call. = FALSE)
  }
  uncovered <- setdiff(unique(as.character(cluster_ids)), names(mapping))
  if (length(uncovered) > 0) {
    stop(sprintf("映射表未确认的簇编号: %s", paste(sort(uncovered), collapse = ", ")),
         call. = FALSE)
  }
  unname(mapping[as.character(cluster_ids)])
}

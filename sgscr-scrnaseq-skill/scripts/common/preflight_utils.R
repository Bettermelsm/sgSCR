# preflight_utils.R — 下游模块输入前置校验（P11）
# Validates: Requirements 13.1, 13.2, 14.1, 14.2, 15.1, 15.2
#
# Property 11: 对任意下游模块与其输入所携带的必需项集合（元数据列/外部依赖/文件），
# 齐备时放行；对任意「缺失必需项子集」，终止并准确报告恰好那组缺失项的名称。

sgs_preflight <- function(module, required, present) {
  required <- as.character(required)
  present <- as.character(present)
  missing_items <- setdiff(required, present)
  if (length(missing_items) > 0) {
    stop(sprintf("[%s] 前置校验失败，缺失必需项: %s",
                 module, paste(sort(missing_items), collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# 产物新鲜度校验（P11 增强，阶段 11 Reviewer 指定优先项）：
# 下游产物 mtime 不得早于任一上游产物；否则说明下游消费了陈旧上游（如 s5 换了
# 映射表重跑而 s6..s10 未跟跑），终止并点名过期文件，不以陈旧状态继续。
sgs_check_freshness <- function(downstream, upstream, tol_sec = 2) {
  downstream <- as.character(downstream)
  upstream <- as.character(upstream)
  missing <- c(downstream, upstream)[!file.exists(c(downstream, upstream))]
  if (length(missing) > 0) {
    stop(sprintf("新鲜度校验失败，文件不存在: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  dn <- file.mtime(downstream)
  up <- file.mtime(upstream)
  stale_upstream <- upstream[up > dn + tol_sec]
  if (length(stale_upstream) > 0) {
    stop(sprintf("新鲜度校验失败：下游产物早于上游（陈旧产物，请重跑下游）— 下游: %s; 过期上游: %s",
                 paste(basename(downstream), collapse = ","),
                 paste(basename(stale_upstream), collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

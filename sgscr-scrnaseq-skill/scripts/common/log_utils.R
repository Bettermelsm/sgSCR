# log_utils.R — 统一日志（Cross_Skill_Convention v1.0 第4节）
# 格式：[<ISO8601时间戳>] [<LEVEL>] [<step_id>] <message>
# 各脚本禁止手写 cat()/print() 当日志，统一用本文件三个函数。

.log <- function(level, step_id, message, con = stderr()) {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  lines <- strsplit(as.character(message), "\n", fixed = TRUE)[[1]]
  for (ln in lines) {
    writeLines(sprintf("[%s] [%s] [%s] %s", ts, level, step_id, ln), con = con)
  }
}

log_info  <- function(step_id, message, con = stderr()) .log("INFO",  step_id, message, con)
log_warn  <- function(step_id, message, con = stderr()) .log("WARN",  step_id, message, con)
log_error <- function(step_id, message, con = stderr()) .log("ERROR", step_id, message, con)

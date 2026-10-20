# manifest_utils.R — Manifest_Record 初始化与只追加写入
# Validates: Requirements 17.1, 17.2, 17.3 (P13, P14)
#
# 设计要点：
# - JSON 格式（程序写/程序读，jsonlite）
# - steps 数组只追加、不改写历史（P14）：每次追加都基于重新读盘的最新状态，
#   条目按完成顺序排列；禁止任何 update/delete 既有条目的 API
# - 字段名与 references/cross_skill_convention.md 登记表一致（P13）
# - run 级决议值（物种/参考集/选定PC数/分辨率/主群映射）由各 Step 经
#   sgs_manifest_set_resolved() 写入 manifest/resolutions，最终汇总时原样导出（P13）

suppressPackageStartupMessages(library(jsonlite))

CONVENTION_VERSION <- "cross_skill_convention v1.0"

# ---- 初始化（S0 开始时调用一次） ----
sgs_manifest_init <- function(manifest_path, run_id, random_seed, species,
                              singleR_reference, config_snapshot_path,
                              input_checksums = list(), excluded_samples = list()) {
  m <- list(
    run_id = run_id,
    spec_feature = "sgscr-scrnaseq-skill",
    created_at = .sgs_iso_now(),
    random_seed = random_seed,
    species = species,
    singleR_reference = singleR_reference,
    config_snapshot_path = config_snapshot_path,
    input_checksums = input_checksums,
    excluded_samples = excluded_samples,
    steps = list(),
    resolutions = list(),          # run 级决议值（P13 回写一致性）
    figures_panel_map = list(),
    convention_version = CONVENTION_VERSION
  )
  .sgs_manifest_write(manifest_path, m)
  invisible(m)
}

# ---- 追加一个 Step 条目（唯一允许的 steps 写操作，只追加） ----
sgs_manifest_append_step <- function(manifest_path, step_id, env_name,
                                     package_versions = list(),
                                     quality_gate,           # list(status=, detail=)
                                     human_confirmation = NULL,
                                     resolved_params = list(),
                                     outputs = list(),
                                     started_at, finished_at) {
  m <- .sgs_manifest_read(manifest_path)
  entry <- list(
    step_id = step_id,
    started_at = started_at,
    finished_at = finished_at,
    duration_sec = as.numeric(difftime(
      .sgs_parse_iso(finished_at), .sgs_parse_iso(started_at), units = "secs")),
    env_name = env_name,
    package_versions = package_versions,
    quality_gate = list(status = quality_gate$status, detail = quality_gate$detail),
    human_confirmation = human_confirmation,
    resolved_params = resolved_params,
    outputs = outputs
  )
  m$steps <- append(m$steps, list(entry))   # 只追加
  .sgs_manifest_write(manifest_path, m)
  invisible(entry)
}

# ---- run 级决议值回写（P13：写入后读出必须与实际使用值完全相等） ----
sgs_manifest_set_resolved <- function(manifest_path, key, value) {
  m <- .sgs_manifest_read(manifest_path)
  if (!is.list(m$resolutions)) m$resolutions <- list()   # 防 JSON 空数组回读为原子向量
  m$resolutions[[key]] <- value
  .sgs_manifest_write(manifest_path, m)
  invisible(value)
}

# 只读窄口：按键读取 resolutions 单项（不导出整体读，避免绕过只追加约束 P14）
sgs_manifest_get_resolved <- function(manifest_path, key) {
  m <- .sgs_manifest_read(manifest_path)
  if (!is.list(m$resolutions)) return(NULL)
  m$resolutions[[key]]
}

# 只读窄口：S11 汇总用——各 step 的 id/gate/产物清单（不暴露整体 manifest）
sgs_manifest_steps_overview <- function(manifest_path) {
  m <- .sgs_manifest_read(manifest_path)
  do.call(rbind, lapply(m$steps, function(s) {
    outs <- s$outputs
    outs <- if (is.null(outs)) character(0) else as.character(unlist(outs))
    data.frame(step_id = s$step_id,
               gate_status = s$quality_gate$status,
               outputs = paste(outs, collapse = ";"),
               stringsAsFactors = FALSE)
  }))
}

sgs_manifest_add_figure_panel <- function(manifest_path, fig_id, panel_id,
                                          file, science_question = NULL) {
  m <- .sgs_manifest_read(manifest_path)
  m$figures_panel_map <- append(m$figures_panel_map, list(list(
    fig_id = fig_id, panel_id = panel_id, file = file,
    science_question = science_question)))
  .sgs_manifest_write(manifest_path, m)
  invisible(TRUE)
}

# 公开 API 即唯一写入口：sgs_manifest_append_step / sgs_manifest_set_resolved /
# sgs_manifest_add_figure_panel。不导出整体读/写函数，防止绕过"只追加"约束（P14）。

# ---- 内部：读/写与时间 ----
.sgs_manifest_read <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("Manifest not found (init first): %s", manifest_path), call. = FALSE)
  }
  jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
}

.sgs_manifest_write <- function(manifest_path, m) {
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  # digits=NA：全精度写出，保证决议值（如选定 PC 数/阈值）JSON 往返后完全一致（P13）
  json <- jsonlite::toJSON(m, auto_unbox = TRUE, null = "null", pretty = TRUE, digits = NA)
  writeLines(json, manifest_path, useBytes = TRUE)
  invisible(manifest_path)
}

.sgs_iso_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
.sgs_parse_iso <- function(s) strptime(s, "%Y-%m-%dT%H:%M:%S%z")

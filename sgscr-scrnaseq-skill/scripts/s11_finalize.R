# s11_finalize.R — 出图与 Manifest 最终汇总（需求 17/18，主环境）
# Validates: Requirements 17.2-17.5, 18.1-18.3 (P16, P17)
#
# 输入：manifest.json（全部前序 Step 条目）+ figures/ 全部产物
# 职责：
#   1. 校验 s0..s10 全部 Step 均有 manifest 条目（缺任一 → FAIL）
#   2. 校验各 Step 声明的 figures 产物均存在（缺 → FAIL 并点名）
#   3. 生成 figures/figure_index.csv（file × step × fig_id × panel_id × science_question）
#   4. 逐图登记 figures_panel_map（P10 图附带目标基因科学问题）
#   5. 追加 s11_finalize 条目：汇总 gate 状态、总耗时、R 版本（P14 只追加）

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s11_finalize.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R", "save_plot.R",
            "qc_gate_utils.R", "preflight_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s11_finalize"
cfg <- sgs_load_config(config_path)
set.seed(cfg$run$random_seed)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

manifest_path <- cfg$paths$manifest_path
fig_dir <- cfg$paths$figures_dir

# ---- P11 前置校验 ----
problems <- c()
if (!file.exists(manifest_path)) problems <- c(problems, "manifest.json")
if (!dir.exists(fig_dir)) problems <- c(problems, "figures 目录")
if (length(problems) > 0) sgs_preflight(STEP, problems, character(0))

ov <- sgs_manifest_steps_overview(manifest_path)
ov$key <- sub("^(s[0-9]+[a-z]?).*", "\\1", ov$step_id)   # s10/s10b/s11 两位编号不能只取前2字符
need <- c(sprintf("s%d", 0:10), "s10b")   # 阶段12：s10b 接入依赖链 s10 -> s10b -> s11
missing_steps <- setdiff(need, ov$key)

# ---- 阶段14（需求 25.2）：样本代表性告警进图注 ----
repr_note <- NA_character_
repr_file <- file.path(dirname(manifest_path), "diagnostics", "sample_representation.csv")
if (file.exists(repr_file)) {
  repr <- read.csv(repr_file, check.names = FALSE)
  drv <- repr[repr$single_sample_driven %in% TRUE, , drop = FALSE]
  if (nrow(drv)) {
    repr_note <- paste(sprintf("cell_group %s (%s) 由单一样本 %s 主导 (%.0f%%)",
                               drv$cell_group, drv$group, drv$dominant_sample,
                               100 * drv$max_fraction), collapse = "; ")
    log_warn(STEP, paste("图注代表性告警:", repr_note))
  }
}

# ---- 各 Step 声明的 figure 产物存在性（按 step 逐一拆分，保留归属） ----
fig_rows <- do.call(rbind, lapply(seq_len(nrow(ov)), function(i) {
  outs <- unique(unlist(strsplit(ov$outputs[i], ";")))
  figs <- outs[nzchar(outs) & grepl(fig_dir, outs, fixed = TRUE)]
  figs <- figs[basename(figs) != "figure_index.csv"]   # 排除本 Step 自身产物，避免自引用
  if (!length(figs)) return(NULL)
  data.frame(step_id = ov$step_id[i], file = figs, stringsAsFactors = FALSE)
}))
fig_outs <- if (is.null(fig_rows)) character(0) else fig_rows$file
missing_figs <- fig_outs[!file.exists(fig_outs)]

# ---- figure_index.csv + figures_panel_map 登记 ----
parse_fp <- sgs_parse_figure_filename   # 已抽到 common/save_plot.R，测试可直接覆盖（Reviewer 隐患项）
sq <- cfg$target_genes$science_questions
sq_text <- if (length(sq)) paste(sprintf("%s: %s", names(sq), unlist(sq)), collapse = "; ") else NA_character_

idx <- do.call(rbind, lapply(seq_len(nrow(fig_rows)), function(i) {
  f <- fig_rows$file[i]
  fp <- parse_fp(f)
  exists <- file.exists(f)
  if (exists && !is.na(fp$fig_id)) {
    sgs_manifest_add_figure_panel(manifest_path, fig_id = fp$fig_id,
      panel_id = fp$panel_id, file = f,
      science_question = if (grepl("^P10", fp$fig_id)) sq_text else NULL,
      source_csv = if (grepl("\\.(pdf|tiff|png)$", f))
        sub("\\.[a-z]+$", ".source.csv", f) else NULL,
      font_family = if (!is.null(cfg$figure$font_family)) cfg$figure$font_family else "sans")
  }
  data.frame(file = basename(f), step_id = fig_rows$step_id[i],
             fig_id = fp$fig_id, panel_id = fp$panel_id, exists = exists,
             science_question = if (identical(fp$fig_id, "P10")) sq_text else NA_character_,
             sample_representation_note = repr_note,
             stringsAsFactors = FALSE)
}))
if (is.null(idx)) idx <- data.frame(file = character(0), step_id = character(0),
  fig_id = character(0), panel_id = character(0), exists = logical(0),
  science_question = character(0), sample_representation_note = character(0))
write.csv(idx, file.path(fig_dir, "figure_index.csv"), row.names = FALSE)

# ---- 最终汇总（追加 s11_finalize 条目，P14 只追加） ----
gate_tab <- paste(sprintf("%s=%s", ov$step_id, ov$gate_status), collapse = ", ")
auto_pass <- length(missing_steps) == 0 && length(missing_figs) == 0 &&
  nrow(fig_rows) > 0
gate <- sgs_qc_gate(
  auto_pass = auto_pass,
  detail = sprintf("steps %d/%d (缺: %s); figures %d 张已登记 (缺: %s); gate: %s",
                   length(intersect(need, ov$key)), length(need),
                   if (length(missing_steps)) paste(missing_steps, collapse = ",") else "无",
                   sum(idx$exists & !is.na(idx$fig_id)),
                   if (length(missing_figs)) paste(basename(missing_figs), collapse = ",") else "无",
                   gate_tab))
log_info(STEP, paste("Quality_Gate:", gate$status, "-", gate$detail))

sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(R = paste(R.version$major, R.version$minor, sep = ".")),
  quality_gate = gate,
  resolved_params = list(
    n_steps = nrow(ov), n_figures = nrow(idx),
    missing_steps = as.list(missing_steps), missing_figures = as.list(basename(missing_figs)),
    gate_status_by_step = as.list(stats::setNames(ov$gate_status, ov$step_id)),
    convention_version = CONVENTION_VERSION),
  outputs = as.list(file.path(fig_dir, "figure_index.csv")),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
if (!auto_pass) quit(status = 1)
log_info(STEP, "S11 complete")

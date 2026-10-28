# s0_env_check.R — S0 环境与数据盘点
# Validates: Requirements 3.3, 3.4, 3.5, 4.3
#
# Quality_Gate（自动）：ssh/frp 任一可达（本脚本在远程执行即视为可达）
#   AND 全部样本三联文件齐备 AND 三个隔离环境均存在
# 产物：{remote_work_dir}/manifest/s0_env_report.json + Manifest 初始化

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s0_env_check.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
common <- file.path(HERE, "common")
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R",
            "triplet_check.R", "env_resolver.R")) source(file.path(common, f))

STEP <- "s0_env_check"
cfg <- sgs_load_config(config_path)

started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# ---- 随机种子贯穿（需求 17.4） ----
set.seed(cfg$run$random_seed)

# ---- 1. 三联文件校验（P1） ----
triplet_res <- sgs_check_all_samples(
  cfg$remote$data_dir, cfg$samples$manifest, cfg$remote$triplet_files,
  excluded = cfg$samples$excluded %||% list())
triplet_gate <- sgs_triplet_gate(triplet_res)
log_info(STEP, paste("triplet check:", triplet_gate$detail))

# ---- 2. 环境存在性（P4 真实路径） ----
# main = R 主环境（系统 R，存在性= Rscript 可用）；三个隔离模块 = conda env
env_status <- list()
for (k in c("main", "cellchat", "trajectory", "scenic")) {
  env_name <- cfg$envs[[k]]
  exists <- if (k == "main") nzchar(Sys.which("Rscript")) else sgs_conda_env_exists(env_name)
  env_status[[k]] <- list(name = env_name, exists = exists)
  if (!exists) {
    log_warn(STEP, sprintf("env missing: %s (envs.%s)", env_name, k))
  }
}
envs_ok <- all(vapply(env_status, function(e) isTRUE(e$exists), logical(1)))

# ---- 3. R 主环境关键包版本盘点（需求 4.3） ----
main_pkgs <- c("Seurat", "SingleR", "clustree", "yaml", "jsonlite")
pkg_versions <- list()
for (p in main_pkgs) {
  pkg_versions[[p]] <- tryCatch(as.character(packageVersion(p)),
                                error = function(e) "NOT_INSTALLED")
  if (pkg_versions[[p]] == "NOT_INSTALLED") {
    log_warn(STEP, paste("main env package not installed:", p))
  }
}

# ---- 4. Quality_Gate 汇总 ----
status <- if (triplet_gate$pass && envs_ok) "PASS" else "FAIL"
detail <- sprintf("triplets=%s; isolated_envs=%s",
                  ifelse(triplet_gate$pass, "ok", triplet_gate$detail),
                  ifelse(envs_ok, "ok",
                         paste("missing:", paste(
                           vapply(Filter(function(e) !e$exists, env_status),
                                  function(e) e$name, character(1)), collapse = ","))))
log_info(STEP, paste("Quality_Gate:", status, "-", detail))

# ---- 5. 产物落盘 ----
finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
manifest_path <- cfg$paths$manifest_path
report <- list(
  step_id = STEP,
  quality_gate = list(status = status, detail = detail),
  env_status = env_status,
  main_env_packages = pkg_versions,
  triplet_check = triplet_res,
  config_snapshot = list(
    data_dir = cfg$remote$data_dir,
    samples = cfg$samples$manifest,
    excluded = cfg$samples$excluded,
    random_seed = cfg$run$random_seed
  )
)
report_path <- file.path(dirname(manifest_path), "s0_env_report.json")
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE, digits = NA), report_path)

# Manifest 初始化（含运行级字段；输入校验和由 S1 导入时补充）
sgs_manifest_init(
  manifest_path, run_id = cfg$run$run_id, random_seed = cfg$run$random_seed,
  species = cfg$species$organism, singleR_reference = cfg$species$singleR_reference,
  config_snapshot_path = normalizePath(config_path),
  excluded_samples = cfg$samples$excluded %||% list())
sgs_manifest_append_step(
  manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = pkg_versions,
  quality_gate = list(status = status, detail = detail),
  started_at = started_at, finished_at = finished_at,
  outputs = list(report_path))

if (status == "FAIL") quit(status = 1)
log_info(STEP, "S0 complete")

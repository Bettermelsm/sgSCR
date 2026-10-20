# env_resolver.R — Step → 隔离环境解析与存在性校验
# Validates: Requirements 4.4, 4.5 (P4)
#
# Feature: sgscr-scrnaseq-skill, Property 4: 对任意已注册模块标识，解析函数返回 Config 中
# 为该模块配置的环境名；对未注册模块或环境不存在/缺依赖的情况，终止并报告缺失的环境名，
# 不静默回退到主环境。
#
# env_exists 可注入（真实实现查 conda env list；测试注入 mock）。

# Step → Config envs 键 的注册表（design.md「C. 环境隔离架构」）
.sgs_STEP_ENV <- c(
  s0 = "main", s1 = "main", s2 = "main", s3 = "main", s4 = "main",
  s5 = "main", s6 = "main",
  s7 = "cellchat", s8 = "trajectory", s9 = "scenic",
  s10 = "main", s11 = "main"
)

# 解析 step_id（如 "s7"）→ conda env 名。未注册 step 报错终止。
sgs_resolve_env <- function(step_id, envs_cfg, step_table = .sgs_STEP_ENV) {
  key <- tolower(sub("_.*$", "", step_id))   # "s7_cellchat" → "s7"
  if (!key %in% names(step_table)) {
    stop(sprintf("Unregistered step: %s (known: %s)", step_id,
                 paste(names(step_table), collapse = ", ")), call. = FALSE)
  }
  env_key <- step_table[[key]]
  if (is.null(envs_cfg[[env_key]])) {
    stop(sprintf("Config envs.%s not set for step %s", env_key, step_id), call. = FALSE)
  }
  envs_cfg[[env_key]]
}

# 环境存在性校验。env_exists: function(env_name) -> logical，缺依赖检查由调用方叠加。
# 校验失败时停止并明确报告环境名——不静默回退主环境。
sgs_require_env <- function(env_name, env_exists) {
  if (!isTRUE(env_exists(env_name))) {
    stop(sprintf("Conda env not found: %s (create it before running this step)",
                 env_name), call. = FALSE)
  }
  invisible(TRUE)
}

# 真实实现：查询 conda env 列表（供 run/生产路径使用）
sgs_conda_env_exists <- function(env_name) {
  conda_bin <- Sys.which("conda")
  if (!nzchar(conda_bin)) {
    # 非交互 ssh 常见：conda 仅在交互 shell PATH 中，尝试常见安装位置
    for (cand in c("~/miniconda3/bin/conda", "~/anaconda3/bin/conda",
                   "/opt/miniconda3/bin/conda", "/opt/anaconda3/bin/conda")) {
      if (file.exists(cand)) { conda_bin <- path.expand(cand); break }
    }
  }
  if (!nzchar(conda_bin)) return(FALSE)   # 无 conda 视作隔离环境不存在，由上层报错
  out <- suppressWarnings(system2(conda_bin, c("env", "list"), stdout = TRUE, stderr = FALSE))
  length(out) > 0 && any(grepl(paste0("(^|/)", env_name, "\\s*$"), trimws(out)))
}

# R 主环境关键包可加载性（S0 环境盘点用；缺依赖即明确报告包名）
sgs_require_packages <- function(pkns) {
  missing <- pkns[!vapply(pkns, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Missing R packages: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

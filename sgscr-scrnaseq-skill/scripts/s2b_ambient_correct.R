# s2b_ambient_correct.R — 环境 RNA 可选校正（阶段13，需求 24.4，主环境）
# 默认 ambient.correction.enabled=false → 跳过（SKIPPED 记入 Manifest），表达矩阵不变。
# 开启时在 QC 之后独立执行，产物 rds/s2b_obj_corrected.rds，不覆盖 s2 产物与原始矩阵。
# 校正后端：decontX（celda 包）。开启而后端缺失 → 明确失败并报告缺失包名，
# 不得静默跳过、不得伪造校正。

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[grep("^--config$", args) + 1]
if (is.na(config_path)) stop("Usage: Rscript s2b_ambient_correct.R --config <config.yaml>", call. = FALSE)

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("log_utils.R", "config_loader.R", "manifest_utils.R",
            "qc_gate_utils.R", "resume_utils.R", "preflight_utils.R")) {
  source(file.path(HERE, "common", f))
}

STEP <- "s2b_ambient_correct"
cfg <- sgs_load_config(config_path)
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

enabled <- isTRUE(cfg$ambient$correction$enabled)
in_rds <- file.path(cfg$paths$rds_dir, "s2_obj_qc.rds")
out_rds <- file.path(cfg$paths$rds_dir, "s2b_obj_corrected.rds")

if (!enabled) {
  # 默认路径：SKIPPED + 表达矩阵不变（24.4；13.6 测试校验逐元素相等）
  gate <- sgs_qc_gate(auto_pass = TRUE, detail = "ambient.correction.enabled=false -> SKIPPED")
  sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
    package_versions = list(),
    quality_gate = gate,
    resolved_params = list(status = "SKIPPED", reason = "correction disabled by default"),
    outputs = list(),
    started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  log_info(STEP, "S2b SKIPPED (correction disabled)")
  quit(status = 0)
}

# ---- 开启路径 ----
pre <- tryCatch(sgs_resume_plan("s2b", cfg$paths$rds_dir), error = function(e) e)
if (inherits(pre, "error")) { log_error(STEP, conditionMessage(pre)); quit(status = 1) }
problems <- c()
if (!file.exists(in_rds)) problems <- c(problems, "s2_obj_qc.rds")
if (length(problems) > 0) sgs_preflight(STEP, problems, character(0))

if (!requireNamespace("celda", quietly = TRUE) ||
    !("decontX" %in% getNamespaceExports("celda"))) {
  # 需求24.4 派工硬约束：开启而后端缺失必须明确失败，报出缺失包名
  log_error(STEP, "ambient.correction.enabled=true 但校正后端缺失：需要 R 包 celda(提供 decontX)。请先安装：BiocManager::install('celda')")
  quit(status = 1)
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(celda))
set.seed(cfg$run$random_seed)

obj <- readRDS(in_rds)
Seurat::DefaultAssay(obj) <- "RNA"
counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
res <- decontX(counts)
decontaminated <- res$decontXcounts
# 写入独立 layer，原始 counts/data 不动；不覆盖 s2 产物
obj[["decontX"]] <- CreateAssayObject(counts = decontaminated)
saveRDS(obj, out_rds)

gate <- sgs_qc_gate(auto_pass = TRUE,
  detail = sprintf("decontX applied on %d cells x %d genes; decontamination %% median=%.2f",
                   ncol(obj), nrow(obj), median(res$decontX$contamination)))
sgs_manifest_append_step(cfg$paths$manifest_path, step_id = STEP, env_name = cfg$envs$main,
  package_versions = list(celda = as.character(packageVersion("celda")),
                          Seurat = as.character(packageVersion("Seurat"))),
  quality_gate = gate,
  resolved_params = list(backend = "decontX/celda", output = out_rds),
  outputs = list(out_rds),
  started_at = started_at, finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
log_info(STEP, "S2b complete")

# Property tests for sgs_manifest_import_steps（P11D-d 任务 2）
#   导入后目标步数核对通过；重复导入不重复追加；
#   产物 MD5 与已导入记录不一致必 FAIL；产物缺失必 FAIL。
# Run: Rscript tests/property/test_manifest_import.R

suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(jsonlite))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "manifest_utils.R"))

tmp <- file.path(tempdir(), sprintf("mfimp_%d", as.integer(Sys.time())))
dir.create(tmp, recursive = TRUE)

mk_artifact <- function(p, content) writeLines(content, p)

src_mp <- file.path(tmp, "src_manifest.json"); dst_mp <- file.path(tmp, "dst_manifest.json")
a1 <- file.path(tmp, "s5_obj_annotated.rds"); a2 <- file.path(tmp, "s2_obj_qc.rds")
mk_artifact(a1, paste(sample(letters, 50, replace = TRUE), collapse = ""))
mk_artifact(a2, paste(sample(letters, 50, replace = TRUE), collapse = ""))

sgs_manifest_init(src_mp, "src_run", 42L, "human", "ref", "cfg.yaml")
sgs_manifest_append_step(src_mp, step_id = "s1_import_integrate", env_name = "main",
  quality_gate = list(status = "PASS", detail = "ok"), outputs = list(a2),
  started_at = "2026-10-19T10:00:00+0800", finished_at = "2026-10-19T10:01:00+0800")
sgs_manifest_append_step(src_mp, step_id = "s5_annotation", env_name = "main",
  quality_gate = list(status = "PASS", detail = "ok"), outputs = list(a1),
  started_at = "2026-10-19T11:00:00+0800", finished_at = "2026-10-19T11:05:00+0800")
sgs_manifest_init(dst_mp, "dst_run", 42L, "human", "ref", "cfg.yaml")
sgs_manifest_append_step(dst_mp, step_id = "s7_cellchat", env_name = "cellchat",
  quality_gate = list(status = "PASS", detail = "ok"), outputs = list(),
  started_at = "2026-10-19T12:00:00+0800", finished_at = "2026-10-19T12:10:00+0800")

test_that("import: steps appended with provenance md5; s11-style step count passes; idempotent re-import", {
  r <- sgs_manifest_import_steps(dst_mp, src_mp, c("s1_import_integrate", "s5_annotation"))
  expect_setequal(r$imported, c("s1_import_integrate", "s5_annotation"))
  m <- jsonlite::fromJSON(dst_mp, simplifyVector = FALSE)
  ids <- vapply(m$steps, function(s) s$step_id, character(1))
  expect_setequal(ids, c("s7_cellchat", "s1_import_integrate", "s5_annotation"))
  # s11 步数核对语义：s0..s10 key 覆盖检查（此处模拟缺 s0/s3/s4 → 导入前应缺，导入 s1/s5 后仅缺 3 个）
  keys <- unique(sub("^(s[0-9]+).*", "\\1", ids))
  expect_length(setdiff(c("s0","s1","s2","s3","s4","s5","s6","s7","s8","s9","s10"), keys), 8)
  # provenance 落上
  s5 <- m$steps[[which(ids == "s5_annotation")]]
  expect_equal(s5$provenance$copied_from_run, "src_run")
  expect_equal(s5$provenance$src_step_started_at, "2026-10-19T11:00:00+0800")
  expect_equal(s5$provenance$artifact_md5[[a1]], digest::digest(file = a1, algo = "md5"))
  # 幂等：重复导入不重复追加，verified 返回
  r2 <- sgs_manifest_import_steps(dst_mp, src_mp, c("s1_import_integrate", "s5_annotation"))
  expect_length(r2$imported, 0)
  m2 <- jsonlite::fromJSON(dst_mp, simplifyVector = FALSE)
  expect_length(m2$steps, 3)
})

test_that("import: tampered artifact md5 mismatch -> FAIL", {
  writeLines("TAMPERED", a1)   # 产物被改动
  expect_error(sgs_manifest_import_steps(dst_mp, src_mp, "s5_annotation"),
               regexp = "MD5 与已导入记录不一致")
})

test_that("import: missing artifact -> FAIL; missing src step -> FAIL", {
  a3 <- file.path(tmp, "ghost.rds")
  sgs_manifest_append_step(src_mp, step_id = "s3_dimreduce", env_name = "main",
    quality_gate = list(status = "PASS", detail = "ok"), outputs = list(a3),
    started_at = "2026-10-19T10:30:00+0800", finished_at = "2026-10-19T10:31:00+0800")
  expect_error(sgs_manifest_import_steps(dst_mp, src_mp, "s3_dimreduce"),
               regexp = "产物缺失")
  expect_error(sgs_manifest_import_steps(dst_mp, src_mp, "s9_scenic"),
               regexp = "缺待导入 step")
})

cat("MANIFEST_IMPORT_TESTS_ALL_PASS\n")

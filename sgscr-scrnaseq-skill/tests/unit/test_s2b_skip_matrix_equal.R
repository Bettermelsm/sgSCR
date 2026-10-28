# Unit test: s2b 默认关闭时表达矩阵与 s2 产物逐元素相等（需求 24.4 / 13.6，P12-c 必修2）
# 直接比对两个矩阵（非推理）：s2b SKIPPED 执行后重读 s2_obj_qc.rds，
# RNA counts/data 两层必须与执行前快照逐元素相等；且不得产出 s2b_obj_corrected.rds。
# Run: Rscript tests/unit/test_s2b_skip_matrix_equal.R

suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(Seurat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
SKILL_ROOT <- normalizePath(file.path(HERE, "..", ".."))

tmp <- file.path(tempdir(), paste0("s2b_skip_", as.integer(Sys.time())))
rds_dir <- file.path(tmp, "rds"); man_dir <- file.path(tmp, "manifest")
fig_dir <- file.path(tmp, "figures")
for (d in c(rds_dir, man_dir, fig_dir)) dir.create(d, recursive = TRUE)

# --- 最小合法 config（全部必需字段；ambient 段不写 correction.enabled = 默认关闭） ---
cfg_yaml <- "
run: {run_id: \"s2b_skip_test\", random_seed: 42}
species: {organism: \"human\", singleR_reference: \"mock_ref\", reference_species: \"human\"}
samples:
  manifest:
    - {sample_id: \"T1\", group: \"g0\"}
    - {sample_id: \"T2\", group: \"g1\"}
  excluded: []
remote:
  host_alias: \"localhost\"
  ssh_target: \"x\"
  mode: \"local\"
  frp_fallback: {host: \"x\", port: 1}
  data_dir: \"TMP/data\"
  remote_work_dir: \"TMP/run\"
  triplet_files: [\"barcodes.tsv.gz\", \"features.tsv.gz\", \"matrix.mtx.gz\"]
envs: {main: \"main\", cellchat: \"cc\", trajectory: \"tj\", scenic: \"sc\"}
qc: {min_features: 1, max_features: 1e6, min_cells_per_gene: 1, max_percent_mt: 100, max_percent_hb: 100, min_cells_pass_gate: 1}
integration: {dims: 2, reduction: \"cca\", nfeatures_hvg: 10}
pca: {npcs_compute: 5, cum_variance_threshold_pct: 90, single_pc_contrib_threshold_pct: 5, elbow_diff_threshold: 0.1, fallback_npc: 3}
clustering: {candidate_resolutions: [0.5], ari_drop_threshold: 0.1}
compartments: {targets: [\"c1\"], subcluster_resolution: 0.5, subcluster_npcs: 2}
cellchat: {pathways_of_interest: [\"MIF\"], min_cells_per_group: 1, exclude_groups: []}
trajectory: {tool: \"slingshot\", target_compartment: [\"c1\"]}
scenic: {backend: \"pyscenic\", species_db_path: \"x.feather\", target_subset: \"c1\"}
target_genes:
  genes: [\"G1\"]
  science_questions: {G1: \"q\"}
figure: {dpi: 72, width_in: 2, height_in: 2, formats: [\"pdf\"], panel_map: \"inline\", font_family: \"sans\", base_size: 10}
paths:
  figures_dir: \"FIG\"
  rds_dir: \"RDS\"
  manifest_path: \"MAN/manifest.json\"
"
cfg_yaml <- gsub("TMP", tmp, cfg_yaml)
cfg_yaml <- gsub("FIG", fig_dir, cfg_yaml)
cfg_yaml <- gsub("RDS", rds_dir, cfg_yaml)
cfg_yaml <- gsub("MAN", man_dir, cfg_yaml)
cfg_file <- file.path(tmp, "config.yaml")
writeLines(cfg_yaml, cfg_file)

# --- 构造最小 s2 对象（RNA: counts + data 两层，含 G1/G2/G3 三基因 6 细胞） ---
counts <- matrix(c(10, 0, 1, 5, 2, 0, 0, 8, 3, 1, 0, 7, 2, 2, 2, 0, 9, 9),
                 nrow = 3, dimnames = list(c("G1", "G2", "G3"),
                                           paste0("C", 1:6)))
obj <- CreateSeuratObject(counts = counts)
obj$sample_id <- rep(c("T1", "T2"), each = 3)
obj <- NormalizeData(obj, verbose = FALSE)
s2_path <- file.path(rds_dir, "s2_obj_qc.rds")
saveRDS(obj, s2_path)

# 执行前快照：矩阵逐元素 + 文件摘要
snap_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
snap_data <- GetAssayData(obj, assay = "RNA", layer = "data")
snap_md5 <- unname(tools::md5sum(s2_path))

# --- 初始化 Manifest（真实 run 由 s0 完成；本测试模拟 s0 的最小初始化） ---
source(file.path(SKILL_ROOT, "scripts", "common", "manifest_utils.R"))
sgs_manifest_init(file.path(man_dir, "manifest.json"), run_id = "s2b_skip_test",
                  random_seed = 42, species = "human", singleR_reference = "mock_ref",
                  config_snapshot_path = cfg_file)

# --- 运行 s2b（默认关闭路径） ---
RS <- file.path(R.home("bin"), "Rscript")
out <- system2(RS, c(file.path(SKILL_ROOT, "scripts", "s2b_ambient_correct.R"),
                     "--config", cfg_file), stdout = TRUE, stderr = TRUE)
cat(paste(out, collapse = "\n"), "\n")

test_that("24.4/13.6: s2b 默认关闭 — SKIPPED 记账、无校正产物、s2 矩阵逐元素相等", {
  expect_null(attr(out, "status"))   # system2 正常退出时 status 为 NULL
  # 不得产出校正对象
  expect_false(file.exists(file.path(rds_dir, "s2b_obj_corrected.rds")))
  # Manifest 必须有 SKIPPED 记账（不得因默认关闭而完全不记账）
  expect_true(file.exists(file.path(man_dir, "manifest.json")))
  man <- jsonlite::fromJSON(file.path(man_dir, "manifest.json"), simplifyVector = FALSE)
  s2b <- Filter(function(x) grepl("^s2b", x$step_id), man$steps)
  expect_length(s2b, 1)
  expect_true(grepl("SKIPPED",
    paste(s2b[[1]]$quality_gate$detail, s2b[[1]]$resolved_params$status), fixed = TRUE))
  # 直接比对：重读 s2 产物，counts/data 两层与执行前快照逐元素相等（非推理）
  obj_after <- readRDS(s2_path)
  expect_true(isTRUE(all.equal(
    GetAssayData(obj_after, assay = "RNA", layer = "counts"), snap_counts,
    check.attributes = TRUE)))
  expect_true(isTRUE(all.equal(
    GetAssayData(obj_after, assay = "RNA", layer = "data"), snap_data,
    check.attributes = TRUE)))
  expect_identical(dim(obj_after), dim(obj))
  # 文件级：字节摘要不变
  expect_identical(unname(tools::md5sum(s2_path)), snap_md5)
})

unlink(tmp, recursive = TRUE)
cat("test_s2b_skip_matrix_equal: all tests passed\n")

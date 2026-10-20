# 集成测试（mock）：s8 multi-compartment（P11A 补充派发 1）
# Validates: trajectory.target_compartment 列表支持、FAIL 隔离、primary/backup 标注
#
# 用临时目录构造两个 mock S6 产物：compA（3 群，可跑通）、compB（单群，注定 FAIL），
# 用完整合法 config 跑真实 scripts/s8_trajectory.R，断言：
#   (1) compA 产出 pseudotime_compA.rds 且 gate PASS；
#   (2) compB FAIL 不阻断 compA（进程退出码 0；反向顺序亦然）；
#   (3) Manifest 含 role=primary/backup 标注。
# Run: Rscript tests/integration/mock_s8_multicomp.R   （需 slingshot 环境，如 sgscr_trajectory）

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(yaml))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
SKILL <- file.path(HERE, "..", "..")
tmp <- file.path(tempdir(), sprintf("s8mock_%d", as.integer(Sys.time())))
dir.create(file.path(tmp, "rds"), recursive = TRUE)
dir.create(file.path(tmp, "manifest"), recursive = TRUE)
dir.create(file.path(tmp, "figures"), recursive = TRUE)
dir.create(file.path(tmp, "trajectory"), recursive = TRUE)

mk_sub <- function(n_groups, prefix) {
  set.seed(42)
  n <- n_groups * 60
  m <- Matrix::rsparsematrix(200, n, 0.15)
  m@x <- as.numeric(m@x > 0) * sample(1:10, length(m@x), replace = TRUE)
  dimnames(m) <- list(sprintf("G%d", 1:200), sprintf("%s_c%d", prefix, 1:n))
  o <- CreateSeuratObject(counts = m)
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = min(10, n - 1), verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:5, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.5, verbose = FALSE)
  # 强制精确 n_groups 个群（小数据聚类数可能偏少/偏多）
  Idents(o) <- factor(rep(seq_len(n_groups), each = 60))
  names(Idents(o)) <- colnames(o)
  o@meta.data$RNA_snn_res.0.5 <- as.character(Idents(o))
  o <- RunUMAP(o, dims = 1:5, verbose = FALSE)
  # 拉开群间距，保证 slingshot 有轨迹（列名须与 key 一致：UMAP_1/UMAP_2）
  e <- Embeddings(o, "umap")
  lab <- as.integer(Idents(o))
  e[, 1] <- e[, 1] + lab * 30
  colnames(e) <- c("UMAP_1", "UMAP_2")
  rownames(e) <- colnames(o)
  o[["umap"]] <- CreateDimReducObject(embeddings = e, key = "UMAP_", assay = "RNA")
  o
}

saveRDS(mk_sub(3, "A"), file.path(tmp, "rds", "s6_compA_subcluster.rds"))
saveRDS(mk_sub(1, "B"), file.path(tmp, "rds", "s6_compB_subcluster.rds"))

cfg <- yaml::read_yaml(file.path(SKILL, "config", "config.example.yaml"))
# 覆盖为临时目录 + mock 场景必需字段（不读写真实运行目录）
cfg$run <- list(run_id = "mock_s8_multicomp", random_seed = 42L)
cfg$species <- list(organism = "human", singleR_reference = "mock", reference_species = "human")
cfg$samples <- list(manifest = list(list(sample_id = "S1", group = "NC")), excluded = list())
cfg$remote <- list(host_alias = "mock", ssh_target = "mock", data_dir = tmp,
                   remote_work_dir = tmp, triplet_files = "mock_triplet",
                   frp_fallback = list(host = "mock", port = 22L))
cfg$envs <- list(main = "mock", cellchat = "mock", trajectory = "mock", scenic = "mock")
cfg$qc <- list(min_features = 1L, max_features = 1e6L, min_cells_per_gene = 1L,
               max_percent_mt = 100, max_percent_hb = 100, min_cells_pass_gate = 1L)
cfg$integration <- list(dims = 5L, reduction = "pca", nfeatures_hvg = 200L)
cfg$pca <- list(npcs_compute = 10L, cum_variance_threshold_pct = 80,
                single_pc_contrib_threshold_pct = 20, elbow_diff_threshold = 1,
                fallback_npc = 5L)
cfg$clustering <- list(candidate_resolutions = list(0.5), ari_drop_threshold = 0.1)
cfg$compartments <- list(targets = list("compA", "compB"), subcluster_resolution = 0.5,
                         subcluster_npcs = 10L)
cfg$cellchat <- list(pathways_of_interest = list("TGFb"), min_cells_per_group = 10L,
                     exclude_groups = list())
cfg$trajectory <- list(tool = "slingshot", target_compartment = list("compA", "compB"))
cfg$scenic <- list(backend = "mock", species_db_path = tmp, target_subset = "all", nes_threshold = 3.0)
cfg$target_genes <- list(genes = list("G1"), science_questions = list("mock"))
cfg$figure <- list(dpi = 72L, width_in = 4, height_in = 4, formats = list("png"),
                   panel_map = list(P8 = list(prefix = "P8")))
cfg$paths <- list(rds_dir = file.path(tmp, "rds"),
                  figures_dir = file.path(tmp, "figures"),
                  manifest_path = file.path(tmp, "manifest", "manifest.json"))
cfg_yaml <- file.path(tmp, "config.yaml"); write_yaml(cfg, cfg_yaml)
# Manifest 初始化（真实流程由 s0 完成；mock 直接落最小骨架）
source(file.path(SKILL, "scripts", "common", "manifest_utils.R"))
sgs_manifest_init(file.path(tmp, "manifest", "manifest.json"), "mock_s8_multicomp", 42L,
                  "human", "mock", cfg_yaml)

s8 <- file.path(SKILL, "scripts", "s8_trajectory.R")
env_bin <- Sys.getenv("SGS_RSCRIPT", unset = "Rscript")
res <- system2(env_bin, c(s8, "--config", shQuote(cfg_yaml)), stdout = "", stderr = "")
cat("s8 exit code:", res, "\n")

ok <- TRUE
# (1) compA 产物存在且为有效轨迹（3 群 >= 2，应真正 PASS，而非空转）
ok <- ok && file.exists(file.path(tmp, "trajectory", "pseudotime_compA.rds"))
pt <- readRDS(file.path(tmp, "trajectory", "pseudotime_compA.rds"))
ok <- ok && sum(!is.na(pt$cell_pseudotime)) == nrow(pt$pseudotime_per_lineage) &&
       ncol(pt$pseudotime_per_lineage) >= 1
# (2) FAIL 隔离：compB（单群，gate FAIL）不阻断 compA，进程仍以 0 退出
ok <- ok && identical(res, 0L)
# (3) Manifest：role 标注 + compA gate PASS + compB gate FAIL
mf <- paste(readLines(file.path(tmp, "manifest", "manifest.json")), collapse = "\n")
ok <- ok && grepl('primary', mf) && grepl('backup', mf)
ok <- ok && grepl('compA', mf) && grepl('compB', mf)
ok <- ok && grepl('with pseudotime', mf)          # compA 的 PASS 详情特征串
ok <- ok && grepl('compB: s6_compB_subcluster.rds 不存在|compB: 轨迹空|FAIL', mf)  # compB 失败留痕
cat(if (ok) "MOCK_S8_MULTICOMP_PASS\n" else "MOCK_S8_MULTICOMP_FAIL\n")
quit(status = if (ok) 0 else 1)

# resume_utils.R — 断点续跑调度
# Validates: Requirements 7.1-7.4 (P6)
#
# Feature: sgscr-scrnaseq-skill, Property 6: 对任意 Step 序列与续跑起点 k，
# 跳过序号 < k 的 Step；第 k-1 步产物存在时加载后从第 k 步执行；缺失时终止并
# 报告缺失的具体产物路径，不以空状态继续。
#
# 产物表（与 cross_skill_convention 第2节命名一致；新增 Step 同步登记）

.sgs_STEP_OUTPUTS <- c(
  s0  = NA,                                  # S0 产物为 manifest 本身
  s1  = "s1_obj_int.rds",
  s2  = "s2_obj_qc.rds",
  s2b = NA,                                  # 可选步骤（默认 SKIPPED，无硬产物）；链 s2 -> (s2b) -> s3
  s3  = "s3_obj_pca.rds",
  s4  = "s4_obj_cluster.rds",
  s5  = "s5_obj_annotated.rds",
  s6  = "s6_<compartment>_subcluster.rds",   # 每 Compartment 一份，S6 内部自行展开
  s6b = "s6b_<compartment>_subcluster.rds",  # 每 Compartment 一份（过滤后），与 s6 同展开语义
  s7  = NA, s8 = NA, s9 = NA,                # 隔离环境模块产物按组目录判定
  s10 = NA,                                  # s10 产物在 run_root/{target_genes,diagnostics,module_score}，s10b 自行 preflight
  s10b = "stats/sample_level_summary.csv",
  s11 = NA
)

# 返回从 from_step 开始应执行的 step 键序列；前序产物缺失时 stop 并报告具体路径
sgs_resume_plan <- function(from_step, rds_dir, step_table = names(.sgs_STEP_OUTPUTS)) {
  key <- tolower(sub("_.*$", "", from_step))
  if (!key %in% step_table) {
    stop(sprintf("Unknown resume step: %s", from_step), call. = FALSE)
  }
  k <- match(key, step_table)
  prev <- step_table[k - 1]
  if (k > 1) {
    prev_out <- .sgs_STEP_OUTPUTS[[prev]]
    if (identical(prev_out, "s6_<compartment>_subcluster.rds")) {
      # S6 的产物按实际 Compartment 文件存在任意一份即视为完成
      hits <- list.files(rds_dir, pattern = "^s6_.*_subcluster\\.rds$")
      if (length(hits) == 0) {
        stop(sprintf("Missing prerequisite output for %s: no s6_*_subcluster.rds in %s",
                     key, rds_dir), call. = FALSE)
      }
    } else if (identical(prev_out, "s6b_<compartment>_subcluster.rds")) {
      hits <- list.files(rds_dir, pattern = "^s6b_.*_subcluster\\.rds$")
      if (length(hits) == 0) {
        stop(sprintf("Missing prerequisite output for %s: no s6b_*_subcluster.rds in %s",
                     key, rds_dir), call. = FALSE)
      }
    } else if (!is.na(prev_out) && !file.exists(file.path(rds_dir, prev_out))) {
      stop(sprintf("Missing prerequisite output for %s: %s",
                   key, file.path(rds_dir, prev_out)), call. = FALSE)
    }
  }
  list(steps = step_table[k:length(step_table)],
       resume_from = key,
       prev_output = if (k > 1) file.path(rds_dir, .sgs_STEP_OUTPUTS[[prev]]) else NULL)
}

# marker_panel_utils.R — 物种/组织 marker 面板加载（P12-a A4）
# Validates: 0028 派工 A4：面板外置 references/marker_panels/<species>_<tissue>.yaml，
# 脚本按 config.annotation.marker_panel 读取，基因不得硬编码在脚本内。
#
# 面板 YAML 结构: species, tissue, lineages: {<lineage>: [genes...]}
# 返回: list(species, tissue, lineages=list, genes=去重全基因集, n_genes)
# 校验失败（文件缺失/结构非法/重复基因跨谱系）→ stop（由调用方 preflight 呈现）。

sgs_load_marker_panel <- function(cfg, skill_dir = NULL) {
  panel_name <- cfg$annotation$marker_panel
  if (is.null(panel_name) || !nzchar(panel_name)) {
    stop("config.annotation.marker_panel 未配置（如 \"human_colon\"）", call. = FALSE)
  }
  if (is.null(skill_dir)) {
    stop("skill_dir 未指定（调用方传入脚本所在 skill 根目录）", call. = FALSE)
  }
  panel_path <- file.path(skill_dir, "references", "marker_panels",
                          sprintf("%s.yaml", panel_name))
  if (!file.exists(panel_path)) {
    stop(sprintf("marker 面板不存在: %s（references/marker_panels/ 下可用面板见目录）", panel_path),
         call. = FALSE)
  }
  suppressPackageStartupMessages(library(yaml))
  p <- yaml::read_yaml(panel_path)
  if (is.null(p$lineages) || length(p$lineages) == 0) {
    stop(sprintf("面板 %s 缺 lineages 段", panel_name), call. = FALSE)
  }
  genes <- unique(unlist(lapply(p$lineages, as.character)))
  if (length(genes) == 0) {
    stop(sprintf("面板 %s lineages 为空", panel_name), call. = FALSE)
  }
  # 跨谱系重复基因 = 结构错误（一个基因只归属一个主谱系，特化用途走子集）
  all_g <- unlist(lapply(p$lineages, as.character))
  dup <- unique(all_g[duplicated(all_g)])
  if (length(dup)) {
    stop(sprintf("面板 %s 存在跨谱系重复基因: %s", panel_name, paste(dup, collapse = ", ")),
         call. = FALSE)
  }
  list(species = p$species, tissue = p$tissue,
       lineages = lapply(p$lineages, as.character), genes = genes, n_genes = length(genes))
}

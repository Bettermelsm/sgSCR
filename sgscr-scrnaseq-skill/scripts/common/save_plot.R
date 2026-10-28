# save_plot.R — Figure_Exporter 统一出口（P16；阶段15：字体配置化 + 每图源数据强制）
# Validates: Requirements 18.1-18.7
# 命名约定（cross_skill_convention 第2节）：Fig<X>_Panel<Y>_<描述>.<格式>
# 阶段15（需求 18.5/18.6/18.7）：
#   - 导出前统一套用 figure.font_family 与 figure.base_size；脚本内不出现字体字面量
#   - 调用方必须传入源数据 df（可复现该图每个数据点的长表，含样本标识列与细胞数列），
#     同名 .source.csv 与图同出；缺源数据视为调用错误，不允许只出图
#   - 源数据路径与字体记入返回值（由调用方随 Manifest figures_panel_map 登记）

suppressPackageStartupMessages(library(ggplot2))

# 校验源数据契约：必须含样本标识列与细胞数列（需求 18.7）
sgs_check_source_df <- function(source_df) {
  if (is.null(source_df) || !is.data.frame(source_df) || !nrow(source_df)) {
    stop("sgs_save_plot: source_df 缺失或为空——每张导出图必须附带可复现源数据（需求18.6）", call. = FALSE)
  }
  need_sample <- intersect(c("sample_id", "sample", "orig.ident"), names(source_df))
  need_ncell  <- intersect(c("n_cells", "n_cell", "cell_count"), names(source_df))
  if (!length(need_sample)) {
    stop("sgs_save_plot: source_df 缺样本标识列（sample_id）", call. = FALSE)
  }
  if (!length(need_ncell)) {
    stop("sgs_save_plot: source_df 缺细胞数列（n_cells）", call. = FALSE)
  }
  invisible(TRUE)
}

sgs_save_plot <- function(cfg, fig_id, panel_id, desc, plt, source_df, step_id = "") {
  sgs_check_source_df(source_df)
  outdir <- cfg$paths$figures_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  base <- sprintf("%s_%s_%s", fig_id, panel_id, desc)
  # 字体配置（需求 18.5）：config 缺省时给设计默认，不允许脚本内字面量散落
  font_family <- if (!is.null(cfg$figure$font_family)) cfg$figure$font_family else "sans"
  base_size <- if (!is.null(cfg$figure$base_size)) cfg$figure$base_size else 10
  plt <- plt + ggplot2::theme_classic(base_size = base_size, base_family = font_family)
  written <- character(0)
  for (fmt in cfg$figure$formats) {
    fpath <- file.path(outdir, paste0(base, ".", fmt))
    ggplot2::ggsave(fpath, plot = plt,
                    width = cfg$figure$width_in, height = cfg$figure$height_in,
                    dpi = cfg$figure$dpi, units = "in")
    written <- c(written, fpath)
  }
  src_path <- file.path(outdir, paste0(base, ".source.csv"))
  write.csv(source_df, src_path, row.names = FALSE)
  list(files = written, source_csv = src_path, font_family = font_family,
       base_size = base_size)
}

sgs_load_plotlib <- function() invisible(NULL)   # 预留：集中 ggplot 主题设置

# ComplexHeatmap 对象专用（netVisual_heatmap 等返回 Heatmap 而非 ggplot）；
# 同样要求源数据（阶段15 无豁免——热图数据矩阵 melt 后即源数据）
sgs_save_heatmap <- function(cfg, fig_id, panel_id, desc, hm, source_df, step_id = "") {
  sgs_check_source_df(source_df)
  outdir <- cfg$paths$figures_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  base <- sprintf("%s_%s_%s", fig_id, panel_id, desc)
  written <- character(0)
  for (fmt in cfg$figure$formats) {
    fpath <- file.path(outdir, paste0(base, ".", fmt))
    if (fmt == "pdf") grDevices::pdf(fpath, width = cfg$figure$width_in,
                                     height = cfg$figure$height_in)
    else grDevices::tiff(fpath, width = cfg$figure$width_in, height = cfg$figure$height_in,
                         units = "in", res = cfg$figure$dpi, compression = "lzw")
    ComplexHeatmap::draw(hm)
    grDevices::dev.off()
    written <- c(written, fpath)
  }
  src_path <- file.path(outdir, paste0(base, ".source.csv"))
  write.csv(source_df, src_path, row.names = FALSE)
  list(files = written, source_csv = src_path,
       font_family = if (!is.null(cfg$figure$font_family)) cfg$figure$font_family else "sans",
       base_size = if (!is.null(cfg$figure$base_size)) cfg$figure$base_size else 10)
}

# figure 产物文件名解析（Fig<X>_Panel<Y>_<描述>.<格式> / P<X>_Panel<Y>_<描述>.<格式>）
# 无前缀的纯数字编号文件按 P 系列归类（约定两种合法前缀之外的风格不属于本 Skill 产物）
sgs_parse_figure_filename <- function(f) {
  b <- basename(f)
  m <- regmatches(b, regexec("^(Fig|P)?([0-9.]+)_Panel([A-Za-z0-9]+)_(.+)\\.[a-z]+$", b))[[1]]
  if (length(m) == 0) {
    return(list(fig_id = NA_character_, panel_id = NA_character_, desc = b, ok = FALSE))
  }
  list(fig_id = sprintf("%s%s", if (nzchar(m[2])) m[2] else "P", m[3]),
       panel_id = paste0("Panel", m[4]), desc = m[5], ok = TRUE)
}

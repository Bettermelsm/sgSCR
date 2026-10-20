# save_plot.R — Figure_Exporter 统一出口（P16 骨架版，阶段 2 起即用，阶段 10 完善）
# Validates: Requirements 18.1-18.3
# 命名约定（cross_skill_convention 第2节）：Fig<X>_Panel<Y>_<描述>.<格式>

suppressPackageStartupMessages(library(ggplot2))

sgs_save_plot <- function(cfg, fig_id, panel_id, desc, plt, step_id = "") {
  outdir <- cfg$paths$figures_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  base <- sprintf("%s_%s_%s", fig_id, panel_id, desc)
  written <- character(0)
  for (fmt in cfg$figure$formats) {
    fpath <- file.path(outdir, paste0(base, ".", fmt))
    ggplot2::ggsave(fpath, plot = plt,
                    width = cfg$figure$width_in, height = cfg$figure$height_in,
                    dpi = cfg$figure$dpi, units = "in")
    written <- c(written, fpath)
  }
  invisible(written)
}

sgs_load_plotlib <- function() invisible(NULL)   # 预留：集中 ggplot 主题设置

# ComplexHeatmap 对象专用（netVisual_heatmap 等返回 Heatmap 而非 ggplot）
sgs_save_heatmap <- function(cfg, fig_id, panel_id, desc, hm, step_id = "") {
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
  invisible(written)
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

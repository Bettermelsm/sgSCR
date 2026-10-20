# Property tests for save_plot.R (P16)
# Feature: sgscr-scrnaseq-skill
#   P16: 对任意 (dpi,宽,高,格式集合,fig_id,panel_id,desc)，导出文件数 == |formats|，
#        实际格式与扩展名对应，文件名包含 fig_id 与 panel_id 且匹配 Fig(X)_Panel(Y)_* 模式
# Run: Rscript tests/property/test_save_plot.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(ggplot2))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "save_plot.R"))

plt <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()

test_that("P16-main: file count == |formats|, extension matches format, name contains ids", {
  property <- forall(
    gen.bind(function(dpi)
      gen.bind(function(w)
        gen.bind(function(h)
          gen.bind(function(nfmt)
            gen.map(function(fid) list(
                dpi = dpi, w = w, h = h,
                formats = c("pdf", "tiff", "png")[seq_len(nfmt)],
                fig_id = fid, panel_id = "PanelA", desc = "mock_plot"),
              gen.element(c("P8", "Fig3", "P10"))),
            gen.element(1:3)),
          gen.element(4:9)),
        gen.element(4:12)),
      gen.element(4:12)),
    function(x) {
      outdir <- file.path(tempdir(), sprintf("p16_%d", as.integer(Sys.time()) + x$dpi))
      cfg <- list(paths = list(figures_dir = outdir),
                  figure = list(dpi = x$dpi, width_in = x$w, height_in = x$h,
                                formats = x$formats))
      written <- sgs_save_plot(cfg, x$fig_id, x$panel_id, x$desc, plt)
      expect_length(written, length(x$formats))                       # 数量 == |formats|
      for (i in seq_along(written)) {
        ext <- tools::file_ext(written[i])
        expect_true(ext %in% x$formats)                               # 扩展名属于指定格式
        expect_true(file.exists(written[i]))
        expect_true(grepl(x$fig_id, basename(written[i]), fixed = TRUE))   # 含 fig_id
        expect_true(grepl(x$panel_id, basename(written[i]), fixed = TRUE)) # 含 panel_id
        expect_match(basename(written[i]),
                     sprintf("^%s_%s_.+\\.%s$", x$fig_id, x$panel_id, ext))  # Fig(X)_Panel(Y)_*
      }
      unlink(outdir, recursive = TRUE)
    }
  )
})

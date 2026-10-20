# Property tests for Cross_Skill_Convention compliance (P17)
# Feature: sgscr-scrnaseq-skill
#   P17: 产物路径命名匹配约定正则；日志行匹配约定格式（时间戳/级别/step_id/消息）
# Run: Rscript tests/property/test_convention.R

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "log_utils.R"))

# 约定正则（cross_skill_convention：Fig<X>_Panel<Y>_<描述>.<格式>；[时间戳] [级别] [step] 消息）
# 注意：匹配需 perl=TRUE——本机 R 的 POSIX 引擎不支持字符类内的 \] 转义
PATH_RX <- "^(Fig|P)?[0-9.]+_Panel[A-Za-z0-9]+_.+\\.(pdf|tiff|png|csv)$"
LOG_RX  <- "^\\[[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9:]{8}\\+?[0-9]*\\] \\[(INFO|WARN|ERROR)\\] \\[[^\\]]+\\] .+$"

test_that("P17-path: constructed figure paths match convention regex", {
  property <- forall(
    gen.bind(function(fid) gen.bind(function(panel) gen.map(
      function(desc) list(fid = fid, panel = panel, desc = desc),
      gen.element(c("umap_major", "MockMajor_3_trajectory", "target_gene_dotplot",
                    "comm_heatmap"))),
      gen.element(c("PanelA", "PanelB", "PanelC"))),
      gen.element(c("P8", "Fig3", "P10", "P4.1"))),
    function(x) {
      for (ext in c("pdf", "tiff")) {
        p <- sprintf("%s_%s_%s.%s", x$fid, x$panel, x$desc, ext)
        expect_match(p, PATH_RX, perl = TRUE)
      }
    }
  )
})

test_that("P17-log: every log line matches convention format (3 levels)", {
  f <- tempfile(fileext = ".log")
  con <- file(f, open = "w")
  log_info("s1_x", "导入完成", con = con)
  log_warn("s10_y", "交集为空", con = con)
  log_error("s9_z", "数据库缺失", con = con)
  close(con)
  lines <- readLines(f)
  expect_length(lines, 3)
  for (ln in lines) expect_match(ln, LOG_RX, perl = TRUE)
  unlink(f)
})

test_that("P17-neg: malformed paths/lines must NOT match", {
  expect_false(grepl(PATH_RX, "heatmap_final.pdf", perl = TRUE))               # 无 Fig/Panel 标识
  expect_false(grepl(PATH_RX, "P8_PanelA.pdf", perl = TRUE))                   # 缺描述段
  expect_false(grepl(LOG_RX, "INFO something", perl = TRUE))                   # 缺时间戳/step 字段
  expect_false(grepl(LOG_RX, "[2026-10-04T10:00:00+0800] [info] [s1] msg", perl = TRUE))  # 级别必须大写
})

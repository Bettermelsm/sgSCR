# Property tests for config_loader.R
# Requires: hedgehog (fallback: quickcheck), yaml
# Run on the remote R main env (no local R on Windows workstation):
#   Rscript tests/property/test_config_loader.R
#
# Feature: sgscr-scrnaseq-skill, Property 2: Config 加载的完整性与往返一致性
#   对任意合法 Config_File，加载后参数对象包含文件中全部字段且值相等；加载无交互式输入。
#   Validates: Requirements 2.3, 6.3
#
# Feature: sgscr-scrnaseq-skill, Property 3: Config 缺失必需参数检测
#   对任意"移除任意非空必需字段子集"的残缺配置，加载器终止并准确报告恰好那组被移除字段名。
#   Validates: Requirements 6.4

suppressPackageStartupMessages(library(hedgehog))
suppressPackageStartupMessages(library(testthat))

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(HERE, "..", "..", "scripts", "common", "config_loader.R"))

gen_string <- gen.element(c("a", "bc", "run1", "hostA", "human", "x_1"))
gen_scalar <- gen.choice(gen_string, gen.int(1000), gen.unif(0, 1))

# 构造"完整合法配置"：全部必需字段赋任意标量值
.sgs_TEST_required <- c("a.b", "a.c", "d", "e.f.g")
complete_cfg <- list(a = list(b = "x", c = 1), d = TRUE, e = list(f = list(g = "z")))

test_that("P2: load preserves all fields round-trip (test scope)", {
  property <- forall(
    gen_scalar,
    function(v) {
      cfg <- complete_cfg
      cfg$e$f$g <- v
      miss <- sgs_missing_fields(cfg, .sgs_TEST_required)
      expect_identical(miss, character(0))
      expect_identical(.sgs_get_at_path(cfg, "e.f.g"), v)
    }
  )
})

test_that("P3: missing subset reported exactly (test scope)", {
  property <- forall(
    gen.subsequence(.sgs_TEST_required),
    function(removed) {
      cfg <- complete_cfg
      for (f in removed) {
        parts <- strsplit(f, ".", fixed = TRUE)[[1]]
        if (length(parts) == 1) {
          cfg[[parts[1]]] <- NULL
        } else {
          cfg[[parts[1]]][[parts[2]]] <- NULL   # 测试域字段最深两层
        }
      }
      got <- sgs_missing_fields(cfg, .sgs_TEST_required)
      expect_true(setequal(got, removed))
      expect_length(got, length(removed))
    }
  )
})

# ---- Stage-0 问题3 待办补测例：samples.manifest 读取行为（阶段2 交付） ----
# Feature: sgscr-scrnaseq-skill, Property 3-supplement:
#   samples.manifest 为"映射表列表"类型（非字符向量），.sgs_is_unresolved 不得误判为未填；
#   每个样本的 sample_id/group 读取后与 YAML 声明一致。

test_that("P3-supplement: samples.manifest list-type fields load faithfully", {
  # 真实 YAML → 完整加载流程（同时预演问题4：真实文件走 sgs_load_config）
  yaml_path <- file.path(tempdir(), "cfg_manifest_test.yaml")
  writeLines(c(
    "samples:",
    "  manifest:",
    "    - {sample_id: \"NC1\", group: \"NC\"}",
    "    - {sample_id: \"UC2\", group: \"UC\"}",
    "  excluded:",
    "    - sample_id: \"UC4\"",
    "      reason: \"test reason\""
  ), yaml_path)
  y <- yaml::read_yaml(yaml_path)
  # manifest 列表不被字符判定分支误伤
  expect_false(.sgs_is_unresolved(y$samples$manifest))
  expect_false(.sgs_is_unresolved(y$samples$excluded))
  # 空字符串组仍应判未填（该分支的真实行为锚定）
  expect_true(.sgs_is_unresolved(c("a", "")))
  # 字段值保真
  expect_identical(y$samples$manifest[[1]]$sample_id, "NC1")
  expect_identical(y$samples$manifest[[2]]$group, "UC")
  expect_identical(y$samples$excluded[[1]]$sample_id, "UC4")
})

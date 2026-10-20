# assay_utils.R — Seurat v5 多 layer 残留清理 + marker 覆盖率门检
# Validates: 阶段 11 Reviewer 派发的 s5 缺陷固化修复
#
# 缺陷背景（2026-10-05 实测）：S1 逐样本 ScaleData 产生 scale.data.<sample> 分样本 layer，
# JoinLayers(assay="RNA") 只合并 counts/data，残留 layer 使批量 FindAllMarkers 对部分簇
# 静默零输出且 Gate 照常 PASS（24,488 行/37 簇，实为 45,359 行/48 簇）。

# 重建 assay 彻底移除 scale.data（含分样本残留）；仅保留 counts/data（验证于本轮诊断，
# 不影响 FindAllMarkers/SingleR 等 data 层消费方）。非 Seurat 对象或 assay 缺失时原样返回。
# 0008 要求：输出必须是 Assay5（v3 Assay 在 SeuratObject 5.2.0 上 JoinLayers 无方法，
# 见 0007 BLOCK），且幂等——对已清理对象重复调用不报错、结果不变。
# 实现：先用 tryCatch 守卫 JoinLayers（v3 分层 Assay 上可用；已是 Assay5 时亦可用），
# 再统一以 CreateAssay5Object 重建，保证无论输入是 v3 还是 v5 分层结构，输出均为 Assay5。
sgs_drop_scale_layers <- function(obj, assay = "RNA") {
  if (!assay %in% names(obj@assays)) return(obj)
  rna <- obj@assays[[assay]]
  if (length(setdiff(Layers(rna), c("counts", "data"))) == 0 &&
      inherits(rna, "Assay5")) return(obj)   # 幂等出口：已是无残留的 Assay5
  obj <- tryCatch(JoinLayers(obj, assay = assay),
                  error = function(e) {
                    # P12-a D2：只吞掉 v3/单层组合"无方法"这一已知错误，其余错误原样上抛
                    msg <- conditionMessage(e)
                    known <- grepl("没有适用于|no applicable method", msg)
                    if (!known) stop(e)
                    obj
                  })
  rna <- obj@assays[[assay]]
  if (length(setdiff(Layers(rna), c("counts", "data"))) == 0 &&
      inherits(rna, "Assay5")) return(obj)          # 重建后二次校验（JoinLayers 生效情形）
  counts <- GetAssayData(rna, layer = "counts")
  data <- GetAssayData(rna, layer = "data")
  rna_new <- CreateAssay5Object(counts = counts, data = data)   # 0008：输出必须 Assay5
  Key(rna_new) <- Key(rna)   # 保留 key：丢失会令下游 JoinLayers 的 Key.Seurat 校验崩溃（0008 实测）
  obj[[assay]] <- rna_new    # 公共赋值接口，避免 @assays 直改破坏对象元数据
  invisible(obj)
}

# 落盘前断言（P11A）：RNA assay 的 layers 只允许 counts/data，任何残留即报错终止——
# 把"静默残留"挡在产物落盘之前，而非靠下游 FindAllMarkers 偶然暴露。
sgs_assert_layers_clean <- function(obj, assay = "RNA") {
  if (!assay %in% names(obj@assays)) return(invisible(NULL))
  lay <- Layers(obj@assays[[assay]])
  bad <- setdiff(lay, c("counts", "data"))
  if (length(bad) > 0) {
    stop(sprintf("%s assay 落盘前仍有残留 layers: %s（sgs_drop_scale_layers 未生效，禁止落盘）",
                 assay, paste(bad, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

# Quality_Gate 新检查项：marker 表实际覆盖的簇数必须等于对象实际簇数。
# 返回 list(ok, missing)；ok=FALSE 时调用方必须判 FAIL 并点名缺失簇号——把同类
# "静默不完整"从靠数据矛盾偶然发现，变成显式报错。
sgs_check_marker_coverage <- function(markers, cluster_ids) {
  cluster_ids <- as.character(cluster_ids)
  covered <- if (is.null(markers) || nrow(markers) == 0) character(0) else {
    as.character(unique(markers$cluster))
  }
  missing <- setdiff(cluster_ids, covered)
  list(ok = length(missing) == 0, missing = missing,
       n_covered = length(covered), n_expected = length(cluster_ids))
}

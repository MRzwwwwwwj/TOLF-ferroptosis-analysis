# ============================================================================
# 47_slc7a11_celltype_expr.R — 2026-11-02 回修（R4：Fig5C 与 Fig6D 细胞定位矛盾）
# 用途：从 sc_annotated.rds 提取 SLC7A11 等关键基因【按细胞类型】的平均表达
#       与表达细胞比例，作为重画 Figure 5C 机制图的硬数据依据。
#       （R4 指 SLC7A11 应在成纤维/内皮而非巨噬，需用真实数据定口径）
# ============================================================================
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

wd <- Sys.getenv("TOLF_WS", ".")  # 工作区根目录（可用环境变量 TOLF_WS 覆盖）
OUT <- file.path(wd, "results")

obj <- readRDS(file.path(wd, "scRNA_out", "sc_annotated.rds"))
cat("对象维度:", nrow(obj), "genes x", ncol(obj), "cells\n")
cat("celltype 分布:\n"); print(table(obj$celltype))

# 兼容 Seurat v4/v5 取 data 层
get_data <- function(o) {
  tryCatch(GetAssayData(o, assay = "RNA", layer = "data"),
           error = function(e) GetAssayData(o, assay = "RNA", slot = "data"))
}
expr <- get_data(obj)
ct  <- as.character(obj$celltype)
grp <- if ("group" %in% colnames(obj@meta.data)) as.character(obj$group) else rep("all", ncol(obj))

genes <- c("SLC7A11","SPP1","HMOX1","STEAP3","CDKN1A","FTH1","FTL","TFRC","GPX4","ACSL4")
genes <- intersect(genes, rownames(expr))
cat("可用基因:", paste(genes, collapse=", "), "\n")

res <- do.call(rbind, lapply(genes, function(g) {
  v <- as.numeric(expr[g, ])
  do.call(rbind, lapply(sort(unique(ct)), function(c) {
    sel <- ct == c
    data.frame(gene = g, celltype = c,
               mean_expr = round(mean(v[sel]), 4),
               pct_expr  = round(mean(v[sel] > 0) * 100, 1),
               n_cells   = sum(sel))
  }))
}))
write.csv(res, file.path(OUT, "SLC7A11_keygenes_by_celltype.csv"), row.names = FALSE)

cat("\n=== 各基因按细胞类型 mean_expr（透视）===\n")
print(tidyr::pivot_wider(res[, c("gene","celltype","mean_expr")],
                         names_from = celltype, values_from = mean_expr))
cat("\n=== SLC7A11 详情 ===\n")
print(res[res$gene == "SLC7A11", c("celltype","mean_expr","pct_expr","n_cells")])
cat("=== 47 done ===\n")

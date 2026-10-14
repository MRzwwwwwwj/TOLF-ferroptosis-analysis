#!/usr/bin/env Rscript
# ============================================================================
# sc_02_annotate.R — 细胞类型注释（标志基因打分法，无参考库依赖）
# 输入：scRNA_out/sc_merged_clustered.rds
# 输出：scRNA_out/sc_annotated.rds, 注释 UMAP, marker 点图, 细胞组成表
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})

OUT <- file.path(getwd(), "scRNA_out")
obj <- readRDS(file.path(OUT, "sc_merged_clustered.rds"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)   # Seurat v5 合并样本层，便于跨群取值

# ---------- 标志基因面板 ----------
markers <- list(
  Fibroblast   = c("COL1A1", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Endothelial  = c("PECAM1", "VWF", "KDR"),
  Pericyte_SMC = c("RGS5", "ACTA2", "MYH11"),
  Macrophage   = c("CD68", "CD14", "CSF1R", "AIF1"),
  SPP1_Macro   = c("SPP1", "CTSB", "TREM2"),            # SPP1+ 巨噬亚群信号
  T_cell       = c("CD3D", "CD3E", "CD8A"),
  NK_cell      = c("NKG7", "GNLY"),
  B_cell       = c("CD79A", "MS4A1"),
  Mast_cell    = c("TPSAB1", "KIT"),
  OsteoChondro = c("SOX9", "RUNX2", "ACAN", "COL2A1", "IBSP"),
  Cycling      = c("MKI67", "TOP2A")
)
all_mk <- unique(unlist(markers))
all_mk <- all_mk[all_mk %in% rownames(obj)]

# 点图：人工检查 + 留档
pdf(file.path(OUT, "Fig_sc_marker_dotplot.pdf"), width = 14, height = 7)
print(DotPlot(obj, features = all_mk, group.by = "seurat_clusters") +
        RotatedAxis() + ggtitle("Canonical markers by cluster"))
dev.off()

# ---------- 每个 cluster 打分 → 最高分的细胞类型作为注释 ----------
avg <- AverageExpression(obj, features = all_mk, group.by = "seurat_clusters",
                         assays = "RNA", slot = "data")$RNA
# 各 marker 面板内取均值表达
panel_score <- sapply(markers, function(g) {
  g <- intersect(g, rownames(avg))
  if (length(g) == 0) return(rep(NA, ncol(avg)))
  colMeans(avg[g, , drop = FALSE])
})
panel_score <- as.data.frame(panel_score)
panel_score$cluster <- rownames(panel_score)
write.csv(panel_score, file.path(OUT, "sc_cluster_marker_scores.csv"), row.names = FALSE)

# 注释规则：每个 cluster 取非 SPP1_Macro 面板中分数最高者（SPP1_Macro 只用于亚群细分）
main_panels <- setdiff(names(markers), "SPP1_Macro")
cl_anno <- apply(panel_score[, main_panels], 1, function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == 0) return("Unknown")
  main_panels[which.max(x)]
})
names(cl_anno) <- sub("^g", "", panel_score$cluster)   # AverageExpression 给列名加了 g 前缀，去掉以对齐 seurat_clusters
cat("=== cluster -> celltype ===\n"); print(cl_anno)
write.csv(data.frame(cluster = names(cl_anno), celltype = cl_anno),
          file.path(OUT, "sc_cluster_annotation.csv"), row.names = FALSE)

obj$celltype <- unname(cl_anno[as.character(obj$seurat_clusters)])

# ---------- SPP1+ 巨噬细胞亚群标记 ----------
# 在注释为 Macrophage 的细胞中，SPP1 表达 > 0 记为 SPP1+ Macro
spp1_expr <- FetchData(obj, vars = "SPP1")[, 1]
is_macro <- obj$celltype == "Macrophage"
obj$celltype_detail <- as.character(obj$celltype)
obj$celltype_detail[is_macro & spp1_expr > 0] <- "SPP1+_Macrophage"
obj$celltype_detail[is_macro & spp1_expr <= 0] <- "SPP1-_Macrophage"

# ---------- 细胞组成 ----------
comp <- obj@meta.data %>%
  count(sample_id, group, dataset, celltype_detail, name = "n") %>%
  group_by(sample_id) %>% mutate(pct = round(100 * n / sum(n), 2)) %>% ungroup()
write.csv(comp, file.path(OUT, "sc_celltype_composition.csv"), row.names = FALSE)
print(comp)

# ---------- 图 ----------
pdf(file.path(OUT, "Fig_sc_UMAP_celltype.pdf"), width = 13, height = 5.5)
p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("Cell types")
p2 <- DimPlot(obj, group.by = "celltype_detail", label = TRUE, repel = TRUE) +
  ggtitle("With SPP1+ macrophage subset")
print(p1 + p2)
dev.off()

saveRDS(obj, file.path(OUT, "sc_annotated.rds"))
cat("=== sc_02 done ===\n")

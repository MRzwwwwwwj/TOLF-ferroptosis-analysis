#!/usr/bin/env Rscript
# ============================================================================
# sc_03_validate.R — 核心验证：铁死亡评分定位 + SPP1+ 巨噬细胞 + 关键基因
# 输入：scRNA_out/sc_annotated.rds, pcd_genesets/gs_ferroptosis.txt
# 输出：Figure 4 候选图 + 评分/统计表
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out")
obj  <- readRDS(file.path(OUT, "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"

# ---------- 基因集 ----------
ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- ferro[ferro != "" & ferro %in% rownames(obj)]
cat("Ferroptosis genes covered:", length(ferro), "\n")
key_genes <- intersect(c("SLC7A11", "STEAP3", "HMOX1", "CDKN1A", "SPP1", "TGFB1",
                         "COL13A1", "CRTAC1"), rownames(obj))

# ---------- 模块评分 ----------
obj <- AddModuleScore(obj, features = list(ferro), name = "Ferroptosis", seed = 42)
colnames(obj@meta.data)[ncol(obj@meta.data)]
score_col <- grep("^Ferroptosis", colnames(obj@meta.data), value = TRUE)[1]
obj$Ferroptosis_score <- obj@meta.data[[score_col]]

# 每细胞评分明细
cell_scores <- obj@meta.data %>%
  select(sample_id, dataset, group, celltype, celltype_detail, Ferroptosis_score)
write.csv(cell_scores, file.path(OUT, "sc_ferroptosis_scores_percell.csv"), row.names = FALSE)

# ---------- 统计 1：细胞类型 × 评分（总体定位） ----------
ct_summary <- cell_scores %>%
  group_by(celltype) %>%
  summarise(n = n(), mean_score = mean(Ferroptosis_score),
            median_score = median(Ferroptosis_score), .groups = "drop") %>%
  arrange(desc(mean_score))
write.csv(ct_summary, file.path(OUT, "sc_ferroptosis_by_celltype.csv"), row.names = FALSE)
print(ct_summary)

# ---------- 统计 2：GSE294458 内 LSS vs LDH（每细胞类型 Wilcoxon） ----------
d294 <- cell_scores %>% filter(dataset == "GSE294458")
wilcox_res <- d294 %>%
  group_by(celltype) %>%
  summarise(n_LSS = sum(group == "LF_hypertrophy"),
            n_LDH = sum(group == "Control_LDH"),
            mean_LSS = mean(Ferroptosis_score[group == "LF_hypertrophy"]),
            mean_LDH = mean(Ferroptosis_score[group == "Control_LDH"]),
            delta = mean_LSS - mean_LDH,
            p_wilcox = tryCatch(
              wilcox.test(Ferroptosis_score[group == "LF_hypertrophy"],
                          Ferroptosis_score[group == "Control_LDH"])$p.value,
              error = function(e) NA_real_), .groups = "drop") %>%
  mutate(p_adj = p.adjust(p_wilcox, "BH"))
write.csv(wilcox_res, file.path(OUT, "sc_ferroptosis_LSS_vs_LDH_wilcox.csv"), row.names = FALSE)
cat("=== LSS vs LDH ferroptosis (GSE294458) ===\n"); print(wilcox_res)

# ---------- 统计 3：SPP1+ vs SPP1- 巨噬细胞铁死亡评分 ----------
macro <- cell_scores %>% filter(celltype == "Macrophage")
if (nrow(macro) > 50 && length(unique(macro$celltype_detail)) == 2) {
  wt <- wilcox.test(Ferroptosis_score ~ celltype_detail, data = macro)
  macro_tab <- macro %>%
    group_by(celltype_detail) %>%
    summarise(n = n(), mean_score = mean(Ferroptosis_score), .groups = "drop")
  macro_tab$p_wilcox <- wt$p.value
  write.csv(macro_tab, file.path(OUT, "sc_SPP1macro_ferroptosis.csv"), row.names = FALSE)
  cat("=== SPP1+ vs SPP1- macrophage ferroptosis ===\n"); print(macro_tab)
}

# ---------- 图：Figure 4 候选 ----------
# A. 铁死亡评分 UMAP
pdf(file.path(OUT, "Fig4A_sc_ferroptosis_UMAP.pdf"), width = 6.5, height = 5.5)
print(FeaturePlot(obj, features = "Ferroptosis_score", order = TRUE) +
        scale_color_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B")) +
        ggtitle("Ferroptosis module score"))
dev.off()

# B. 评分 by 细胞类型 小提琴
pdf(file.path(OUT, "Fig4B_sc_ferroptosis_by_celltype.pdf"), width = 9, height = 5)
print(VlnPlot(obj, features = "Ferroptosis_score", group.by = "celltype",
              pt.size = 0) + RotatedAxis() +
        ggtitle("Ferroptosis score by cell type"))
dev.off()

# C. LSS vs LDH 分组小提琴（GSE294458）
obj294 <- subset(obj, subset = dataset == "GSE294458")
pdf(file.path(OUT, "Fig4C_sc_ferroptosis_LSS_vs_LDH.pdf"), width = 10, height = 5)
print(VlnPlot(obj294, features = "Ferroptosis_score",
              group.by = "celltype", split.by = "group", pt.size = 0) +
        RotatedAxis() + ggtitle("Ferroptosis score: LSS vs LDH (GSE294458)"))
dev.off()

# D. 关键基因表达
pdf(file.path(OUT, "Fig4D_sc_keygenes_featureplot.pdf"), width = 13, height = 9)
print(FeaturePlot(obj, features = key_genes, ncol = 3, order = TRUE) &
        scale_color_gradientn(colors = c("#4393C3", "#F7F7F7", "#D73027")))
dev.off()

# E. SPP1+/- 巨噬细胞铁死亡对比
if (sum(obj$celltype == "Macrophage") > 50) {
  pdf(file.path(OUT, "Fig4E_sc_SPP1macro_ferroptosis.pdf"), width = 6, height = 5)
  print(VlnPlot(subset(obj, subset = celltype == "Macrophage"),
                features = "Ferroptosis_score", group.by = "celltype_detail",
                pt.size = 0) + ggtitle("Ferroptosis: SPP1+ vs SPP1- macrophages"))
  dev.off()
}

cat("=== sc_03 done ===\n")

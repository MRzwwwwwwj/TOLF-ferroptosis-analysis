#!/usr/bin/env Rscript
# ============================================================================
# sc_06d_monocle3_nopartition.R — 巨噬细胞拟时序【全细胞连通图重跑版】
# 背景：sc_06 原版 learn_graph(use_partition=TRUE 默认) 按分区分别学图，
#       导致 1,310/3,367 (38.9%) 细胞拟时序=Inf（不在 root 所在分区）。
#       全量 Spearman rho=0.312 主要靠 Inf 细胞并列最高秩撑起；
#       轨迹内 2,057 细胞 rho 仅 0.047。审稿人可复现该问题。
# 改动：learn_graph(use_partition = FALSE) —— 全细胞学单张连通图，
#       所有细胞获得有限拟时序，重新计算全部统计量。
# 输入：scRNA_out/sc_annotated.rds + pcd_genesets/gs_ferroptosis.txt（与 sc_06 相同）
# 输出：scRNA_out/monocle_v2/（新目录，不覆盖旧 monocle/）
# 运行：cd 到含 scRNA_out/ 与 pcd_genesets/ 的工作目录后：
#       Rscript sc_06d_monocle3_nopartition.R 2>&1 | tee monocle_v2.log
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(ggplot2)
})
set.seed(42)

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out", "monocle_v2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 取巨噬细胞亚群（与 sc_06 完全一致） ----------
obj <- readRDS(file.path(BASE, "scRNA_out", "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"
try(obj <- JoinLayers(obj), silent = TRUE)

macro <- subset(obj, celltype_detail %in% c("SPP1+_Macrophage", "SPP1-_Macrophage"))
macro$celltype_detail <- droplevels(factor(macro$celltype_detail))
cat("=== 巨噬细胞亚群构成 ===\n"); print(table(macro$celltype_detail))

# ---------- 铁死亡模块评分（与 sc_06 完全一致） ----------
ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- ferro[ferro != "" & ferro %in% rownames(macro)]
macro <- AddModuleScore(macro, features = list(ferro), name = "FerroScore", seed = 42)
score_col <- grep("^FerroScore", colnames(macro@meta.data), value = TRUE)[1]
macro$Ferroptosis_score <- macro@meta.data[[score_col]]

# ---------- 构建 monocle3 对象 ----------
counts <- tryCatch(GetAssayData(macro, layer = "counts"),
                   error = function(e) GetAssayData(macro, slot = "counts"))
cds <- new_cell_data_set(counts,
                         cell_metadata = macro@meta.data,
                         gene_metadata = data.frame(gene_short_name = rownames(macro),
                                                    row.names = rownames(macro)))
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds <- cluster_cells(cds)

# 分区数（诊断输出：原版按分区学图 → 跨分区细胞 Inf）
np <- tryCatch(length(unique(partitions(cds))), error = function(e) NA)
cat("=== cluster_cells 分区数:", np, "===\n")

for (cn in c("celltype_detail", "group", "sample_id", "Ferroptosis_score")) {
  colData(cds)[[cn]] <- macro@meta.data[colnames(cds), cn]
}

# ★ 核心改动：全细胞学单张连通图（忽略分区）
cds <- learn_graph(cds, use_partition = FALSE)

# ---------- 根细胞：铁死亡评分最低（与 sc_06 相同逻辑） ----------
root_cell <- rownames(colData(cds))[which.min(colData(cds)$Ferroptosis_score)]
cat("Root cell:", root_cell, "(lowest ferroptosis score)\n")
cds <- order_cells(cds, root_cells = root_cell)

# ---------- 图 ----------
g1 <- plot_cells(cds, color_cells_by = "celltype_detail", label_branch_points = TRUE,
                 label_leaves = FALSE, cell_size = 0.6) + ggtitle("Celltype")
g2 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = TRUE,
                 cell_size = 0.6) + ggtitle("Pseudotime")
g3 <- plot_cells(cds, color_cells_by = "Ferroptosis_score", label_branch_points = TRUE,
                 cell_size = 0.6) + ggtitle("Ferroptosis score")
ggsave(file.path(OUT, "Fig_mono_A_trajectory.png"), g1, width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "Fig_mono_B_pseudotime.png"), g2, width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "Fig_mono_C_ferroscore.png"), g3, width = 6.5, height = 5.5, dpi = 300)

# ---------- 拟时序 × 铁死亡评分相关性 ----------
pt <- as.numeric(pseudotime(cds))
fs <- colData(cds)$Ferroptosis_score
n_inf <- sum(!is.finite(pt))
cat(sprintf("\n=== Inf 拟时序细胞数: %d / %d (%.1f%%) ===\n",
            n_inf, length(pt), 100 * n_inf / length(pt)))
sp <- cor.test(pt, fs, method = "spearman")
cat(sprintf("Spearman(pseudotime, ferroptosis score): rho=%.4f, p=%.3g\n",
            sp$estimate, sp$p.value))

pd <- data.frame(pseudotime = pt, ferroptosis = fs,
                 celltype = colData(cds)$celltype_detail)
write.csv(pd, file.path(OUT, "pseudotime_ferroscore_percell.csv"), row.names = FALSE)

p4 <- ggplot(pd, aes(pseudotime, ferroptosis, color = celltype)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  labs(x = "Pseudotime", y = "Ferroptosis module score",
       subtitle = sprintf("Spearman rho=%.2f, P=%.1g", sp$estimate, sp$p.value)) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_mono_D_pt_vs_ferro.png"), p4, width = 7, height = 5, dpi = 300)

p5 <- ggplot(pd, aes(celltype, pseudotime, fill = celltype)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  labs(y = "Pseudotime", x = NULL) + theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_mono_E_pt_by_celltype.png"), p5, width = 5.5, height = 5, dpi = 300)

wt <- wilcox.test(pseudotime ~ celltype, data = pd)
cat(sprintf("Pseudotime SPP1+ vs SPP1- Wilcoxon P=%.3g\n", wt$p.value))

saveRDS(cds, file.path(OUT, "cds_macro.rds"))

# ---------- Moran's I 驱动基因 ----------
pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
gd <- as.data.frame(pr_test_res)
mi_col <- if ("moranI" %in% colnames(gd)) "moranI" else "morans_I"
q_col  <- if ("moranI_q" %in% colnames(gd)) "moranI_q" else "q_value"
gd$moranI_use <- suppressWarnings(as.numeric(gd[[mi_col]]))
gd <- gd[!is.na(gd$moranI_use), ]
gd$gene <- rownames(gd)
top <- head(gd[order(-gd$moranI_use), ], 100)
write.csv(top[, c("gene", "gene_short_name", mi_col, "p_value", q_col)],
          file.path(OUT, "genes_moransI_top100.csv"), row.names = FALSE)
cat("\n=== Moran's I Top5 ===\n")
print(head(top[, c("gene_short_name", mi_col)], 5))

key_g <- intersect(c("SPP1", "SLC7A11", "STEAP3", "HMOX1", "CDKN1A", "FTL", "TFRC"),
                   rownames(cds))
if (length(key_g) >= 3) {
  g6 <- plot_genes_in_pseudotime(cds[key_g, ], color_cells_by = "celltype_detail",
                                 min_expr = 0.5, ncol = 2)
  ggsave(file.path(OUT, "Fig_mono_F_keygenes.png"), g6, width = 9, height = 10, dpi = 300)
}

cat("\n=== sc_06d done, 输出目录:", OUT, "===\n")

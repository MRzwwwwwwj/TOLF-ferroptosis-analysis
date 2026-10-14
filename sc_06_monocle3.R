#!/usr/bin/env Rscript
# ============================================================================
# sc_06_monocle3.R — 巨噬细胞拟时序分析（第一队列 GSE294458/267819）
# 目标：SPP1- -> SPP1+ 巨噬细胞的铁死亡表型演化轨迹
# 输入：scRNA_out/sc_annotated.rds + pcd_genesets/gs_ferroptosis.txt
# 输出：scRNA_out/monocle/ 目录（轨迹图 + 拟时序×铁死亡评分 + 驱动基因表）
#
# ★ 依赖：monocle3 包（需手动安装，见 加分项分析操作指南.md）
# ★ 若 monocle3 装不上，改跑 sc_06b_slingshot.R（结果等价可写入论文）
# ★ 人工复核点见脚本内 [复核] 标记
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(ggplot2)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out", "monocle")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 取巨噬细胞亚群 ----------
obj <- readRDS(file.path(BASE, "scRNA_out", "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"
try(obj <- JoinLayers(obj), silent = TRUE)

macro <- subset(obj, celltype_detail %in% c("SPP1+_Macrophage", "SPP1-_Macrophage"))
macro$celltype_detail <- droplevels(factor(macro$celltype_detail))
cat("=== 巨噬细胞亚群构成 ===\n"); print(table(macro$celltype_detail))
# [复核] SPP1+ 约 343 / SPP1- 约 3,024（总计 3,367），做拟时序足够

# ---------- 铁死亡模块评分（亚群内重新计算） ----------
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

# 标签转移到 cds
for (cn in c("celltype_detail", "group", "sample_id", "Ferroptosis_score")) {
  colData(cds)[[cn]] <- macro@meta.data[colnames(cds), cn]
}

cds <- learn_graph(cds)

# ---------- 根细胞选择：铁死亡评分最低的细胞（假设 SPP1- 低铁死亡为起点） ----------
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
sp <- cor.test(pt, fs, method = "spearman")
cat(sprintf("\nSpearman(pseudotime, ferroptosis score): rho=%.3f, p=%.2g\n",
            sp$estimate, sp$p.value))
# [复核] 预期正相关（SPP1+ 高铁死亡位于轨迹末端）；若为负，把 root 改为评分最高细胞重跑

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

# 拟时序在两亚群的分布
p5 <- ggplot(pd, aes(celltype, pseudotime, fill = celltype)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  labs(y = "Pseudotime", x = NULL) + theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_mono_E_pt_by_celltype.png"), p5, width = 5.5, height = 5, dpi = 300)

wt <- wilcox.test(pseudotime ~ celltype, data = pd)
cat(sprintf("Pseudotime SPP1+ vs SPP1- Wilcoxon P=%.2g\n", wt$p.value))

# 保存 cds 对象（便于复查/改图，且避免后续步骤失败时重跑全程）
saveRDS(cds, file.path(OUT, "cds_macro.rds"))

# ---------- 沿轨迹的驱动基因（Moran's I） ----------
# 注意：graph_test 返回逐基因统计表（DataFrame），不能赋回 cds
# ★ 本机 monocle3 版本输出列为 status/p_value/morans_test_statistic/morans_I/
#   gene_short_name/q_value（新版叫 moranI），此处两种列名兼容
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

# 关键基因轨迹图（铁死亡 + SPP1 相关）
key_g <- intersect(c("SPP1", "SLC7A11", "STEAP3", "HMOX1", "CDKN1A", "FTL", "TFRC"),
                   rownames(cds))
if (length(key_g) >= 3) {
  g6 <- plot_genes_in_pseudotime(cds[key_g, ], color_cells_by = "celltype_detail",
                                 min_expr = 0.5, ncol = 2)
  ggsave(file.path(OUT, "Fig_mono_F_keygenes.png"), g6, width = 9, height = 10, dpi = 300)
}

cat("\n=== sc_06 Monocle3 done ===\n")
cat("输出目录:", OUT, "\n")

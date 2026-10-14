#!/usr/bin/env Rscript
# ============================================================================
# sc_06b_slingshot.R — 巨噬细胞拟时序分析（slingshot 备选方案）
# 用途：monocle3 在 Windows 上安装失败时的等价替代（Bioconductor 二进制包，易装）
# 结果口径与 sc_06_monocle3.R 一致，可直接替换写入论文（方法学写 slingshot）
# 输入：scRNA_out/sc_annotated.rds + pcd_genesets/gs_ferroptosis.txt
# 输出：scRNA_out/slingshot/ 目录
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment)
  library(ggplot2)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out", "slingshot")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(BASE, "scRNA_out", "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"
try(obj <- JoinLayers(obj), silent = TRUE)

macro <- subset(obj, celltype_detail %in% c("SPP1+_Macrophage", "SPP1-_Macrophage"))
cat("=== 巨噬细胞亚群构成 ===\n"); print(table(macro$celltype_detail))

# ---------- 铁死亡评分 ----------
ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- ferro[ferro != "" & ferro %in% rownames(macro)]
macro <- AddModuleScore(macro, features = list(ferro), name = "FerroScore", seed = 42)
score_col <- grep("^FerroScore", colnames(macro@meta.data), value = TRUE)[1]
macro$Ferroptosis_score <- macro@meta.data[[score_col]]

# ---------- 亚群内独立降维（保证轨迹在同一亚群结构内） ----------
macro <- RunPCA(macro, npcs = 30, verbose = FALSE)
macro <- RunUMAP(macro, dims = 1:20, verbose = FALSE)

logdat <- tryCatch(GetAssayData(macro, layer = "data"),
                   error = function(e) GetAssayData(macro, slot = "data"))
sce <- SingleCellExperiment(
  assays = list(logcounts = logdat),
  colData = macro@meta.data[colnames(macro), ],
  reducedDims = list(UMAP = Embeddings(macro, "umap"),
                     PCA  = Embeddings(macro, "pca")[, 1:20])
)

# ---------- slingshot 轨迹（起点锁定 SPP1- 巨噬细胞） ----------
sce <- slingshot(sce, clusterLabels = "celltype_detail", reducedDim = "UMAP",
                 start.clus = "SPP1-_Macrophage")

pt <- slingPseudotime(sce)[, 1]
fs <- colData(sce)$Ferroptosis_score
ok <- !is.na(pt)
sp <- suppressWarnings(cor.test(pt[ok], fs[ok], method = "spearman"))
cat(sprintf("\nSpearman(pseudotime, ferroptosis score): rho=%.3f, p=%.2g\n",
            sp$estimate, sp$p.value))

pd <- data.frame(UMAP1 = reducedDims(sce)$UMAP[, 1],
                 UMAP2 = reducedDims(sce)$UMAP[, 2],
                 pseudotime = pt, ferroptosis = fs,
                 celltype = colData(sce)$celltype_detail)
write.csv(pd, file.path(OUT, "pseudotime_ferroscore_percell.csv"), row.names = FALSE)

# ---------- 图 ----------
p1 <- ggplot(pd, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(size = 0.4, alpha = 0.6) + labs(color = NULL) + theme_bw(base_size = 12)
p2 <- ggplot(pd, aes(UMAP1, UMAP2, color = pseudotime)) +
  geom_point(size = 0.4, alpha = 0.6) + scale_color_viridis_c() + theme_bw(base_size = 12)
p3 <- ggplot(pd, aes(UMAP1, UMAP2, color = ferroptosis)) +
  geom_point(size = 0.4, alpha = 0.6) + scale_color_viridis_c() + theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_sl_A_celltype.png"),  p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUT, "Fig_sl_B_pseudotime.png"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUT, "Fig_sl_C_ferroscore.png"), p3, width = 6, height = 5, dpi = 300)

p4 <- ggplot(pd[ok, ], aes(pseudotime, ferroptosis, color = celltype)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_smooth(method = "loess", color = "black") +
  labs(x = "Pseudotime", y = "Ferroptosis module score",
       subtitle = sprintf("Spearman rho=%.2f, P=%.1g", sp$estimate, sp$p.value)) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_sl_D_pt_vs_ferro.png"), p4, width = 7, height = 5, dpi = 300)

p5 <- ggplot(pd[ok, ], aes(celltype, pseudotime, fill = celltype)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  labs(y = "Pseudotime", x = NULL) + theme_bw(base_size = 12)
ggsave(file.path(OUT, "Fig_sl_E_pt_by_celltype.png"), p5, width = 5.5, height = 5, dpi = 300)

wt <- wilcox.test(pseudotime ~ celltype, data = pd[ok, ])
cat(sprintf("Pseudotime SPP1+ vs SPP1- Wilcoxon P=%.2g\n", wt$p.value))

cat("\n=== sc_06b slingshot done ===\n")
cat("输出目录:", OUT, "\n")

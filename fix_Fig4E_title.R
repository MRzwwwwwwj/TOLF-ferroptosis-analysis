#!/usr/bin/env Rscript
# ============================================================================
# fix_Fig4E_title.R — 只重画 Fig4E（标题左缘 "F" 被裁）
# 根因：sc_03_validate.R 中 pdf(width=6) 放不下左对齐长标题
# 做法：width 6 -> 7.5，其余参数与 sc_03 完全一致（AddModuleScore seed=42）
# 输出：scRNA_out/Fig4E_sc_SPP1macro_ferroptosis.pdf（覆盖）
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out")
obj  <- readRDS(file.path(OUT, "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"

ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- ferro[ferro != "" & ferro %in% rownames(obj)]
obj <- AddModuleScore(obj, features = list(ferro), name = "Ferroptosis", seed = 42)
score_col <- grep("^Ferroptosis", colnames(obj@meta.data), value = TRUE)[1]
obj$Ferroptosis_score <- obj@meta.data[[score_col]]

pdf(file.path(OUT, "Fig4E_sc_SPP1macro_ferroptosis.pdf"), width = 7.5, height = 5)
print(VlnPlot(subset(obj, subset = celltype == "Macrophage"),
              features = "Ferroptosis_score", group.by = "celltype_detail",
              pt.size = 0) + ggtitle("Ferroptosis: SPP1+ vs SPP1- macrophages"))
dev.off()
cat("=== fix_Fig4E done ===\n")

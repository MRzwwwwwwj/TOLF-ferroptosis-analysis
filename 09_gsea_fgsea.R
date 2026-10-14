#!/usr/bin/env Rscript
# ============================================================================
# 09_gsea_fgsea.R — 加分项③：经典 GSEA（fgsea，基于全基因 t 统计量排序）
# 基因集：5 套 PCD + m6A 调控酶
# 输出：results/GSEA_fgsea_PCD.csv + 铁死亡/焦亡 enrichment 曲线图 + NES 点图
# ============================================================================
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(ggplot2)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")

# ---------- 排序统计量（limma t） ----------
deg <- fread(file.path(R, "DEG_limma_full.csv"), data.table = FALSE)
stats <- setNames(deg$t, deg$Gene)
stats <- sort(stats[!is.na(stats)], decreasing = TRUE)
cat("Ranked genes:", length(stats), "\n")

# ---------- 基因集 ----------
gs_files <- list.files(file.path(BASE, "pcd_genesets"), pattern = "\\.txt$", full.names = TRUE)
pathways <- lapply(gs_files, function(f) {
  g <- trimws(readLines(f)); g <- g[g != "" & g != "NA"]
  g
})
names(pathways) <- sub("^gs_|\\.txt$", "", basename(f <- gs_files))
# m6A 调控酶作为第 6 套
pathways[["m6A_regulators"]] <- c("METTL3","METTL14","METTL16","WTAP","VIRMA","ZC3H13",
                                  "RBM15","RBM15B","CBLL1","FTO","ALKBH5","ALKBH3",
                                  "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2",
                                  "IGF2BP1","IGF2BP2","IGF2BP3","HNRNPA2B1","LRPPRC","ELAVL1")
pathways <- lapply(pathways, function(g) intersect(g, names(stats)))
cat("Gene set sizes:\n"); print(sapply(pathways, length))

# ---------- fgsea ----------
set.seed(42)
fg <- fgsea(pathways = pathways, stats = stats, minSize = 5, maxSize = 1000)
fg <- as.data.frame(fg)[order(fg$pval), ]
fg$leadingEdge <- sapply(fg$leadingEdge, function(x) paste(x, collapse = ";"))
fwrite(fg, file.path(R, "GSEA_fgsea_PCD.csv"))
cat("=== fgsea results ===\n")
print(fg[, c("pathway", "size", "ES", "NES", "pval", "padj")], digits = 3, row.names = FALSE)

# ---------- enrichment 曲线（显著或铁死亡必画） ----------
plot_sets <- unique(c("ferroptosis", fg$pathway[fg$padj < 0.10]))
for (pw in plot_sets) {
  p <- plotEnrichment(pathways[[pw]], stats) +
    labs(title = paste0("GSEA: ", pw),
         subtitle = sprintf("NES=%.2f, padj=%.3g",
                            fg$NES[fg$pathway == pw], fg$padj[fg$pathway == pw])) +
    theme_bw(base_size = 12)
  ggsave(file.path(F, sprintf("Fig_GSEA_%s.png", pw)), p, width = 6.5, height = 5, dpi = 300)
}

# ---------- NES 汇总点图 ----------
fg$lab <- sprintf("%s (n=%d)", fg$pathway, fg$size)
p <- ggplot(fg, aes(x = NES, y = reorder(lab, NES),
                    color = padj < 0.05, size = -log10(pmax(padj, 1e-5)))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"),
                     name = "padj<0.05") +
  labs(x = "Normalized enrichment score (NES)", y = NULL,
       size = expression(-log[10](padj)),
       title = "GSEA of PCD gene sets (ranked by limma t, GSE113212)") +
  theme_bw(base_size = 12)
ggsave(file.path(F, "Fig_GSEA_NES_dotplot.png"), p, width = 8, height = 4.5, dpi = 300)

cat("=== 09 done ===\n")

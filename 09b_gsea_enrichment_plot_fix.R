#!/usr/bin/env Rscript
# ============================================================================
# 09b_gsea_enrichment_plot_fix.R — 手工重绘铁死亡 GSEA 富集曲线
# 背景：fgsea::plotEnrichment 在本机 fgsea 版本下产出空图（tops/bottoms 警告），
#       此处按标准 GSEA running-sum 定义手绘，ES/NES 与 GSEA_fgsea_PCD.csv 对齐。
# 输入：results/DEG_limma_full.csv, pcd_genesets/gs_ferroptosis.txt,
#       results/GSEA_fgsea_PCD.csv（取 ES/NES/pval/padj 用于标注）
# 输出：figures/Fig_GSEA_ferroptosis.png（覆盖空图）
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- getwd(); R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")

deg <- fread(file.path(R, "DEG_limma_full.csv"), data.table = FALSE)
stats <- deg$t; names(stats) <- deg$Gene
stats <- sort(stats[!is.na(stats)], decreasing = TRUE)

gs <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
gs <- intersect(gs[gs != "" & gs != "NA"], names(stats))

# fgsea 官方结果（校验 + 标注）
fg <- fread(file.path(R, "GSEA_fgsea_PCD.csv"), data.table = FALSE)
row <- fg[fg$pathway == "ferroptosis.txt", ]

# ---------- 标准 GSEA running enrichment score（gseaParam = 1，与 fgsea 一致） ----------
hits <- names(stats) %in% gs
w <- abs(stats)                # 权重 |t|^1
N <- length(stats); Nh <- sum(hits)
Phit <- ifelse(hits, w / sum(w[hits]), 0)
Pmis <- ifelse(hits, 0, 1 / (N - Nh))
es <- cumsum(Phit - Pmis)      # running ES
peak <- which.max(abs(es))
cat(sprintf("手绘峰值 ES = %.4f @ rank %d ; fgsea 官方 ES = %.4f\n",
            es[peak], peak, row$ES))

df <- data.frame(rank = seq_len(N), es = es)
tick_df <- data.frame(rank = which(hits))

p <- ggplot(df, aes(x = rank, y = es)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_line(color = "#B2182B", linewidth = 0.9) +
  geom_segment(data = tick_df,
               aes(x = rank, xend = rank, y = min(df$es), yend = min(df$es) + 0.04 * diff(range(df$es))),
               color = "grey20", linewidth = 0.25, inherit.aes = FALSE) +
  labs(title = "GSEA: ferroptosis",
       subtitle = sprintf("NES = %.2f, P = %.3g, adj. P = %.3g (fgsea; leading edge at rank %d)",
                          row$NES, row$pval, row$padj, peak),
       x = "rank in ordered gene list (limma t, decreasing)", y = "running enrichment score") +
  theme_bw(base_size = 12)

ggsave(file.path(F, "Fig_GSEA_ferroptosis.png"), p, width = 6.5, height = 5, dpi = 300)
cat("=== 09b done: Fig_GSEA_ferroptosis.png 已重绘 ===\n")

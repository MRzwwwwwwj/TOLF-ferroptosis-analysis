# ============================================================
# TOLF Step 02b: Fig1A 火山图重标（2026-10-13 复核遗留项）
# 不重跑 limma，直接读 results/DEG_limma_full.csv
# 按 Figure 1 图注标注: TGFB1, BMPR1B, COL13A1, COL11A1
# 输出覆盖 figures/Fig_volcano_GSE113212.png/pdf，并重拼 Fig1
# ============================================================
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(magick)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
deg <- read.csv(file.path(wd, "results", "DEG_limma_full.csv"), stringsAsFactors = FALSE)

deg$sig <- "NS"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC > 1]  <- "Up (adj.P<0.05)"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC < -1] <- "Down (adj.P<0.05)"

lab_genes <- c("TGFB1", "BMPR1B", "COL13A1", "COL11A1")
lab_df <- subset(deg, Gene %in% lab_genes)
cat("标注基因核对:\n")
print(lab_df[, c("Gene", "logFC", "P.Value", "adj.P.Val")])
stopifnot(nrow(lab_df) == 4)

set.seed(42)
p <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Up (adj.P<0.05)" = "#D85A30",
                                "Down (adj.P<0.05)" = "#185FA5", "NS" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3) +
  geom_text_repel(data = lab_df, aes(label = Gene),
                  color = "black", size = 3.6, fontface = "italic",
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0, segment.size = 0.3,
                  max.overlaps = Inf, seed = 42) +
  labs(title = "GSE113212: elderly hypertrophied vs young control",
       x = "log2 fold change", y = "-log10 adjusted P") +
  theme_bw(base_size = 12)

ggsave(file.path(wd, "figures", "Fig_volcano_GSE113212.pdf"), p, width = 6, height = 5)
ggsave(file.path(wd, "figures", "Fig_volcano_GSE113212.png"), p, width = 6, height = 5, dpi = 300)
cat("volcano relabeled & saved.\n")

# ---------- 仅重拼 Fig1（不动其它拼版，14b/14c 的修复保持原样） ----------
FIG <- file.path(wd, "figures")
OUT <- file.path(FIG, "panels")
load_h <- function(path, h = 1800) image_scale(image_read(path), paste0("x", h))
tag <- function(img, label) image_annotate(img, label, size = 90, weight = 700,
                                           gravity = "northwest", location = "+30+20", color = "black")
imgs <- list(load_h(file.path(FIG, "Fig_volcano_GSE113212.png")),
             load_h(file.path(FIG, "Fig_PCD_scores_boxplot.png")),
             load_h(file.path(FIG, "Fig_immune_boxplot.png")),
             load_h(file.path(FIG, "Fig_PCD_immune_corheatmap.png")))
tagged <- Map(tag, imgs, c("A", "B", "C", "D"))
row1 <- image_append(c(tagged[[1]], tagged[[2]]))
row2 <- image_append(c(tagged[[3]], tagged[[4]]))
final <- image_border(image_append(c(row1, row2), stack = TRUE), "white", "40x40")
image_write(final, file.path(OUT, "Fig1_discovery_overview.png"), format = "png", density = 300)
cat("Fig1 reassembled:", image_info(final)$width, "x", image_info(final)$height, "\n")

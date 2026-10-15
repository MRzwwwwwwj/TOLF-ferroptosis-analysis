# ============================================================
# TOLF Step 02: limma 差异分析 (GSE113212)
# 对比: elderly_hypertrophied vs young_control
# 注意: n=4v4, 功效有限 -> 主阈值 adj.P<0.05 & |log2FC|>1;
#       若显著基因过少, 降级报告 P<0.01 & |log2FC|>1 (探索性, 如实标注)
# ============================================================

.libPaths("D:/ARLARS/R-4.5.3/library")
suppressPackageStartupMessages({library(limma); library(ggplot2)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
dir.create(file.path(wd, "figures"), showWarnings = FALSE)

expr <- read.csv(file.path(wd, "results", "expr_gene_GSE113212.csv"),
                 row.names = 1, check.names = FALSE)
group <- read.csv(file.path(wd, "results", "group_GSE113212.csv"),
                  stringsAsFactors = FALSE)
expr <- expr[, group$sample]
g <- factor(group$group, levels = c("young_control", "elderly_hypertrophied"))

design <- model.matrix(~0 + g)
colnames(design) <- levels(g)
fit <- lmFit(expr, design)
cm <- makeContrasts(DEG = elderly_hypertrophied - young_control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
deg <- topTable(fit2, number = Inf, sort.by = "P")
deg$Gene <- rownames(deg)

write.csv(deg, file.path(wd, "results", "DEG_limma_full.csv"), row.names = FALSE)

sig_strict <- subset(deg, adj.P.Val < 0.05 & abs(logFC) > 1)
sig_loose  <- subset(deg, P.Value < 0.01 & abs(logFC) > 1)
cat("严格阈值 (adj.P<0.05 & |log2FC|>1):", nrow(sig_strict),
    " 上调", sum(sig_strict$logFC > 0), " 下调", sum(sig_strict$logFC < 0), "\n")
cat("宽松阈值 (P<0.01 & |log2FC|>1):", nrow(sig_loose),
    " 上调", sum(sig_loose$logFC > 0), " 下调", sum(sig_loose$logFC < 0), "\n")
write.csv(sig_strict, file.path(wd, "results", "DEG_sig_adjP.csv"), row.names = FALSE)
write.csv(sig_loose,  file.path(wd, "results", "DEG_sig_exploratory.csv"), row.names = FALSE)

# ---------- 火山图 ----------
deg$sig <- "NS"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC > 1]  <- "Up (adj.P<0.05)"
deg$sig[deg$adj.P.Val < 0.05 & deg$logFC < -1] <- "Down (adj.P<0.05)"
p <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Up (adj.P<0.05)" = "#D85A30",
                                "Down (adj.P<0.05)" = "#185FA5", "NS" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3) +
  labs(title = "GSE113212: elderly hypertrophied vs young control",
       x = "log2 fold change", y = "-log10 adjusted P") +
  theme_bw(base_size = 12)
ggsave(file.path(wd, "figures", "Fig_volcano_GSE113212.pdf"), p, width = 6, height = 5)
ggsave(file.path(wd, "figures", "Fig_volcano_GSE113212.png"), p, width = 6, height = 5, dpi = 300)

cat("\nTop 15 差异基因:\n")
print(head(deg[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 15))

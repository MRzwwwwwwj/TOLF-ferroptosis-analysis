# ============================================================
# TOLF Step 03: 多模式 PCD 评分 (GSVA + ssGSEA)  【核心创新图】
# 输入: results/expr_gene_GSE113212.csv + pcd_genesets/gs_*.txt
# 输出: results/PCD_scores_GSE113212.csv (GSVA 与 ssGSEA 两套)
#       figures/Fig_PCD_scores_boxplot.pdf/png
# 方法学: 5 个模块分数分别报告, 不合成单一总分(避免直接相加);
#         n=4v4 用 wilcox 检验, 如实标注探索性
# ============================================================

.libPaths("D:/ARLARS/R-4.5.3/library")
suppressPackageStartupMessages({library(GSVA); library(ggplot2); library(reshape2)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
dir.create(file.path(wd, "figures"), showWarnings = FALSE)

expr <- as.matrix(read.csv(file.path(wd, "results", "expr_gene_GSE113212.csv"),
                           row.names = 1, check.names = FALSE))
group <- read.csv(file.path(wd, "results", "group_GSE113212.csv"), stringsAsFactors = FALSE)
expr <- expr[, group$sample]

# ---------- 读 5 套 PCD 基因集 ----------
gs_dir <- file.path(wd, "pcd_genesets")
sets <- list()
for (f in list.files(gs_dir, pattern = "^gs_.*\\.txt$", full.names = TRUE)) {
  nm <- sub("^gs_|\\.txt$", "", basename(f))
  sets[[nm]] <- trimws(readLines(f))
  sets[[nm]] <- sets[[nm]][sets[[nm]] != ""]
}
cat("基因集:", sapply(sets, length), "\n")

# 覆盖度检查
for (nm in names(sets)) {
  ov <- sum(sets[[nm]] %in% rownames(expr))
  cat(sprintf("%-15s 覆盖 %d/%d\n", nm, ov, length(sets[[nm]])))
}

# ---------- GSVA 与 ssGSEA (GSVA>=1.50 新 API) ----------
gp <- gsvaParam(expr, sets, kcdf = "Gaussian")
gsva_sc <- gsva(gp, verbose = FALSE)
sp <- ssgseaParam(expr, sets)
ssgsea_sc <- gsva(sp, verbose = FALSE)

write.csv(data.frame(Module = rownames(gsva_sc), gsva_sc, check.names = FALSE),
          file.path(wd, "results", "PCD_scores_GSVA_GSE113212.csv"), row.names = FALSE)
write.csv(data.frame(Module = rownames(ssgsea_sc), ssgsea_sc, check.names = FALSE),
          file.path(wd, "results", "PCD_scores_ssGSEA_GSE113212.csv"), row.names = FALSE)

# ---------- 组间比较 (wilcox, 探索性) ----------
cmp <- function(sc, tag) {
  out <- data.frame()
  for (m in rownames(sc)) {
    y <- sc[m, group$group == "young_control"]
    e <- sc[m, group$group == "elderly_hypertrophied"]
    wt <- suppressWarnings(wilcox.test(e, y))
    out <- rbind(out, data.frame(Module = m, method = tag,
                                 mean_young = mean(y), mean_elderly = mean(e),
                                 direction = ifelse(mean(e) > mean(y), "up_in_elderly", "down_in_elderly"),
                                 wilcox_P = wt$p.value))
  }
  out
}
res <- rbind(cmp(gsva_sc, "GSVA"), cmp(ssgsea_sc, "ssGSEA"))
write.csv(res, file.path(wd, "results", "PCD_scores_groupcompare.csv"), row.names = FALSE)
print(res)

# ---------- 箱线图 (ssGSEA 主图) ----------
df <- melt(ssgsea_sc)
colnames(df) <- c("Module", "sample", "score")
df$group <- group$group[match(df$sample, group$sample)]
p <- ggplot(df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2) +
  scale_fill_manual(values = c("young_control" = "#185FA5", "elderly_hypertrophied" = "#D85A30")) +
  facet_wrap(~Module, scales = "free_y", nrow = 2) +
  labs(title = "Multi-mode PCD scores in GSE113212 (ssGSEA)", y = "ssGSEA score", x = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(wd, "figures", "Fig_PCD_scores_boxplot.pdf"), p, width = 9, height = 5.5)
ggsave(file.path(wd, "figures", "Fig_PCD_scores_boxplot.png"), p, width = 9, height = 5.5, dpi = 300)
cat("\n=== 完成 ===\n")

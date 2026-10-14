#!/usr/bin/env Rscript
# ============================================================================
# 13_ckpt_M1M2.R — 一档加分项②：铁死亡评分 × 免疫检查点 / M1-M2 极化
#   A. 5 套 PCD 模块评分 × 免疫检查点基因 Spearman 相关热图（铁死亡特异性）
#   B. 铁死亡评分 × 巨噬细胞极化（ssGSEA M1/M2 评分 + CIBERSORT M0/M1/M2 浸润）
# 输入：results/expr_gene_GSE113212.csv, PCD_scores_ssGSEA_GSE113212.csv,
#       immune_cibersort_GSE113212.csv
# 输出：figures/Fig_ckpt_A/B, results/checkpoint_cor.csv, M1M2_*.csv
# ★ 包：GSVA (ssGSEA), pheatmap, stats::cor.test（包函数）
# ★ [复核] 标记为人工复核点
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(GSVA); library(pheatmap); library(ggplot2)
})

BASE <- getwd(); R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")

expr <- fread(file.path(R, "expr_gene_GSE113212.csv"), data.table = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
mat <- as.matrix(expr)

pcd <- fread(file.path(R, "PCD_scores_ssGSEA_GSE113212.csv"), data.table = FALSE)
pcd_scores <- as.matrix(pcd[, -1]); rownames(pcd_scores) <- sub("\\.txt$", "", pcd$Module)
samples <- colnames(pcd_scores)

## ---------- 1. 免疫检查点基因集（文献常用 24 个，含抑制/共刺激） ----------
ck_genes <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","VSIR",
              "ICOS","ICOSLG","CD276","VTCN1","IDO1","IDO2","CD80","CD86",
              "TNFRSF9","TNFRSF4","TNFRSF18","CD40","CD70","NRP1","LGALS9","LILRB2")
ck_use <- intersect(ck_genes, rownames(mat))
cat("检查点基因覆盖:", length(ck_use), "/", length(ck_genes), "\n")
cat("缺失:", paste(setdiff(ck_genes, ck_use), collapse = ", "), "\n")

## ---------- 2. Spearman 相关：5 PCD 模块 × 检查点 ----------
mods <- c("ferroptosis","cuproptosis","disulfidptosis","panoptosis","pyroptosis")
mods <- intersect(mods, rownames(pcd_scores))
cor_mat <- matrix(NA_real_, nrow = length(mods), ncol = length(ck_use),
                  dimnames = list(mods, ck_use))
p_mat   <- cor_mat
for (m in mods) for (g in ck_use) {
  ct <- suppressWarnings(cor.test(pcd_scores[m, ], mat[g, ], method = "spearman"))
  cor_mat[m, g] <- ct$estimate; p_mat[m, g] <- ct$p.value
}
fwrite(data.frame(module = rep(mods, length(ck_use)),
                  checkpoint = rep(ck_use, each = length(mods)),
                  rho = as.vector(cor_mat), P = as.vector(p_mat)),
       file.path(R, "checkpoint_cor.csv"))

# 热图（显示 r 值 + 显著性星号）
lab <- matrix("", nrow(cor_mat), ncol(cor_mat), dimnames = dimnames(cor_mat))
lab[] <- sprintf("%.2f%s", cor_mat,
                 ifelse(p_mat < 0.001, "***", ifelse(p_mat < 0.01, "**",
                        ifelse(p_mat < 0.05, "*", ""))))
rownames(cor_mat) <- paste0(rownames(cor_mat), " (", 1:nrow(cor_mat), ")")
png(file.path(F, "Fig_ckpt_A_ckpt_corheatmap.png"), width = 3600, height = 1500, res = 300)
pheatmap(cor_mat, display_numbers = lab, cluster_cols = TRUE, cluster_rows = FALSE,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color = NA, fontsize_number = 8,
         main = "Spearman correlation: PCD module scores x immune checkpoints")
dev.off()
# [复核] 重点看 ferroptosis 行的检查点相关性是否强于其他 PCD 行（铁死亡特异性）

## ---------- 3. M1/M2 极化 ----------
m1_set <- c("CD80","CD86","NOS2","IL1B","IL6","TNF","CXCL9","CXCL10","CXCL11",
            "SOCS3","CCR7","TLR2","TLR4")
m2_set <- c("CD163","MRC1","MSR1","CD209","CCL18","CCL22","IL10","ARG1",
            "ALOX15","MMP9","TGFB1")
m1_use <- intersect(m1_set, rownames(mat)); m2_use <- intersect(m2_set, rownames(mat))
cat("M1 覆盖:", length(m1_use), "/", length(m1_set),
    "; M2 覆盖:", length(m2_use), "/", length(m2_set), "\n")

par <- ssgseaParam(exprData = mat[, samples],
                   geneSets = list(M1 = m1_use, M2 = m2_use))
mm <- gsva(par, verbose = FALSE)
fwrite(data.frame(module = rownames(mm), t(mm)),
       file.path(R, "M1M2_ssgsea_GSE113212.csv"))

# 与 5 个 PCD 模块的相关
rows <- list(`M1 score` = mm["M1", ], `M2 score` = mm["M2", ],
             `M1/M2 ratio` = mm["M1", ] - mm["M2", ])
# CIBERSORT 巨噬细胞亚群
cib <- fread(file.path(R, "immune_cibersort_GSE113212.csv"), data.table = FALSE)
mac_cols <- grep("Macrophage", colnames(cib), value = TRUE)
for (cc in mac_cols) rows[[cc]] <- cib[[cc]][match(samples, cib$Sample)]

res <- do.call(rbind, lapply(names(rows), function(nm) {
  do.call(rbind, lapply(mods, function(m) {
    ct <- suppressWarnings(cor.test(rows[[nm]], pcd_scores[m, ], method = "spearman"))
    data.frame(marker = nm, module = m, rho = unname(ct$estimate), P = ct$p.value)
  }))
}))
fwrite(res, file.path(R, "M1M2_PCD_cor.csv"))

# 图：分组条形图（按 module 分面）
res$sig <- ifelse(res$P < 0.001, "***", ifelse(res$P < 0.01, "**",
              ifelse(res$P < 0.05, "*", "")))
p <- ggplot(res, aes(marker, rho, fill = rho > 0)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sig), vjust = ifelse(res$rho >= 0, -0.3, 1.1), size = 3.2) +
  facet_wrap(~ module, nrow = 1) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
  labs(x = NULL, y = "Spearman rho (vs PCD module score)",
       title = "Macrophage polarization x PCD module scores") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(F, "Fig_ckpt_B_M1M2_cor.png"), p, width = 11, height = 4.5, dpi = 300)
# [复核] 预期 ferroptosis 面板中 M2/CIBERSORT M2 呈正相关（SPP1+ 属 M2 样）

# 附：铁死亡 vs M2 散点（叙述用数字）
ct <- suppressWarnings(cor.test(pcd_scores["ferroptosis", ], mm["M2", ], method = "spearman"))
cat(sprintf("\nferroptosis ssGSEA x M2 score: rho=%.3f, P=%.3g\n", ct$estimate, ct$p.value))
ct2 <- suppressWarnings(cor.test(pcd_scores["ferroptosis", ], mm["M1", ], method = "spearman"))
cat(sprintf("ferroptosis ssGSEA x M1 score: rho=%.3f, P=%.3g\n", ct2$estimate, ct2$p.value))

cat("\n=== 13 checkpoint/M1M2 done ===\n")

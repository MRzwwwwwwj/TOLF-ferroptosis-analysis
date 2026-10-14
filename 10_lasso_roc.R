#!/usr/bin/env Rscript
# ============================================================================
# 10_lasso_roc.R — 加分项①：LASSO-logistic 诊断模型 + ROC（LOOCV）
# 特征：探索性 DEG（720 个）；响应：老年肥厚 vs 年轻对照
# 全部用标准包：glmnet（LASSO + LOOCV）、pROC（ROC/AUC/CI）
# 输出：results/LASSO_signature_genes.csv, results/ROC_AUC.csv,
#       figures/Fig_LASSO_cvcurve.png, Fig_LASSO_signature.png, Fig_ROC.png
# 复核：install.packages 不可省——本脚本依赖 glmnet + pROC
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")
dir.create(F, showWarnings = FALSE, recursive = TRUE)

# ---------- 数据 ----------
expr <- fread(file.path(R, "expr_gene_GSE113212.csv"), data.table = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
grp  <- fread(file.path(R, "group_GSE113212.csv"), data.table = FALSE)
y    <- as.integer(grp$group == "elderly_hypertrophied")
names(y) <- grp$sample
y <- y[colnames(expr)]

deg_exp <- fread(file.path(R, "DEG_sig_exploratory.csv"), data.table = FALSE)
sig_genes <- intersect(deg_exp$Gene, rownames(expr))   # 探索性阈值 DEG
cat("Features (exploratory DEG):", length(sig_genes), "\n")

X <- t(as.matrix(expr[sig_genes, , drop = FALSE]))     # 样本 x 基因

# ---------- LASSO：LOOCV（nfolds = n） ----------
set.seed(42)
cvfit <- cv.glmnet(X, y, family = "binomial",
                   nfolds = nrow(X),           # Leave-One-Out
                   type.measure = "deviance",  # LOOCV 下 AUC 未定义，deviance 为 glmnet 标准 CV 指标
                   alpha = 1, standardize = TRUE, keep = TRUE)

# CV 曲线（deviance；ROC 由下方 fit.preval 样本外预测经 pROC 计算）
p_cv <- ggplot(data.frame(lambda = log(cvfit$lambda),
                          dev = cvfit$cvm,
                          lo = cvfit$cvlo, hi = cvfit$cvup),
               aes(x = lambda, y = dev)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), color = "grey70", width = 0.05) +
  geom_line(color = "#2166AC") + geom_point(size = 1.5, color = "#2166AC") +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed", color = "#B2182B") +
  annotate("text", x = log(cvfit$lambda.min), y = max(cvfit$cvm),
           label = "lambda.min", hjust = -0.1, color = "#B2182B") +
  labs(x = "log(lambda)", y = "Binomial deviance (LOOCV)",
       title = "LASSO logistic regression (LOOCV), GSE113212") +
  theme_bw(base_size = 12)
ggsave(file.path(F, "Fig_LASSO_cvcurve.png"), p_cv, width = 7, height = 5, dpi = 300)

cat("lambda.min =", cvfit$lambda.min,
    " LOOCV deviance =", cvfit$cvm[cvfit$lambda == cvfit$lambda.min], "\n")

# ---------- 签名基因（lambda.min 非零系数） ----------
co <- coef(cvfit, s = "lambda.min")
sig <- summary(co)
sig_genes_sel <- rownames(co)[sig[-1, 1]]          # 排除截距
beta <- sig[-1, ncol(sig)]
lasso_tab <- data.frame(Gene = sig_genes_sel,
                        beta = round(as.numeric(beta), 4),
                        direction = ifelse(beta > 0, "up_in_hypertrophy", "down_in_hypertrophy"))
lasso_tab <- lasso_tab[order(-abs(lasso_tab$beta)), ]
fwrite(lasso_tab, file.path(R, "LASSO_signature_genes.csv"))
cat("=== LASSO signature (", nrow(lasso_tab), "genes) ===\n"); print(lasso_tab)

# 系数条形图
p_sig <- ggplot(lasso_tab, aes(x = beta, y = reorder(Gene, beta),
                               fill = beta > 0)) +
  geom_col() + scale_fill_manual(values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"),
                                 guide = "none") +
  labs(x = "LASSO coefficient (lambda.min)", y = NULL,
       title = "LASSO signature genes") +
  theme_bw(base_size = 12)
ggsave(file.path(F, "Fig_LASSO_signature.png"), p_sig,
       width = 7, height = max(3, 0.35 * nrow(lasso_tab) + 1.5), dpi = 300)

# ---------- ROC：三个模型对比 ----------
# 1) LASSO 全模型（LOOCV 样本外预测，cv.glmnet keep=TRUE 的 fit.preval）
prev <- cvfit$fit.preval[, which(cvfit$lambda == cvfit$lambda.min)]
roc1 <- roc(y, as.numeric(prev), quiet = TRUE)

# 2) 单指标：铁死亡 GSVA 评分
pcd <- fread(file.path(R, "PCD_scores_GSVA_GSE113212.csv"), data.table = FALSE)
ferr <- as.numeric(pcd[grep("ferroptosis", pcd$Module, ignore.case = TRUE)[1], -1])
roc2 <- roc(y, ferr, quiet = TRUE)

# 3) 单指标：巨噬细胞 ssGSEA
imm <- fread(file.path(R, "immune_ssgsea_GSE113212.csv"), data.table = FALSE)
mac <- as.numeric(imm[grep("^Macrophages$", imm$CellType)[1], -1])
roc3 <- roc(y, mac, quiet = TRUE)

auc_tab <- data.frame(
  Model = c("LASSO signature (LOOCV)", "Ferroptosis GSVA score", "Macrophage ssGSEA score"),
  AUC = c(roc1$auc, roc2$auc, roc3$auc),
  CI95 = sprintf("%.3f-%.3f",
                 ci.auc(roc1, quiet = TRUE)[1], ci.auc(roc1, quiet = TRUE)[3]))
fwrite(auc_tab, file.path(R, "ROC_AUC.csv"))
cat("=== ROC AUC ===\n"); print(auc_tab)

roc_df1 <- data.frame(FPR = 1 - roc1$specificities, TPR = roc1$sensitivities,
                      Model = sprintf("LASSO signature (AUC=%.2f)", roc1$auc))
roc_df2 <- data.frame(FPR = 1 - roc2$specificities, TPR = roc2$sensitivities,
                      Model = sprintf("Ferroptosis score (AUC=%.2f)", roc2$auc))
roc_df3 <- data.frame(FPR = 1 - roc3$specificities, TPR = roc3$sensitivities,
                      Model = sprintf("Macrophage score (AUC=%.2f)", roc3$auc))
p_roc <- ggplot(rbind(roc_df1, roc_df2, roc_df3),
                aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "ROC: hypertrophied vs control LF (GSE113212, n=8)") +
  coord_equal() + theme_bw(base_size = 12) +
  theme(legend.position = c(0.72, 0.22),
        legend.background = element_rect(fill = "white", color = "grey80"))
# 2026-10-14 修复：长标题右缘被裁——缩短标题（去掉 "discriminating"），恢复 width 6
ggsave(file.path(F, "Fig_ROC_LASSO_PCD.png"), p_roc, width = 6, height = 6, dpi = 300)

cat("=== 10 done ===\n")

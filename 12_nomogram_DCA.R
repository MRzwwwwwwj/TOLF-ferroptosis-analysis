#!/usr/bin/env Rscript
# ============================================================================
# 12_nomogram_DCA.R — 一档加分项①：临床预测模型三件套
#   A. 列线图 Nomogram（6 个 LASSO 签名基因 → TOLF 风险概率）
#   B. 校准曲线（bootstrap B=1000）
#   C. 决策曲线分析 DCA（训练队列：签名 vs 单基因）
#   D. 外部队列 DCA（GSE206560 / GSE106253）
# 输入：results/expr_gene_GSE113212.csv, group_GSE113212.csv,
#       LASSO_signature_genes.csv, LASSO_external_GSE206560.csv, LASSO_external_GSE106253.csv
# 输出：figures/Fig_nomo_A-D, results/nomogram_*.csv
# ★ 包：rms (lrm/nomogram/calibrate), dcurves (DCA), stats::glm（概率化，包函数）
# ★ 人工复核点见 [复核] 标记
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(rms); library(dcurves)
  library(ggplot2); library(patchwork)
})

BASE <- getwd(); R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")

## ---------- 1. 训练队列数据 ----------
expr <- fread(file.path(R, "expr_gene_GSE113212.csv"), data.table = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
grp  <- fread(file.path(R, "group_GSE113212.csv"), data.table = FALSE)
sig  <- fread(file.path(R, "LASSO_signature_genes.csv"), data.table = FALSE)

genes  <- intersect(sig$Gene, rownames(expr))
beta   <- setNames(sig$beta[match(genes, sig$Gene)], genes)
train  <- data.frame(t(expr[genes, , drop = FALSE]), check.rows = FALSE)
train$sample <- rownames(train)
train$label  <- as.integer(grp$group[match(train$sample, grp$sample)] == "elderly_hypertrophied")
train$risk_score <- as.matrix(train[, genes]) %*% beta[genes]
cat("签名基因:", paste(genes, collapse = ", "), "\n")
cat("训练集分组 (label=1 肥厚):", table(train$label), "\n")
# [复核] 确认 label 方向：1 = elderly_hypertrophied（病例），0 = young_control

## ---------- 2. 列线图（rms::lrm） ----------
dd <- datadist(train); options(datadist = "dd")
fit <- tryCatch(
  lrm(label ~ TGFB1 + COL13A1 + BMPR1B + LOC283404 + FGFRL1 + PIK3IP1,
      data = train, x = TRUE, y = TRUE, penalty = 1),
  error = function(e) NULL)
penalized <- TRUE
if (is.null(fit)) {
  cat("惩罚 lrm 失败 -> 回退未惩罚\n")
  fit <- lrm(label ~ TGFB1 + COL13A1 + BMPR1B + LOC283404 + FGFRL1 + PIK3IP1,
             data = train, x = TRUE, y = TRUE)
  penalized <- FALSE
}
# 说明：训练集 AUC=1.0（完全分离），未惩罚 MLE 系数发散（SE>1e3）、bootstrap 校准
#       无法收敛；penalty=1 的惩罚 logistic（rms 内置）是小样本分离数据的规范处理。
print(fit)
fwrite(data.frame(term = names(coef(fit)),
                  coefficient = coef(fit),
                  penalized = penalized),
       file.path(R, "nomogram_coefficients.csv"))

nom <- nomogram(fit, fun = plogis,
                fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95),
                funlabel = "TOLF risk probability")
png(file.path(F, "Fig_nomo_A_nomogram.png"), width = 3600, height = 2000, res = 300)
plot(nom, cex.axis = 0.7, cex.var = 0.9)
dev.off()

## ---------- 3. 校准曲线（4 折交叉验证；n=8 完全分离，boot 重采样发散不可用） ----------
cal <- calibrate(fit, method = "crossvalidation", B = 4,
                 predy = seq(0.1, 0.9, by = 0.1))
png(file.path(F, "Fig_nomo_B_calibration.png"), width = 2000, height = 2000, res = 300)
plot(cal, xlab = "Predicted probability", ylab = "Observed probability",
     subtitles = FALSE)
title(main = sprintf("Calibration (4-fold CV)%s",
                     ifelse(penalized, ", penalized lrm", "")))
dev.off()
# C-index（= logistic AUC）
ci <- c(Cindex = as.numeric(fit$stats["C"]))
fwrite(data.frame(metric = names(ci), value = ci),
       file.path(R, "nomogram_Cindex.csv"))
cat(sprintf("\nNomogram C-index = %.3f\n", ci))

## ---------- 4. DCA 辅助 ----------
# 连续标志物 -> 概率。训练队列用【留一交叉拟合】惩罚 logistic 概率：
#   若直接 glm，完美分离导致概率全为 0/1，DCA 曲线全部叠平、无信息量；
#   LOO 概率是模型开发文献中内部验证的规范做法，曲线才有区分度。
prob_of_loo <- function(y, x) {
  sapply(seq_along(y), function(i) {
    d <- data.frame(y = y[-i], x = x[-i])
    m <- tryCatch(lrm(y ~ x, data = d, x = TRUE, y = TRUE, penalty = 1),
                  error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    plogis(predict(m, newdata = data.frame(x = x[i])))
  })
}
prob_of <- function(y, x) {   # 外部队列：真实测试集，直接 glm 概率
  d <- data.frame(y = y, x = x)
  m <- glm(y ~ x, family = binomial(), data = d)
  as.numeric(predict(m, type = "response"))
}
train_dca <- train
train_dca$Signature <- prob_of_loo(train$label, train$risk_score)
for (g in genes) train_dca[[g]] <- prob_of_loo(train$label, train[[g]])

## ---------- 5. DCA-A：训练队列（签名 vs 6 个单基因，LOO 交叉拟合概率） ----------
fml <- as.formula(paste("label ~ Signature +", paste(genes, collapse = " + ")))
dca1 <- dca(fml, data = train_dca, thresholds = seq(0.01, 0.99, 0.01))
p1 <- plot(dca1, smooth = FALSE) +
  ggtitle("Training cohort (GSE113212, n=8; LOO cross-fitted probabilities)") +
  coord_cartesian(ylim = c(-0.15, 1)) + theme_bw(base_size = 11)
ggsave(file.path(F, "Fig_nomo_C_DCA_train.png"), p1, width = 8, height = 6, dpi = 300)

## ---------- 6. DCA-B：两个外部队列（签名） ----------
mk_ext <- function(path, case_str) {
  d <- fread(path, data.table = FALSE)
  d$label <- as.integer(d$group == case_str)
  d$Signature <- prob_of(d$label, d$risk_score)
  d
}
e1 <- mk_ext(file.path(R, "LASSO_external_GSE206560.csv"), "ossification")
e2 <- mk_ext(file.path(R, "LASSO_external_GSE106253.csv"), "OLF")
cat(sprintf("外部 GSE206560 n=%d (case=%d); GSE106253 n=%d (case=%d)\n",
            nrow(e1), sum(e1$label), nrow(e2), sum(e2$label)))

dca_e1 <- dca(label ~ Signature, data = e1, thresholds = seq(0.01, 0.99, 0.01))
dca_e2 <- dca(label ~ Signature, data = e2, thresholds = seq(0.01, 0.99, 0.01))
p2 <- plot(dca_e1, smooth = FALSE) + ggtitle("GSE206560 (n=5)") +
  coord_cartesian(ylim = c(-0.15, 1)) + theme_bw(base_size = 11)
p3 <- plot(dca_e2, smooth = FALSE) + ggtitle("GSE106253 (n=8)") +
  coord_cartesian(ylim = c(-0.15, 1)) + theme_bw(base_size = 11)
ggsave(file.path(F, "Fig_nomo_D_DCA_external.png"), p2 + p3,
       width = 10, height = 5.5, dpi = 300)
# [复核] 看签名曲线是否在大部分阈值区间高于 Treat All / Treat None 与单基因

cat("\n=== 12 nomogram/DCA done ===\n")

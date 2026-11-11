#!/usr/bin/env Rscript
# ============================================================================
# 83_auc_ci_round2.R — 二审意见 2：为所有 AUC 估计补 95% CI
# 输入（rerun_2026-10-13/results/）：expr_gene_GSE113212.csv, group_GSE113212.csv,
#   DEG_sig_exploratory.csv, PCD_scores_GSVA_GSE113212.csv,
#   immune_ssgsea_GSE113212.csv, LASSO_external_GSE206560.csv,
#   LASSO_external_GSE106253.csv
# 输出：results/AUC_with_CI_round2.csv + 控制台完整日志
# 方法：DeLong（pROC::ci.auc）+ 2,000 次分层 bootstrap 百分位 CI + 精确 Wilcoxon P
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC)
})
set.seed(42)
R <- Sys.getenv("TOLF_RES", "rerun_2026-10-13/results")  # 主管线结果目录

rd <- function(f) {
  x <- fread(file.path(R, f), data.table = FALSE)
  colnames(x)[1] <- "ID"
  x
}

# ---------- 1. 发现队列：LASSO LOOCV 样本外预测 ----------
expr <- rd("expr_gene_GSE113212.csv"); rownames(expr) <- expr$ID; expr$ID <- NULL
grp  <- rd("group_GSE113212.csv")
y    <- as.integer(grp$group == "elderly_hypertrophied"); names(y) <- grp$ID
y    <- y[colnames(expr)]
deg  <- rd("DEG_sig_exploratory.csv")          # 基因名在列 "Gene"（首列是 logFC）
sig  <- intersect(deg$Gene, rownames(expr))
X    <- t(as.matrix(expr[sig, , drop = FALSE]))
cvfit <- cv.glmnet(X, y, family = "binomial", nfolds = nrow(X),
                   type.measure = "deviance", alpha = 1, standardize = TRUE, keep = TRUE)
p_oof <- as.numeric(cvfit$fit.preval[, which(cvfit$lambda == cvfit$lambda.min)])
p_oof <- 1 / (1 + exp(-p_oof))          # link -> probability

# ---------- 2. 模块分数（GSVA ferroptosis / ssGSEA macrophage） ----------
gsva <- rd("PCD_scores_GSVA_GSE113212.csv"); rownames(gsva) <- gsva$ID; gsva$ID <- NULL
fer  <- as.numeric(gsva["ferroptosis.txt", names(y)])
imm  <- rd("immune_ssgsea_GSE113212.csv"); rownames(imm) <- imm$ID; imm$ID <- NULL
mac  <- as.numeric(imm["Macrophages", names(y)])

# ---------- 3. 外部队列 ----------
e1 <- rd("LASSO_external_GSE206560.csv")
e2 <- rd("LASSO_external_GSE106253.csv")

# ---------- 4. 统一计算 ----------
summ <- function(label, resp, pred, n_boot = 2000) {
  r  <- roc(resp, pred, quiet = TRUE, direction = "<")
  ci_d <- as.numeric(ci.auc(r, method = "delong"))
  ci_b <- tryCatch(as.numeric(ci.auc(r, method = "bootstrap", boot.n = n_boot,
                                     boot.stratified = TRUE)), error = function(e) c(NA, NA, NA))
  pw <- suppressWarnings(wilcox.test(pred[resp == 1], pred[resp == 0], exact = TRUE)$p.value)
  data.frame(Model = label, n_case = sum(resp == 1), n_ctrl = sum(resp == 0),
             AUC = round(as.numeric(auc(r)), 3),
             CI95_DeLong = sprintf("%.3f-%.3f", ci_d[1], ci_d[3]),
             CI95_bootstrap = ifelse(any(is.na(ci_b)), "NA",
                                     sprintf("%.3f-%.3f", ci_b[1], ci_b[3])),
             P_exact_Wilcoxon = signif(pw, 3),
             stringsAsFactors = FALSE)
}

res <- rbind(
  summ("LASSO signature, LOOCV out-of-fold (GSE113212, 4 vs 4)", y, p_oof),
  summ("Ferroptosis GSVA score (GSE113212, 4 vs 4)",             y, fer),
  summ("Macrophage ssGSEA score (GSE113212, 4 vs 4)",             y, mac),
  summ("LASSO signature, external (GSE206560 Input, 3 vs 2)",
       as.integer(e1$group != "non_ossification"), e1$risk_score),
  summ("LASSO signature, external (GSE106253, 4 vs 4)",
       as.integer(e2$group != "Control"), e2$risk_score)
)

# 训练集内部 LOOCV AUC 的“nested”复核（沿用正文已有的置换 P）
cat("\n=== AUC with 95% CI (round 2) ===\n")
print(res, row.names = FALSE)
fwrite(res, file.path(R, "AUC_with_CI_round2.csv"))

cat("\n--- 关键校验 ---\n")
cat("完全分离队列的 DeLong CI 是否退化(1.000-1.000):",
    all(res$CI95_DeLong == "1.000-1.000"), "\n")
cat("bootstrap CI 是否同样退化:",
    all(res$CI95_bootstrap %in% c("1.000-1.000", "NA")), "\n")
cat("LOOCV 预测值范围:", sprintf("%.4f - %.4f", min(p_oof), max(p_oof)), "\n")

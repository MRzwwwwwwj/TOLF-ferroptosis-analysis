#!/usr/bin/env Rscript
# ============================================================================
# 39_nested_cv_signature.R — E 包：严格嵌套交叉验证 + permutation test
# 回应 R1/R3/R4："特征选择在 CV 外、AUC=1.0 系过拟合、lambda.min 在边界"
# 设计：
#   外层 LOOCV（8 折）：每折在 7 个训练样本上【重做 limma 筛选 + cv.glmnet】，
#     预测留出样本 → 8 个真正样本外预测 → AUC + bootstrap 95% CI
#   permutation：200 次标签打乱，每次完整重跑 nested CV → null 分布 + 精确 p
#   基因稳定性：记录每折 lambda.min 非零系数基因 → 入选频率
#   lambda 边界检查：全模型 + 各折 lambda.min 是否落在网格端点
# 输出：results/nested_CV_*.csv, figures/Fig_nested_CV_permutation.png
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(limma)
  library(ggplot2); library(parallel)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")

expr <- read.csv(file.path(R, "expr_gene_GSE113212.csv"), row.names = 1, check.names = FALSE)
group <- read.csv(file.path(R, "group_GSE113212.csv"), stringsAsFactors = FALSE)
expr <- expr[, group$sample]
y <- as.integer(group$group == "elderly_hypertrophied")
n <- length(y)
cat("n =", n, " labels:", paste(y, collapse = ","), "\n")

emat <- as.matrix(expr)

# ---------- 一轮完整 nested CV ----------
nested_once <- function(ylabels, seed = 42, verbose = FALSE) {
  preds <- numeric(n); genes_used <- vector("list", n)
  lam_pos <- integer(n); lam_len <- integer(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    g <- factor(ylabels[tr], levels = c(0, 1))
    design <- model.matrix(~ 0 + g); colnames(design) <- c("ctrl", "case")
    fit <- suppressWarnings(lmFit(emat[, tr], design))
    cm <- makeContrasts(DEG = case - ctrl, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm))
    tt <- topTable(fit2, number = Inf, sort.by = "none")
    sel <- rownames(tt)[tt$P.Value < 0.01 & abs(tt$logFC) > 1]
    if (length(sel) < 5) {
      sel <- rownames(tt)[order(tt$P.Value)][seq_len(min(50, nrow(tt)))]
    }
    Xtr <- t(emat[sel, tr, drop = FALSE])
    set.seed(seed + i)
    cv <- suppressWarnings(cv.glmnet(Xtr, ylabels[tr], family = "binomial",
                                     nfolds = length(tr), type.measure = "deviance",
                                     alpha = 1, standardize = TRUE))
    preds[i] <- as.numeric(predict(cv, newx = t(emat[sel, i, drop = FALSE]),
                                   s = "lambda.min", type = "link"))
    co <- coef(cv, s = "lambda.min")
    genes_used[[i]] <- rownames(co)[which(as.numeric(co) != 0)][-1]
    lam_pos[i] <- which(cv$lambda == cv$lambda.min)[1]
    lam_len[i] <- length(cv$lambda)
    if (verbose) cat(sprintf("  fold %d: %d DEGs, %d genes selected, lambda.min %d/%d\n",
                             i, length(sel), length(genes_used[[i]]), lam_pos[i], lam_len[i]))
  }
  list(pred = preds, genes = genes_used, lam_pos = lam_pos, lam_len = lam_len)
}

# ---------- [1] 全模型 lambda 边界检查（回应"lambda.min 在边界"） ----------
cat("\n[1] 全模型 lambda.min 位置\n")
deg_exp <- fread(file.path(R, "DEG_sig_exploratory.csv"), data.table = FALSE)
sel_full <- intersect(deg_exp$Gene, rownames(emat))
set.seed(42)
cv_full <- cv.glmnet(t(emat[sel_full, ]), y, family = "binomial", nfolds = n,
                     type.measure = "deviance", alpha = 1, standardize = TRUE)
pos_full <- which(cv_full$lambda == cv_full$lambda.min)[1]
cat(sprintf("  全模型 lambda.min 位于网格 %d/%d（1=最强正则端；%d=最弱端）\n",
            pos_full, length(cv_full$lambda), length(cv_full$lambda)))

# ---------- [2] 观测数据的 nested CV ----------
cat("\n[2] 观测数据 nested CV（外层 LOOCV，折内重做 limma+cv.glmnet）\n")
obs <- nested_once(y, seed = 42, verbose = TRUE)
pred_df <- data.frame(sample = group$sample, label = y,
                      oof_pred = round(obs$pred, 4),
                      fold_genes = sapply(obs$genes, function(g) paste(g, collapse = ";")),
                      lambda_pos = obs$lam_pos, lambda_grid = obs$lam_len)
fwrite(pred_df, file.path(R, "nested_CV_predictions.csv"))

roc_obs <- roc(y, obs$pred, quiet = TRUE)
auc_obs <- as.numeric(roc_obs$auc)
cat("\n>>> nested-CV out-of-fold AUC =", round(auc_obs, 4), "\n")
print(pred_df[, c("sample", "label", "oof_pred")])

# bootstrap 95% CI（固定样本外预测，重抽样样本对）
set.seed(123)
B_boot <- 2000
boot_auc <- replicate(B_boot, {
  idx <- sample(seq_len(n), replace = TRUE)
  if (length(unique(y[idx])) < 2) return(NA)
  as.numeric(auc(roc(y[idx], obs$pred[idx], quiet = TRUE)))
})
boot_auc <- boot_auc[!is.na(boot_auc)]
ci_lo <- quantile(boot_auc, 0.025); ci_hi <- quantile(boot_auc, 0.975)
cat(sprintf(">>> bootstrap 95%% CI: %.3f - %.3f (B=%d)\n", ci_lo, ci_hi, B_boot))

# ---------- [3] 基因稳定性 ----------
gene_freq <- sort(table(unlist(obs$genes)) / n, decreasing = TRUE)
stab_df <- data.frame(gene = names(gene_freq), fold_freq = as.numeric(gene_freq))
fwrite(stab_df, file.path(R, "nested_CV_gene_stability.csv"))
cat("\n[3] 折内入选基因频率（top 15）:\n"); print(head(stab_df, 15))
orig6 <- c("TGFB1", "COL13A1", "BMPR1B", "LOC283404", "FGFRL1", "PIK3IP1")
cat("原六基因在 nested 各折的入选频率:\n")
for (g6 in orig6) cat(sprintf("  %s: %.2f\n", g6, ifelse(g6 %in% names(gene_freq), gene_freq[[g6]], 0)))

# ---------- [4] permutation test（并行） ----------
B_perm <- 200
cat(sprintf("\n[4] permutation test：%d 次标签打乱完整重跑 nested CV（并行）\n", B_perm))
run_one_perm <- function(b, emat, n, y) {
  suppressPackageStartupMessages({library(glmnet); library(pROC); library(limma)})
  set.seed(1000 + b)
  yb <- sample(y)
  preds <- numeric(n)
  for (i in seq_len(n)) {
    tr <- setdiff(seq_len(n), i)
    g <- factor(yb[tr], levels = c(0, 1))
    design <- model.matrix(~ 0 + g); colnames(design) <- c("ctrl", "case")
    fit <- suppressWarnings(lmFit(emat[, tr], design))
    cm <- makeContrasts(DEG = case - ctrl, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm))
    tt <- topTable(fit2, number = Inf, sort.by = "none")
    sel <- rownames(tt)[tt$P.Value < 0.01 & abs(tt$logFC) > 1]
    if (length(sel) < 5) sel <- rownames(tt)[order(tt$P.Value)][seq_len(min(50, nrow(tt)))]
    Xtr <- t(emat[sel, tr, drop = FALSE])
    set.seed(2000 + b * 10 + i)
    cv <- suppressWarnings(cv.glmnet(Xtr, yb[tr], family = "binomial",
                                     nfolds = length(tr), type.measure = "deviance",
                                     alpha = 1, standardize = TRUE))
    preds[i] <- as.numeric(predict(cv, newx = t(emat[sel, i, drop = FALSE]),
                                   s = "lambda.min", type = "link"))
  }
  as.numeric(auc(roc(yb, preds, quiet = TRUE)))
}

ncores <- max(1, detectCores() - 1)
cat("使用", ncores, "核并行；预计耗时 =", B_perm, "轮 x ~35s /", ncores, "\n")
t0 <- Sys.time()
cl <- makeCluster(ncores)
clusterExport(cl, c("emat", "n", "y"), envir = environment())
null_auc <- parLapply(cl, seq_len(B_perm), run_one_perm, emat = emat, n = n, y = y)
stopCluster(cl)
null_auc <- as.numeric(null_auc)
cat("permutation 耗时:", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n")

p_perm <- (1 + sum(null_auc >= auc_obs)) / (B_perm + 1)
cat(sprintf(">>> permutation p = %.4f（观测 AUC=%.3f；null 均值=%.3f, 95%% 上限=%.3f）\n",
            p_perm, auc_obs, mean(null_auc), quantile(null_auc, 0.95)))

fwrite(data.frame(b = seq_len(B_perm), null_auc = null_auc),
       file.path(R, "nested_CV_permutation_null.csv"))

# ---------- [5] 汇总 + 图 ----------
summary_df <- data.frame(
  metric = c("Full-model LOOCV AUC (特征选择在 CV 外, 原投稿口径)",
             "Full-model lambda.min grid position",
             "Nested-CV out-of-fold AUC (特征选择在 CV 内, 严格无泄漏)",
             "Nested-CV AUC bootstrap 95% CI",
             sprintf("Permutation p (B=%d)", B_perm),
             "Permutation null AUC mean",
             "Permutation null AUC 95% upper"),
  value = c("1.000", sprintf("%d/%d", pos_full, length(cv_full$lambda)),
            sprintf("%.3f", auc_obs),
            sprintf("%.3f-%.3f", ci_lo, ci_hi),
            sprintf("%.4f", p_perm),
            sprintf("%.3f", mean(null_auc)),
            sprintf("%.3f", quantile(null_auc, 0.95))))
fwrite(summary_df, file.path(R, "nested_CV_summary.csv"))
print(summary_df)

p_null <- ggplot(data.frame(null_auc = null_auc), aes(null_auc)) +
  geom_histogram(bins = 30, fill = "#92C5DE", color = "white") +
  geom_vline(xintercept = auc_obs, color = "#B2182B", linewidth = 1.2) +
  annotate("text", x = auc_obs, y = Inf, vjust = 2, hjust = -0.1, color = "#B2182B",
           label = sprintf("observed nested-CV AUC = %.2f\npermutation p = %.4f", auc_obs, p_perm)) +
  labs(x = "Permuted-label nested-CV AUC (null distribution)",
       y = "Count",
       title = sprintf("Nested CV with in-fold feature selection, GSE113212 (n=%d)", n)) +
  theme_bw(base_size = 12)
ggsave(file.path(F, "Fig_nested_CV_permutation.png"), p_null, width = 7, height = 5, dpi = 300)
cat("\n图已存 figures/Fig_nested_CV_permutation.png\n")
cat("=== 39 done ===\n")

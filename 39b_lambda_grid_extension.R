#!/usr/bin/env Rscript
# ============================================================================
# 39b_lambda_grid_extension.R — 回应 R4 "lambda.min 在网格边界"
# 39 已证实：全模型与 nested 各折 lambda.min 均在 100/100（最弱正则端）。
# 本脚本：lambda.min.ratio 从默认 1e-2 扩到 1e-5（网格延到更弱正则），
# 检验 lambda.min 是否仍贴边、选出基因与样本外预测是否改变。
# 输出：results/lambda_grid_extension.csv
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(glmnet); library(pROC)})

BASE <- getwd(); R <- file.path(BASE, "results")
expr <- read.csv(file.path(R, "expr_gene_GSE113212.csv"), row.names = 1, check.names = FALSE)
group <- read.csv(file.path(R, "group_GSE113212.csv"), stringsAsFactors = FALSE)
expr <- expr[, group$sample]
y <- as.integer(group$group == "elderly_hypertrophied")
emat <- as.matrix(expr)
deg_exp <- fread(file.path(R, "DEG_sig_exploratory.csv"), data.table = FALSE)
sel <- intersect(deg_exp$Gene, rownames(emat))
X <- t(emat[sel, ])

res <- list()
for (ratio in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  set.seed(42)
  cv <- suppressWarnings(cv.glmnet(X, y, family = "binomial", nfolds = length(y),
                                   type.measure = "deviance", alpha = 1,
                                   standardize = TRUE, keep = TRUE,
                                   lambda.min.ratio = ratio, nlambda = 200))
  pos <- which(cv$lambda == cv$lambda.min)[1]
  co <- coef(cv, s = "lambda.min")
  genes <- rownames(co)[which(as.numeric(co) != 0)][-1]
  prev <- cv$fit.preval[, which(cv$lambda == cv$lambda.min)]
  aucv <- as.numeric(auc(roc(y, as.numeric(prev), quiet = TRUE)))
  res[[length(res) + 1]] <- data.frame(
    lambda.min.ratio = ratio, grid_pos = pos, grid_len = length(cv$lambda),
    at_edge = pos == length(cv$lambda),
    n_genes = length(genes), genes = paste(genes, collapse = ";"),
    loocv_auc = aucv)
}
out <- do.call(rbind, res)
fwrite(out, file.path(R, "lambda_grid_extension.csv"))
print(out)
cat("=== 39b done ===\n")

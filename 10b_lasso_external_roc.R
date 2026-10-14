#!/usr/bin/env Rscript
# ============================================================================
# 10b_lasso_external_roc.R — LASSO 6 基因签名在 GSE206560（骨化 vs 非骨化，2v3）的外部验证
# 风险评分 = sum(beta_g * expr_g)；ROC 用 pROC
# 输入：GSE206560/GSE206560_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
#       + Human.GRCh38.p13.annot.tsv.gz（ENSEMBL->Symbol）
#       + results/LASSO_signature_genes.csv
# 输出：results/LASSO_external_GSE206560.csv, figures/Fig_ROC_LASSO_external.png
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(pROC); library(ggplot2)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")
D <- file.path(BASE, "GSE206560")

# ---------- TPM 矩阵 ----------
tpm <- fread(file.path(D, "GSE206560_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
             data.table = FALSE)
annot <- fread(file.path(D, "Human.GRCh38.p13.annot.tsv.gz"), data.table = FALSE)

# 样本名 GSM 对应（07 脚本已确认：6256886/87/88/89/90，其中 886 缺失 counts；
# TPM 文件实际列以文件头为准）
cat("TPM cols:", paste(colnames(tpm)[1:8], collapse = ","), "\n")

# Entrez GeneID -> Symbol（用 match，不能用命名向量整数索引——会按位置错位查找）
m <- match(tpm[[1]], annot$GeneID)
sym <- annot$Symbol[m]
ok <- !is.na(sym) & sym != ""
tpm <- tpm[ok, , drop = FALSE]
tpm$Gene <- sym[ok]
# 按表达量排序后再去重（同 symbol 多 GeneID 保留高表达那条，与 07 一致）
o <- order(rowSums(as.matrix(tpm[, grep("GSM", colnames(tpm)), drop = FALSE])),
           decreasing = TRUE)
tpm <- tpm[o, , drop = FALSE]
tpm <- tpm[!duplicated(tpm$Gene), ]
mat <- as.matrix(tpm[, grep("GSM", colnames(tpm)), drop = FALSE])
rownames(mat) <- tpm$Gene
exprv <- log2(mat + 1)

# ---------- 分组（TPM 实际列：885,887 非骨化 Input；888,889,890 骨化 Input） ----------
non_oss <- intersect(colnames(exprv), c("GSM6256885", "GSM6256886", "GSM6256887"))
oss     <- intersect(colnames(exprv), c("GSM6256888", "GSM6256889", "GSM6256890"))
cat("Non-ossification:", paste(non_oss, collapse = ","), "\n")
cat("Ossification:", paste(oss, collapse = ","), "\n")
stopifnot(length(non_oss) >= 2, length(oss) >= 2)

# ---------- LASSO 签名 ----------
lasso <- fread(file.path(R, "LASSO_signature_genes.csv"), data.table = FALSE)
cat("Signature genes:", paste(lasso$Gene, collapse = ", "), "\n")
avail <- lasso$Gene[lasso$Gene %in% rownames(exprv)]
miss  <- setdiff(lasso$Gene, avail)
if (length(miss)) cat("Missing in GSE206560:", paste(miss, collapse = ", "), "\n")

# 风险评分（可得基因子集，系数用 LASSO 原值）
beta <- setNames(lasso$beta[match(avail, lasso$Gene)], avail)
score <- as.numeric(t(exprv[avail, , drop = FALSE]) %*% beta)
names(score) <- colnames(exprv)

grp <- ifelse(names(score) %in% oss, 1, 0)     # 1 = 骨化
roc_obj <- roc(grp, score, quiet = TRUE)

res_tab <- data.frame(
  sample = names(score), group = ifelse(grp == 1, "ossification", "non_ossification"),
  risk_score = round(score, 4))
fwrite(res_tab, file.path(R, "LASSO_external_GSE206560.csv"))

auc_ext <- data.frame(Model = "LASSO signature (external, GSE206560 3v2)",
                      AUC = roc_obj$auc,
                      CI95 = sprintf("%.3f-%.3f", ci.auc(roc_obj, quiet = TRUE)[1],
                                     ci.auc(roc_obj, quiet = TRUE)[3]))
auc_ext$p_value <- wilcox.test(score ~ factor(grp))$p.value
fwrite(auc_ext, file.path(R, "LASSO_external_AUC.csv"))
cat("=== External ROC ===\n"); print(auc_ext)
cat("Scores by group:\n"); print(tapply(score, grp, function(x) round(x, 3)))

p <- ggplot(data.frame(FPR = 1 - roc_obj$specificities,
                       TPR = roc_obj$sensitivities),
            aes(FPR, TPR)) +
  geom_line(color = "#B2182B", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = sprintf("LASSO signature, external validation (GSE206560, %dv%d)\nAUC = %.2f",
                       length(oss), length(non_oss), roc_obj$auc)) +
  coord_equal() + theme_bw(base_size = 12)
ggsave(file.path(F, "Fig_ROC_LASSO_external.png"), p, width = 6, height = 6, dpi = 300)

cat("=== 10b done ===\n")

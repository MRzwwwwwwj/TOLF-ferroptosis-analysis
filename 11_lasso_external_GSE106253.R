#!/usr/bin/env Rscript
# ============================================================================
# 11_lasso_external_GSE106253.R — 第三队列 bulk 外部验证（GSE106253, 4 OLF vs 4 健康对照）
# 平台：GPL21827 Agilent-079487 Arraystar Human LncRNA microarray V4（含 7,506 mRNA 探针）
# 内容：(a) LASSO 签名风险评分 + ROC（可得基因子集）
#       (b) 铁死亡基因集 ssGSEA 组间比较（防止签名基因探针覆盖不足）
# 输入：GSE106253/ 目录（需先下载数据，见 加分项分析操作指南.md）
#       + results/LASSO_signature_genes.csv
# 输出：results/LASSO_external_GSE106253.csv, LASSO_external_AUC_GSE106253.csv,
#       results/PCD_scores_ssGSEA_GSE106253.csv, figures/Fig_ROC_LASSO_GSE106253.png
#
# ★ 人工复核点见 [复核] 标记
# ============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(data.table); library(pROC)
  library(ggplot2)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")
D <- file.path(BASE, "GSE106253")
dir.create(D, showWarnings = FALSE, recursive = TRUE)

# ---------- 读取数据（优先本地文件；缺失则自动从 GEO 下载） ----------
cat("读取 GSE106253（本地有 series matrix 则直接用，否则联网下载）...\n")
gse <- getGEO("GSE106253", GSEMatrix = TRUE, destdir = D)   # getGPL=TRUE 默认带注释
eset <- gse[[1]]
expr <- exprs(eset)
pd   <- pData(eset)
fd   <- fData(eset)
cat("维度:", nrow(expr), "探针 x", ncol(expr), "样本\n")

# ---------- 分组识别（characteristics 里找 OLF vs healthy/normal） ----------
char_cols <- grep("characteristics", colnames(pd), value = TRUE)
grp_raw <- NULL
for (cc in char_cols) {
  v <- as.character(pd[[cc]])
  if (any(grepl("ossif|OLF", v, ignore.case = TRUE)) &&
      any(grepl("normal|healthy|control", v, ignore.case = TRUE))) {
    grp_raw <- v; break
  }
}
if (is.null(grp_raw)) {   # 兜底：按 GSM 顺序前 4 后 4（GSE106253 = 4 OLF + 4 healthy）
  grp_raw <- rep(c("OLF patient", "healthy control"), each = 4)
  cat("[警告] characteristics 未识别分组，按前4=OLF/后4=健康对照处理，请人工核对！\n")
}
grp_raw <- gsub("^.*:\\s*", "", grp_raw)
case  <- grepl("ossif|OLF", grp_raw, ignore.case = TRUE)
ctrl  <- grepl("normal|healthy|control", grp_raw, ignore.case = TRUE)
stopifnot(sum(case) == 4, sum(ctrl) == 4)
grp <- ifelse(case, "OLF", "Control")
names(grp) <- rownames(pd)
cat("分组:"); print(table(grp))
# [复核] 对照打印的 pd 样本名与 GEO 页面 GSM 顺序核对分组正确性

# ---------- 探针 -> 基因符号 ----------
# GPL21827 平台在 GEO 无基因注释（2026-10-10 已核实），使用 11a 脚本产出的
# 定向序列比对注释（results/probe2gene_GPL21827_targeted.csv，全部 0 错配）
p2g_file <- file.path(R, "probe2gene_GPL21827_targeted.csv")
stopifnot(file.exists(p2g_file))   # 没有就先跑 scripts/11a_prep_annot_GPL21827.R
p2g <- fread(p2g_file, data.table = FALSE)
cat("使用定向比对注释:", nrow(p2g), "对探针->基因\n")

sym <- p2g$gene; names(sym) <- p2g$probe_id
ok  <- p2g$probe_id %in% rownames(expr)
stopifnot(all(ok))                                  # 44 探针须全部在矩阵中
expr <- expr[p2g$probe_id[ok], , drop = FALSE]      # 字符索引（勿用逻辑索引：会回收错位）
sym <- sym[ok]
# 同 symbol 多探针取方差最大者
v <- apply(expr, 1, var)
o <- order(v, decreasing = TRUE)
expr <- expr[o, , drop = FALSE]; sym <- sym[o]
keep <- !duplicated(sym)
mat <- expr[keep, , drop = FALSE]
rownames(mat) <- sym[keep]
stopifnot(!any(is.na(rownames(mat))))               # 行名必须全为基因符号

# log2 检测（Agilent 原始值为线性刻度）
if (max(mat, na.rm = TRUE) > 100) {
  mat <- log2(mat + 1)
  cat("数据为线性刻度，已 log2 转换\n")
}
mat[is.na(mat)] <- 0

# ---------- (0) 基因级差异表落盘（2026-10-13 复核补充） ----------
# 之前三基因 logFC 靠人工现算，现固化为 csv，供审稿复核直接引用
gl <- data.frame(
  gene = rownames(mat),
  mean_OLF = round(rowMeans(mat[, grp[colnames(mat)] == "OLF", drop = FALSE]), 4),
  mean_Control = round(rowMeans(mat[, grp[colnames(mat)] == "Control", drop = FALSE]), 4)
)
gl$logFC <- round(gl$mean_OLF - gl$mean_Control, 4)
gl$wilcox_p <- signif(apply(mat, 1, function(x)
  wilcox.test(x ~ factor(grp[colnames(mat)]))$p.value), 4)
gl <- gl[order(gl$wilcox_p), ]
fwrite(gl, file.path(R, "DEG_GSE106253_genelevel.csv"))
cat("基因级差异表已落盘: results/DEG_GSE106253_genelevel.csv (", nrow(gl), "个基因)\n")

# ---------- (a) LASSO 签名风险评分 ----------
lasso <- fread(file.path(R, "LASSO_signature_genes.csv"), data.table = FALSE)
avail <- lasso$Gene[lasso$Gene %in% rownames(mat)]
miss  <- setdiff(lasso$Gene, avail)
cat("\n签名基因可用:", paste(avail, collapse = ", "), "\n")
if (length(miss)) cat("签名基因缺失(平台无探针):", paste(miss, collapse = ", "), "\n")

res_tab <- NULL; roc_obj <- NULL
if (length(avail) >= 2) {
  beta <- setNames(lasso$beta[match(avail, lasso$Gene)], avail)
  # 与 10b 一致：不跨数据集 z-score，直接 log 表达 × 系数（方向由 beta 符号保证）
  score <- as.numeric(t(mat[avail, , drop = FALSE]) %*% beta)
  names(score) <- colnames(mat)
  roc_obj <- roc(ifelse(grp[colnames(mat)] == "OLF", 1, 0), score, quiet = TRUE)

  res_tab <- data.frame(sample = names(score),
                        group = grp[colnames(mat)],
                        risk_score = round(score, 4))
  fwrite(res_tab, file.path(R, "LASSO_external_GSE106253.csv"))

  auc_tab <- data.frame(
    Model = sprintf("LASSO signature (external, GSE106253 %dv%d, %d/%d genes)",
                    sum(grp == "OLF"), sum(grp == "Control"), length(avail), nrow(lasso)),
    AUC = roc_obj$auc,
    CI95 = sprintf("%.3f-%.3f", ci.auc(roc_obj, quiet = TRUE)[1],
                   ci.auc(roc_obj, quiet = TRUE)[3]))
  auc_tab$p_value <- wilcox.test(score ~ factor(grp[colnames(mat)]))$p.value
  fwrite(auc_tab, file.path(R, "LASSO_external_AUC_GSE106253.csv"))
  cat("=== GSE106253 LASSO ROC ===\n"); print(auc_tab)

  p <- ggplot(data.frame(FPR = 1 - roc_obj$specificities,
                         TPR = roc_obj$sensitivities),
              aes(FPR, TPR)) +
    geom_line(color = "#B2182B", linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = sprintf("LASSO signature, third cohort (GSE106253, 4v4)\nAUC = %.2f",
                         roc_obj$auc)) +
    coord_equal() + theme_bw(base_size = 12)
  ggsave(file.path(F, "Fig_ROC_LASSO_GSE106253.png"), p, width = 6, height = 6, dpi = 300)
} else {
  cat("[警告] 可用签名基因 < 2，跳过 ROC，仅做基因集打分\n")
}

# ---------- (b) 铁死亡基因集 ssGSEA 组间比较 ----------
ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- intersect(ferro, rownames(mat))
cat("\n铁死亡基因覆盖:", length(ferro), "/", 36, "\n")

# [方法学说明] GPL21827 无全平台基因注释（仅定向注释 42 基因），ssGSEA 需要
# 全转录组秩背景，故本队列改用透明替代：每基因跨样本 z-score 后取均值
# （ferroptosis gene-set aggregate score）。论文方法学中如实描述。
if (length(ferro) >= 10) {
  fz <- t(scale(t(mat[ferro, , drop = FALSE])))     # 行内 z-score
  agg <- as.numeric(colMeans(fz))
  ssdf <- data.frame(sample = colnames(mat),
                     group = grp[colnames(mat)],
                     ferroptosis_zmean = round(agg, 4))
  fwrite(ssdf, file.path(R, "PCD_scores_zmean_GSE106253.csv"))

  wt <- wilcox.test(ferroptosis_zmean ~ group, data = ssdf)
  cat(sprintf("铁死亡聚合评分(z-mean) OLF vs Control: Wilcoxon P=%.3g\n", wt$p.value))
  cat("均值: Control =", round(mean(ssdf$ferroptosis_zmean[ssdf$group == "Control"]), 3),
      " OLF =", round(mean(ssdf$ferroptosis_zmean[ssdf$group == "OLF"]), 3), "\n")
  p2 <- ggplot(ssdf, aes(group, ferroptosis_zmean, fill = group)) +
    geom_boxplot(width = 0.5, outlier.size = 2) +
    geom_jitter(width = 0.08, size = 2.5, alpha = 0.7) +
    labs(y = "Ferroptosis gene-set score (z-mean)", x = NULL,
         subtitle = sprintf("Wilcoxon P = %.2g", wt$p.value)) +
    theme_bw(base_size = 12)
  ggsave(file.path(F, "Fig_ferroptosis_zmean_GSE106253.png"), p2,
         width = 5, height = 5, dpi = 300)
  # [复核] 本队列为健康人对照 + lncRNA 芯片（mRNA 单探针），基因集层面方向
  # 若与发现队列不一致，属平台/对照异质性，论文中如实报告并讨论
}

cat("\n=== 11 done ===\n")

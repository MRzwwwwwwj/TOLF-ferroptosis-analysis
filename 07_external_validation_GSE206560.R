# =====================================================================
# 07_external_validation_GSE206560.R  (2026-10-09 按实际文件结构重写)
# GSE206560 (m6A MeRIP-seq, n=12) 骨化 vs 非骨化 外部验证
#   NCBI counts 实际仅含 5 个 Input 样本(非骨化2 vs 骨化3, GSM6256886缺失):
#   差异表达 + PCD 评分 + 免疫浸润(ssGSEA) + 与主数据集一致性
# 依赖: edgeR, limma, GSVA, data.table, ggplot2
# 运行(PowerShell):
#   & "D:\ARLARS\R-4.5.3\bin\Rscript.exe" scripts/07_external_validation_GSE206560.R
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(GSVA); library(ggplot2)
})

setwd("C:/Users/10215/WorkBuddy/2026-10-03-00-10-07")
DATA_DIR <- "GSE206560"; OUT_DIR <- "results"; FIG_DIR <- "figures"

# ---------------- 1. 样本元数据 ----------------
sm <- readLines(gzfile(file.path(DATA_DIR, "GSE206560_series_matrix.txt.gz")))
getrow <- function(tag) {
  l <- sm[grepl(tag, sm, fixed = TRUE)][1]
  strsplit(gsub('"', "", l), "\t", fixed = TRUE)[[1]][-1]
}
ttl <- getrow("!Sample_title")
geo <- getrow("!Sample_geo_accession")
abl <- sm[grepl("!Sample_characteristics_ch1", sm, fixed = TRUE) & grepl("antibody", sm)][1]
ab  <- sub("antibody:\\s*", "", strsplit(gsub('"', "", abl), "\t", fixed = TRUE)[[1]][-1])
meta <- data.frame(geo = geo, title = ttl, antibody = ab, stringsAsFactors = FALSE)
meta$assay <- ifelse(grepl("IP", meta$title, fixed = TRUE), "IP", "Input")
meta$oss   <- ifelse(grepl("non.?ossified", meta$title, ignore.case = TRUE),
                     "NonOssified", "Ossified")
cat("=== GSE206560 12 样本元数据 ===\n"); print(meta)

# ---------------- 2. counts + 注释 → 基因矩阵 ----------------
cnt <- fread(file.path(DATA_DIR, "GSE206560_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
             header = TRUE, data.table = FALSE)
annot <- fread(file.path(DATA_DIR, "Human.GRCh38.p13.annot.tsv.gz"),
               header = TRUE, data.table = FALSE)
cat("counts 文件样本:", paste(colnames(cnt)[-1], collapse = ","), "\n")

map_sym <- function(mat, annot) {
  m <- match(mat[[1]], annot$GeneID)
  sym <- annot$Symbol[m]
  ok <- !is.na(sym) & sym != ""
  mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(rowSums(as.matrix(mat[, -1, drop = FALSE])), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  keep <- !duplicated(sym)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- sym[keep]
  mat
}
cnt <- map_sym(cnt, annot)
cat("counts:", nrow(cnt), "genes x", ncol(cnt) - 1, "samples\n")

sid  <- intersect(meta$geo[meta$assay == "Input"], colnames(cnt)[-1])
grp  <- factor(setNames(meta$oss, meta$geo)[sid], levels = c("NonOssified", "Ossified"))
cnt_in <- cnt[, sid, drop = FALSE]
cat("Input 验证样本:", paste(sid, grp, sep = "=", collapse = ", "), "\n")
if (length(sid) < 4) stop("Input 样本不足 4 个, 中止")

# ---------------- 3. 差异分析: edgeR-TMM + voom + limma ----------------
mat <- as.matrix(cnt_in); storage.mode(mat) <- "integer"
dge <- DGEList(counts = mat, group = grp)
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ grp)
v <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(v, design))
deg <- topTable(fit, coef = "grpOssified", number = Inf, sort.by = "P")
deg$gene <- rownames(deg)
fwrite(deg, file.path(OUT_DIR, "DEG_GSE206560_full.csv"))
sig <- subset(deg, adj.P.Val < 0.05)
cat("GSE206560 显著基因(adj.P<0.05):", nrow(sig), "\n")
fwrite(sig, file.path(OUT_DIR, "DEG_GSE206560_sig.csv"))

# ---------------- 4. 与主数据集(GSE113212)方向一致性 ----------------
main_deg <- fread(file.path(OUT_DIR, "DEG_limma_full.csv"), data.table = FALSE)
names(main_deg)[names(main_deg) == "Gene"] <- "gene"   # 主表列名为大写 Gene
common <- intersect(deg$gene, main_deg$gene)
d1 <- deg[match(common, deg$gene), ]
d2 <- main_deg[match(common, main_deg$gene), ]
ct <- suppressWarnings(cor(d1$logFC, d2$logFC, use = "complete.obs"))
ct_p <- suppressWarnings(cor.test(d1$logFC, d2$logFC, use = "complete.obs")$p.value)
cat(sprintf("\n两数据集 logFC 一致性: Pearson r = %.3f (P = %.2e, n = %d)\n", ct, ct_p, length(common)))
keyg <- c("TGFB1", "BMPR1B", "COL13A1", "COL11A1", "STEAP3", "HMOX1", "CDKN1A")
tab <- deg[deg$gene %in% keyg, c("gene", "logFC", "P.Value", "adj.P.Val")]
cat("=== 关键基因在 GSE206560 的表现 ===\n"); print(tab, row.names = FALSE)
fwrite(tab, file.path(OUT_DIR, "keygenes_GSE206560.csv"))
writeLines(sprintf("pearson_r\t%.4f\tp_value\t%.3e\tn_common\t%d", ct, ct_p, length(common)),
           file.path(OUT_DIR, "consistency_GSE113212_GSE206560.txt"))

# ---------------- 5. PCD 评分 (TPM log2+1, GSVA + ssGSEA) ----------------
tpm <- fread(file.path(DATA_DIR, "GSE206560_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
             header = TRUE, data.table = FALSE)
tpm <- map_sym(tpm, annot)
tpm_in <- log2(tpm[, sid, drop = FALSE] + 1)
cat("TPM 矩阵:", nrow(tpm_in), "genes x", ncol(tpm_in), "samples\n")

gs <- list()
for (f in list.files("pcd_genesets", pattern = "^gs_.*\\.txt$", full.names = TRUE)) {
  nm <- sub("\\.txt$", "", sub("^gs_", "", basename(f)))  # sub 只替换首个匹配, 需两步
  g <- trimws(readLines(f)); g <- intersect(g[nzchar(g)], rownames(tpm_in))
  gs[[nm]] <- g
}
gs <- gs[sapply(gs, length) >= 5]
cat("PCD 基因集覆盖:", paste(sapply(gs, length), names(gs), sep = "/", collapse = ", "), "\n")

pcd <- gsva(gsvaParam(as.matrix(tpm_in), gs, kcdf = "Gaussian"), verbose = FALSE)
pcd_ss <- gsva(ssgseaParam(as.matrix(tpm_in), gs), verbose = FALSE)
fwrite(data.table(PCD = rownames(pcd), as.data.table(pcd)), file.path(OUT_DIR, "PCD_scores_GSVA_GSE206560.csv"))
fwrite(data.table(PCD = rownames(pcd_ss), as.data.table(pcd_ss)), file.path(OUT_DIR, "PCD_scores_ssGSEA_GSE206560.csv"))

cmp <- do.call(rbind, lapply(rownames(pcd), function(g) {
  x <- pcd[g, grp == "Ossified"]; y <- pcd[g, grp == "NonOssified"]
  data.frame(PCD = g, mean_Ossified = round(mean(x), 3), mean_NonOssified = round(mean(y), 3),
             direction = ifelse(mean(x) > mean(y), "Up", "Down"),
             P.wilcox = tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA))
}))
cat("\n=== PCD 组间比较 (GSE206560, 2v3 不平衡, 仅报方向) ===\n"); print(cmp, row.names = FALSE)
fwrite(cmp, file.path(OUT_DIR, "PCD_groupcompare_GSE206560.csv"))

df <- data.frame(PCD = rep(rownames(pcd), ncol(pcd)),
                 score = as.vector(pcd),
                 group = rep(as.character(grp), each = nrow(pcd)))
p <- ggplot(df, aes(PCD, score, fill = group)) +
  geom_boxplot(outlier.size = 1, width = 0.6) +
  geom_jitter(aes(color = group), width = 0.12, size = 2, show.legend = FALSE) +
  scale_fill_manual(values = c(NonOssified = "#4C72B0", Ossified = "#C44E52")) +
  scale_color_manual(values = c(NonOssified = "#2a4a7f", Ossified = "#8f3033")) +
  labs(x = NULL, y = "GSVA score", fill = NULL,
       title = "PCD scores in GSE206560 (Input samples, n=2 vs 3)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "Fig_PCD_scores_GSE206560.png"), p, width = 8, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig_PCD_scores_GSE206560.pdf"), p, width = 8, height = 5)

# ---------------- 6. 免疫浸润 (ssGSEA-Bindea) ----------------
sigb <- read.csv("immune_signatures_bindea.csv", stringsAsFactors = FALSE)
sigb <- sigb[!sigb$celltype %in% c("Blood_vessels", "Lymph_vessels", "Normal_mucosa", "SW480_cancer_cells"), ]
sig28 <- lapply(split(sigb$gene, sigb$celltype), function(g) intersect(na.omit(g), rownames(tpm_in)))
sig28 <- sig28[sapply(sig28, length) >= 3]
imm <- gsva(ssgseaParam(as.matrix(tpm_in), sig28), verbose = FALSE)
fwrite(data.table(cell = rownames(imm), as.data.table(imm)), file.path(OUT_DIR, "immune_ssgsea_GSE206560.csv"))

# ---------------- 7. PCD × 免疫 相关性 ----------------
cor_mat <- cor(t(pcd), t(imm), method = "spearman")
pv <- outer(rownames(pcd), rownames(imm),
            Vectorize(function(a, b) tryCatch(cor.test(pcd[a, ], imm[b, ], method = "spearman")$p.value,
                                              error = function(e) NA)))
dimnames(pv) <- list(rownames(pcd), rownames(imm))
fwrite(data.table(PCD = rep(rownames(pcd), ncol(cor_mat)),
                  immune = rep(colnames(cor_mat), each = nrow(cor_mat)),
                  rho = round(as.vector(cor_mat), 3), P = round(as.vector(pv), 4)),
       file.path(OUT_DIR, "PCD_immune_cor_GSE206560.csv"))
core <- grep("Macro|Monocyte", colnames(cor_mat), value = TRUE)
fer_row <- grep("^ferroptosis$", rownames(cor_mat), value = TRUE, ignore.case = TRUE)
core_tab <- data.frame(immune = core,
                       rho = round(cor_mat[fer_row, core], 3),
                       P = round(pv[fer_row, core], 4))
cat("\n=== 核心发现验证: Ferroptosis × 巨噬细胞 (GSE206560, n=5) ===\n"); print(core_tab, row.names = FALSE)
fwrite(core_tab, file.path(OUT_DIR, "corefinding_validation_GSE206560.csv"))

cat("\n[07] GSE206560 外部验证完成。结果在 results/, 图在 figures/。\n")

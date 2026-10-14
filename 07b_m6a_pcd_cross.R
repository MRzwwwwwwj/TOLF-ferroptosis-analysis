# =====================================================================
# 07b_m6a_pcd_cross.R
# GSE206560 m6A MeRIP-seq × PCD 交叉分析
#   IP 样本 peak → 基因级 m6A 状态 → 骨化组高/低甲基化基因
#   → PCD 基因集富集(Fisher) → 与主数据集/验证集表达方向联动
# 前置: scripts/07b_prep_m6a.py 已生成 results/m6a_sites_all.csv
# 运行: Rterm.exe --no-restore --file=scripts/07b_m6a_pcd_cross.R
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

setwd("C:/Users/10215/WorkBuddy/2026-10-03-00-10-07")
OUT_DIR <- "results"; FIG_DIR <- "figures"

# ---------------- 1. 读入位点表, 统计各样本 peak 数 ----------------
sites <- fread(file.path(OUT_DIR, "m6a_sites_all.csv"))
samples <- c(A01 = "NonOssified", A02 = "NonOssified", A03 = "NonOssified",
             B01 = "Ossified", B02 = "Ossified", B03 = "Ossified")
cnt <- sites[, .(n_sites = .N), by = .(sample, group)]
cnt$group <- factor(cnt$group, levels = c("NonOssified", "Ossified"))
cat("=== 各样本 m6A 位点数 ===\n"); print(cnt)
w <- tryCatch(wilcox.test(n_sites ~ group, cnt)$p.value, error = function(e) NA)
cat(sprintf("骨化 vs 非骨化 peak 数 Wilcoxon P = %.2f (3v3, 最小可达 0.1)\n", w))

p1 <- ggplot(cnt, aes(sample, n_sites, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_sites), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(NonOssified = "#4C72B0", Ossified = "#C44E52")) +
  labs(x = NULL, y = "Number of m6A peak sites", fill = NULL,
       title = "m6A sites per sample (GSE206560, IP)") +
  theme_bw(base_size = 12)
ggsave(file.path(FIG_DIR, "Fig_m6a_peakcounts_GSE206560.png"), p1, width = 7, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig_m6a_peakcounts_GSE206560.pdf"), p1, width = 7, height = 5)
fwrite(cnt, file.path(OUT_DIR, "m6a_peakcounts_GSE206560.csv"))

# ---------------- 2. 基因级 m6A 状态矩阵 ----------------
gm <- dcast(sites, gene ~ sample, fun.aggregate = length, value.var = "peakID")
gmat <- as.matrix(gm[, -1]); rownames(gmat) <- gm$gene
storage.mode(gmat) <- "integer"
oss   <- c("B01", "B02", "B03"); non <- c("A01", "A02", "A03")
present <- gmat > 0
n_oss <- rowSums(present[, oss]); n_non <- rowSums(present[, non])
bg <- rownames(gmat)                       # 背景: 任一样本检出 m6A 的基因
cat("m6A 检出基因数(背景):", length(bg), "\n")

hyper <- bg[n_oss >= 2 & n_non == 0]       # 骨化组稳定甲基化, 非骨化组全无
hypo  <- bg[n_non >= 2 & n_oss == 0]       # 反向
both  <- bg[n_oss >= 2 & n_non >= 2]       # 两组共有
cat("骨化组高甲基化(hyper):", length(hyper), "| 非骨化组特异(hypo):", length(hypo),
    "| 两组共有:", length(both), "\n")
fwrite(data.table(gene = bg, n_ossified = n_oss[bg], n_nonossified = n_non[bg],
                  category = ifelse(bg %in% hyper, "Ossified_hyper",
                             ifelse(bg %in% hypo, "NonOssified_specific", "Common"))),
       file.path(OUT_DIR, "m6a_gene_status_GSE206560.csv"))

# ---------------- 3. PCD 基因集 m6A 富集 (Fisher) ----------------
gs <- list()
for (f in list.files("pcd_genesets", pattern = "^gs_.*\\.txt$", full.names = TRUE)) {
  nm <- sub("\\.txt$", "", sub("^gs_", "", basename(f)))
  gs[[nm]] <- trimws(readLines(f))
}
enr <- do.call(rbind, lapply(names(gs), function(nm) {
  gl <- intersect(gs[[nm]], bg)                       # 在背景内可检的 PCD 基因
  hit <- intersect(gl, hyper)
  a <- length(hit); b <- length(gl) - a
  c <- length(hyper) - a; d <- length(bg) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")
  data.frame(PCD = nm, pcd_in_bg = length(gl), hyper_hits = a,
             expected = round(length(gl) * length(hyper) / length(bg), 1),
             odds_ratio = round(unname(ft$estimate), 2), P = signif(ft$p.value, 3))
}))
cat("\n=== PCD 基因集在骨化组高甲基化基因中的富集 ===\n"); print(enr, row.names = FALSE)
fwrite(enr, file.path(OUT_DIR, "m6a_PCD_enrichment_GSE206560.csv"))

# ---------------- 4. 具体命中基因表 ----------------
hits <- do.call(rbind, lapply(names(gs), function(nm) {
  hit <- intersect(intersect(gs[[nm]], bg), hyper)
  if (length(hit)) data.frame(PCD = nm, gene = hit) else NULL
}))
if (!is.null(hits)) {
  cat("\n=== 骨化组高甲基化的 PCD 基因 ===\n"); print(hits, row.names = FALSE)
  fwrite(hits, file.path(OUT_DIR, "m6a_PCD_hyper_genes_GSE206560.csv"))
}

# ---------------- 5. 与表达方向联动(主数据集 + 验证集 DE) ----------------
main_deg <- fread(file.path(OUT_DIR, "DEG_limma_full.csv"), data.table = FALSE)
names(main_deg)[names(main_deg) == "Gene"] <- "gene"
val_deg <- fread(file.path(OUT_DIR, "DEG_GSE206560_full.csv"), data.table = FALSE)
val_deg <- val_deg[!duplicated(val_deg$gene), ]
lk <- function(genes) {
  m <- main_deg[match(genes, main_deg$gene), c("logFC", "adj.P.Val")]
  v <- val_deg[match(genes, val_deg$gene), c("logFC", "adj.P.Val")]
  data.frame(gene = genes,
             main_logFC = round(m$logFC, 2), main_adjP = signif(m$adj.P.Val, 2),
             val_logFC = round(v$logFC, 2), val_adjP = signif(v$adj.P.Val, 2))
}
if (!is.null(hits)) {
  cat("\n=== 高甲基化 PCD 基因的表达联动(主=GSE113212 老年肥厚↑; 验证=GSE206560 骨化↑) ===\n")
  print(lk(hits$gene), row.names = FALSE)
  fwrite(lk(hits$gene), file.path(OUT_DIR, "m6a_PCD_hyper_exprlink.csv"))
}

# ---------------- 6. m6A 调控酶在验证集的表达 ----------------
regs <- c("METTL3","METTL14","METTL16","WTAP","VIRMA","ZC3H13","RBM15","RBM15B",
          "CBLL1","KIAA1429","ALKBH5","FTO","YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3",
          "IGF2BP1","IGF2BP2","IGF2BP3","HNRNPA2B1","RBMX","ELAVL1","FMR1","LRPPRC")
rt <- lk(intersect(regs, val_deg$gene))
rt$m6a_status <- ifelse(rt$gene %in% hyper, "Ossified_hyper",
                 ifelse(rt$gene %in% bg, "methylated", "not_methylated"))
cat("\n=== m6A 调控酶在 GSE206560 的表达与自身 m6A 状态 ===\n"); print(rt, row.names = FALSE)
fwrite(rt, file.path(OUT_DIR, "m6a_regulators_GSE206560.csv"))

# ---------------- 7. 关键基因 m6A 状态 ----------------
keyg <- c("TGFB1","BMPR1B","COL13A1","COL11A1","STEAP3","HMOX1","CDKN1A","BMP2","RUNX2")
kt <- data.frame(gene = keyg,
                 n_oss = n_oss[keyg], n_non = n_non[keyg],
                 category = ifelse(keyg %in% hyper, "Ossified_hyper",
                            ifelse(keyg %in% bg, "methylated", "not_methylated")))
kt <- cbind(kt, lk(keyg)[, -1])
cat("\n=== 关键基因的 m6A 状态与表达 ===\n"); print(kt, row.names = FALSE)
fwrite(kt, file.path(OUT_DIR, "keygenes_m6a_status_GSE206560.csv"))

cat("\n[07b] m6A × PCD 交叉分析完成。\n")

#!/usr/bin/env Rscript
# ============================================================================
# 11a_prep_annot_GPL21827.R — GSE106253 平台探针→基因 定向序列比对注释
#
# 背景：GPL21827 (Arraystar lncRNA V4) 在 GEO 无基因符号注释（原数据集作者也是
#       自己比对序列注释的，见 Kong 2019 方法）。本脚本按 MDC Berlin 团队发表的
#       同平台再注释流程（补充材料）实现：
#       探针 60mer 序列 → 比对目标基因的转录本序列（Biostrings，容错 0-4 错配）
#
# 转录本来源：NCBI RefSeq mRNA（rentrez E-utilities；2026-10-10 实测 Ensembl
#            全部镜像在国内不可达，NCBI 可达）
#
# 只注释分析需要的基因（铁死亡 36 + LASSO 4 + 内参 2 + 相关基因），靶向比对，
# 计算量小（61k 探针 × ~45 基因转录本，分钟级）。
#
# 输入：GSE106253/GPL21827.soft.gz（纯文本 SOFT，getGEO 已下载）
# 输出：results/probe2gene_GPL21827_targeted.csv（probe_id, gene, mismatches）
# ============================================================================
suppressPackageStartupMessages({
  library(Biostrings); library(data.table); library(rentrez)
})

BASE <- getwd()
D <- file.path(BASE, "GSE106253")
R <- file.path(BASE, "results")

## ---------- 目标基因 ----------
ferro <- trimws(readLines(file.path(BASE, "pcd_genesets", "gs_ferroptosis.txt")))
ferro <- ferro[nzchar(ferro)]
lasso <- fread(file.path(R, "LASSO_signature_genes.csv"), data.table = FALSE)$Gene
target_genes <- sort(unique(c(ferro, lasso, "ACTB", "GAPDH", "SPP1", "TGFB1")))
cat("目标基因数:", length(target_genes), "\n")

## ---------- 1. 读取平台探针表 ----------
cat("读取 GPL21827 平台表...\n")
ln <- readLines(file.path(D, "GPL21827.soft.gz"), warn = FALSE)
i0 <- grep("^!platform_table_begin", ln); i1 <- grep("^!platform_table_end", ln)
stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0)
plat <- fread(text = paste(ln[(i0 + 1):(i1 - 1)], collapse = "\n"),
              sep = "\t", header = TRUE, data.table = FALSE, quote = "")
cat("平台探针数:", nrow(plat), "列:", paste(colnames(plat), collapse = ","), "\n")
plat <- plat[, c("ID", "SEQUENCE", "TRANSCRIPT_TYPE")]
plat <- plat[plat$SEQUENCE != "" & !is.na(plat$SEQUENCE), ]
# 只留纯 ACGT 的 60mer（PDict 要求）
is_acgt <- grepl("^[ACGT]+$", plat$SEQUENCE) & nchar(plat$SEQUENCE) == 60
plat <- plat[is_acgt, ]
cat("可比对探针:", nrow(plat), "\n")

## ---------- 2. 获取目标基因 RefSeq mRNA 序列（NCBI rentrez，带缓存） ----------
cache <- file.path(D, "target_tx_refseq.csv")
if (file.exists(cache)) {
  tx <- fread(cache, data.table = FALSE)
  cat("使用缓存转录本:", nrow(tx), "条\n")
} else {
  cat("经 NCBI E-utilities 获取 RefSeq mRNA...\n")
  tx <- data.frame(hgnc_symbol = character(), cdna = character(), acc = character())
  for (g in target_genes) {
    term <- sprintf('%s[Gene Name] AND "Homo sapiens"[Organism] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]', g)
    ids <- try(entrez_search(db = "nucleotide", term = term, retmax = 60), silent = TRUE)
    if (inherits(ids, "try-error") || length(ids$ids) == 0) {
      cat(sprintf("  %-10s: 无 RefSeq 记录\n", g)); next
    }
    fa <- try(entrez_fetch(db = "nucleotide", id = ids$ids, rettype = "fasta"), silent = TRUE)
    if (inherits(fa, "try-error")) { cat(sprintf("  %-10s: efetch 失败\n", g)); next }
    # 解析 FASTA：只保留 NM_（已审查转录本），排除 XM_ 预测版本
    recs <- strsplit(sub("^>", "", fa), ">")[[1]]
    for (rec in recs) {
      lns <- strsplit(rec, "\n")[[1]]
      hd <- lns[1]
      if (!grepl("^NM_", hd)) next
      seq <- paste(lns[-1], collapse = "")
    seq <- gsub("[^ACGTNacgtn]", "", seq)
      if (nchar(seq) >= 60) {
        acc <- sub("^ +(NM_[0-9]+\\.[0-9]+).*", "\\1", hd)
        tx <- rbind(tx, data.frame(hgnc_symbol = g, cdna = toupper(seq), acc = acc))
      }
    }
    cat(sprintf("  %-10s: %d 条 NM_ 转录本\n", g,
                sum(tx$hgnc_symbol == g)))
    Sys.sleep(0.4)   # NCBI 限速礼貌间隔
  }
  stopifnot(nrow(tx) > 0)
  fwrite(tx, cache)
}

genes_found <- unique(tx$hgnc_symbol)
cat("有转录本序列的基因:", length(genes_found), "/", length(target_genes), "\n")
miss <- setdiff(target_genes, genes_found)
if (length(miss)) cat("无序列基因:", paste(miss, collapse = ", "), "\n")

## ---------- 3. 定向比对（0-4 错配逐级放宽） ----------
cat("构建探针字典...\n")
pdict <- PDict(plat$SEQUENCE, tb.start = 1, tb.end = 40)   # 信任带支持错配

hits <- list()
for (g in genes_found) {
  seqs <- tx$cdna[tx$hgnc_symbol == g]
  sub <- DNAStringSet(seqs)
  best_mm <- NA
  for (mm in 0:4) {
    cnt <- vcountPDict(pdict, sub, max.mismatch = mm)
    tot <- rowSums(cnt)
    if (any(tot > 0)) { best_mm <- mm; break }
  }
  if (is.na(best_mm)) next
  tot <- rowSums(vcountPDict(pdict, sub, max.mismatch = best_mm))
  idx <- which(tot > 0)
  hits[[g]] <- data.frame(probe_id = plat$ID[idx], gene = g, mismatches = best_mm)
  cat(sprintf("  %-10s: %d 探针命中 (mm<=%d)\n", g, length(idx), best_mm))
}

hit_tab <- do.call(rbind, hits)
# 跨基因歧义探针剔除（同一探针命中多个基因）
amb <- hit_tab$probe_id[duplicated(hit_tab$probe_id) |
                          hit_tab$probe_id %in% hit_tab$probe_id[duplicated(hit_tab$probe_id)]]
hit_tab <- hit_tab[!(hit_tab$probe_id %in% amb), ]
cat("\n剔除歧义探针:", length(unique(amb)), "\n")
cat("最终 注释探针对:", nrow(hit_tab), "覆盖基因:", length(unique(hit_tab$gene)), "\n")

fwrite(hit_tab, file.path(R, "probe2gene_GPL21827_targeted.csv"))

# 覆盖情况速览
key <- intersect(c(lasso, "SPP1"), unique(hit_tab$gene))
cat("\nLASSO/SPP1 覆盖:", paste(key, collapse = ", "), "\n")
key_miss <- setdiff(c(lasso, "SPP1"), unique(hit_tab$gene))
if (length(key_miss)) cat("LASSO/SPP1 未覆盖:", paste(key_miss, collapse = ", "), "\n")
cat("\n=== 11a done ===\n")

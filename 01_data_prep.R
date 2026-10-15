# ============================================================
# TOLF Step 01: GSE113212 数据整理
# 输入: GSE113212_series_matrix.txt.gz (已下载, log2 标准化)
#       GSE113212_RAW/ (Agilent FE txt, 提供探针->基因映射)
# 输出: results/expr_gene_GSE113212.csv  (基因 x 样本, 已按基因聚合)
#       results/group_GSE113212.csv      (分组: young_control vs elderly_hypertrophied)
#       results/probe2gene_GPL17077.csv
# ============================================================

.libPaths("D:/ARLARS/R-4.5.3/library")
suppressPackageStartupMessages(library(limma))

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
dir.create(file.path(wd, "results"), showWarnings = FALSE)

sm <- file.path(wd, "GSE113212_series_matrix.txt.gz")
con <- gzfile(sm, "rt")
lines <- readLines(con)
close(con)

# ---------- 1. 样本元数据 ----------
get_meta <- function(key) {
  ln <- lines[grepl(paste0("^!", key), lines)][1]
  vals <- regmatches(ln, gregexpr('"[^"]*"', ln))[[1]]
  gsub('"', "", vals)
}
titles <- get_meta("Sample_title")
gsm    <- get_meta("Sample_geo_accession")
cat("样本:\n"); print(data.frame(gsm, title = titles))

# ---------- 2. 表达矩阵 ----------
beg <- which(lines == "!series_matrix_table_begin")
end <- which(lines == "!series_matrix_table_end")
tab <- read.delim(text = lines[(beg + 1):(end - 1)], check.names = FALSE)
probe_ids <- tab[[1]]
expr <- as.matrix(tab[, -1])
rownames(expr) <- probe_ids
storage.mode(expr) <- "double"
colnames(expr) <- gsm
cat("探针数:", nrow(expr), " 样本数:", ncol(expr), "\n")

# ---------- 3. 探针 -> 基因注释 ----------
# 注: Agilent FE 原始文件不含 GeneName 列, 故用 GPL17077 平台文本表
gpl_tab <- file.path(wd, "GPL17077_table.txt")
stopifnot("GPL17077_table.txt 未下载完成" = file.exists(gpl_tab))
pl_lines <- readLines(gpl_tab)
pb <- which(pl_lines == "!platform_table_begin")
pe <- which(pl_lines == "!platform_table_end")
pl <- read.delim(text = pl_lines[(pb + 1):(pe - 1)], check.names = FALSE,
                 quote = "", stringsAsFactors = FALSE)
cat("平台表列名:", paste(colnames(pl), collapse = ", "), "\n")
sym_col <- intersect(c("GENE_SYMBOL", "GENE", "Gene Symbol", "GB_LIST", "GENE_NAME"),
                     colnames(pl))[1]
ctrl_col <- intersect(c("CONTROL_TYPE", "ControlType"), colnames(pl))[1]
stopifnot("平台表中未找到基因符号列" = !is.na(sym_col))
annot <- data.frame(
  ProbeID = pl$ID,
  GeneSymbol = pl[[sym_col]],
  ControlType = if (!is.na(ctrl_col)) pl[[ctrl_col]] else 0,
  stringsAsFactors = FALSE
)
write.csv(annot, file.path(wd, "results", "probe2gene_GPL17077.csv"), row.names = FALSE)
cat("注释探针数:", nrow(annot), "\n")

# ---------- 4. 过滤 + 基因级聚合 ----------
annot$GeneSymbol <- trimws(annot$GeneSymbol)
ctl <- tolower(as.character(annot$ControlType))
keep <- annot$GeneSymbol != "" & !is.na(annot$GeneSymbol) & ctl %in% c("0", "false", "")
cat("对照探针类型取值:", paste(unique(ctl), collapse = ", "), "\n")
annot_k <- annot[keep, ]
# 一探针对多基因(以 /// 分隔)的取第一个, 并在文中说明
annot_k$GeneSymbol <- sub("///.*", "", annot_k$GeneSymbol)

idx <- match(annot_k$ProbeID, rownames(expr))
idx <- idx[!is.na(idx)]
expr_k <- expr[idx, ]
rownames(expr_k) <- annot_k$GeneSymbol[match(rownames(expr_k), annot_k$ProbeID)]

# 多探针同基因 -> limma::avereps 取均值
expr_g <- avereps(expr_k, ID = rownames(expr_k))
cat("基因级矩阵:", nrow(expr_g), "基因 x", ncol(expr_g), "样本\n")

# ---------- 5. 分组 ----------
grp <- ifelse(grepl("young|Young", titles), "young_control",
       ifelse(grepl("elderly|Elderly|old", titles), "elderly_hypertrophied", NA))
if (any(is.na(grp))) {
  cat("!! 标题未含 young/elderly 关键词, 打印标题人工核对:\n"); print(titles)
}
group <- data.frame(sample = gsm, title = titles, group = grp)

# ---------- 6. 输出 ----------
write.csv(data.frame(Gene = rownames(expr_g), expr_g, check.names = FALSE),
          file.path(wd, "results", "expr_gene_GSE113212.csv"), row.names = FALSE)
write.csv(group, file.path(wd, "results", "group_GSE113212.csv"), row.names = FALSE)

cat("\n=== 完成 ===\n")
print(table(group$group))
cat("表达范围:", range(expr_g), "\n")

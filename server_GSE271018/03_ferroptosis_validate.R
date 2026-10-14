#!/usr/bin/env Rscript
# 03_ferroptosis_validate.R — GSE271018 第二重单细胞验证（主分析）
#
# 设计：逐样本处理（不合并大对象），每个样本打分后立即落盘并释放内存。
#       兼容 Seurat v3 / v4 / v5（自动分支），不依赖 harmony / patchwork / dplyr。
#
# 流程：
#   每个 RDS -> 找注释列（优先用原文现成注释；无则标志基因面板兜底）
#            -> AddModuleScore 铁死亡打分（36 基因集）
#            -> 提取关键基因表达（SLC7A11/STEAP3/HMOX1/CDKN1A/SPP1/CRTAC1）
#            -> SPP1+/- 巨噬细胞标记 -> 每细胞一行写盘 -> 释放内存
#   汇总：细胞类型排序 / 退变 vs 对照 Wilcoxon / SPP1+/- 比较 / 4 张图
#
# 用法：Rscript 03_ferroptosis_validate.R > log03.txt 2>&1
# 长任务建议：nohup Rscript 03_ferroptosis_validate.R > log03.txt 2>&1 &
#
# ★ 运行前先跑 02_probe_rds.R：若 probe 报告里有注释列名，填到下面
#   CELLTYPE_COL_OVERRIDE（留空 = 自动识别）。

CELLTYPE_COL_OVERRIDE <- ""     # 例如 "cell_type" 或 "seurat_clusters"
DATA_DIR  <- "GSE271018"
GENE_FILE <- "gs_ferroptosis.txt"
OUT_SCORE <- "GSE271018_scores"
OUT_RES   <- "GSE271018_results"
dir.create(OUT_SCORE, showWarnings = FALSE)
dir.create(OUT_RES,   showWarnings = FALSE)

suppressMessages({ library(Seurat); library(ggplot2) })

## ---------- 分组映射（按文件名；请对照 series_matrix 核实） ----------
# Tendon_C1/C3/C5 -> Control；Tendon_N0/N1/N2/N3/N5/N10 -> Degenerative
infer_group <- function(fname) {
  if (grepl("Tendon_C", fname, ignore.case = TRUE)) return("Control")
  if (grepl("Tendon_N", fname, ignore.case = TRUE)) return("Degenerative")
  return(NA)
}

## ---------- 版本兼容辅助 ----------
get_counts <- function(obj) {
  m <- tryCatch(GetAssayData(obj, layer = "counts"),   # Seurat v5
                error = function(e) NULL)
  if (is.null(m))
    m <- GetAssayData(obj, slot = "counts")            # Seurat v3/v4
  m
}
get_data <- function(obj) {
  m <- tryCatch(GetAssayData(obj, layer = "data"),
                error = function(e) NULL)
  if (is.null(m))
    m <- GetAssayData(obj, slot = "data")
  m
}

## ---------- 标志基因面板（无现成注释时兜底） ----------
marker_panel <- list(
  Fibroblast     = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  Macrophage     = c("CD68", "CD163", "LYZ", "C1QA", "C1QB", "CD14"),
  Endothelial    = c("PECAM1", "VWF", "KDR", "CLDN5"),
  Pericyte_SMC   = c("ACTA2", "RGS5", "MYH11", "TAGLN"),
  T_NK           = c("CD3D", "CD3E", "NKG7", "GNLY", "IL7R"),
  B_cell         = c("CD79A", "MS4A1", "CD19"),
  Mast_cell      = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
  Chondrocyte    = c("ACAN", "SOX9", "COL2A1", "CRTAC1")
)

## ---------- 自动识别注释列 ----------
pick_celltype_col <- function(md) {
  if (nzchar(CELLTYPE_COL_OVERRIDE) && CELLTYPE_COL_OVERRIDE %in% colnames(md))
    return(CELLTYPE_COL_OVERRIDE)
  cand <- grep("cell.?type|annotation|celltype|sub.?type|major|label",
               colnames(md), ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, c("nCount_RNA", "nFeature_RNA"))
  for (cn in cand) {
    v <- as.character(md[[cn]])
    lv <- unique(v)
    if (length(lv) >= 2 && length(lv) <= 40 &&
        any(grepl("macro|fibro|mono", lv, ignore.case = TRUE)))
      return(cn)
  }
  # 次级候选：任何含 macro/fibro 的低基数字符列
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.character(v) || is.factor(v)) {
      lv <- unique(as.character(v))
      if (length(lv) >= 2 && length(lv) <= 40 &&
          any(grepl("macro|fibro", lv, ignore.case = TRUE)))
        return(cn)
    }
  }
  return(NA)
}

## ---------- 无注释兜底：聚类 + 标志基因面板打分 ----------
annotate_by_markers <- function(obj) {
  message("    no annotation column -> marker-panel fallback (clustering)")
  obj <- tryCatch(NormalizeData(obj, verbose = FALSE), error = function(e) obj)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
  avg <- tryCatch(AverageExpression(obj, assays = "RNA")$RNA,
                  error = function(e) NULL)
  if (is.null(avg))
    avg <- AggregateExpression(obj, assays = "RNA")$RNA   # v5 命名
  cl <- colnames(avg); cl_clean <- sub("^g", "", cl)
  score <- sapply(marker_panel, function(g) {
    g <- intersect(g, rownames(avg))
    if (length(g) == 0) return(rep(0, length(cl)))
    colMeans(as.matrix(avg[g, , drop = FALSE]))
  })
  best <- apply(score, 1, function(x) {
    if (max(x) == 0) return("Unknown")
    colnames(score)[which.max(x)]
  })
  names(best) <- cl_clean   # v5 的 AverageExpression 给列名加 g 前缀，须去掉再映射
  ct <- unname(best[as.character(obj$seurat_clusters)])
  as.character(ct)
}

## ---------- 主循环 ----------
ferro_genes <- readLines(GENE_FILE)
ferro_genes <- trimws(ferro_genes[nzchar(trimws(ferro_genes))])
key_genes <- c("SLC7A11", "STEAP3", "HMOX1", "CDKN1A", "SPP1", "CRTAC1")

rds_files <- sort(list.files(DATA_DIR, pattern = "\\.rds(\\.gz)?$",
                             full.names = TRUE, ignore.case = TRUE))
stopifnot(length(rds_files) > 0)

cat("=== sample -> group mapping (请核对!) ===\n")
for (f in rds_files)
  cat(sprintf("  %-40s -> %s\n", basename(f), infer_group(basename(f))))
cat("\n")

for (f in rds_files) {
  sid <- sub("\\.rds(\\.gz)?$", "", basename(f), ignore.case = TRUE)
  grp <- infer_group(basename(f))
  out_csv <- file.path(OUT_SCORE, paste0(sid, ".csv"))
  if (file.exists(out_csv)) { cat("[skip]", sid, "already done\n"); next }
  cat("=====", sid, "(", grp, ") =====\n")

  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error") || !inherits(obj, "Seurat")) {
    cat("  !! read failed / not Seurat — skipped\n"); next
  }
  cat("  cells:", ncol(obj), " genes:", nrow(obj), "\n")

  # 老对象兼容：过滤细胞后残留的 graph 与当前细胞集不一致，
  # Seurat5 任何 [[<- 写入都会 validObject 报错 —— 直接丢弃旧 graph/降维
  #（后续聚类会重新计算）
  try(obj@graphs <- list(), silent = TRUE)
  ok <- tryCatch({ validObject(obj); TRUE }, error = function(e) FALSE)
  if (!ok) {
    cat("  object still invalid after dropping graphs -> rebuild from counts\n")
    md0 <- obj@meta.data
    cm  <- get_counts(obj)
    obj <- CreateSeuratObject(counts = cm, meta.data = md0)
    rm(md0, cm); gc(verbose = FALSE)
  }

  # 注释
  md <- obj@meta.data
  ct_col <- pick_celltype_col(md)
  if (!is.na(ct_col)) {
    cat("  using annotation column:", ct_col, "\n")
    celltype <- as.character(md[[ct_col]])
  } else {
    celltype <- annotate_by_markers(obj)
  }
  celltype[is.na(celltype)] <- "Unknown"

  # 铁死亡打分
  gs <- intersect(ferro_genes, rownames(obj))
  cat("  ferroptosis genes found:", length(gs), "/", length(ferro_genes), "\n")
  obj <- tryCatch(AddModuleScore(obj, features = list(gs), name = "FERRO",
                                 seed = 42),
                  error = function(e) {
                    message("  AddModuleScore failed (", conditionMessage(e),
                            ") — trying JoinLayers then retry")
                    obj2 <- tryCatch(JoinLayers(obj), error = function(e2) obj)
                    AddModuleScore(obj2, features = list(gs), name = "FERRO",
                                   seed = 42)
                  })
  ferro <- obj@meta.data[["FERRO1"]]

  # 关键基因表达
  dat <- get_data(obj)
  kg <- intersect(key_genes, rownames(dat))
  kg_mat <- if (length(kg) > 0)
    t(as.matrix(dat[kg, , drop = FALSE])) else NULL

  # SPP1 亚群（巨噬细胞内）
  spp1_expr <- if ("SPP1" %in% rownames(dat)) as.numeric(dat["SPP1", ]) else rep(NA, ncol(obj))
  is_macro <- grepl("macro|mono", celltype, ignore.case = TRUE)
  spp1_status <- ifelse(!is_macro, NA,
                 ifelse(is.na(spp1_expr), NA,
                 ifelse(spp1_expr > 0, "SPP1+_Macrophage", "SPP1-_Macrophage")))

  res <- data.frame(cell = colnames(obj), sample = sid, group = grp,
                    cluster = as.character(obj$seurat_clusters),
                    celltype = celltype, ferro_score = ferro,
                    spp1_status = spp1_status, stringsAsFactors = FALSE)
  if (!is.null(kg_mat)) res <- cbind(res, as.data.frame(kg_mat))
  write.csv(res, out_csv, row.names = FALSE)
  cat("  saved:", out_csv, "\n")

  rm(obj, dat, kg_mat, res); gc(verbose = FALSE)
}

## ---------- 汇总分析 ----------
cat("\n===== combining =====\n")
csvs <- list.files(OUT_SCORE, pattern = "\\.csv$", full.names = TRUE)
dfs <- lapply(csvs, function(x) read.csv(x, stringsAsFactors = FALSE))
# 各样本关键基因检出数不同导致列数不一致 -> 按列名对齐，缺失补 NA
allcols <- unique(unlist(lapply(dfs, colnames)))
dfs <- lapply(dfs, function(d) {
  miss <- setdiff(allcols, colnames(d))
  for (m in miss) d[[m]] <- NA
  d[, allcols, drop = FALSE]
})
all <- do.call(rbind, dfs)
write.csv(all, file.path(OUT_RES, "GSE271018_percell_scores.csv"),
          row.names = FALSE)
cat("total cells:", nrow(all), "\n")

# 1) 按细胞类型排序
ct_tab <- aggregate(ferro_score ~ celltype, all,
                    function(x) c(n = length(x), mean = mean(x),
                                  median = median(x)))
ct_out <- data.frame(celltype = ct_tab$celltype, ct_tab$ferro_score)
ct_out <- ct_out[order(-ct_out$mean), ]
write.csv(ct_out, file.path(OUT_RES, "GSE271018_ferroptosis_by_celltype.csv"),
          row.names = FALSE)
print(ct_out)

# 2) 退变 vs 对照（每细胞类型 Wilcoxon + BH）
grp_stat <- do.call(rbind, lapply(split(all, all$celltype), function(d) {
  d1 <- d$ferro_score[d$group == "Degenerative"]
  d0 <- d$ferro_score[d$group == "Control"]
  if (length(d1) < 10 || length(d0) < 10) return(NULL)
  p <- tryCatch(wilcox.test(d1, d0)$p.value, error = function(e) NA)
  data.frame(celltype = d$celltype[1], n_deg = length(d1),
             n_ctrl = length(d0), mean_deg = mean(d1),
             mean_ctrl = mean(d0), delta = mean(d1) - mean(d0),
             p_wilcox = p)
}))
if (!is.null(grp_stat)) {
  grp_stat$p_adj <- p.adjust(grp_stat$p_wilcox, method = "BH")
  grp_stat <- grp_stat[order(grp_stat$p_wilcox), ]
  write.csv(grp_stat, file.path(OUT_RES, "GSE271018_deg_vs_ctrl_wilcox.csv"),
            row.names = FALSE)
  print(grp_stat)
}

# 3) SPP1+ vs SPP1- 巨噬细胞
sp <- all[!is.na(all$spp1_status), ]
if (nrow(sp) > 50) {
  p_sp <- wilcox.test(ferro_score ~ spp1_status, data = sp)$p.value
  sp_tab <- aggregate(ferro_score ~ spp1_status, sp,
                      function(x) c(n = length(x), mean = mean(x)))
  sp_out <- data.frame(spp1_status = sp_tab$spp1_status,
                       sp_tab$ferro_score, p_wilcox = p_sp)
  write.csv(sp_out, file.path(OUT_RES, "GSE271018_SPP1macro_ferroptosis.csv"),
            row.names = FALSE)
  print(sp_out)
}

## ---------- 图（ggplot2 only） ----------
plot_df <- all
if (nrow(plot_df) > 100000)
  plot_df <- plot_df[sample(nrow(plot_df), 100000), ]   # 画图降采样

p1 <- ggplot(ct_out, aes(reorder(celltype, mean), mean)) +
  geom_col(fill = "#4393C3") + coord_flip() +
  labs(x = NULL, y = "Mean ferroptosis module score",
       title = "Ferroptosis activity by cell type (GSE271018)") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_RES, "FigS5A_ferroptosis_by_celltype.pdf"),
       p1, width = 6.5, height = 4.5)

if (!is.null(grp_stat)) {
  p2 <- ggplot(plot_df[plot_df$group %in% c("Control", "Degenerative"), ],
               aes(group, ferro_score, fill = group)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_wrap(~celltype, scales = "free_y") +
    scale_fill_manual(values = c(Control = "#74ADD1", Degenerative = "#D73027")) +
    labs(x = NULL, y = "Ferroptosis score",
         title = "Degenerative vs Control per cell type") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_RES, "FigS5B_deg_vs_ctrl.pdf"), p2,
         width = 10, height = 6)
}

if (nrow(sp) > 50) {
  p3 <- ggplot(sp, aes(spp1_status, ferro_score, fill = spp1_status)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("SPP1+_Macrophage" = "#D73027",
                                 "SPP1-_Macrophage" = "#74ADD1")) +
    labs(x = NULL, y = "Ferroptosis score",
         title = sprintf("SPP1+ vs SPP1- macrophages (P = %.2e)", p_sp)) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(OUT_RES, "FigS5C_SPP1macro.pdf"), p3,
         width = 4.5, height = 5)
}

kg_cols <- intersect(key_genes, colnames(all))
if (length(kg_cols) > 0) {
  kg_long <- do.call(rbind, lapply(kg_cols, function(g) {
    data.frame(celltype = all$celltype, gene = g, expr = all[[g]])
  }))
  kg_sum <- aggregate(expr ~ celltype + gene, kg_long, mean)
  p4 <- ggplot(kg_sum, aes(gene, celltype, fill = expr)) +
    geom_tile() +
    scale_fill_gradient(low = "#F7F7F7", high = "#B2182B") +
    labs(title = "Key gene mean expression by cell type",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_RES, "FigS5D_keygenes_heatmap.pdf"), p4,
         width = 6, height = 5)
}

cat("\n[ALL DONE] 结果在 ", OUT_RES, "/ — 打包下载后发给本地\n")

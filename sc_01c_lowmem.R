#!/usr/bin/env Rscript
# ============================================================================
# sc_01c_lowmem.R — 低内存版：逐样本 QC → 稀疏矩阵直接 cbind → 一次建对象
# 规避 Seurat merge() 内存峰值。输出与 sc_01 相同。
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(Matrix)
})
BASE_DIR <- getwd()
OUT <- file.path(BASE_DIR, "scRNA_out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

samples <- data.frame(
  sample_id = c("LSS", "LDH", "LFH1", "LFH2", "LFH3"),
  dataset   = c("GSE294458", "GSE294458", "GSE267819", "GSE267819", "GSE267819"),
  group     = c("LF_hypertrophy", "Control_LDH",
                "LF_hypertrophy", "LF_hypertrophy", "LF_hypertrophy"),
  prefix    = c("GSM8906515_LSS", "GSM8906516_LDH",
                "GSM8279256_LFH1", "GSM8279257_LFH2", "GSM8279258_LFH3"),
  stringsAsFactors = FALSE)

read_10x_sample <- function(dataset, prefix) {
  d <- file.path(BASE_DIR, dataset)
  tmp <- file.path(d, paste0("mtx_", prefix))
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(d, paste0(prefix, "_barcodes.tsv.gz")), file.path(tmp, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(prefix, "_features.tsv.gz")), file.path(tmp, "features.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(prefix, "_matrix.mtx.gz")),    file.path(tmp, "matrix.mtx.gz"),    overwrite = TRUE)
  Read10X(data.dir = tmp)
}

# ---------- Pass A: 逐样本读取 + QC，矩阵与元数据分别存盘后立即释放 ----------
qc_stats <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  cat(">>> [Pass A]", s$sample_id, "\n")
  mat <- read_10x_sample(s$dataset, s$prefix)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  obj <- CreateSeuratObject(counts = mat, project = s$sample_id, min.cells = 3, min.features = 200)
  pct_mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  keep <- obj$nFeature_RNA >= 200 & obj$nFeature_RNA <= 6000 & pct_mt < 15
  m <- obj[["RNA"]]$counts[, keep, drop = FALSE]
  meta <- data.frame(
    cell_id = paste(s$sample_id, colnames(m), sep = "_"),
    sample_id = s$sample_id, dataset = s$dataset, group = s$group,
    nFeature_RNA = obj$nFeature_RNA[keep], nCount_RNA = obj$nCount_RNA[keep],
    percent.mt = pct_mt[keep], stringsAsFactors = FALSE)
  colnames(m) <- meta$cell_id
  saveRDS(m, file.path(OUT, paste0("tmp_counts_", s$sample_id, ".rds")))
  saveRDS(meta, file.path(OUT, paste0("tmp_meta_", s$sample_id, ".rds")))
  qc_stats[[s$sample_id]] <- data.frame(sample_id = s$sample_id, group = s$group,
                                        n_cells_raw = ncol(obj), n_cells_post = ncol(m),
                                        median_nFeature = median(obj$nFeature_RNA[keep]),
                                        median_nCount = median(obj$nCount_RNA[keep]),
                                        median_pctMT = round(median(pct_mt[keep]), 2))
  rm(mat, obj, m, meta, pct_mt, keep); gc(verbose = FALSE)
}
qc_df <- bind_rows(qc_stats)
write.csv(qc_df, file.path(OUT, "sc_QC_summary.csv"), row.names = FALSE)
print(qc_df)

# ---------- Pass B: 公共基因（交集）逐样本 cbind ----------
cat(">>> [Pass B] combining matrices\n")
mats <- list()
for (sid in samples$sample_id) mats[[sid]] <- readRDS(file.path(OUT, paste0("tmp_counts_", sid, ".rds")))
common_genes <- Reduce(intersect, lapply(mats, rownames))
cat(">>> common genes:", length(common_genes), "\n")
big <- NULL
for (sid in samples$sample_id) {
  m <- mats[[sid]][common_genes, , drop = FALSE]
  big <- if (is.null(big)) m else cbind(big, m)
  mats[[sid]] <- NA; rm(m); gc(verbose = FALSE)
  cat(">>> cbind", sid, "done, total cells:", ncol(big), "\n")
}
rm(mats); gc(verbose = FALSE)

metas <- bind_rows(lapply(samples$sample_id,
                          function(sid) readRDS(file.path(OUT, paste0("tmp_meta_", sid, ".rds")))))
rownames(metas) <- metas$cell_id
metas <- metas[colnames(big), ]

merged <- CreateSeuratObject(counts = big, meta.data = metas, project = "LF_scRNA")
rm(big, metas); gc(verbose = FALSE)
cat(">>> merged object cells:", ncol(merged), "\n")

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = "percent.mt", verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

if (requireNamespace("harmony", quietly = TRUE)) {
  cat(">>> Running Harmony\n")
  merged <- harmony::RunHarmony(merged, group.by.vars = "sample_id", verbose = FALSE)
  red <- "harmony"
} else { red <- "pca" }

merged <- RunUMAP(merged, reduction = red, dims = 1:30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = red, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)

saveRDS(merged, file.path(OUT, "sc_merged_clustered.rds"))

pdf(file.path(OUT, "Fig_sc_UMAP_clusters.pdf"), width = 12, height = 5)
p1 <- DimPlot(merged, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(merged, group.by = "sample_id") + ggtitle("Samples")
print(p1 + p2)
dev.off()

pdf(file.path(OUT, "Fig_sc_UMAP_group.pdf"), width = 6, height = 5)
print(DimPlot(merged, group.by = "group") + ggtitle("Group"))
dev.off()

# 清理临时矩阵文件
invisible(file.remove(list.files(OUT, pattern = "^tmp_(counts|meta)_.*\\.rds$", full.names = TRUE)))
cat("=== sc_01c done. Cells total:", ncol(merged), "===\n")

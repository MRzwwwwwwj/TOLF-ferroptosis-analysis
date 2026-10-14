# ============================================================
# 06b 富集分析: DEG x PCD 交集 + g:Profiler (GO/KEGG/Reactome)
# 绕开 GO.db/org.Hs.eg.db (本地安装受阻), 使用 g:Profiler 官方 API
# 参考文献: Kolberg L, et al. Nucleic Acids Res. 2023 (g:Profiler2)
# ============================================================
wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
setwd(wd)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(ggplot2)
})

# ---------- 1. 读入 DEG 与 PCD 基因集 ----------
deg <- read.csv("results/DEG_sig_exploratory.csv", stringsAsFactors = FALSE)
deg_up <- deg$Gene[deg$logFC > 0]; deg_dn <- deg$Gene[deg$logFC < 0]
cat("DEG up:", length(deg_up), " down:", length(deg_dn), "\n")

pcs <- c("ferroptosis","pyroptosis","cuproptosis","disulfidptosis","panoptosis")
sets <- lapply(pcs, function(p) {
  g <- readLines(file.path("pcd_genesets", paste0("gs_", p, ".txt")))
  g <- trimws(g[nzchar(trimws(g)) & !startsWith(trimws(g), "#")])
  g
}); names(sets) <- pcs

# ---------- 2. 交集: OPCDEGs ----------
inter_rows <- list()
opcd_up <- c(); opcd_dn <- c()
for (p in pcs) {
  u <- intersect(deg_up, sets[[p]]); d <- intersect(deg_dn, sets[[p]])
  opcd_up <- c(opcd_up, u); opcd_dn <- c(opcd_dn, d)
  if (length(u)) for (g in u) inter_rows[[length(inter_rows)+1]] <- data.frame(Gene=g, PCD=p, Direction="Up")
  if (length(d)) for (g in d) inter_rows[[length(inter_rows)+1]] <- data.frame(Gene=g, PCD=p, Direction="Down")
  cat(sprintf("%-14s 集大小%3d | up交集%2d | down交集%2d\n", p, length(sets[[p]]), length(u), length(d)))
}
opcd <- unique(c(opcd_up, opcd_dn))
cat("OPCDEG 总数:", length(opcd), "\n")
if (length(inter_rows)) {
  inter_df <- do.call(rbind, inter_rows)
  write.csv(inter_df, "results/OPCDEGs_intersection.csv", row.names = FALSE)
}

# ---------- 3. g:Profiler API 封装 ----------
gp_enrich <- function(query, background, sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC")) {
  body <- list(organism = "hsapiens",
               query = query,
               sources = sources,
               domain_scope = "custom",
               background = background,
               user_threshold = 0.05,
               significance_threshold_method = "g_SCS")
  r <- POST("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            body = toJSON(body, auto_unbox = TRUE), encode = "json",
            timeout(120))
  stop_for_status(r)
  res <- fromJSON(content(r, "text"))$result
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res
}

bg <- read.csv("results/expr_gene_GSE113212.csv", stringsAsFactors = FALSE)
bg_genes <- unique(bg$Gene)
cat("背景基因(芯片检出):", length(bg_genes), "\n")

queries <- list(DEG_up = deg_up,
                DEG_dn = deg_dn,
                DEG_all = c(deg_up, deg_dn),
                OPCDEG_ferroptosis = intersect(c(opcd_up, opcd_dn), sets$ferroptosis))
queries <- queries[sapply(queries, length) >= 5]   # API 需要足够基因

all_res <- list()
for (nm in names(queries)) {
  cat("富集查询:", nm, "(", length(queries[[nm]]), "genes ) ... ")
  res <- try(gp_enrich(queries[[nm]], bg_genes), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res)) { cat("无显著条目或失败\n"); next }
  res$query_set <- nm
  all_res[[length(all_res)+1]] <- res
  cat(nrow(res), "条显著条目\n")
  Sys.sleep(2)
}

if (length(all_res)) {
  enr <- do.call(rbind, all_res)
  keep <- c("query_set","source","native","name","p_value","significant","term_size",
            "query_size","intersection_size","effective_domain_size","intersection")
  enr <- enr[, intersect(keep, colnames(enr))]
  write.csv(enr, "results/enrichment_gprofiler_all.csv", row.names = FALSE)

  # ---------- 4. 点图: 每个来源 top 条目 ----------
  enr$ GeneRatio <- enr$intersection_size / enr$query_size
  topn <- function(df, n = 8) {
    df <- df[order(df$p_value), ]
    df <- df[!duplicated(df$name), ]
    head(df, n)
  }
  for (qs in unique(enr$query_set)) {
    sub <- enr[enr$query_set == qs, ]
    for (src in c("GO:BP","KEGG","REAC")) {
      s2 <- sub[sub$source == src, ]
      if (nrow(s2) == 0) next
      plt <- topn(s2)
      plt$name <- factor(plt$name, levels = rev(plt$name))
      p <- ggplot(plt, aes(x = GeneRatio, y = name, color = p_value, size = intersection_size)) +
        geom_point() +
        scale_color_gradient(low = "#d73027", high = "#4575b4", name = "p-value") +
        scale_size_continuous(name = "Count") +
        labs(x = "Gene ratio", y = NULL,
             title = paste0("Enrichment (", src, ") - ", qs)) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(size = 11, face = "bold"))
      ggsave(file.path("figures", paste0("Fig_enrich_", gsub("[:]", "_", src), "_", qs, ".png")),
             p, width = 9, height = 4.2, dpi = 300)
    }
  }
  cat("富集结果与点图已保存\n")
} else {
  cat("!! 无富集结果\n")
}

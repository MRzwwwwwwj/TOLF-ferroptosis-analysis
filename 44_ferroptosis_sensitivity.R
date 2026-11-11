# ============================================================
# D2-G1: 铁死亡评分三套敏感性分析（回应 R4 评分前提质疑）
#   三套: S1 促/抑拆分; S2 去五套 PCD 重叠基因; S3 去铁代谢基因
#   每套都重跑 GSVA + ssGSEA + 组间 Wilcoxon + 与巨噬细胞相关
#   输出全部写入 results/，绝不覆盖原 PCD_scores_*
# ============================================================
suppressPackageStartupMessages({library(GSVA)})

wd <- Sys.getenv("TOLF_WS", ".")  # 工作区根目录（可用环境变量 TOLF_WS 覆盖）
res <- file.path(wd, "results")
gs_dir <- file.path(wd, "pcd_genesets")

expr <- as.matrix(read.csv(file.path(res, "expr_gene_GSE113212.csv"),
                           row.names = 1, check.names = FALSE))
group <- read.csv(file.path(res, "group_GSE113212.csv"), stringsAsFactors = FALSE)
expr <- expr[, group$sample]

sets <- list()
for (f in list.files(gs_dir, pattern = "^gs_.*\\.txt$", full.names = TRUE)) {
  nm <- sub("\\.txt$", "", sub("^gs_", "", basename(f)))  # 防御性：sub 只替换首个匹配，.txt 必须单独去掉
  sets[[nm]] <- trimws(readLines(f)); sets[[nm]] <- sets[[nm]][sets[[nm]] != ""]
}
fer <- sets$ferroptosis
cat("原始铁死亡基因集:", length(fer), "个; 表达矩阵覆盖",
    sum(fer %in% rownames(expr)), "个\n")

# ---- 构造三套敏感性基因集 ----
# S1: FerrDb 促/抑标注（driver vs suppressor），仅保留表达矩阵覆盖的
driver <- c("ACSL4", "LPCAT3", "ALOX15", "ALOX12", "TFRC", "NCOA4",
            "HMOX1", "GLS2", "GLS", "SAT1", "TP53", "CDKN1A", "AIFM2",
            "DHODH", "ATG5", "ATG7", "BECN1", "SLC40A1", "STEAP3",
            "ACSL3", "POR", "CYBB", "ACACA")
suppressor <- c("GPX4", "SLC7A11", "SLC3A2", "FTH1", "FTL", "NFE2L2",
                "KEAP1", "GCH1", "HSPB1", "CARS", "FADS2", "SCD")
# S2: 去掉与其他四套 PCD 重叠的基因
other <- unique(unlist(sets[names(sets) != "ferroptosis"]))
# S3: 去铁代谢/铁稳态基因
iron <- c("TFRC", "FTH1", "FTL", "NCOA4", "SLC40A1", "CP", "STEAP3",
          "HMOX1", "HSPB1")

inv <- function(x) intersect(x, rownames(expr))
variants <- list(
  original   = inv(fer),
  S1_driver  = inv(intersect(fer, driver)),
  S1_suppressor = inv(intersect(fer, suppressor)),
  S2_noOverlap = inv(setdiff(fer, other)),
  S3_noIron  = inv(setdiff(fer, iron))
)
for (nm in names(variants)) cat(sprintf("  %-16s %2d 基因\n", nm, length(variants[[nm]])))

# ---- 组间比较函数 ----
cmp <- function(sc, tag, setname) {
  out <- data.frame()
  for (m in rownames(sc)) {
    y <- sc[m, group$group == "young_control"]
    e <- sc[m, group$group == "elderly_hypertrophied"]
    wt <- suppressWarnings(wilcox.test(e, y))
    out <- rbind(out, data.frame(set = setname, method = tag, Module = m,
                                mean_young = mean(y), mean_elderly = mean(e),
                                delta = mean(e) - mean(y),
                                direction = ifelse(mean(e) > mean(y), "up_in_hypertrophy", "down_in_hypertrophy"),
                                wilcox_P = wt$p.value))
  }
  out
}

all_res <- data.frame()
cor_res <- data.frame()
imm <- as.matrix(read.csv(file.path(res, "immune_ssgsea_GSE113212.csv"),
                          row.names = 1, check.names = FALSE)[, group$sample])

for (nm in names(variants)) {
  g <- list(ferroptosis = variants[[nm]])
  if (length(variants[[nm]]) < 2) { cat("跳过(基因过少):", nm, "\n"); next }
  gp <- tryCatch(gsvaParam(expr, g, kcdf = "Gaussian", minSize = 1), error = function(e) NULL)
  sp <- tryCatch(ssgseaParam(expr, g, minSize = 1), error = function(e) NULL)
  if (is.null(gp) || is.null(sp)) { cat("gsvaParam 失败:", nm, "\n"); next }
  gv <- suppressWarnings(tryCatch(gsva(gp, verbose = FALSE), error = function(e) NULL))
  ss <- suppressWarnings(tryCatch(gsva(sp, verbose = FALSE), error = function(e) NULL))
  if (is.null(gv) || is.null(ss)) { cat("gsva 运行失败:", nm, "\n"); next }
  all_res <- rbind(all_res, cmp(gv, "GSVA", nm), cmp(ss, "ssGSEA", nm))

  # 与巨噬细胞（Bindea 28 细胞集里的 Macrophages）相关
  for (meth in c("GSVA", "ssGSEA")) {
    sc <- if (meth == "GSVA") gv else ss
    for (cell in rownames(imm)) {
      ct <- suppressWarnings(cor.test(as.numeric(sc["ferroptosis", ]),
                                       as.numeric(imm[cell, ]), method = "pearson"))
      cor_res <- rbind(cor_res, data.frame(set = nm, method = meth, immune_cell = cell,
                                           r = unname(ct$estimate), P = ct$p.value))
    }
  }
}
write.csv(all_res, file.path(res, "ferroptosis_sensitivity_groupcompare.csv"), row.names = FALSE)
write.csv(cor_res, file.path(res, "ferroptosis_sensitivity_immune_cor.csv"), row.names = FALSE)

cat("\n=== 组间比较（全部敏感性集）===\n")
print(all_res[all_res$Module == "ferroptosis", ], row.names = FALSE)
cat("\n=== 与 Macrophages 相关（各敏感性集）===\n")
print(cor_res[grepl("Macroph", cor_res$immune_cell, ignore.case = TRUE), ], row.names = FALSE)
cat("\n=== 完成 ===\n")

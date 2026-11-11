# TOLF Analysis Code

Analysis code for the manuscript:

**"Landscape of Multi-Modal Programmed Cell Death and Immune Microenvironment Remodeling in Thoracic Ossification of the Ligamentum Flavum: Integrative Bioinformatic Analysis with Epitranscriptomic and Single-Cell Validation"**

*(submitted to the International Journal of Molecular Sciences (IJMS), MDPI)*

**Permanent archive (Figshare)**: https://doi.org/10.6084/m9.figshare.34489533

**Source-code repository (GitHub)**: https://github.com/MRzwwwwwwj/TOLF-ferroptosis-analysis

## Overview

All analyses are purely computational, based on publicly available GEO datasets:

| Dataset | Type | Use |
|---|---|---|
| GSE113212 | Bulk microarray (GPL17077) | Discovery cohort (hypertrophied vs non-hypertrophied LF) |
| GSE206560 | m6A MeRIP-seq (OLF) | External validation + gene-level m6A layer |
| GSE106253 | Bulk microarray (GPL21827, OLF) | Third bulk cohort (LASSO signature validation) |
| GSE294458 | scRNA-seq (OLF) | Single-cell discovery |
| GSE267819 | scRNA-seq (OLF) | Single-cell validation |
| GSE271018 | scRNA-seq (degenerative spinal ligaments) | Independent single-cell replication cohort |

## System Requirements

- **R 4.5.3** (Windows, local bulk + single-cell pipelines); server-side scripts (GSE271018) run on Linux with **R 4.4.3 / Seurat 5.5.1**
- Key packages: `limma`, `GSVA`, `data.table`, `ggplot2`, `pheatmap`, `e1071`, `preprocessCore`, `glmnet`, `pROC`, `fgsea`, `Seurat` (5.x), `CellChat`, `monocle3`, `slingshot`, `harmony`
- Full per-script run order and parameters are documented in the manuscript's Section 2 and in `rerun_logs_2026-10-13/`.

## Pipeline Structure

### Bulk pipeline (local, run in numeric order)
`01_data_prep.R` → `02_limma.R` (DEGs) → `03_gsva_pcd.R` (ssGSEA scores of 5 PCD modalities) → `04_immune.R` / `04b_cibersort.R` / `04c_cibersort_post.R` (immune deconvolution) → `05_cor_heatmap.R` → `06b_enrichment.R` (g:Profiler) → `07_external_validation_GSE206560.R` + `07b_m6a_pcd_cross.R` (m6A layer; `07b_prep_m6a.py` is a small Python helper) → `08_m6a_regulators_GSE113212.R` → `09_gsea_fgsea.R` → `10_lasso_roc.R` / `10b_lasso_external_roc.R` / `11_lasso_external_GSE106253.R` (signature + external ROC) → `12_nomogram_DCA.R` → `13_ckpt_M1M2.R` → `14_figure_panels.R` (+ `14b`–`14e` panel reassembly fixes).

### Single-cell pipeline (local)
`sc_01c_lowmem.R` (integration, low-memory settings) → `sc_02_annotate.R` → `sc_03_validate.R` → `sc_04_cross_validation.R` → `sc_05_cellchat.R` → `sc_06_monocle3.R` / `sc_06b_slingshot.R` / `sc_06c_finish.R` / `sc_06d_monocle3_nopartition.R` → `sc_07_cellchat_cross.R`.

### Server-side replication cohort (GSE271018)
- `server_GSE271018/`: ferroptosis-score validation (see `README_服务器操作指令.md`)
- `server_cellchat_GSE271018/`: CellChat replication (see `README_服务器CellChat复现指令.md`)

### Reproducibility evidence
`rerun_logs_2026-10-13/` and `server_rerun_2026-10-13/` contain logs from a full independent re-run of the entire pipeline (byte-identical CSV outputs; figures verified pixel-level). See the manuscript's reproducibility statement.

## Reproducibility Notes

- Random components are seeded: CIBERSORT uses its built-in seed; `AddModuleScore` uses `seed = 42`; fgsea sets seeds explicitly. Byte-identical reproduction of all CSV outputs is expected on the same R/package versions.
- Raw input data are NOT included; download the GEO datasets above from NCBI GEO.
- Long-running scripts (sc_05, sc_06) support checkpoint resume.

## CIBERSORT License Notice

`04b_cibersort.R` / `04c_cibersort_post.R` rely on the **CIBERSORT source code (`CIBERSORT.R`) and the LM22 signature matrix**, developed by Newman et al. (Stanford University; Nat Methods 2015, PMID 25822800). These two files are **NOT open-source** and are shared by their authors for **non-commercial academic research only**. They are **not redistributed in this repository**; to reproduce the deconvolution step, request them from the authors at https://cibersort.stanford.edu and place them next to `04b_cibersort.R`. All other scripts in this repository are provided under the MIT license.

## Contact

Corresponding author: Dingsheng Zha — see the manuscript for contact details.

## 2026-10-14 补充：图形修复脚本（不影响数据）

- `fix_Fig4E_title.R`：只重画 Fig4E（sc_03 的 pdf width=6 裁掉标题左缘 "F"），width=7.5，AddModuleScore seed=42 不变
- `server_fix_FigS5C_title.R`：服务器端只重画 FigS5C（width=4.5 裁掉 P 值指数），读 GSE271018_percell_scores.csv，width=5.5；用法 `Rscript fix_FigS5C_title.R`（OUT_RES 指向 server_GSE271018/GSE271018_results）
- `10_lasso_roc.R` 已修正：ROC 标题缩短（去 "discriminating"）防右缘裁切；重跑验证 LASSO/ROC csv 逐字节零漂移
- `14_figure_panels.R` 已修正：FigS7 段路径 monocle → monocle_v2（旧路径会把生产图覆盖回 ρ=0.31 归档版）
- 矢量拼版：11 张图另有 17cm 宽 PDF 版（figures/panels_pdf/，由 figure-merge 技能生成），投稿优先使用


## 2026-10-15 补充：标签清洗与投稿图组装

- `05_cor_heatmap.R` 已更新：主图**固定使用 ssGSEA28**（论文报告的 115 对 PCD × 免疫相关即出自该矩阵）；CIBERSORT (LM22) 相关矩阵由 `04b_cibersort.R` / `04c_cibersort_post.R` 单独产出，不进本脚本，以免 `pick` 逻辑在不同环境下切换主图数据集。读取 PCD 评分后**防御性去除模块名的 `.txt` 后缀**（旧版结果 CSV 曾残留文件名后缀并污染 Figure 1B/1D 与 Table S2）；图中模块名首字母大写、细胞名下划线转空格——**仅显示层处理，写出 CSV 的取值与统计结果零改动**（重跑后 115 对相关的 r 与 P 与论文正文完全一致）。
- `34_rebuild_fig1_pdf.py`：把面板 PNG 按固定坐标组装为**单页投稿 PDF**（Figure 1：页面 481.89 × 270.66 pt，A–D 四幅面板按原坐标嵌入，标签用 Helvetica 15pt）。原始投稿 PDF 的组装步骤脚本已遗失，此脚本补齐该环节，保证 panels_pdf 可由脚本重建。
- `make_suppl_tables.py`：由结果 CSV 汇总生成 Table S1–S9（不做任何重新计算）。路径常量支持环境变量 `TOLF_WS` / `TOLF_RES` 覆盖。
- `Rrun.sh` 及全部 R 脚本的 `.libPaths` 已统一到本机当前 R 安装路径（本机迁移所致，与他人环境无关；使用时按自身 R 路径调整即可）。

## Version 2 (2026-10-16)

Submitted alongside the IJMS manuscript. Changes relative to Version 1 (2026-10-14):

- `05_cor_heatmap.R`: the main PCD x immune correlation figure (Figure 1D) is pinned to the
  ssGSEA28 matrix; CIBERSORT (LM22) correlations are produced separately by
  `04b_cibersort.R` / `04c_cibersort_post.R`. Module names are defensively stripped of a
  residual `.txt` filename suffix inherited from legacy result CSVs, which had leaked into
  Figure 1B/1D and Supplementary Table S2. Module and cell-type names are prettified for
  display only. **All statistics are unchanged**: r and P for the 115 correlation pairs
  reproduce exactly.
- `32_clean_pcd_txt_labels.py` (new): the exact cleanup script used to purge the `.txt`
  suffix from the result CSVs, retained for provenance. Backups are written before any file
  is modified.
- `34_rebuild_fig1_pdf.py` (new): assembles panel PNGs into the single-page submission PDF
  for Figure 1 (481.89 x 270.66 pt; panels A-D at the original coordinates; Helvetica 15pt
  labels). This step was missing from Version 1.
- `make_suppl_tables.py` (new): builds Supplementary Tables S1-S9 from the result CSVs with
  no recomputation. Paths overridable via `TOLF_WS` / `TOLF_RES`.
- `make_graphical_abstract_v2.py` (new): generates the graphical abstract used for
  submission (1100 x 620 px). Paths overridable via `TOLF_WS` / `TOLF_SUB`.

The Figshare DOI is unchanged across versions.

## Version 3 (2026-11-11)

Released upon manuscript acceptance. Adds the robustness / sensitivity analyses written
during peer review (revision rounds 1 and 2) that produce the **accepted** version of the
Supplementary Materials (Tables S10-S15, Figures S8 and S10). All core results (bulk,
GSVA, immune, m6A, GSEA, LASSO signature, single-cell, CellChat) are unchanged from
Version 2 and reproduce exactly; the scripts below are *additional* reviewer-requested
validations, not modifications to the main pipeline.

| Script | Produces | Reviewer point addressed |
|---|---|---|
| `39_nested_cv_signature.R` | Table S10 (panels A-C) + Fig S8 data | strict nested CV + permutation test against the "AUC = 1.0 = overfitting" concern (R1/R3/R4) |
| `39b_lambda_grid_extension.R` | Table S10 (lambda-grid extension) | "lambda.min sits at the grid boundary" (R4) |
| `83_auc_ci_round2.R` | Table S10 (panel D: AUC 95% CI) | AUC confidence intervals, DeLong + 2,000-sample stratified bootstrap + exact Clopper-Pearson (round-2 comment 2) |
| `44_ferroptosis_sensitivity.R` | Table S11 | ferroptosis-score sensitivity across gene-set definitions (R4) |
| `43_pseudobulk_BH.py` | Table S12 | sample-level pseudobulk + Benjamini-Hochberg correction (R4 "pseudo-replication / multiple testing") |
| `45_m6A_depth_and_H2_trace.py` | Table S13 | m6A peak-count sensitivity to sequencing depth (R4) |
| `47_slc7a11_celltype_expr.R` | Table S14 | key-gene expression by single-cell type (R4 cell-localisation query) |
| `54_add_tableS15_parameters.py` | Table S15 | full single-cell / CellChat / monocle3 / robustness parameters (R1) |
| `49_make_figS8.py` | Figure S8 | rendering of the nested-CV null + sensitivity forest plot |
| `53_make_figS10_study_design.py` | Figure S10 | study-design flowchart (R4 / R2) |
| `48_build_suppl_new_tables.py` | (assembles Tables S10-S14) | pipeline glue that collects the above outputs |

**Run prerequisites.** These scripts continue from intermediate results produced by the
main pipeline (Steps `01`-`14`, `sc_*`) and the GEO raw data listed above; they are not
standalone. Set the working directory via `TOLF_WS` (workspace root, default `.`) and, for
`83_auc_ci_round2.R`, the main-pipeline results directory via `TOLF_RES` (default
`rerun_2026-10-13/results`). The machine-specific `.libPaths` lines have been removed for
portability; use your own R library paths. Raw GEO data are NOT bundled - download them
from NCBI GEO.

The Figshare DOI (10.6084/m9.figshare.34489533) is a concept DOI and resolves to the
latest version; it is unchanged across Versions 1-3.

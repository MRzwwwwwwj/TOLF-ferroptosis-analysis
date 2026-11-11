# -*- coding: utf-8 -*-
import os
r"""
53_make_figS10_study_design.py — 回修补充图 Figure S10
回应 R4 Minor 4.1 / R2 2.2：Study design + GEO 检索策略 + 纳入/排除 + 数据集流程
输出：figures/Fig_S10_study_design.png
纯英文、白色背景、300 dpi，可直接作为 IJMS 补充图。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.environ.get("TOLF_WS", "."), "figures", "Fig_S10_study_design.png")

BLUE = "#185FA5"   # discovery bulk
ORANGE = "#D85A30" # validation bulk
GREEN = "#2F7A52"  # single-cell
GREY = "#5A6570"   # search / output
LIGHT = "#F4F7FA"
EDGE = "#C9D3DE"

fig, ax = plt.subplots(figsize=(16, 9.6), dpi=300)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
fig.patch.set_facecolor("white")

def box(x, y, w, h, text, fc, ec, fontsize=9.2, bold=False, align="center", lw=1.2, alpha=1.0):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.008",
                       linewidth=lw, edgecolor=ec, facecolor=fc, alpha=alpha)
    ax.add_patch(p)
    ha = {"center": "center", "left": "left"}[align]
    tx = x + w / 2 if align == "center" else x + 0.012
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fontsize,
            fontweight="bold" if bold else "normal", color="#1B2733", linespacing=1.35)
    return p

def arrow(x1, y1, x2, y2, color=GREY, lw=1.5, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=13, linewidth=lw, color=color,
                                 linestyle=ls, shrinkA=2, shrinkB=2))

# ---------- 标题 ----------
ax.text(0.5, 0.972, "Figure S10. Study design and dataset flow", ha="center", va="center",
        fontsize=15, fontweight="bold", color="#1B2733")
ax.text(0.5, 0.943, "Re-analysis of six public GEO datasets; no new human specimens or wet-lab experiments",
        ha="center", va="center", fontsize=9.5, color=GREY, style="italic")

# ---------- 顶部：GEO 检索条 ----------
box(0.025, 0.855, 0.95, 0.072,
    "GEO search (accessed 2026):  'ligamentum flavum'  |  'ossification of the ligamentum flavum'  |  'hypertrophy'  |  'lumbar spinal stenosis'\n"
    "Inclusion: human surgical ligament tissue; case-control or ossified vs non-ossified design; >=2 samples per group; gene-level expression or MeRIP-seq peak tables; documented group labels.\n"
    "Exclusion: non-human, cell-line-only, single-sample or unlabeled series.",
    LIGHT, EDGE, fontsize=8.5, bold=False, align="left", lw=1.0)

# ---------- 列标题 ----------
ax.text(0.135, 0.833, "DATASET", ha="center", fontsize=9.5, fontweight="bold", color=GREY)
ax.text(0.405, 0.833, "ROLE / n / SEGMENT", ha="center", fontsize=9.5, fontweight="bold", color=GREY)
ax.text(0.755, 0.833, "ANALYSIS OUTPUT", ha="center", fontsize=9.5, fontweight="bold", color=GREY)

rows = [
    # (y, color, dataset, role, output)
    (0.700, BLUE, "GSE113212\nAgilent GPL17077; bulk",
     "DISCOVERY (n = 8)\n4 elderly hypertrophic (LSS) vs 4 young\nnon-hypertrophic (LDH)  |  LUMBAR",
     "DEGs (720 exploratory); 5 PCD module scores;\nimmune infiltration; PCD-immune correlations;\nm6A-regulator expression; LASSO signature TRAINING"),
    (0.560, ORANGE, "GSE206560\nHiSeq 4000 MeRIP-seq; bulk",
     "VALIDATION 1 (n = 3 vs 2 Input;\n3 vs 3 IP)  |  OLF vs non-ossified  |  THORACIC",
     "External PCD direction; m6A peak counts;\nSLC7A11 hypermethylation/downregulation;\nexternal ROC; depth sensitivity"),
    (0.420, ORANGE, "GSE106253\nArraystar LncRNA V4.0; bulk",
     "VALIDATION 2 (n = 4 vs 4)\nOLF vs non-ossified trauma controls  |  THORACIC*",
     "Second external ROC (5/6 genes);\nSLC7A11 downregulation; signature portability"),
    (0.280, GREEN, "GSE294458 + GSE267819\nscRNA-seq (2 + 3 samples)",
     "SINGLE-CELL DISCOVERY\n33,817 cells; 7 cell types  |  LUMBAR",
     "Cell-type ferroptosis ranking; SPP1+ subset;\nSLC7A11 localization; CellChat (SPP1-CD44/\nintegrin); monocle3 pseudotime"),
    (0.140, GREEN, "GSE271018\nscRNA-seq (9 samples)",
     "SINGLE-CELL REPLICATION\n51,686 cells; 6 degenerative vs 3 traumatic  |  LUMBAR",
     "Independent replication; sample-level\npseudobulk (6 vs 3); BH correction of\nall correlation matrices"),
]
for y, c, ds, role, out in rows:
    box(0.025, y, 0.215, 0.118, ds, "white", c, fontsize=9.3, bold=True, lw=1.8)
    box(0.258, y, 0.295, 0.118, role, LIGHT, c, fontsize=8.9, lw=1.2)
    box(0.585, y, 0.390, 0.118, out, "white", c, fontsize=8.7, lw=1.2)
    arrow(0.243, y + 0.059, 0.256, y + 0.059, c, lw=1.6)
    arrow(0.555, y + 0.059, 0.583, y + 0.059, c, lw=1.6)

# ---------- 底部：证据层级说明 ----------
box(0.025, 0.026, 0.95, 0.088,
    "EVIDENCE TIER (stated explicitly in the revised manuscript)\n"
    "Computational observations:  cross-cohort score differences  |  m6A peak counts  |  expression differences  |  single-cell localization  |  signature ROC\n"
    "Untested hypotheses (candidate / association only):  m6A -> SLC7A11 regulation  |  direction of SPP1-mediated signalling  |  any diagnostic or therapeutic application\n"
    "Robustness added in revision:  leakage-free nested CV + permutation  |  pseudobulk  |  BH correction  |  scoring sensitivity  |  m6A depth sensitivity  (Tables S10-S14, Fig. S8)\n"
    "* Thoracic origin is documented for GSE206560 and probable but not explicitly documented for GSE106253; segment-stratified inference is not possible with these public data.",
    LIGHT, EDGE, fontsize=8.1, align="left", lw=1.0)

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
plt.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.12)
print("saved:", OUT)

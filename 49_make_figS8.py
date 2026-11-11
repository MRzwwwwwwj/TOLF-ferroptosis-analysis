# -*- coding: utf-8 -*-
r"""
49_make_figS8.py — 2026-11-02 回修：新增补充图 Figure S8
(A) Nested-CV permutation null 分布 vs observed AUC（回应过拟合质疑）
(B) 铁死亡评分敏感性：5 套基因集 × 2 评分法与巨噬细胞浸润的相关（森林图）
输出 figures/Fig_S8_robustness_sensitivity.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WS = os.environ.get("TOLF_WS", ".")
RES = os.path.join(WS, "results")
FIG = os.path.join(WS, "figures")
OUT = os.path.join(FIG, "Fig_S8_robustness_sensitivity.png")

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False})

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4),
                         gridspec_kw={"width_ratios": [1, 1.15]})

# ---------------- Panel A: permutation null ----------------
ax = axes[0]
null = pd.read_csv(os.path.join(RES, "nested_CV_permutation_null.csv"))["null_auc"].values
obs = 1.0
perm_p = 0.234
ax.hist(null, bins=18, color="#9bb7d4", edgecolor="white", alpha=0.9)
ax.axvline(obs, color="#c0392b", lw=2, ls="--")
ax.text(obs - 0.03, ax.get_ylim()[1] * 0.55, "observed\nAUC = 1.00",
        ha="right", va="center", color="#c0392b", fontsize=9, fontweight="bold")
ax.text(0.03, 0.96,
        f"null mean = {null.mean():.2f}\npermutation p = {perm_p:.3f} (B=200)",
        transform=ax.transAxes, ha="left", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4f6f8", ec="#cccccc"))
ax.set_xlabel("Nested-CV AUC under permuted labels")
ax.set_ylabel("Frequency")
ax.set_title("A  Permutation null distribution", loc="left", fontweight="bold")
ax.set_xlim(0, 1.05)

# ---------------- Panel B: ferroptosis sensitivity forest ----------------
ax = axes[1]
d = pd.read_csv(os.path.join(RES, "ferroptosis_sensitivity_immune_cor.csv"))
d = d[d["immune_cell"] == "Macrophages"].copy()
order = ["original", "S1_driver", "S1_suppressor", "S2_noOverlap", "S3_noIron"]
labels = {"original": "Full 36-gene set", "S1_driver": "Pro-ferroptosis only",
          "S1_suppressor": "Suppressors only", "S2_noOverlap": "No PCD overlap",
          "S3_noIron": "Iron-metabolism removed"}
d["lab"] = d["set"].map(labels)
d["y"] = d["set"].map({s: i for i, s in enumerate(order[::-1])})
colors = {"GSVA": "#185FA5", "ssGSEA": "#D85A30"}
offsets = {"GSVA": 0.18, "ssGSEA": -0.18}
for m in ["GSVA", "ssGSEA"]:
    sub = d[d["method"] == m]
    ax.scatter(sub["r"], sub["y"] + offsets[m], s=55, color=colors[m],
               label=m, zorder=3,
               marker="o" if m == "GSVA" else "s",
               edgecolor="white", linewidth=0.6)
    for _, rr in sub.iterrows():
        star = "*" if rr["P"] < 0.05 else "ns"
        ax.text(rr["r"] + 0.02, rr["y"] + offsets[m], f"{rr['r']:.2f}{star}",
                va="center", fontsize=7.5, color=colors[m])
ax.axvline(0, color="#888888", lw=1)
ax.set_yticks(sorted(d["y"].unique()))
ax.set_yticklabels([labels[s] for s in order[::-1]])
ax.set_xlabel("Pearson r  (ferroptosis score vs macrophage infiltration)")
ax.set_title("B  Scoring sensitivity (all directions positive)", loc="left",
             fontweight="bold")
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.set_xlim(-0.1, 1.15)
ax.grid(axis="x", ls=":", alpha=0.4)

fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print("已生成:", OUT, "%.1f KB" % (os.path.getsize(OUT) / 1024))

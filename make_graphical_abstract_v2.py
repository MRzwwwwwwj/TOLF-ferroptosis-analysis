# -*- coding: utf-8 -*-
"""
GA v2：按导师意见重做的简洁版图形摘要（2026-10-15）
针对导师反馈："颜色一看就不对劲"、"图摘要包括排版有问题"
设计原则（对标 Cell/NEJM 系图形摘要惯例 + IJMS 合规）：
  - 单色系（蓝灰）+ 唯一强调色（橙，仅用于铁死亡节点）
  - 严格左→右三段式（DATA → ANALYSIS → FINDING），无斜跨箭头
  - 底部一条机制链（m6A–SLC7A11–ferroptosis–SPP1+ Mac 轴）
  - 与正文 Figure 5C（FigDraw 细胞示意图）构图/元素/配色完全不同 —— 满足 IJMS
    "GA 不得与正文任何图相同" 的硬规定
输出：1100 x 620 px PNG。MDPI IJMS 官方要求 "minimum required size ... 560 x 1100
  pixels (height x width)"，即宽 1100 x 高 560；本图 1100 x 620 满足要求。
路径可用环境变量覆盖：TOLF_WS（工作区）/ TOLF_SUB（投稿包根目录，可省略）
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT1 = os.path.join(os.environ.get("TOLF_WS", r"C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"), "figures", "Graphical_Abstract_v2.png")
OUT2 = os.path.join(os.environ.get("TOLF_SUB", os.environ.get("TOLF_WS", ".")),
                    "figures_按正文编号", "Graphical_Abstract_v2.png")
FONT = "Arial"

# 配色（单色系蓝灰 + 橙色强调）
INK      = "#22303E"   # 主文字
MUTE     = "#5A6B7B"   # 次要文字
PANEL_BG = "#F3F6FA"   # 面板底
PANEL_BD = "#C7D3E0"   # 面板边
TEAL     = "#2C6E8F"   # 主色
TEAL_LT  = "#DCE9F2"   # 主色浅底
ORANGE   = "#D95F02"   # 强调色（仅铁死亡）
ORANGE_LT= "#FBEADB"

W_IN, H_IN, DPI = 11.0, 6.2, 100
fig, ax = plt.subplots(figsize=(W_IN, H_IN), dpi=DPI)
ax.set_xlim(0, 11.0); ax.set_ylim(0, 6.2); ax.axis("off")
fig.patch.set_facecolor("white")


def panel(x, y, w, h, title):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                       linewidth=1.2, edgecolor=PANEL_BD,
                       facecolor=PANEL_BG, mutation_scale=0.5, zorder=1)
    ax.add_patch(p)
    ax.text(x + 0.22, y + h - 0.34, title, ha="left", va="center",
            fontsize=12, fontweight="bold", color=TEAL, fontname=FONT, zorder=4)


def item(x, y, w, h, text, fc="white", tc=INK, fs=9.5, bold=False, ec=PANEL_BD):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                       linewidth=1.0, edgecolor=ec, facecolor=fc,
                       mutation_scale=0.5, zorder=3)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight="bold" if bold else "normal",
            color=tc, fontname=FONT, zorder=5, linespacing=1.3)


def arrow(x1, y1, x2, y2, color=MUTE, lw=2.2, ms=15):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                        mutation_scale=ms, linewidth=lw, color=color,
                        zorder=2, shrinkA=3, shrinkB=3)
    ax.add_patch(a)


# ---------- 上部三段面板 ----------
PY, PH = 2.55, 3.30        # 面板 y / 高
PW = 3.10                  # 面板宽
X1, X2, X3 = 0.25, 3.95, 7.65

panel(X1, PY, PW, PH, "DATA")
item(X1 + 0.18, PY + 1.95, PW - 0.36, 0.82,
     "Bulk transcriptome\nGSE113212 (n = 8)\nLFH vs. control", fs=9.5)
item(X1 + 0.18, PY + 1.00, PW - 0.36, 0.82,
     "m6A MeRIP-seq\nGSE206560", fs=9.5)
item(X1 + 0.18, PY + 0.05, PW - 0.36, 0.82,
     "scRNA-seq, 85,503 cells\n3 independent cohorts", fs=9.5)

panel(X2, PY, PW, PH, "ANALYSIS")
item(X2 + 0.18, PY + 1.95, PW - 0.36, 0.82,
     "Five PCD modalities scored\n(GSVA / ssGSEA / GSEA)", fs=9.5)
item(X2 + 0.18, PY + 1.00, PW - 0.36, 0.82,
     "m6A regulator landscape\n+ MeRIP differential peaks", fs=9.5)
item(X2 + 0.18, PY + 0.05, PW - 0.36, 0.82,
     "CellChat + monocle3\npseudotime", fs=9.5)

panel(X3, PY, PW, PH, "FINDING")
item(X3 + 0.18, PY + 1.95, PW - 0.36, 0.82,
     "Ferroptosis most activated;\nlocalized to SPP1+ macrophages", fs=9.5, bold=True)
item(X3 + 0.18, PY + 1.00, PW - 0.36, 0.82,
     "SPP1\u2013CD44 / integrin signaling\nto fibroblasts", fs=9.5)
item(X3 + 0.18, PY + 0.05, PW - 0.36, 0.82,
     "SLC7A11 hypermethylated\nand silenced (logFC = \u22122.81)", fs=9.5)

arrow(X1 + PW + 0.06, PY + PH / 2, X2 - 0.06, PY + PH / 2)
arrow(X2 + PW + 0.06, PY + PH / 2, X3 - 0.06, PY + PH / 2)

# ---------- 底部机制链 ----------
ax.text(0.25, 2.18, "Convergent mechanism", ha="left", va="center",
        fontsize=11, fontweight="bold", color=INK, fontname=FONT)

chain = [
    ("FTO \u2193", False),
    ("m6A hyper-\nmethylation\nof SLC7A11", False),
    ("SLC7A11 \u2193", False),
    ("Ferroptosis \u2191\nSPP1+ macrophage", True),   # 唯一强调节点
    ("SPP1\u2013CD44 /\nintegrin", False),
    ("Fibroblast\nactivation,\nECM remodeling", False),
    ("Ligamentum\nflavum\nossification", False),
]
n = len(chain)
gap = 0.16
cw = (10.5 - gap * (n - 1)) / n
cy, ch = 0.55, 1.30
x = 0.25
for i, (txt, hot) in enumerate(chain):
    if hot:
        item(x, cy, cw, ch, txt, fc=ORANGE, tc="white", fs=8, bold=True, ec="none")
    else:
        item(x, cy, cw, ch, txt, fc=TEAL_LT, tc=INK, fs=9, ec=PANEL_BD)
    if i < n - 1:
        arrow(x + cw + 0.015, cy + ch / 2, x + cw + gap - 0.015, cy + ch / 2,
              color=MUTE, lw=1.8, ms=11)
    x += cw + gap

# 面板到底部链的引导箭头（一根，居中）
arrow(5.5, PY - 0.02, 5.5, cy + ch + 0.06, color=PANEL_BD, lw=1.6, ms=12)

fig.savefig(OUT1, dpi=DPI, bbox_inches=None, facecolor="white")
fig.savefig(OUT2, dpi=DPI, bbox_inches=None, facecolor="white")
plt.close(fig)

# 校验尺寸
import struct
with open(OUT1, "rb") as fh:
    d = fh.read(33)
w, h = struct.unpack(">II", d[16:24])
print("GA v2 已生成:", OUT1)
print("尺寸:", w, "x", h, "px（IJMS 最小 1100 x 560）", "OK" if (w >= 1100 and h >= 560) else "FAIL")
print("副本:", OUT2)

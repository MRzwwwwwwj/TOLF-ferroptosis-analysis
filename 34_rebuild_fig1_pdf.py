# -*- coding: utf-8 -*-
r"""
34_rebuild_fig1_pdf.py — 重建 panels_pdf/Fig1_discovery_overview.pdf（2026-10-15）
精确复刻旧版布局（1 页 481.89x270.66 pt，4 面板位图 + Helvetica 15pt 标签），
其中 B（箱线图）与 D（相关热图）换用去 ".txt" 后的新 PNG；A/C 不变。
"""
import os, shutil
import pymupdf
from PIL import Image

WS = r"C:\Users\10215\WorkBuddy\2026-10-03-00-10-07"
FIG = os.path.join(WS, "figures")
OUT = os.path.join(FIG, "panels_pdf", "Fig1_discovery_overview.pdf")
BAK = os.path.join(WS, "_backup", "pcd_txt_cleanup_2026-10-15", "Fig1_discovery_overview_old.pdf")
os.makedirs(os.path.dirname(BAK), exist_ok=True)
if not os.path.exists(BAK):
    shutil.copy2(OUT, BAK)
    print("旧版已备份:", BAK)

W, H = 481.8897705078125, 270.6602783203125
PANELS = [
    ("A", os.path.join(FIG, "Fig_volcano_GSE113212.png"),       (31.8, 30.0, 201.4, 171.3), (16.0, 28.0)),
    ("B", os.path.join(FIG, "Fig_PCD_scores_boxplot.png"),      (232.9, 30.0, 464.2, 171.3), (213.1841278076172, 28.0)),
    ("C", os.path.join(FIG, "Fig_immune_boxplot.png"),          (51.0, 195.3, 262.1, 258.7), (16.0, 193.32009887695312)),
    ("D", os.path.join(FIG, "Fig_PCD_immune_corheatmap.png"),   (318.1, 195.3, 458.9, 258.7), (293.13385009765625, 193.32009887695312)),
]

doc = pymupdf.open()
pg = doc.new_page(width=W, height=H)
for label, png, bbox, origin in PANELS:
    with Image.open(png) as im:
        print(f"  {label}: {os.path.basename(png)} {im.size}")
    pg.insert_image(pymupdf.Rect(*bbox), filename=png)
for label, png, bbox, origin in PANELS:
    pg.insert_text(pymupdf.Point(*origin), label, fontname="helv", fontsize=15, color=(0, 0, 0))
doc.save(OUT, deflate=True)
doc.close()
print("已重建:", OUT, "%.1f KB" % (os.path.getsize(OUT) / 1024))

# 渲染校验
doc = pymupdf.open(OUT)
pm = doc[0].get_pixmap(dpi=150)
chk = os.path.join(FIG, "_check_txt", "Fig1_new_pdf_render.png")
pm.save(chk)
doc.close()
print("渲染校验图:", chk)

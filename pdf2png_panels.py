# -*- coding: utf-8 -*-
"""PDF -> PNG 转换（供拼版用，300dpi）"""
import fitz, os

BASE = r"C:\Users\10215\WorkBuddy\2026-10-03-00-10-07"
OUT = os.path.join(BASE, "figures", "panels_src")
os.makedirs(OUT, exist_ok=True)

jobs = {
    "scRNA_out/Fig4A_sc_ferroptosis_UMAP.pdf": "Fig4A.png",
    "scRNA_out/Fig4B_sc_ferroptosis_by_celltype.pdf": "Fig4B.png",
    "scRNA_out/Fig4C_sc_ferroptosis_LSS_vs_LDH.pdf": "Fig4C.png",
    "scRNA_out/Fig4D_sc_keygenes_featureplot.pdf": "Fig4D.png",
    "scRNA_out/Fig4E_sc_SPP1macro_ferroptosis.pdf": "Fig4E.png",
    "GSE271018_results/FigS5A_ferroptosis_by_celltype.pdf": "FigS5A.png",
    "GSE271018_results/FigS5B_deg_vs_ctrl.pdf": "FigS5B.png",
    "GSE271018_results/FigS5C_SPP1macro.pdf": "FigS5C.png",
    "GSE271018_results/FigS5D_keygenes_heatmap.pdf": "FigS5D.png",
    "scRNA_out/cellchat/Fig_cc_E_SPP1_network.pdf": "Fig_cc_E.png",
    "scRNA_out/cellchat/Fig_cc_D_signalingRole_heatmap.pdf": "Fig_cc_D.png",
    "GSE271018_cellchat/Fig_cc_E_SPP1_network.pdf": "Fig_cc_E_srv.png",
}

for src, dst in jobs.items():
    p = os.path.join(BASE, src)
    if not os.path.exists(p):
        print("MISSING:", src); continue
    doc = fitz.open(p)
    page = doc[0]
    pix = page.get_pixmap(dpi=300)
    pix.save(os.path.join(OUT, dst))
    print("OK:", dst, pix.width, "x", pix.height)
    doc.close()
print("=== convert done ===")

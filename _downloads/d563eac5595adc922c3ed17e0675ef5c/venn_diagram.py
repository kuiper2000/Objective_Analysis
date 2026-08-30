"""
venn_diagram.py
---------------
Draws a Venn diagram using matplotlib patches.

Layout
------
A rectangle represents the sample space (total probability = 1).
Two partially overlapping circles represent events E1 and E2.
Regions are labelled with P(E1), P(E2), P(E1∩E2), and P(E1∪E2).
A text box explains the addition rule.

Probabilities used (chosen so they are easy to read off):
    P(E1) = 0.45
    P(E2) = 0.40
    P(E1 ∩ E2) = 0.15
    P(E1 ∪ E2) = P(E1) + P(E2) − P(E1 ∩ E2) = 0.70
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# ── Probabilities ─────────────────────────────────────────────────────────────
p_e1        = 0.45
p_e2        = 0.40
p_inter     = 0.15
p_union     = p_e1 + p_e2 - p_inter   # = 0.70

# ── Circle geometry ───────────────────────────────────────────────────────────
# Place both circles on a unit canvas [0,1] x [0,1]
cx1, cy = 0.37, 0.50   # centre of E1
cx2     = 0.63          # centre of E2 (same y)
r       = 0.22          # radius of both circles

fig, ax = plt.subplots(figsize=(9, 6))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.axis("off")

# ── Rectangle (sample space) ──────────────────────────────────────────────────
rect = mpatches.FancyBboxPatch((0.03, 0.08), 0.94, 0.84,
                               boxstyle="round,pad=0.01",
                               linewidth=2, edgecolor="black",
                               facecolor="#F5F5F5")
ax.add_patch(rect)

# ── Circles ───────────────────────────────────────────────────────────────────
circle1 = plt.Circle((cx1, cy), r, color="#2196F3", alpha=0.35, zorder=2)
circle2 = plt.Circle((cx2, cy), r, color="#FF9800", alpha=0.35, zorder=2)
# Redraw outlines on top for clarity
outline1 = plt.Circle((cx1, cy), r, fill=False, edgecolor="#1565C0",
                       linewidth=2.2, zorder=3)
outline2 = plt.Circle((cx2, cy), r, fill=False, edgecolor="#E65100",
                       linewidth=2.2, zorder=3)

for patch in [circle1, circle2, outline1, outline2]:
    ax.add_patch(patch)

# ── Labels inside regions ─────────────────────────────────────────────────────
kw_label = dict(ha="center", va="center", fontsize=12, fontweight="bold",
                zorder=5)

# E1 only region (left lune)
ax.text(cx1 - r*0.60, cy + 0.04,
        f"$E_1$ only\n$P = {p_e1 - p_inter:.2f}$",
        color="#1565C0", **kw_label)

# E2 only region (right lune)
ax.text(cx2 + r*0.60, cy + 0.04,
        f"$E_2$ only\n$P = {p_e2 - p_inter:.2f}$",
        color="#E65100", **kw_label)

# Intersection
ax.text((cx1 + cx2) / 2, cy,
        f"$E_1 \\cap E_2$\n$P = {p_inter:.2f}$",
        color="#4A148C", fontsize=11, fontweight="bold",
        ha="center", va="center", zorder=5)

# Outside (complement of union)
ax.text(0.88, 0.88,
        f"$P(\\overline{{E_1 \\cup E_2}})$\n$= {1 - p_union:.2f}$",
        fontsize=10, ha="center", va="center", color="#555555")

# ── Circle labels (event names) ───────────────────────────────────────────────
ax.text(cx1, cy + r + 0.04, f"$E_1$\n$P(E_1)={p_e1}$",
        ha="center", va="bottom", fontsize=12, color="#1565C0", fontweight="bold")
ax.text(cx2, cy + r + 0.04, f"$E_2$\n$P(E_2)={p_e2}$",
        ha="center", va="bottom", fontsize=12, color="#E65100", fontweight="bold")

# ── Rectangle label (sample space) ───────────────────────────────────────────
ax.text(0.06, 0.93, "$\\Omega$ (sample space, $P=1$)",
        fontsize=11, color="black", va="top")

# ── Union brace/annotation ────────────────────────────────────────────────────
ax.annotate("",
            xy=(cx2 + r, cy - r - 0.08),
            xytext=(cx1 - r, cy - r - 0.08),
            arrowprops=dict(arrowstyle="<->", color="#4A148C", lw=2))
ax.text((cx1 + cx2) / 2, cy - r - 0.12,
        f"$P(E_1 \\cup E_2) = {p_union:.2f}$",
        ha="center", va="top", fontsize=12, color="#4A148C", fontweight="bold")

# ── Text box: addition rule ───────────────────────────────────────────────────
formula = (r"$P(E_1 \cup E_2) = P(E_1) + P(E_2) - P(E_1 \cap E_2)$"
           "\n"
           rf"$= {p_e1} + {p_e2} - {p_inter} = {p_union}$")
ax.text(0.50, 0.11, formula,
        ha="center", va="center", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", fc="#E8F5E9", ec="#388E3C", lw=1.5))

# ── Title ─────────────────────────────────────────────────────────────────────
ax.set_title("Venn Diagram — Addition Rule of Probability",
             fontsize=14, fontweight="bold", pad=10)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "venn_diagram.png", dpi=150, bbox_inches="tight")
plt.show()

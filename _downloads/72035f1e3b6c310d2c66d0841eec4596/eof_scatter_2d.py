"""
eof_scatter_2d.py
-----------------
The simplest possible EOF analysis: two grid points, so every observation is a
dot on a plane and the EOFs are literally the axes of the scatter cloud.

Data: eight observations of anomalies at two stations, X = [8 times x 2 points].
The covariance matrix is C = [[6, 4], [4, 6]], whose eigenvectors are
(1,1)/sqrt(2) with lambda_1 = 10 and (1,-1)/sqrt(2) with lambda_2 = 2.

Left panel : the dots in data space (station A, station B), with the two EOF
             axes drawn with half-length sqrt(lambda_i) so their tips land
             exactly on the 1-sigma covariance ellipse.  The least-squares
             regression line of y on x is shown for contrast -- it is NOT
             EOF 1 (slope 2/3 versus 1), because regression minimizes vertical
             distance while EOF 1 minimizes perpendicular distance.

Right panel: the same dots after projecting onto the EOFs, Z = X E.  The cloud
             is now axis-aligned: the rotation has removed the correlation, and
             the spread along each axis is sqrt(lambda_i).

Both panels use an equal aspect ratio -- without it the 45-degree tilt and the
perpendicularity of the two EOFs are not faithfully drawn.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# ── Colors ───────────────────────────────────────────────────────────────────
C_EOF1 = "#2a78d6"   # blue
C_EOF2 = "#eb6834"   # orange
C_REG  = "#1baf7a"   # aqua
C_DOT  = "#0b0b0b"   # primary ink
C_MUTE = "#8c8b86"   # recessive grid / ellipse

# ── The data ─────────────────────────────────────────────────────────────────
X = np.array([
    [-4, -2],
    [-2, -4],
    [-2,  0],
    [ 0, -2],
    [ 0,  2],
    [ 2,  0],
    [ 2,  4],
    [ 4,  2],
], dtype=float)

M = X.shape[0]
X -= X.mean(axis=0)                 # already zero-mean, but be explicit

# ── Step 2: covariance matrix ────────────────────────────────────────────────
C = (X.T @ X) / M                   # -> [[6, 4], [4, 6]]

# ── Step 3: eigenanalysis (eigh: C is symmetric) ─────────────────────────────
lam, E = np.linalg.eigh(C)
order = np.argsort(lam)[::-1]       # eigh returns ascending; we want descending
lam, E = lam[order], E[:, order]

# Fix the arbitrary sign so EOF 1 points into the first quadrant (cosmetic only)
for j in range(E.shape[1]):
    if E[0, j] < 0:
        E[:, j] *= -1

# ── Step 5: the PCs ──────────────────────────────────────────────────────────
Z = X @ E
frac = lam / lam.sum()

print(f"C =\n{C}")
print(f"lambda      = {lam}")
print(f"EOF 1       = {E[:, 0]}   explains {100 * frac[0]:.1f}%")
print(f"EOF 2       = {E[:, 1]}   explains {100 * frac[1]:.1f}%")
print(f"var(PC)     = {Z.var(axis=0)}   (equals lambda)")
print(f"PC1 . PC2   = {Z[:, 0] @ Z[:, 1]:.1e}   (uncorrelated)")

# ── Figure ───────────────────────────────────────────────────────────────────
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.5, 5.8))
LIM = 5.6

for ax in (ax0, ax1):
    ax.set_aspect("equal")                      # essential: angles must be true
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.axhline(0, color=C_MUTE, lw=0.8, zorder=0)
    ax.axvline(0, color=C_MUTE, lw=0.8, zorder=0)
    ax.grid(True, color=C_MUTE, alpha=0.22, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

# --- Left panel: data space --------------------------------------------------
# 1-sigma covariance ellipse; its semi-axes are exactly sqrt(lambda_i)
angle = np.degrees(np.arctan2(E[1, 0], E[0, 0]))
ax0.add_patch(Ellipse((0, 0),
                      width=2 * np.sqrt(lam[0]), height=2 * np.sqrt(lam[1]),
                      angle=angle, facecolor=C_MUTE, alpha=0.12,
                      edgecolor=C_MUTE, lw=1.2, ls="--", zorder=1,
                      label=r"1-$\sigma$ covariance ellipse"))

# Least-squares regression of y on x -- deliberately NOT the same as EOF 1
slope = C[0, 1] / C[0, 0]
xx = np.array([-LIM, LIM])
ax0.plot(xx, slope * xx, color=C_REG, lw=2, ls=(0, (5, 3)), zorder=2,
         label=f"regression of $y$ on $x$ (slope {slope:.2f})")

# The two EOF axes, half-length sqrt(lambda) -- drawn double-ended because the
# sign of an eigenvector is arbitrary: an axis has a direction, not an arrowhead
tip_label_offset = [(0.62, -0.60), (0.52, -0.30)]   # placed in clear space
for j, (color, name) in enumerate([(C_EOF1, "EOF 1"), (C_EOF2, "EOF 2")]):
    v = E[:, j] * np.sqrt(lam[j])
    ax0.plot([-v[0], v[0]], [-v[1], v[1]], color=color, lw=2.6,
             solid_capstyle="round", zorder=4,
             label=f"{name}:  $\\lambda_{j+1}={lam[j]:.0f}$  "
                   f"({100 * frac[j]:.1f}% of variance)")
    ax0.scatter(*v, s=42, color=color, zorder=5)
    ax0.scatter(*(-v), s=42, color=color, zorder=5)
    ax0.text(v[0] + tip_label_offset[j][0], v[1] + tip_label_offset[j][1],
             name, color=color, fontsize=10.5, fontweight="bold",
             ha="left", va="top", zorder=7)

ax0.scatter(X[:, 0], X[:, 1], s=85, color=C_DOT, zorder=6,
            edgecolor="white", linewidth=1.6, label="observation ($M=8$)")

ax0.legend(loc="upper left", fontsize=8.6, framealpha=0.92,
           borderpad=0.6, labelspacing=0.55, handlelength=1.8)

ax0.set_xlabel("anomaly at station A,  $x$", fontsize=12)
ax0.set_ylabel("anomaly at station B,  $y$", fontsize=12)
ax0.set_title("Data space: the EOFs are the axes of the cloud\n"
              r"$\mathbf{C}=\frac{1}{M}\mathbf{X}^T\mathbf{X}$ = [[6, 4], [4, 6]]",
              fontsize=12)

# --- Right panel: EOF space --------------------------------------------------
for j, color in enumerate([C_EOF1, C_EOF2]):
    s = np.sqrt(lam[j])
    if j == 0:
        ax1.plot([-s, s], [0, 0], color=color, lw=2.6,
                 solid_capstyle="round", zorder=4)
    else:
        ax1.plot([0, 0], [-s, s], color=color, lw=2.6,
                 solid_capstyle="round", zorder=4)

ax1.scatter(Z[:, 0], Z[:, 1], s=85, color=C_DOT, zorder=6,
            edgecolor="white", linewidth=1.6)

ax1.text(np.sqrt(lam[0]) + 0.25, 0.45,
         f"$\\sqrt{{\\lambda_1}}={np.sqrt(lam[0]):.2f}$",
         color=C_EOF1, fontsize=10.5, fontweight="bold", va="bottom")
ax1.text(0.3, np.sqrt(lam[1]) + 0.2,
         f"$\\sqrt{{\\lambda_2}}={np.sqrt(lam[1]):.2f}$",
         color=C_EOF2, fontsize=10.5, fontweight="bold", va="bottom")

ax1.text(0.03, 0.03,
         "the same eight dots,\nrotated into their own frame:\n"
         r"$\mathrm{corr}(z_1,z_2)=0$",
         transform=ax1.transAxes, fontsize=10, va="bottom", color="#52514e")

ax1.set_xlabel(r"PC 1,  $z_1 = (x+y)/\sqrt{2}$", fontsize=12)
ax1.set_ylabel(r"PC 2,  $z_2 = (x-y)/\sqrt{2}$", fontsize=12)
ax1.set_title("EOF space: a rigid rotation\n"
              r"$\mathbf{Z}=\mathbf{X}\mathbf{E}$,  "
              r"$\mathrm{var}(z_i)=\lambda_i$", fontsize=12)

fig.tight_layout()

outfile = Path(__file__).resolve().parent / "eof_scatter_2d.png"
fig.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
print(f"\nwrote {outfile}")
plt.show()

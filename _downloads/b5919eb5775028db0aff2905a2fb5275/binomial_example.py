"""
binomial_example.py
-------------------
Three panels illustrating binomial distributions in meteorological contexts.

Left   : N=20, p=1/6 (die rolling).  Mark k=4 with a different colour.
Middle : N=30, p=0.5 (true/false test).  Shade k ≥ 12 in red; label P(X ≥ 12).
Right  : N=48, p=0.5 (CMIP5 models agreeing with observation).
         Shade k ≥ 31 in red; overlay normal approximation.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.subplots_adjust(wspace=0.38)

# ─── Left: die rolling (N=20, p=1/6) ─────────────────────────────────────────
ax = axes[0]
N1, p1 = 20, 1/6
k1 = np.arange(0, N1 + 1)
pmf1 = stats.binom.pmf(k1, N1, p1)

colors1 = ["#E53935" if k == 4 else "#1976D2" for k in k1]
ax.bar(k1, pmf1, color=colors1, edgecolor="white", linewidth=0.5)

ax.set_xlabel("$k$ (number of successes)", fontsize=12)
ax.set_ylabel("$P(X = k)$", fontsize=12)
ax.set_title(f"Die Rolling\n$N={N1}$, $p=1/6$", fontsize=12, fontweight="bold")

# Legend
from matplotlib.patches import Patch
handles = [Patch(color="#1976D2", label="PMF"),
           Patch(color="#E53935", label="$k=4$ (highlighted)")]
ax.legend(handles=handles, fontsize=9.5)

# Annotate k=4
p4 = stats.binom.pmf(4, N1, p1)
ax.text(4, p4 + 0.005, f"$P(X=4)$\n$={p4:.4f}$",
        ha="center", fontsize=8.5, color="#E53935")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Middle: true/false test (N=30, p=0.5) ────────────────────────────────────
ax = axes[1]
N2, p2, k_crit2 = 30, 0.5, 12
k2 = np.arange(0, N2 + 1)
pmf2 = stats.binom.pmf(k2, N2, p2)
prob_tail2 = stats.binom.sf(k_crit2 - 1, N2, p2)   # P(X >= 12)

colors2 = ["#E53935" if k >= k_crit2 else "#1976D2" for k in k2]
ax.bar(k2, pmf2, color=colors2, edgecolor="white", linewidth=0.5)

ax.set_xlabel("$k$ (number correct)", fontsize=12)
ax.set_ylabel("$P(X = k)$", fontsize=12)
ax.set_title(f"True/False Test\n$N={N2}$, $p=0.5$", fontsize=12, fontweight="bold")

handles = [Patch(color="#1976D2", label="$P(X < 12)$"),
           Patch(color="#E53935",
                 label=f"$P(X \\geq 12) = {prob_tail2:.4f}$")]
ax.legend(handles=handles, fontsize=9.5)

ax.axvline(k_crit2 - 0.5, color="black", linestyle="--", linewidth=1.2)
ax.text(k_crit2 + 1, pmf2.max() * 0.85,
        f"$P(X \\geq {k_crit2})$\n$= {prob_tail2:.4f}$",
        ha="left", fontsize=9, color="#E53935",
        bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: CMIP5 models (N=48, p=0.5) ────────────────────────────────────────
ax = axes[2]
N3, p3, k_crit3 = 48, 0.5, 31
k3 = np.arange(0, N3 + 1)
pmf3 = stats.binom.pmf(k3, N3, p3)
prob_tail3 = stats.binom.sf(k_crit3 - 1, N3, p3)   # P(X >= 31)

colors3 = ["#E53935" if k >= k_crit3 else "#1976D2" for k in k3]
ax.bar(k3, pmf3, color=colors3, edgecolor="white", linewidth=0.5,
       label="Binomial PMF")

# Normal approximation
mu_norm  = N3 * p3
sig_norm = np.sqrt(N3 * p3 * (1 - p3))
x_cont   = np.linspace(0, N3, 500)
pdf_norm = stats.norm.pdf(x_cont, mu_norm, sig_norm)
ax.plot(x_cont, pdf_norm, "k-", linewidth=2.0,
        label=rf"$\mathcal{{N}}({mu_norm:.0f},\,{sig_norm:.2f}^2)$ approx.")

ax.set_xlabel("$k$ (models agreeing)", fontsize=12)
ax.set_ylabel("$P(X = k)$", fontsize=12)
ax.set_title(f"CMIP5 Models\n$N={N3}$, $p=0.5$", fontsize=12, fontweight="bold")

handles_extra = [Patch(color="#1976D2", label="$P(X < 31)$"),
                 Patch(color="#E53935",
                       label=f"$P(X \\geq 31) = {prob_tail3:.4f}$")]
ax.legend(handles=ax.get_legend_handles_labels()[0][:2] + handles_extra,
          fontsize=8.5)

ax.axvline(k_crit3 - 0.5, color="black", linestyle="--", linewidth=1.2)
ax.text(k_crit3 + 1, pmf3.max() * 0.75,
        f"$P(X \\geq {k_crit3})$\n$= {prob_tail3:.4f}$",
        ha="left", fontsize=9, color="#E53935",
        bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.8))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "binomial_example.png", dpi=150, bbox_inches="tight")
plt.show()

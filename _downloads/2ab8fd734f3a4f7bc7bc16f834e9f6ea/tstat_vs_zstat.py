"""
tstat_vs_zstat.py
-----------------
Left  : t-distributions for nu = 1, 3, 10, 30 and the standard normal,
        all on the same axes.  Heavier tails for smaller nu are highlighted.
Right : critical value t_c (two-tailed, alpha=0.05) vs degrees of freedom
        nu = 1..50.  Show z_c = 1.96 as a horizontal dashed line.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

x = np.linspace(-5, 5, 800)

nus    = [1, 3, 10, 30]
colors = ["#D32F2F", "#E64A19", "#1976D2", "#388E3C"]

# ── Right panel data ──────────────────────────────────────────────────────────
nu_range  = np.arange(1, 51)
t_crits   = stats.t.ppf(0.975, df=nu_range)   # two-tailed 5%
z_crit    = 1.96

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.subplots_adjust(wspace=0.38)

# ─── Left: PDF comparison ─────────────────────────────────────────────────────
ax = axes[0]

for nu, col in zip(nus, colors):
    lbl = rf"$t(\nu={nu})$"
    ax.plot(x, stats.t.pdf(x, df=nu), color=col, linewidth=2.0, label=lbl)

# Standard normal
ax.plot(x, stats.norm.pdf(x), "k--", linewidth=2.2,
        label=r"$\mathcal{N}(0,1)$  ($\nu\to\infty$)")

# Shade tail region to highlight differences
x_tail = x[x > 2.5]
for nu, col in zip(nus, colors):
    ax.fill_between(x_tail, stats.t.pdf(x_tail, df=nu),
                    stats.norm.pdf(x_tail), alpha=0.15, color=col)

ax.set_xlim(-5, 5)
ax.set_ylim(0, 0.42)
ax.set_xlabel("$t$", fontsize=13)
ax.set_ylabel("Probability density $f(t)$", fontsize=12)
ax.set_title("$t$-distributions vs. Standard Normal\n"
             "(heavier tails for smaller $\\nu$)", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper right")

# Inset annotation arrow to tail
ax.annotate("Heavier tails\n(small $\\nu$)",
            xy=(3.5, 0.022), xytext=(4.0, 0.10),
            fontsize=9, ha="center", color="#D32F2F",
            arrowprops=dict(arrowstyle="->", color="#D32F2F"))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: critical value vs nu ─────────────────────────────────────────────
ax = axes[1]
ax.plot(nu_range, t_crits, "o-", color="#1976D2", markersize=4,
        linewidth=1.8, label=r"$t_c(\nu)$ — two-tailed $\alpha=0.05$")
ax.axhline(z_crit, color="#E53935", linestyle="--", linewidth=1.8,
           label=f"$z_c = {z_crit}$ (normal limit)")

# Annotate convergence
ax.annotate(f"Converges to $z_c={z_crit}$",
            xy=(40, z_crit + 0.01), xytext=(25, 2.5),
            fontsize=9, color="#E53935",
            arrowprops=dict(arrowstyle="->", color="#E53935"))

# Annotate some specific values
for nu_ann in [1, 3, 5, 10, 30]:
    tc = stats.t.ppf(0.975, df=nu_ann)
    ax.annotate(f"$t_c={tc:.2f}$",
                xy=(nu_ann, tc), xytext=(nu_ann + 1.5, tc + 0.18),
                fontsize=7.5, color="#1976D2")

ax.set_xlabel("Degrees of freedom $\\nu$", fontsize=12)
ax.set_ylabel("Critical value $t_c$ (two-tailed, $\\alpha=0.05$)", fontsize=11)
ax.set_title("Critical Value $t_c$ vs. Degrees of Freedom\n"
             "Convergence to $z_c = 1.96$", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(0, 51)
ax.set_ylim(1.5, 14)
ax.set_yscale("log")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "tstat_vs_zstat.png", dpi=150, bbox_inches="tight")
plt.show()

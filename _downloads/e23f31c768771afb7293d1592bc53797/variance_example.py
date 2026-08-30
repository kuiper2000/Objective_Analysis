"""
variance_example.py
-------------------
Illustrates sample variance for a graduate atmospheric science statistics course.

Left panel : scatter plot of 50 random N(0,1) points.  The sample mean is shown
             as a horizontal dashed line; vertical lines depict the individual
             deviations (x_i - x_bar).

Right panel: how the sample-variance estimate converges to the true variance
             (sigma^2 = 1) as sample size N grows from 5 to 500, comparing
             the biased (1/N) and unbiased (1/(N-1)) estimators.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

rng = np.random.default_rng(42)

# ── Left panel data ──────────────────────────────────────────────────────────
N_scatter = 50
x = rng.standard_normal(N_scatter)
x_bar = x.mean()

# ── Right panel data ─────────────────────────────────────────────────────────
N_values = np.arange(5, 501)
n_trials = 2000          # average over many realisations to get a smooth curve
var_biased   = np.zeros(len(N_values))
var_unbiased = np.zeros(len(N_values))

for j, n in enumerate(N_values):
    samples = rng.standard_normal((n_trials, n))
    var_biased[j]   = np.mean(samples.var(axis=1, ddof=0))
    var_unbiased[j] = np.mean(samples.var(axis=1, ddof=1))

# ── Figure ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 5))
gs  = gridspec.GridSpec(1, 2, wspace=0.35)

# --- Left panel --------------------------------------------------------------
ax0 = fig.add_subplot(gs[0])

ax0.axhline(x_bar, color="black", linestyle="--", linewidth=1.4,
            label=rf"$\bar{{x}}$ = {x_bar:.3f}")

# Shade deviations as vertical line segments
for i, xi in enumerate(x):
    color = "#2196F3" if xi >= x_bar else "#F44336"
    ax0.plot([i, i], [x_bar, xi], color=color, linewidth=0.9, alpha=0.7)

ax0.scatter(range(N_scatter), x, s=30, color="black", zorder=5)

ax0.set_xlabel("Sample index $i$", fontsize=12)
ax0.set_ylabel("Value", fontsize=12)
ax0.set_title("Deviations from the sample mean\n"
              r"$N=50$, $X_i \sim \mathcal{N}(0,1)$", fontsize=12)
ax0.legend(fontsize=11)
ax0.set_xlim(-1, N_scatter)

# Annotate variance
s2 = np.var(x, ddof=1)
ax0.text(0.97, 0.05,
         rf"$s^2$ (unbiased) $= {s2:.3f}$",
         transform=ax0.transAxes, ha="right", va="bottom",
         fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="wheat", alpha=0.6))

# --- Right panel -------------------------------------------------------------
ax1 = fig.add_subplot(gs[1])

ax1.plot(N_values, var_biased,   color="#F44336", linewidth=1.8,
         label=r"Biased $\hat{\sigma}^2 = \frac{1}{N}\sum(x_i-\bar{x})^2$")
ax1.plot(N_values, var_unbiased, color="#2196F3", linewidth=1.8,
         label=r"Unbiased $s^2 = \frac{1}{N-1}\sum(x_i-\bar{x})^2$")
ax1.axhline(1.0, color="black", linestyle="--", linewidth=1.4,
            label=r"True $\sigma^2 = 1$")

ax1.set_xlabel("Sample size $N$", fontsize=12)
ax1.set_ylabel(r"Expected variance estimate $\langle\hat{\sigma}^2\rangle$",
               fontsize=12)
ax1.set_title("Convergence of variance estimators\n"
              r"(averaged over 2000 realisations)", fontsize=12)
ax1.legend(fontsize=9, loc="upper right")
ax1.set_xlim(5, 500)
ax1.set_ylim(0.7, 1.15)

# ── Save & show ──────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "variance_example.png", dpi=150, bbox_inches="tight")
plt.show()

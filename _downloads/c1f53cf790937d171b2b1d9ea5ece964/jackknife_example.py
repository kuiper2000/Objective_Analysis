"""
jackknife_example.py
--------------------
Demonstrates jackknife variance estimation of the mean.

Left  : strip plot of N=20 samples from N(5, 2^2) with the sample mean.
Right : all N jackknife leave-one-out mean estimates as dots; the jackknife
        standard error shown as error bars; compared with the analytical
        standard error sigma/sqrt(N).
"""

import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(21)

# ── Data ─────────────────────────────────────────────────────────────────────
N       = 20
mu_true = 5.0
sigma   = 2.0

data   = rng.normal(mu_true, sigma, N)
x_bar  = data.mean()

# ── Jackknife ─────────────────────────────────────────────────────────────────
jack_means = np.array([np.delete(data, i).mean() for i in range(N)])

# Jackknife estimate of the mean
jack_mean_est = jack_means.mean()

# Jackknife variance of the estimator (of the mean)
jack_var  = (N - 1) / N * np.sum((jack_means - jack_mean_est) ** 2)
jack_se   = np.sqrt(jack_var)

# Analytical standard error
analytical_se = sigma / np.sqrt(N)
sample_se     = data.std(ddof=1) / np.sqrt(N)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(wspace=0.40)

# ─── Left: strip plot ─────────────────────────────────────────────────────────
ax = axes[0]

# Jitter x positions for visibility
rng2   = np.random.default_rng(99)
jitter = rng2.uniform(-0.12, 0.12, N)
ax.scatter(jitter, data, s=60, color="#1976D2", zorder=3,
           edgecolors="white", linewidth=0.6, label="Observations")

ax.axhline(x_bar, color="#E53935", linewidth=2.0, linestyle="--",
           label=f"Sample mean $\\bar{{x}}={x_bar:.3f}$")
ax.axhline(mu_true, color="#388E3C", linewidth=1.5, linestyle=":",
           label=f"True $\\mu={mu_true}$")

# Deviation markers
for i in range(N):
    ax.plot([jitter[i], jitter[i]], [x_bar, data[i]],
            color="#90CAF9", linewidth=0.8, zorder=2)

ax.set_xlim(-0.5, 0.5)
ax.set_xticks([])
ax.set_xlabel("")
ax.set_ylabel("Value", fontsize=12)
ax.set_title(f"Original Data\n$N={N}$, $\\mathcal{{N}}({mu_true},{sigma}^2)$",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9.5, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

# ─── Right: jackknife estimates ────────────────────────────────────────────────
ax = axes[1]

# Each jackknife estimate as a dot
ax.plot(range(1, N + 1), jack_means, "o", color="#7B1FA2",
        markersize=7, zorder=3, label="Jackknife leave-one-out means")

# Jackknife mean ± SE band
ax.axhline(jack_mean_est, color="#7B1FA2", linestyle="--", linewidth=1.5,
           label=f"Jackknife mean = {jack_mean_est:.3f}")
ax.axhspan(jack_mean_est - jack_se, jack_mean_est + jack_se,
           alpha=0.15, color="#7B1FA2",
           label=f"Jackknife SE = {jack_se:.4f}")

# Sample mean
ax.axhline(x_bar, color="#E53935", linestyle="--", linewidth=1.5,
           label=f"Sample mean $\\bar{{x}} = {x_bar:.3f}$")

# True mean
ax.axhline(mu_true, color="#388E3C", linestyle=":", linewidth=1.5,
           label=f"True $\\mu = {mu_true}$")

# Annotate comparison of SEs
text = (f"Jackknife SE = {jack_se:.4f}\n"
        f"Analytical SE = {analytical_se:.4f}\n"
        f"Sample SE = {sample_se:.4f}")
ax.text(0.97, 0.05, text,
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=9, bbox=dict(boxstyle="round", fc="lightyellow", ec="goldenrod",
                              alpha=0.8))

ax.set_xlabel("Left-out observation index $i$", fontsize=12)
ax.set_ylabel("Jackknife mean estimate", fontsize=12)
ax.set_title("Jackknife Leave-One-Out Estimates\nof the Mean",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")
ax.set_xlim(0, N + 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "jackknife_example.png", dpi=150, bbox_inches="tight")
plt.show()

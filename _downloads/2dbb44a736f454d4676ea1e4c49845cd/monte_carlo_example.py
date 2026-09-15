"""
monte_carlo_example.py
----------------------
Monte Carlo test for the maximum of N=31 samples from N(0,1).

- Run 10 000 experiments: draw 31 samples from N(0,1), record the maximum.
- Plot histogram of the maximum values.
- Mark the 95th percentile of the empirical distribution.
- Mark 2.2 with a vertical line and report its p-value (fraction with max > 2.2).
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

rng = np.random.default_rng(55)

# ── Monte Carlo ───────────────────────────────────────────────────────────────
N_samp   = 31
n_exp    = 10_000
max_vals = rng.standard_normal((n_exp, N_samp)).max(axis=1)

pct_95   = np.percentile(max_vals, 95)
threshold = 2.2
p_value  = np.mean(max_vals > threshold)

print(f"Empirical 95th percentile of max: {pct_95:.4f}")
print(f"P(max > {threshold}): {p_value:.4f}")

# ── Theoretical distribution of max of N iid N(0,1) ──────────────────────────
# CDF: F_max(x) = [Phi(x)]^N;  PDF: f_max(x) = N * phi(x) * [Phi(x)]^(N-1)
x = np.linspace(0.5, 4.5, 600)
cdf_max = stats.norm.cdf(x) ** N_samp
pdf_max = N_samp * stats.norm.pdf(x) * stats.norm.cdf(x) ** (N_samp - 1)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))

# Histogram
counts, bin_edges, patches = ax.hist(
    max_vals, bins=60, density=True,
    color="#90CAF9", edgecolor="white", linewidth=0.4,
    label=f"MC max values ({n_exp:,} experiments)")

# Shade region > 95th percentile
for patch, left in zip(patches, bin_edges[:-1]):
    if left >= pct_95:
        patch.set_facecolor("#E53935")
        patch.set_alpha(0.80)

# Theoretical PDF
ax.plot(x, pdf_max, "k-", linewidth=2.0,
        label=r"Theoretical PDF of $\max(X_1,\ldots,X_{31})$")

# 95th percentile line
ax.axvline(pct_95, color="#E53935", linestyle="--", linewidth=2.0,
           label=f"95th pct = {pct_95:.3f}")

# 2.2 sigma line
ax.axvline(threshold, color="#7B1FA2", linestyle="-", linewidth=2.0,
           label=f"$z = {threshold}$  ($p = {p_value:.4f}$)")

# Annotations
ax.text(pct_95 + 0.05, ax.get_ylim()[1] * 0.72,
        f"5% rejection\nregion",
        color="#E53935", fontsize=9.5, va="center")

ax.text(threshold - 0.07, ax.get_ylim()[1] * 0.55,
        f"$p(\\max > {threshold})$\n$= {p_value:.4f}$",
        color="#7B1FA2", fontsize=10, ha="right",
        bbox=dict(boxstyle="round", fc="lavender", alpha=0.8))

ax.set_xlabel(r"Maximum of $N=31$ draws from $\mathcal{N}(0,1)$", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title(f"Monte Carlo Test: Distribution of Sample Maximum\n"
             f"$N=31$, {n_exp:,} experiments",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper left")
ax.set_xlim(0.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "monte_carlo_example.png", dpi=150, bbox_inches="tight")
plt.show()

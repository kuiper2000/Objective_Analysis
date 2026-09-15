"""
central_limit_theorem.py
------------------------
Demonstrates the Central Limit Theorem using a lognormal parent distribution.

For N = 1, 5, and 30:
  - Draw 5000 samples of size N from Lognormal(mu_ln, sigma_ln).
  - Compute the sample mean for each draw.
  - Plot the histogram of sample means together with the theoretical N(mu, sigma/sqrt(N))
    predicted by the CLT.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

rng = np.random.default_rng(13)

# ── Lognormal parameters ──────────────────────────────────────────────────────
# Lognormal(0, 1): mu_log=0, sigma_log=1
# E[X] = exp(mu + sigma^2/2)  = exp(0.5)
# Var[X]= (exp(sigma^2)-1)*exp(2*mu+sigma^2) = (e-1)*e
mu_log    = 0.0
sigma_log = 1.0
mu_true   = np.exp(mu_log + sigma_log**2 / 2)
var_true  = (np.exp(sigma_log**2) - 1) * np.exp(2 * mu_log + sigma_log**2)
sigma_true= np.sqrt(var_true)

Ns       = [1, 5, 30]
n_boot   = 5000

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.subplots_adjust(wspace=0.38)

colors = ["#7B1FA2", "#1976D2", "#388E3C"]

for ax, N, col in zip(axes, Ns, colors):
    # Sample means
    draws    = rng.lognormal(mu_log, sigma_log, (n_boot, N))
    means    = draws.mean(axis=1)

    # CLT prediction
    mu_clt    = mu_true
    sigma_clt = sigma_true / np.sqrt(N)

    # Plot histogram
    counts, bin_edges, _ = ax.hist(
        means, bins=50, density=True,
        color=col, alpha=0.55, edgecolor="white", linewidth=0.4,
        label=f"5000 sample means")

    # Theoretical normal
    x_range = np.linspace(means.min() - 0.5, means.max() + 0.5, 500)
    pdf_clt = stats.norm.pdf(x_range, mu_clt, sigma_clt)
    ax.plot(x_range, pdf_clt, "k-", linewidth=2.2,
            label=r"$\mathcal{N}(\mu,\,\sigma/\sqrt{N})$" + f"\n$\\mu={mu_clt:.2f}$, "
                  f"$\\sigma/\\sqrt{{N}}={sigma_clt:.3f}$")

    # Actual parent distribution for N=1
    if N == 1:
        x_ln = np.linspace(1e-3, 10, 600)
        ax.plot(x_ln, stats.lognorm.pdf(x_ln, s=sigma_log),
                linestyle="--", color="#E65100", linewidth=1.8,
                label="Lognormal parent PDF")

    ax.set_xlabel(r"$\bar{X}$", fontsize=13)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(f"$N = {N}$", fontsize=14, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # KS test p-value to quantify normality
    ks_stat, ks_p = stats.kstest(
        (means - mu_clt) / sigma_clt, "norm")
    ax.text(0.97, 0.60,
            f"KS $p$-value\n$= {ks_p:.3f}$",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))

fig.suptitle("Central Limit Theorem — Lognormal parent distribution",
             fontsize=13, fontweight="bold", y=1.01)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "central_limit_theorem.png", dpi=150, bbox_inches="tight")
plt.show()

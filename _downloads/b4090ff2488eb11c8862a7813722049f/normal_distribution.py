"""
normal_distribution.py
----------------------
Three panels illustrating the standard normal distribution.

Left   : histogram of 1000 N(0,1) samples overlaid with the theoretical PDF.
Middle : two-tailed 95% CI — shade |z| > 1.96 in red.
Right  : one-tailed 5% significance level — shade z > 1.645 in red.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

rng = np.random.default_rng(0)

# ── Data ─────────────────────────────────────────────────────────────────────
N = 1000
samples = rng.standard_normal(N)
x = np.linspace(-4.5, 4.5, 800)
pdf = stats.norm.pdf(x)

mu_hat  = samples.mean()
std_hat = samples.std(ddof=1)

z_two   = 1.96      # two-tailed 5%
z_one   = 1.645     # one-tailed 5%

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.subplots_adjust(wspace=0.38)

# ─── Left: histogram + PDF ────────────────────────────────────────────────────
ax = axes[0]
ax.hist(samples, bins=35, density=True, color="#90CAF9", edgecolor="white",
        linewidth=0.5, label="1000 samples")
ax.plot(x, pdf, "k-", linewidth=2.0, label="Theoretical $\\mathcal{N}(0,1)$")
ax.axvline(mu_hat,  color="#E53935", linestyle="--", linewidth=1.5,
           label=f"Sample mean = {mu_hat:.3f}")
ax.axvline(mu_hat + std_hat, color="#43A047", linestyle=":", linewidth=1.5,
           label=f"$\\pm$ sample std = {std_hat:.3f}")
ax.axvline(mu_hat - std_hat, color="#43A047", linestyle=":", linewidth=1.5)
ax.set_xlabel("$z$", fontsize=13)
ax.set_ylabel("Density", fontsize=13)
ax.set_title("Histogram vs. Theoretical PDF\n"
             r"$N=1000$, $\mathcal{N}(0,1)$", fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")
ax.set_xlim(-4.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Middle: two-tailed ───────────────────────────────────────────────────────
ax = axes[1]
ax.plot(x, pdf, "k-", linewidth=2.0)

# Shade tails
for x_fill, side in [(x[x <= -z_two], "left"), (x[x >= z_two], "right")]:
    ax.fill_between(x_fill, stats.norm.pdf(x_fill),
                    alpha=0.65, color="#E53935",
                    label="Rejection region (2.5% each)" if side == "left" else "_")

# Shade centre
xc = x[(x >= -z_two) & (x <= z_two)]
ax.fill_between(xc, stats.norm.pdf(xc), alpha=0.25, color="#1976D2",
                label="95% in middle")

# Critical value lines
ax.axvline(-z_two, color="#E53935", linestyle="--", linewidth=1.5)
ax.axvline( z_two, color="#E53935", linestyle="--", linewidth=1.5)
ax.text(-z_two - 0.05, 0.03, f"$-z_{{c}}$\n$={-z_two}$",
        ha="right", fontsize=9, color="#E53935")
ax.text( z_two + 0.05, 0.03, f"$z_{{c}}$\n$={z_two}$",
        ha="left",  fontsize=9, color="#E53935")
ax.text(0, 0.22, "95%", ha="center", fontsize=14, color="#1976D2",
        fontweight="bold")
ax.text(-3.3, 0.065, "2.5%", ha="center", fontsize=10, color="#E53935",
        fontweight="bold")
ax.text( 3.3, 0.065, "2.5%", ha="center", fontsize=10, color="#E53935",
        fontweight="bold")

ax.set_xlabel("$z$", fontsize=13)
ax.set_ylabel("$f(z)$", fontsize=13)
ax.set_title("Two-tailed test\n"
             r"$\alpha = 0.05$, $z_c = \pm 1.96$", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(-4.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: one-tailed ────────────────────────────────────────────────────────
ax = axes[2]
ax.plot(x, pdf, "k-", linewidth=2.0)

# Shade right tail
x_tail = x[x >= z_one]
ax.fill_between(x_tail, stats.norm.pdf(x_tail),
                alpha=0.65, color="#E53935", label="Rejection region (5%)")

# Shade left of z_one
x_main = x[x <= z_one]
ax.fill_between(x_main, stats.norm.pdf(x_main),
                alpha=0.20, color="#1976D2", label="95% (one-sided)")

ax.axvline(z_one, color="#E53935", linestyle="--", linewidth=1.5)
ax.text(z_one + 0.08, 0.03, f"$z_c = {z_one}$",
        ha="left", fontsize=9.5, color="#E53935")
ax.text(3.0, 0.07, "5%", ha="center", fontsize=11,
        color="#E53935", fontweight="bold")
ax.text(-0.5, 0.22, "95%", ha="center", fontsize=14,
        color="#1976D2", fontweight="bold")

ax.set_xlabel("$z$", fontsize=13)
ax.set_ylabel("$f(z)$", fontsize=13)
ax.set_title("One-tailed test (right)\n"
             r"$\alpha = 0.05$, $z_c = 1.645$", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper left")
ax.set_xlim(-4.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "normal_distribution.png", dpi=150, bbox_inches="tight")
plt.show()

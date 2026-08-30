"""
bootstrap_example.py
--------------------
Bootstrap hypothesis test for 500 hPa geopotential height.

Setup
-----
- 200 days of Z500 drawn from N(5886, 40^2) m.
- 20 "aerosol days" randomly selected; their observed mean is ~5906 m.
- 2500 bootstrap experiments: randomly draw 20 days from the 200-day record
  and compute the mean each time.
- Plot the histogram of bootstrap means, mark the 95th percentile, the observed
  value, and shade the rejection region.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

rng = np.random.default_rng(42)

# ── Simulate data ─────────────────────────────────────────────────────────────
mu_bg, sigma_bg = 5886.0, 40.0
N_total = 200
N_aerosol = 20

Z500 = rng.normal(mu_bg, sigma_bg, N_total)

# The "aerosol days" have a prescribed observed mean of ~5906 m.
# In a real study these would come from independent observations; here we
# fix the mean directly to keep the pedagogical example reproducible.
obs_mean = 5906.0   # observed mean of the 20 aerosol days (m)

print(f"Prescribed aerosol-day mean: {obs_mean:.2f} m")

# ── Bootstrap ─────────────────────────────────────────────────────────────────
# Each bootstrap experiment draws 20 days randomly (with replacement) from the
# 200-day background record to build the null distribution of the mean.
n_boot = 2500
boot_means = np.array([
    rng.choice(Z500, size=N_aerosol, replace=True).mean()
    for _ in range(n_boot)
])

pct_95 = np.percentile(boot_means, 95)
p_value = np.mean(boot_means >= obs_mean)

print(f"Bootstrap 95th percentile: {pct_95:.2f} m")
print(f"Bootstrap p-value (one-tailed): {p_value:.4f}")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))

# Histogram
counts, bin_edges, patches = ax.hist(
    boot_means, bins=50, density=True,
    color="#90CAF9", edgecolor="white", linewidth=0.5,
    label=f"Bootstrap means ($n={n_boot}$)")

# Shade rejection region (> 95th percentile)
for patch, left in zip(patches, bin_edges[:-1]):
    if left >= pct_95:
        patch.set_facecolor("#E53935")
        patch.set_alpha(0.80)

# Mark 95th percentile
ax.axvline(pct_95, color="#E53935", linestyle="--", linewidth=2.0,
           label=f"95th percentile = {pct_95:.1f} m")

# Mark observed value
ax.axvline(obs_mean, color="#4A148C", linestyle="-", linewidth=2.2,
           label=f"Observed mean = {obs_mean:.1f} m")

# Overlay theoretical null distribution (mean of 20 from the background)
x_range = np.linspace(boot_means.min() - 20, boot_means.max() + 20, 500)
sigma_null = sigma_bg / np.sqrt(N_aerosol)
pdf_null   = stats.norm.pdf(x_range, mu_bg, sigma_null)
ax.plot(x_range, pdf_null, "k:", linewidth=1.8,
        label=r"Theoretical $\mathcal{N}(\mu_{bg},\,\sigma_{bg}/\sqrt{N})$")

# Annotations
ax.text(pct_95 + 3, ax.get_ylim()[1] * 0.55,
        f"Rejection\nregion\n($\\alpha=0.05$)",
        color="#E53935", fontsize=9.5, va="center")

verdict = "Reject $H_0$" if obs_mean > pct_95 else "Fail to reject $H_0$"
ax.text(0.02, 0.96,
        f"$p$-value = {p_value:.4f}\n{verdict}",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=10, bbox=dict(boxstyle="round", fc="lightyellow",
                               ec="goldenrod", alpha=0.85))

ax.set_xlabel("500 hPa geopotential height mean (m)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("Bootstrap Test: 500 hPa Geopotential Height\n"
             f"$N_{{bg}}={N_total}$ days, $N_{{aerosol}}={N_aerosol}$ days, "
             f"{n_boot} bootstrap samples",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "bootstrap_example.png", dpi=150, bbox_inches="tight")
plt.show()

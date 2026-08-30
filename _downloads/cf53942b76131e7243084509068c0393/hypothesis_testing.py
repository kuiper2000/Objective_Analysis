"""
hypothesis_testing.py
---------------------
Left  : two-tailed hypothesis test.  Standard normal PDF; rejection regions
        (|z| > 1.96) shaded red; example test statistic z_obs = 2.3 shown;
        p-value region annotated.
Right : one-tailed (right tail) test.  Rejection region z > 1.645 shaded red;
        same z_obs shown.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

x = np.linspace(-4.5, 4.5, 900)
pdf = stats.norm.pdf(x)

z_obs      = 2.3
z_two      = 1.96
z_one      = 1.645
p_two      = 2 * (1 - stats.norm.cdf(abs(z_obs)))   # two-tailed p-value
p_one      = 1 - stats.norm.cdf(z_obs)               # one-tailed p-value

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.subplots_adjust(wspace=0.40)

# ─── Left: two-tailed ─────────────────────────────────────────────────────────
ax = axes[0]
ax.plot(x, pdf, "k-", linewidth=2.0)

# Rejection regions
x_left  = x[x <= -z_two]
x_right = x[x >= z_two]
ax.fill_between(x_left,  stats.norm.pdf(x_left),  alpha=0.70,
                color="#E53935", label=r"Rejection region ($\alpha=0.05$)")
ax.fill_between(x_right, stats.norm.pdf(x_right), alpha=0.70,
                color="#E53935")

# p-value region (beyond |z_obs| in both tails)
x_pval_r = x[x >= z_obs]
x_pval_l = x[x <= -z_obs]
ax.fill_between(x_pval_r, stats.norm.pdf(x_pval_r), alpha=0.55,
                color="#AB47BC", label=f"$p$-value region: $p={p_two:.4f}$")
ax.fill_between(x_pval_l, stats.norm.pdf(x_pval_l), alpha=0.55,
                color="#AB47BC")

# Non-rejection region
x_mid = x[(x >= -z_two) & (x <= z_two)]
ax.fill_between(x_mid, stats.norm.pdf(x_mid), alpha=0.15,
                color="#1976D2", label="Non-rejection region (95%)")

# Critical value lines
ax.axvline(-z_two, color="#E53935", linestyle="--", linewidth=1.5)
ax.axvline( z_two, color="#E53935", linestyle="--", linewidth=1.5)
ax.text(-z_two - 0.06, 0.025, f"$-z_c={-z_two}$",
        ha="right", fontsize=9, color="#E53935")
ax.text( z_two + 0.06, 0.025, f"$z_c={z_two}$",
        ha="left", fontsize=9, color="#E53935")

# Test statistic
ax.axvline(z_obs, color="#4A148C", linewidth=2.2, linestyle="-",
           label=f"Test statistic $z_{{obs}}={z_obs}$")
ax.text(z_obs + 0.08, 0.19, f"$z_{{obs}}={z_obs}$\n$p={p_two:.4f}$",
        ha="left", fontsize=9, color="#4A148C",
        bbox=dict(boxstyle="round", fc="lavender", alpha=0.8))

ax.set_xlabel("$z$", fontsize=13)
ax.set_ylabel("$f(z)$", fontsize=13)
ax.set_title("Two-tailed Hypothesis Test\n"
             r"$H_0: \mu=0$, $\alpha=0.05$, $z_c = \pm1.96$",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper left")
ax.set_xlim(-4.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: one-tailed ────────────────────────────────────────────────────────
ax = axes[1]
ax.plot(x, pdf, "k-", linewidth=2.0)

# Rejection region
x_rej = x[x >= z_one]
ax.fill_between(x_rej, stats.norm.pdf(x_rej), alpha=0.70,
                color="#E53935", label=r"Rejection region ($\alpha=0.05$)")

# p-value region
x_pval = x[x >= z_obs]
ax.fill_between(x_pval, stats.norm.pdf(x_pval), alpha=0.55,
                color="#AB47BC", label=f"$p$-value region: $p={p_one:.4f}$")

# Non-rejection
x_main = x[x <= z_one]
ax.fill_between(x_main, stats.norm.pdf(x_main), alpha=0.15,
                color="#1976D2", label="Non-rejection region (95%)")

# Critical value line
ax.axvline(z_one, color="#E53935", linestyle="--", linewidth=1.5)
ax.text(z_one + 0.08, 0.025, f"$z_c={z_one}$",
        ha="left", fontsize=9, color="#E53935")

# Test statistic
ax.axvline(z_obs, color="#4A148C", linewidth=2.2, linestyle="-",
           label=f"Test statistic $z_{{obs}}={z_obs}$")
ax.text(z_obs + 0.08, 0.19, f"$z_{{obs}}={z_obs}$\n$p={p_one:.4f}$",
        ha="left", fontsize=9, color="#4A148C",
        bbox=dict(boxstyle="round", fc="lavender", alpha=0.8))

ax.text(3.5, 0.06, "5%", ha="center", fontsize=10,
        color="#E53935", fontweight="bold")
ax.text(0, 0.22, "95%", ha="center", fontsize=13,
        color="#1976D2", fontweight="bold")

ax.set_xlabel("$z$", fontsize=13)
ax.set_ylabel("$f(z)$", fontsize=13)
ax.set_title("One-tailed Hypothesis Test (right)\n"
             r"$H_0: \mu\leq0$, $\alpha=0.05$, $z_c = 1.645$",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper left")
ax.set_xlim(-4.5, 4.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "hypothesis_testing.png", dpi=150, bbox_inches="tight")
plt.show()

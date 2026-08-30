"""
skewness_kurtosis.py
--------------------
Illustrates skewness and kurtosis using probability density functions.

Left panel  : standard normal (skewness=0), a right-skewed lognormal, and a
              left-skewed (reflected lognormal) distribution.
Right panel : standard normal (excess kurtosis=0), a leptokurtic t(3), and a
              platykurtic uniform distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

x = np.linspace(-5, 5, 1000)

# ── Skewness distributions ───────────────────────────────────────────────────
# Standard normal  (skewness = 0)
pdf_norm = stats.norm.pdf(x)

# Right-skewed: lognormal shifted so that its mode is near 0
# lognormal with sigma=0.8 -> skewness = (e^s^2+2)*sqrt(e^s^2-1) > 0
ln_sigma = 0.8
# Plot over positive axis; shift to centre the bulk near 0
x_ln = np.linspace(0.001, 6, 1000)
x_ln_shifted = x_ln - 2.0        # shift for visual alignment
pdf_rskew = stats.lognorm.pdf(x_ln, s=ln_sigma)

# Left-skewed: mirror of lognormal
x_lskew = -x_ln + 2.0
pdf_lskew = stats.lognorm.pdf(x_ln, s=ln_sigma)

# Compute moments for annotations
skew_norm  = stats.norm.stats(moments="s")
skew_right = stats.lognorm.stats(s=ln_sigma, moments="s")
skew_left  = -float(skew_right)   # mirror

# ── Kurtosis distributions ───────────────────────────────────────────────────
# Standard normal excess kurtosis = 0
pdf_norm_k = stats.norm.pdf(x)

# Leptokurtic: t-distribution with nu=3
nu = 3
pdf_lepto = stats.t.pdf(x, df=nu)
kurt_lepto = stats.t.stats(df=nu, moments="k")   # excess kurtosis

# Platykurtic: uniform on [-sqrt(3), sqrt(3)] so variance=1
a_u = -np.sqrt(3)
b_u =  np.sqrt(3)
x_unif = np.linspace(-4, 4, 1000)
pdf_platy = stats.uniform.pdf(x_unif, loc=a_u, scale=b_u - a_u)
kurt_unif = stats.uniform.stats(moments="k")     # excess kurtosis = -1.2

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.subplots_adjust(wspace=0.35)

# --- Left: Skewness ----------------------------------------------------------
ax = axes[0]
ax.plot(x, pdf_norm, "k-",    linewidth=2.2,
        label=f"Normal  (skewness = 0)")
ax.plot(x_ln_shifted, pdf_rskew, color="#E65100", linewidth=2.2,
        label=f"Lognormal right-skewed\n(skewness $\\approx$ {float(skew_right):.2f})")
ax.plot(x_lskew, pdf_lskew, color="#1565C0", linewidth=2.2,
        label=f"Reflected lognormal left-skewed\n(skewness $\\approx$ {skew_left:.2f})")

ax.set_xlim(-4.5, 4.5)
ax.set_ylim(0, 0.75)
ax.set_xlabel("$x$", fontsize=13)
ax.set_ylabel("Probability density $f(x)$", fontsize=12)
ax.set_title("Skewness", fontsize=14, fontweight="bold")
ax.legend(fontsize=9.5, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Annotate tail direction
ax.annotate("heavy right tail", xy=(3.2, 0.06), fontsize=8.5,
            color="#E65100", ha="left")
ax.annotate("heavy left tail",  xy=(-4.0, 0.06), fontsize=8.5,
            color="#1565C0", ha="left")

# --- Right: Kurtosis ---------------------------------------------------------
ax = axes[1]

# Normal
ax.plot(x, pdf_norm_k, "k-", linewidth=2.2,
        label=r"Normal — excess kurtosis = 0")

# Leptokurtic
ax.plot(x, pdf_lepto, color="#C62828", linewidth=2.2,
        label=rf"$t$-dist ($\nu={nu}$) — leptokurtic"
              f"\nexcess kurtosis = {float(kurt_lepto):.1f} (if $\\nu>4$) / $\\infty$ else")

# Platykurtic
ax.plot(x_unif, pdf_platy, color="#1B5E20", linewidth=2.2,
        label=f"Uniform$[-\\sqrt{{3}},\\sqrt{{3}}]$ — platykurtic\n"
              f"excess kurtosis = {float(kurt_unif):.2f}")

ax.set_xlim(-4.5, 4.5)
ax.set_ylim(0, 0.52)
ax.set_xlabel("$x$", fontsize=13)
ax.set_ylabel("Probability density $f(x)$", fontsize=12)
ax.set_title("Kurtosis", fontsize=14, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Shade tails of t to highlight heavy tails
x_tail = np.linspace(-4.5, -2.5, 200)
ax.fill_between(x_tail, stats.t.pdf(x_tail, df=nu),
                stats.norm.pdf(x_tail), alpha=0.25, color="#C62828",
                label="_nolegend_")
x_tail2 = np.linspace(2.5, 4.5, 200)
ax.fill_between(x_tail2, stats.t.pdf(x_tail2, df=nu),
                stats.norm.pdf(x_tail2), alpha=0.25, color="#C62828")

# ── Save & show ──────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "skewness_kurtosis.png", dpi=150, bbox_inches="tight")
plt.show()

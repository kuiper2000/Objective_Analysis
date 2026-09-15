"""
what_is_pdf.py
--------------
Introduces the probability density function (PDF) and CDF.

Left   : Gaussian PDF f(x), mu=0, sigma=1.  Shade P(-1 ≤ x ≤ 1) = 68.27%.
Middle : CDF F(x) for the same distribution.  Mark F(1) and F(-1) with dashed
         lines and show F(1) - F(-1) = 68.27%.
Right  : Uniform PDF on [-2, 2].  Shade area between 0 and 1 and label P.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

x = np.linspace(-4, 4, 1000)
mu, sigma = 0.0, 1.0

pdf = stats.norm.pdf(x, mu, sigma)
cdf = stats.norm.cdf(x, mu, sigma)

# Shaded region [-1, 1]
xa  = np.linspace(-1, 1, 400)
pxa = stats.norm.pdf(xa, mu, sigma)
prob_inner = stats.norm.cdf(1, mu, sigma) - stats.norm.cdf(-1, mu, sigma)

# Uniform distribution on [-2, 2]
a_u, b_u = -2.0, 2.0
pdf_u = np.where((x >= a_u) & (x <= b_u), 1.0 / (b_u - a_u), 0.0)
x_shade_u = np.linspace(0, 1, 300)
prob_u = (1.0 - 0.0) / (b_u - a_u)    # = 0.25

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.38)

# ─── Left: Gaussian PDF ──────────────────────────────────────────────────────
ax = axes[0]
ax.plot(x, pdf, "k-", linewidth=2)
ax.fill_between(xa, pxa, alpha=0.45, color="#1976D2",
                label=f"$P(-1 \\leq x \\leq 1)={prob_inner*100:.2f}\\%$")
ax.axvline(-1, color="#1976D2", linestyle="--", linewidth=1.2)
ax.axvline( 1, color="#1976D2", linestyle="--", linewidth=1.2)
ax.annotate("$-1\\sigma$", xy=(-1, 0), xytext=(-1.5, 0.08),
            fontsize=10, color="#1976D2",
            arrowprops=dict(arrowstyle="->", color="#1976D2"))
ax.annotate("$+1\\sigma$", xy=(1, 0), xytext=(1.2, 0.08),
            fontsize=10, color="#1976D2",
            arrowprops=dict(arrowstyle="->", color="#1976D2"))
ax.set_xlabel("$x$", fontsize=13)
ax.set_ylabel("$f(x)$", fontsize=13)
ax.set_title("Gaussian PDF\n"
             r"$\mu=0,\ \sigma=1$", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(-4, 4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Middle: CDF ─────────────────────────────────────────────────────────────
ax = axes[1]
ax.plot(x, cdf, "k-", linewidth=2, label="$F(x)$")

f_m1 = stats.norm.cdf(-1, mu, sigma)
f_p1 = stats.norm.cdf( 1, mu, sigma)

# Dashed guide lines
for xv, fv, col, lbl in [(-1, f_m1, "#E53935", f"$F(-1)={f_m1:.4f}$"),
                           ( 1, f_p1, "#43A047", f"$F(1)={f_p1:.4f}$")]:
    ax.plot([xv, xv], [0, fv], linestyle="--", color=col, linewidth=1.4)
    ax.plot([-4, xv], [fv, fv], linestyle="--", color=col, linewidth=1.4,
            label=lbl)

# Brace annotation for the difference
ax.annotate("", xy=(-4.3, f_p1), xytext=(-4.3, f_m1),
            arrowprops=dict(arrowstyle="<->", color="#1976D2", lw=1.8))
ax.text(-3.9, (f_m1 + f_p1) / 2,
        f"$F(1)-F(-1)$\n$={prob_inner*100:.2f}\\%$",
        fontsize=9, color="#1976D2", va="center")

ax.set_xlabel("$x$", fontsize=13)
ax.set_ylabel("$F(x) = P(X \\leq x)$", fontsize=12)
ax.set_title("Cumulative Distribution Function (CDF)\n"
             r"$\mu=0,\ \sigma=1$", fontsize=12, fontweight="bold")
ax.legend(fontsize=9.5, loc="upper left")
ax.set_xlim(-4, 4)
ax.set_ylim(-0.05, 1.05)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: Uniform PDF ──────────────────────────────────────────────────────
ax = axes[2]
ax.plot(x, pdf_u, "k-", linewidth=2, label="Uniform$[-2,2]$")
pdf_u_shade = np.full_like(x_shade_u, 1.0 / (b_u - a_u))
ax.fill_between(x_shade_u, pdf_u_shade, alpha=0.5, color="#FB8C00",
                label=f"$P(0 \\leq x \\leq 1) = {prob_u:.2f}$")
ax.axvline(0, color="#FB8C00", linestyle="--", linewidth=1.2)
ax.axvline(1, color="#FB8C00", linestyle="--", linewidth=1.2)
ax.text(0.5, 0.22,
        f"$P={prob_u:.2f}$",
        ha="center", fontsize=11, color="#E65100", fontweight="bold")
ax.set_xlabel("$x$", fontsize=13)
ax.set_ylabel("$f(x)$", fontsize=13)
ax.set_title("Uniform PDF on $[-2, 2]$", fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(-3, 3)
ax.set_ylim(0, 0.45)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "what_is_pdf.png", dpi=150, bbox_inches="tight")
plt.show()

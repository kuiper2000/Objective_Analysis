"""
bayesian_sst_example.py
-----------------------
Bayesian SST estimation example.

Model
-----
  X ~ Uniform(-2, 2)          prior on true SST anomaly
  W ~ N(0, 0.75^2)            observation noise
  Y = X + W                   observed value

Three panels
------------
Left   : PDFs of X (prior), W (noise), and Y (marginal, via convolution).
Middle : Posterior PDF f(X | Y = 2.1) computed numerically via Bayes theorem.
         Mark the posterior mean E[X | Y=2.1].
Right  : Conditional expectation E[X | Y=y] as a function of y in [-4, 4],
         showing saturation near ±2.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# ── Grid setup ────────────────────────────────────────────────────────────────
a, b      = -2.0, 2.0          # uniform support
sigma_w   = 0.75               # noise std

dx        = 0.005
x_grid    = np.arange(-4, 4 + dx, dx)   # fine grid for computations

# ── Prior: f_X(x) ─────────────────────────────────────────────────────────────
def f_X(x):
    return np.where((x >= a) & (x <= b), 1.0 / (b - a), 0.0)

# ── Noise: f_W(w) ──────────────────────────────────────────────────────────────
def f_W(w):
    return stats.norm.pdf(w, 0, sigma_w)

# ── Marginal of Y via convolution: f_Y(y) = integral f_X(x) f_W(y-x) dx ─────
# Computed numerically on the fine grid
y_grid = x_grid.copy()
fX_vals = f_X(x_grid)
fW_vals = f_W(x_grid)
# Convolution: (f_X * f_W)(y) — use full convolution then trim
conv_full = np.convolve(fX_vals, fW_vals, mode="full") * dx
# The output grid spans [x_grid[0]+x_grid[0], x_grid[-1]+x_grid[-1]]
x_start = x_grid[0] + x_grid[0]
y_conv_grid = np.linspace(x_start, -x_start, len(conv_full))

def f_Y_at(y_val):
    """Evaluate f_Y at a scalar or array y_val by interpolation."""
    return np.interp(y_val, y_conv_grid, conv_full, left=0.0, right=0.0)

# ── Likelihood: f_{Y|X}(y | x) = f_W(y - x) ─────────────────────────────────
def likelihood(y_obs, x):
    return stats.norm.pdf(y_obs - x, 0, sigma_w)

# ── Posterior: f_{X|Y}(x | y_obs) via Bayes theorem ──────────────────────────
def posterior_pdf(y_obs):
    """Returns (x_values, posterior_pdf_values) for x in [a,b] ± 2*sigma."""
    x_range = np.linspace(a - 2 * sigma_w, b + 2 * sigma_w, 2000)
    lik      = likelihood(y_obs, x_range)
    prior    = f_X(x_range)
    unnorm   = lik * prior
    norm_const = np.trapezoid(unnorm, x_range)
    return x_range, unnorm / norm_const

# ── Conditional expectation E[X | Y=y] ────────────────────────────────────────
y_range = np.linspace(-4, 4, 300)
E_X_given_Y = np.zeros_like(y_range)
for j, yv in enumerate(y_range):
    x_r, post = posterior_pdf(yv)
    E_X_given_Y[j] = np.trapezoid(x_r * post, x_r)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.subplots_adjust(wspace=0.40)

# ─── Left: Marginal PDFs ──────────────────────────────────────────────────────
ax = axes[0]
x_plot = np.linspace(-4, 4, 800)

ax.plot(x_plot, f_X(x_plot), color="#1976D2", linewidth=2.2,
        label=r"$f_X(x)$ — Uniform$[-2,2]$ (prior)")
ax.plot(x_plot, f_W(x_plot), color="#388E3C", linewidth=2.2,
        label=r"$f_W(w)$ — $\mathcal{N}(0,0.75^2)$ (noise)")
ax.plot(x_plot, f_Y_at(x_plot), color="#E65100", linewidth=2.2,
        label=r"$f_Y(y)$ — marginal (convolution)")

ax.set_xlabel("Value", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("Prior, Noise, and Marginal PDFs\n"
             r"$Y = X + W$", fontsize=12, fontweight="bold")
ax.legend(fontsize=9.5, loc="upper right")
ax.set_xlim(-4, 4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Middle: Posterior ────────────────────────────────────────────────────────
ax = axes[1]
y_obs = 2.1
x_post, post = posterior_pdf(y_obs)
E_post = np.trapezoid(x_post * post, x_post)

ax.plot(x_post, post, color="#7B1FA2", linewidth=2.2,
        label=r"$f_{X|Y}(x\,|\,Y=2.1)$")
ax.axvline(E_post, color="#E53935", linestyle="--", linewidth=1.8,
           label=f"Posterior mean $E[X|Y=2.1]={E_post:.3f}$")
ax.fill_between(x_post, post, alpha=0.20, color="#7B1FA2")

# Mark prior support edges
ax.axvline(a, color="#1976D2", linestyle=":", linewidth=1.2,
           label=f"Prior bounds $x={a},{b}$")
ax.axvline(b, color="#1976D2", linestyle=":", linewidth=1.2)

ax.set_xlabel("$x$ (true SST anomaly)", fontsize=12)
ax.set_ylabel(r"$f_{X|Y}(x\,|\,Y=2.1)$", fontsize=12)
ax.set_title(r"Posterior PDF $f(X|Y=2.1)$" + "\n(Bayes theorem)",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9.5, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ─── Right: E[X|Y=y] vs y ────────────────────────────────────────────────────
ax = axes[2]
ax.plot(y_range, E_X_given_Y, color="#E65100", linewidth=2.2,
        label=r"$E[X\,|\,Y=y]$")
ax.axhline(a, color="#1976D2", linestyle="--", linewidth=1.2,
           label=f"Prior bounds $x=\\pm{abs(a)}$")
ax.axhline(b, color="#1976D2", linestyle="--", linewidth=1.2)
ax.plot([0, 0], [E_X_given_Y.min(), E_X_given_Y.max()],
        color="gray", linestyle=":", linewidth=1.0)

# Identity line for reference
ax.plot(y_range, y_range, "k:", linewidth=1.0, alpha=0.4,
        label="$E[X|Y]=y$ (identity, no noise)")

ax.set_xlabel("Observed value $y$", fontsize=12)
ax.set_ylabel(r"$E[X\,|\,Y=y]$", fontsize=12)
ax.set_title("Conditional Expectation\n"
             r"$E[X\,|\,Y=y]$ — saturation near $\pm 2$",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=10, loc="upper left")
ax.set_xlim(-4, 4)
ax.set_ylim(-2.2, 2.2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Annotate saturation
ax.annotate("Saturation at $+2$",
            xy=(3.5, 1.92), xytext=(1.5, 2.1),
            fontsize=8.5, color="#E65100",
            arrowprops=dict(arrowstyle="->", color="#E65100"))
ax.annotate("Saturation at $-2$",
            xy=(-3.5, -1.92), xytext=(-3.8, -1.5),
            fontsize=8.5, color="#E65100",
            arrowprops=dict(arrowstyle="->", color="#E65100"))

# ── Save & show ───────────────────────────────────────────────────────────────
outdir = "/Users/kaichiht/kuiper2000@gmail.com - Google Drive/My Drive/website-hugo/Objective_Analysis/General_statistical_test/"
fig.savefig(outdir + "bayesian_sst_example.png", dpi=150, bbox_inches="tight")
plt.show()

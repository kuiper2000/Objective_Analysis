(Regression)=
# Week 5-8: Regression & AR1

## Linear Regression

In linear regression, the goal is to determine
- the linear fit of $X$ and $Y$ (the regression coefficient)
- the robustness of the fit (the correlation coefficient)

Regression is simple but powerful, however, this also makes it easily misused. In a sense, the entire class (EOFs; Fourier analysis) is all based on regression.

How do we find the slope and y-intercept of the line that best fits the observed data? Let's assume $x(t)$ and $y(t)$ are time series sampled at $N$ time steps (so that each point represents a time step).

First, we have to define what "best-fit" means. For now, we will use the conventional definition which means that we want to reduce the sum of the squared errors of $y$.

Using the method of least squares:

$$\hat{y}(t) = a_1 x(t) + a_0$$

where
- $\hat{y}(t)$ denotes the estimate of $y(t)$ based on the linear relationship with $x(t)$
- $a_1$ denotes the slope, a.k.a. the regression coefficient
- $a_0$ denotes the y-intercept

Define the error of the fit as the sum of squares of the $y(\text{estimate}) - y(\text{actual})$:

$$Q = \sum_{i=1}^{N}(\hat{y}_i - y_i)^2 = \sum_{i=1}^{N}(a_1 x_i + a_0 - y_i)^2$$

where the subscript $i$ denotes the time step.

The error is squared so that
- the error is positive definite (don't want positive and negative errors canceling out)
- the minimization of $Q$ (the derivative of $Q$) is a linear problem

Note: the square causes larger errors to be more heavily weighted.

We now follow steps from our college Calculus I course and find the $a_1$ and $a_0$ that minimize $Q$ (sometimes called the cost function):

$$\begin{aligned}
\frac{dQ}{da_0} &= 0 \\
0 &= 2\sum_{i=1}^{N}(a_1 x_i + a_0 - y_i) \\
0 &= a_1\sum_{i=1}^{N} x_i + a_0 N - \sum_{i=1}^{N} y_i
\end{aligned}$$

$$\begin{aligned}
\frac{dQ}{da_1} &= 0 \\
0 &= 2\sum_{i=1}^{N}(a_1 x_i + a_0 - y_i)x_i \\
0 &= a_1\sum_{i=1}^{N} x_i^2 + a_0\sum_{i=1}^{N} x_i - \sum_{i=1}^{N} x_i y_i
\end{aligned}$$

Divide through by $N$ and move the $y$ terms to the left-hand side, where overbars denote the mean and primes denote departures from the mean:

$$\begin{aligned}
\overline{y} &= a_1 \overline{x} + a_0 \\
\overline{xy} &= a_1 \overline{x^2} + a_0 \overline{x}
\end{aligned}$$

Two equations, two unknowns. The solutions are:

$$a_1 = \frac{\overline{xy} - \overline{x}\cdot\overline{y}}{\overline{x^2} - \overline{x}^2}$$

Note that,

$$\overline{xy} = \overline{x}\cdot\overline{y} + \overline{x'y'} \qquad \text{and} \qquad \overline{x^2} = \overline{x}^2 + \overline{x'^2}$$

Hence,

$$\boxed{a_1 = \frac{\overline{x'y'}}{\overline{x'^2}}, \qquad a_0 = \overline{y} - a_1\overline{x}}$$

**$a_1$** (regression coefficient / slope):
- slope of the best-fit line
- equal to the covariance of $x$ and $y$ divided by the variance of $x$

**$a_0$** (y-intercept):
- note that if the means of the time series are 0 (they are anomalies), then $a_0 = 0$

One can put confidence limits on the slope $a_1$ in a number of ways. For example, a jackknife approach can be used to determine the sensitivity to removing a single point. Alternatively, the standard error $\sigma_{a_1}$ of the slope is given by

$$\sigma_{a_1}^2 = \frac{\dfrac{1}{N-2}\sum_{i=1}^{N}(y_i - \hat{y}_i)^2}{\sum_{i=1}^{N}(x_i - \overline{x})^2}$$

where we assume that $x$ is known exactly. Intuitively, this is the error in our $y$ estimate divided by the variance in our $x$ values.

The confidence interval for the true slope $b$ is then:

$$a_1 - t_{N-2,\alpha}\cdot\sigma_{a_1} < b < a_1 + t_{N-2,\alpha}\cdot\sigma_{a_1}$$

The $N-2$ comes from the fact that two degrees of freedom were used to estimate $a_1$ and $a_0$.

### How good is the fit?

How much we "believe" the regression coefficient ($a_1$) depends on the spread of the dots about the best-fit line. If the dots are closely packed about the regression line, then the fit is good. The spread of the dots is given by the correlation coefficient $r$.

By definition, the total variance of $y(t)$ is $\frac{1}{N}\sum_{i=1}^{N}(y_i - \overline{y})^2$, and the total variance of the fit $\hat{y}(t)$ is $\frac{1}{N}\sum_{i=1}^{N}(\hat{y}_i - \overline{y})^2$, where we use the fact that $\overline{\hat{y}} = a_1\overline{x} + a_0 = \overline{y}$.

The percent of the total variance in $y(t)$ explained by the fit $\hat{y}(t)$:

$$\begin{aligned}
r^2 &= \frac{\text{explained variance}}{\text{total variance}}
= \frac{\sum_{i=1}^{N}(\hat{y}_i - \overline{y})^2}{\sum_{i=1}^{N}(y_i - \overline{y})^2}
= \frac{\sum_{i=1}^{N}(a_1 x_i')^2}{\sum_{i=1}^{N}(y_i'^2)} \\
&= \left(\frac{\overline{x'y'}}{\overline{x'^2}}\right)^2 \frac{\sum_{i=1}^{N}(x_i')^2}{\sum_{i=1}^{N}(y_i'^2)}
= \frac{(\overline{x'y'})^2}{\overline{x'^2}\cdot\overline{y'^2}}
\end{aligned}$$

Hence:

$$\boxed{r = \frac{\overline{x'y'}}{\sigma_x\sigma_y}}$$

where $\sigma_x = (\overline{x'^2})^{1/2}$.

- $r^2$ is the fraction of variance explained by the linear least-squares fit; it always lies between 0 and 1
- $r$ varies between $-1$ and $1$

Relationships between $r$ and $r^2$:

| $r$ | $r^2$ |
|-----|-------|
| 0.99 | 0.98 |
| 0.90 | 0.81 |
| 0.70 | 0.49 |
| 0.50 | 0.25 |
| 0.25 | 0.06 |

### Relationship between the slope and the correlation coefficient

Since $a_1 = \overline{x'y'}/\overline{x'^2}$, it follows that:

$$a_1 = r\frac{\sigma_y}{\sigma_x}$$

- the regression coefficient can be thought of as the correlation coefficient multiplied by the ratio of the standard deviations of $y$ and $x$
- regression coefficients give information about the correlation coefficient and the relative amplitudes of variations of $y$ and $x$
- in the special case where $x$ and $y$ are standardized, the correlation coefficient and the regression coefficient are equal

### General comments on linear regression

- only works for linear relationships
- does not reveal relationships that are lagged or out of phase
- need to be careful about estimating the true sample size (more on this later)
- correlation does **NOT** reveal cause and effect
- flipping $x$ and $y$ will not give the same results — it is very important to physically justify your choice of $x$ and $y$!

### Orthogonal Least Squares Regression

If you can't justify which of your data is dependent and which is independent, **orthogonal least squares** may actually be what you want. In this case, you minimize the orthogonal (perpendicular) distances to the fit line, rather than the vertical distances.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: Orthogonal Least Squares vs Ordinary Least Squares* — a diagram comparing two fits: (left) **ordinary least squares**, minimizing *vertical offsets*; (right) **orthogonal least squares**, minimizing *perpendicular offsets*. The two approaches can give different answers.
:::

It just so happens that in 2-dimensions, EOF analysis (to be discussed later) gives you the orthogonal least squares fit. So, in a few weeks, you will be capable of calculating this too.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: `LSQ_OLS.py`*
:::

### Filtering with linear regression

Consider the decomposition of a variable $y$ into a fraction that is linearly congruent with $x$ and the fraction uncorrelated with $x$:

$$y(t) = y(t)_{\text{fitted}} + y(t)_{\text{residual}}$$

The fit of $y(t)$ from $x(t)$ is:

$$y(t)_{\text{fitted}} = a_1 x(t) + a_0, \qquad a_1 = \frac{\overline{x'y'}}{\overline{x'^2}} = r\frac{\sigma_y}{\sigma_x}$$

If the means of $y(t)$ and $x(t)$ are zero:

$$y(t)_{\text{residual}} = y(t) - \frac{\overline{x'y'}}{\overline{x'^2}}\cdot x(t) = y(t) - r\frac{\sigma_y}{\sigma_x}\cdot x(t)$$

- $y(t)_{\text{fitted}}$ represents the LSQ fit of $x(t)$ to $y(t)$
- by construction, $y(t)_{\text{residual}}$ is uncorrelated with $x(t)$
- the fraction of variance of $y(t)$ explained by $x(t)$ is $r^2$
- the fraction of variance of $y(t)$ not explained by $x(t)$ is $1 - r^2$

### Signal-to-Noise Ratio and Nonlinear Trend Detection

#### Signal and noise in the regression framework

Linear regression provides a natural decomposition of a time series into a **signal** (the fitted trend) and **noise** (the residuals):

$$y(t) = \underbrace{f(t)}_{\text{signal}} + \underbrace{\eta(t)}_{\text{noise}}$$

For a linear trend $f(t) = a_1 t + a_0$, the signal and noise variances are:

$$\sigma_{\text{signal}}^2 = \frac{1}{N}\sum_{i=1}^{N}(\hat{y}_i - \overline{y})^2 = r^2\,\sigma_y^2, \qquad \sigma_{\text{noise}}^2 = (1-r^2)\,\sigma_y^2$$

so the **signal-to-noise ratio** is:

$$\boxed{\text{SNR} = \frac{\sigma_{\text{signal}}^2}{\sigma_{\text{noise}}^2} = \frac{r^2}{1-r^2}}$$

This directly links the familiar correlation coefficient $r$ to the detectability of a trend:

| $r$ | SNR |
|-----|-----|
| 0.99 | 49 |
| 0.70 | 0.96 |
| 0.50 | 0.33 |
| 0.25 | 0.07 |

When SNR $\gg 1$, the forced signal dominates and is easy to detect. When SNR $\ll 1$, internal variability overwhelms the signal.

#### Testing whether a linear trend is significant

The significance of the slope $a_1$ is a direct application of our earlier t-test. Under $H_0: a_1 = 0$, the test statistic is:

$$t = \frac{a_1}{\sigma_{a_1}} \sim t_{N^*-2}$$

where $\sigma_{a_1}$ is the standard error of the slope (see above) and **$N^*$ is the effective sample size** (Leith formula), not $N$, to account for autocorrelation in the residuals. This is a common mistake — if the residuals are red noise with lag-1 autocorrelation $\alpha$, using $N$ instead of $N^*$ inflates the t-statistic and leads to spurious detections.

#### Why a linear fit can fail: the nonlinear forced response problem

In climate science, the forced response to greenhouse warming is often **not linear in time** — it can accelerate or decelerate depending on emission scenarios. If we force a linear fit onto a nonlinear trend:

- part of the signal leaks into the residuals
- $\sigma_{\text{noise}}^2$ is overestimated
- SNR is underestimated, making detection harder than it really is

The remedy is to replace the linear model with a **polynomial regression**:

$$f(t) = a_0 + a_1 t + a_2 t^2 + \dots + a_p t^p$$

A 2nd-order (quadratic) polynomial is often sufficient for capturing a time-varying forced trend. This is just multiple regression with predictors $t, t^2, \dots, t^p$ — the normal equations still apply:

$$\mathbf{a} = \mathbf{C}_{xx}^{-1}\,\mathbf{C}_{xy}$$

where the predictor matrix now contains powers of time.

:::{admonition} Example / deeper dive
:class: note
**Why does a non-linear trend matter for SNR?**

Suppose the true forced response is $f(t) = 0.01\,t^2$ (accelerating trend), but we fit a linear model. The linear fit will explain only part of the variance in $f(t)$, absorbing the remainder into the residuals. This inflates $\sigma_\eta^2$, suppresses SNR, and causes us to underestimate our confidence in the forced response.

A quadratic fit, by contrast, captures the curvature exactly, so the residuals contain only internal variability — giving a faithful estimate of noise and a correct SNR.
:::

#### Confidence intervals for a nonlinear trend with AR1 noise

Once the polynomial trend $f(t)$ has been removed, the residuals $\eta(t) = y(t) - f(t)$ are typically modeled as an AR1 process:

$$\eta(t) = \alpha\,\eta(t-\Delta t) + \epsilon(t), \qquad \epsilon \sim \mathcal{N}(0,\,\sigma_\epsilon^2)$$

The procedure for constructing confidence intervals is then:

1. **Fit** $f(t)$ by polynomial regression of order $p$ (e.g., $p = 2$).
2. **Estimate** the lag-1 autocorrelation $\alpha = \rho(\Delta t)$ of the residuals $\eta(t)$.
3. **Compute** the effective sample size using the Leith formula:
   $$N^* \approx N\frac{1-\alpha}{1+\alpha}$$
4. **Compute** the standard deviation of the residuals $\sigma_\eta$ and the standard error of the fit:
   $$\sigma_f = \frac{\sigma_\eta}{\sqrt{N^*}}$$
5. **Construct** the pointwise confidence interval:
   $$\text{CI}(t) = f(t) \pm t_{N^*-p,\,\alpha/2}\cdot\sigma_f$$

where $t_{N^*-p,\,\alpha/2}$ is the t-critical value with $N^*-p$ degrees of freedom (we lose $p$ degrees of freedom for the $p$ fitted polynomial coefficients).

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure: `nonlinear_snr.py`* — (left) linear fit with CI using $N$ (orange) vs. $N^*$ (blue); (right) quadratic fit with the same CI comparison. Note how the quadratic fit captures the true forced response (green dashed) much better and yields a higher SNR.

```{figure} nonlinear_snr.png
:width: 100%
:align: center
95% CI (gray shading) for a quadratic forced trend $f(t)=a_2 t^2$ embedded in AR1 noise ($\alpha=0.70$). Each panel uses the same single realization (gray dots), with the blue dots indicating the data window used for fitting. The fitted quadratic trend (red) is compared to the true forced response (green dashed), which is exactly zero at $t=0$. The CI is wide and overlaps zero for short records ($T=10, 30$ yr — signal undetected), but narrows progressively until the forced response clearly emerges from noise ($T=60, 100$ yr — signal detected). The CI width uses the hat-matrix leverage inflated by $\sqrt{T/N^*}$ to account for AR1 autocorrelation.
```
:::

This framework was applied to detecting the forced response of atmospheric rivers under greenhouse warming in {cite:t}`Tseng2021`, who showed that a second-order polynomial fit combined with an AR1 noise model can identify the time of emergence of a forced signal even from a single model realization — without needing a large ensemble.

:::{admonition} Example / deeper dive
:class: note
**Connecting back to the linear case**

For $p = 1$ (linear trend), the framework above reduces exactly to the standard t-test for the regression slope:

$$t = \frac{a_1}{\sigma_{a_1}}, \qquad \text{with } \nu = N^* - 2$$

The nonlinear case is simply the natural generalization: fit a richer model for $f(t)$, then use the residual AR1 structure to correct the degrees of freedom. The key insight is that **SNR, polynomial fitting, AR1 noise, and effective sample size are all part of one unified regression framework**.
:::

## Theory of Correlation (Pearson's Correlation)

### Statistical significance of correlations

The correlation $r$ between two time series $x(t)$ and $y(t)$ gives a measure of how well the two time series vary linearly together. $-1 \leq r \leq 1$, with numbers closer to $\pm1$ implying a stronger linear relationship.

We denote the sample correlation as $r$ and the theoretical true value as $\rho$.

If $\rho = 0$, we can use the t-statistic:

$$t = r\frac{\sqrt{N-2}}{\sqrt{1-r^2}}$$

:::{admonition} Example: testing the hypothesis that $\rho = 0$
:class: note
We have two time series, each of length 20, correlated at $r = 0.6$. Does this exceed the 95% confidence interval under $H_0: \rho = 0$?

We had no prior knowledge of the sign of the correlation, so we use a two-tailed t-test. For $\nu = N-2 = 18$, the critical value is $t_c = 2.1$.

$$t = 0.6\frac{\sqrt{20-2}}{\sqrt{1-0.6^2}} = 3.18$$

Since $t = 3.18 > t_c = 2.1$, we can reject the null hypothesis.
:::

:::{admonition} Example: confidence limits on the true correlation
:class: note
What are the 95% confidence limits on the true correlation if you drew 21 samples and obtained $r = 0.8$?

$$\begin{aligned}
Z &= \frac{1}{2}\ln\left(\frac{1+0.8}{1-0.8}\right) = 1.099 \\
\sigma_Z &= \frac{1}{\sqrt{21-3}} = 0.235
\end{aligned}$$

With $t_{0.025} = 2.1$ (for $\nu = 21-3 = 18$):

$$Z - 2.1\sigma_Z \leq \mu_Z \leq Z + 2.1\sigma_Z \quad \Rightarrow \quad 0.61 \leq \mu_Z \leq 1.59$$

Converting back to correlation via $\rho = \tanh(\mu_Z)$:

$$0.54 \leq \rho \leq 0.92$$
:::

The above statistic only works if the underlying distributions are normal, or if $N$ is large enough for the CLT to apply (roughly $N > 20$).

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: `testing_normality_of_correlations.py`*
:::

If $\rho \neq 0$, we must use the **Fisher-Z Transformation**. When the true correlation is not zero, the distribution of $r$ is not symmetric, so we cannot directly use the normal/t distribution. The Fisher-Z transformation converts $r$ into a quantity that is approximately normally distributed:

$$Z = \frac{1}{2}\ln\left(\frac{1+r}{1-r}\right) = \tanh^{-1}(r)$$

The Fisher-Z statistic is normally distributed with:

$$\begin{aligned}
\mu_Z &= \frac{1}{2}\ln\left(\frac{1+\rho}{1-\rho}\right) \\
\sigma_Z &= \frac{1}{\sqrt{N-3}}
\end{aligned}$$

The confidence bounds for $Z$ are:

$$Z - t_c\sigma_Z \leq \mu_Z \leq Z + t_c\sigma_Z$$

To convert back from $\mu_Z$ to the actual correlation $\rho$:

$$\rho = \frac{e^{2\mu_Z} - 1}{e^{2\mu_Z}+1} = \tanh(\mu_Z)$$

### Comparing two non-zero sample correlations

To test whether two correlations $r_1$ (from sample $N_1$) and $r_2$ (from sample $N_2$) are significantly different, apply the Fisher-Z to each:

$$\begin{aligned}
Z_1 &= \frac{1}{2}\ln\left(\frac{1+r_1}{1-r_1}\right) \\
Z_2 &= \frac{1}{2}\ln\left(\frac{1+r_2}{1-r_2}\right)
\end{aligned}$$

Then use the z-score for the difference of means:

$$z = \frac{Z_1 - Z_2 - \delta_{1,2}}{\sigma_{1,2}}, \qquad \sigma_{1,2} = \sqrt{\frac{1}{N_1-3} + \frac{1}{N_2-3}}$$

where $\delta_{1,2} = \mu_1 - \mu_2$ is the hypothesized difference (typically 0 if $H_0: \rho_1 = \rho_2$).

### Spearman's rank correlation

Spearman's rank correlation is a **nonparametric** test for whether paired data monotonically co-vary. No normality assumption is needed.

The original data $x_i$ and $y_i$ are converted into ranks $X_i$ and $Y_i$, and the correlation is computed on the ranks:

$$\rho_s = \frac{\sum_i (X_i - \overline{X})(Y_i - \overline{Y})}{\sqrt{\sum_i (X_i - \overline{X})^2 \sum_i (Y_i - \overline{Y})^2}}$$

When there are duplicate values, ranks are set to the average position. The standard error is:

$$\sigma_{\rho_s} = \frac{0.6325}{\sqrt{N-1}}$$

Significance can be tested using the Fisher-Z test or the t-test (for $H_0: \rho = 0$), as for Pearson's $r$.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: see slides 08\_correlation.pdf*
:::

Note: a second nonparametric method is Kendall's Tau Rank Correlation — not covered here.

## Autocorrelation & Estimating the Number of Independent Samples

Thus far, we have assumed that our time series have no intrinsic memory. Now, we will discuss these assumptions and how to determine the true number of degrees of freedom in an autocorrelated data set.

### Stationarity

Stationarity implies that the statistics of a time series (mean and higher-order moments) are independent of time — unchanging in time. In general, we will assume this is the case. This means one should **remove any trend** in the data before performing the analysis, using the linear regression method discussed above.

### Autocorrelation

The **autocovariance function** $\gamma(\tau)$ is the covariance of a time series with itself at lag $\tau$:

$$\gamma(\tau) = \frac{1}{(t_N - \tau) - t_1}\sum_{t=t_1}^{N-\tau}\left[x'(t)\cdot x'(t+\tau)\right]$$

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: Draw out example of how autocovariance works.*
:::

For a time series with integer positions $k = 1, 2, \dots, N$:

$$\gamma(\tau) = \overline{x'(t)\,x'(t+\tau)}$$

At $\tau = 0$: $\gamma(0) = \overline{x'^2} = \text{variance}$.

The **autocorrelation** $\rho(\tau)$ is $\gamma(\tau)$ normalized by $\gamma(0)$ — simply the correlation of a time series with itself at another time.

Notes:
- $\gamma$ is symmetric about $\tau = 0$
- $-1 \leq \rho(\tau) \leq 1$
- $\rho(0) = 1$
- if the time series is not periodic, $\rho(\tau) \rightarrow 0$ as $\tau \rightarrow \infty$

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: see slides 08\_correlation.pdf*
:::

### The first-order autoregressive model (AR1 / red noise)

Also referred to as a "first order Markov process" or "red noise."

**Red noise:** "today is like yesterday plus noise"

$$x(t) = a\cdot x(t-\Delta t) + b\cdot \epsilon(t)$$

where:
- $x$ is a standardized variable (zero mean, unit variance)
- $\Delta t$ is the (constant) time interval between data points
- $a \in [0,1]$ measures the memory of the previous state
- $\epsilon(t) \sim \mathcal{N}(0,1)$ is white noise

**Deriving $a$:** Multiply both sides by $x(t-\Delta t)$ and time-average:

$$\overline{x(t)\,x(t-\Delta t)} = a\underbrace{\overline{x^2(t-\Delta t)}}_{=1} + b\underbrace{\overline{\epsilon(t)\,x(t-\Delta t)}}_{=0}$$

$$\Rightarrow \quad a = \overline{x(t)\,x(t-\Delta t)} = \rho(\Delta t) = \rho(1)$$

**Deriving $b$:** Square both sides and time-average:

$$\begin{aligned}
\overline{x^2(t)} &= a^2\overline{x^2(t-\Delta t)} + b^2\overline{\epsilon^2(t)} \\
1 &= a^2 + b^2 \quad \Rightarrow \quad b = \sqrt{1-a^2}
\end{aligned}$$

**Autocorrelation of red noise:** Multiplying the recursion two steps forward by $x(t)$ and averaging shows that $\rho(2\Delta t) = \rho^2(\Delta t)$, and more generally:

$$\rho(n\Delta t) = \rho^n(\Delta t) = e^{-n\Delta t/T_e}$$

The autocorrelation decays **exponentially** with an e-folding time:

$$T_e = \frac{-\Delta t}{\ln(a)}$$

The e-folding time $T_e$ is the lag at which $\rho$ drops to $1/e \approx 0.368$. For example, if $\Delta t = 1$ day and $a = \rho(1) = 0.6$, then $T_e = 2$ days.

### White noise

White noise is the special case of AR1 with $a = 0$ (i.e. $\rho(\tau > 0) = 0$). It has equal power at all frequencies and zero autocorrelation — no memory of previous time steps. In geophysics, white noise is generally assumed to be normally distributed.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: `correlation_with_memory_examples.py`*
:::

:::{admonition} Side note: diffusion models in machine learning and their connection to AR processes
:class: note

**Diffusion models** (also called Denoising Diffusion Probabilistic Models, DDPMs) have become the dominant generative model architecture in machine learning — used for image synthesis, weather downscaling, and bias correction. At their core, they are built on a process that is mathematically identical to AR(1).

**The forward process — adding noise step by step**

Given a data sample $\mathbf{x}_0$ (e.g., a high-resolution precipitation field), a diffusion model defines a sequence of increasingly noisy versions $\mathbf{x}_1, \mathbf{x}_2, \dots, \mathbf{x}_T$:

$$\mathbf{x}_t = \sqrt{1-\beta_t}\,\mathbf{x}_{t-1} + \sqrt{\beta_t}\,\boldsymbol{\epsilon}_t, \qquad \boldsymbol{\epsilon}_t \sim \mathcal{N}(\mathbf{0}, \mathbf{I})$$

where $\beta_t \in (0,1)$ is a **noise schedule** (a small, pre-defined sequence that increases with $t$). This is **exactly an AR(1) process** with time-varying coefficient $a_t = \sqrt{1-\beta_t}$ and noise amplitude $\sqrt{\beta_t}$.

By the end of the forward chain ($t = T$, typically $T = 1000$ steps), $\mathbf{x}_T \approx \mathcal{N}(\mathbf{0}, \mathbf{I})$ — pure Gaussian noise, regardless of what $\mathbf{x}_0$ was.

Using the telescoping property of AR(1), one can skip directly to any step:

$$\mathbf{x}_t = \sqrt{\bar{\alpha}_t}\,\mathbf{x}_0 + \sqrt{1-\bar{\alpha}_t}\,\boldsymbol{\epsilon}, \qquad \bar{\alpha}_t = \prod_{s=1}^{t}(1-\beta_s)$$

This is the closed-form solution for the AR(1) recursion: $\bar{\alpha}_t$ plays the role of $a^t$ (the t-step autocorrelation) in our notation.

**The reverse process — learning to denoise**

The model learns the reverse: given $\mathbf{x}_t$, predict $\mathbf{x}_{t-1}$ (i.e., remove one step of noise). This reverse distribution is intractable analytically, so a neural network (usually a **U-Net**) $\boldsymbol{\epsilon}_\theta(\mathbf{x}_t, t)$ is trained to predict the noise $\boldsymbol{\epsilon}$ that was added at step $t$. Sampling then iterates:

$$\mathbf{x}_{t-1} = \frac{1}{\sqrt{1-\beta_t}}\!\left(\mathbf{x}_t - \frac{\beta_t}{\sqrt{1-\bar{\alpha}_t}}\,\boldsymbol{\epsilon}_\theta(\mathbf{x}_t, t)\right) + \sqrt{\beta_t}\,\mathbf{z}, \qquad \mathbf{z} \sim \mathcal{N}(\mathbf{0},\mathbf{I})$$

starting from $\mathbf{x}_T \sim \mathcal{N}(\mathbf{0},\mathbf{I})$ and working backwards to $\mathbf{x}_0$.

**Application to downscaling**

For climate downscaling, $\mathbf{x}_0$ is a high-resolution field (e.g., 1 km precipitation) and the U-Net is **conditioned** on a low-resolution field $\mathbf{y}$ (e.g., 25 km GCM output). The model learns the conditional distribution $p(\mathbf{x}_0 \mid \mathbf{y})$, allowing it to generate statistically realistic high-resolution fields consistent with the coarse model output.

**Summary of the AR1–diffusion connection**

| AR(1) concept | Diffusion model equivalent |
|---|---|
| $a = \sqrt{1-\beta}$ (memory coefficient) | Noise schedule $\beta_t$ |
| $b = \sqrt{1-a^2}$ (noise amplitude) | $\sqrt{\beta_t}$ |
| $a^n$ (n-step autocorrelation) | $\sqrt{\bar{\alpha}_t}$ |
| White noise limit ($a \to 0$) | $\mathbf{x}_T \to \mathcal{N}(\mathbf{0},\mathbf{I})$ as $T\to\infty$ |
| Stationary variance = 1 | Unit variance of the noise prior |

The key innovation of diffusion models is not the forward AR(1) process — that is just Gaussian noise addition — but the learned *reverse* process, which turns pure noise back into structured data. The statistical machinery of AR1 you learned here is the exact mathematical foundation.
:::

### Effective sample size $N^*$

Persistence in a data set leads to **overestimation of the sample size**, because each data point is not independent of its neighbors. If persistence is ignored, the standard error of the mean is underestimated and the t-statistic is inflated.

The solution is to introduce an **effective sample size** $N^* \leq N$ and substitute it for $N$ in formulas.

For a first-order AR process, $N^*$ can be estimated with the **Bretherton/Wilks approximation** (Wilks, p. 127):

$$\frac{N^*}{N} \approx \frac{1-\rho(\Delta t)}{1+\rho(\Delta t)}$$

- if $\rho(1) = 0$ (white noise): $N^* = N$
- as $\rho(1)$ increases, $N^*$ decreases

An equivalent form from Leith (*J. Appl. Meteor.*, 1973):

$$N^* \approx \frac{N\Delta t}{2T_e} = \frac{\text{total record length}}{2 \times \text{e-folding time}}$$

The factor of 2 reflects the fact that any point in red noise can be predicted by points both before and after it. The Leith formula can also be written as:

$$\frac{N^*}{N} \approx \frac{\ln a}{-2}$$

The table below shows $N^*/N$ as a function of lag-1 autocorrelation:

| $\rho(\Delta t)$ | $<0.1$ | $0.3$ | $0.5$ | $0.7$ | $0.9$ |
|---|---|---|---|---|---|
| $N^*/N$ | $\approx 1$ | $0.60$ | $0.35$ | $0.18$ | $0.053$ |

Bretherton et al. (*J. Climate*, 1999) proposed a less conservative approximation (use for variance/higher-order moments):

$$\frac{N^*}{N} \approx \frac{1-\rho^2(\Delta t)}{1+\rho^2(\Delta t)}$$

This yields nearly twice as many degrees of freedom as the Leith formula. For testing the **mean**, use the Leith formula; for testing **variance**, the Bretherton formula may be used.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: `effective_sample_size.py`* — plot of $N^*/N$ vs $\rho(\Delta t)$ comparing the Leith and Bretherton approximations.
:::

## Multiple Regression (Multi-linear Regression)

*Basic idea:* Generalize the regression coefficient derivation to multiple **linear** predictors:

$$\hat{y} = a_0 + a_1 x_1 + a_2 x_2 + \dots + a_n x_n$$

The fit now lives in a multi-dimensional predictor space.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: 3-D scatter diagram* — points $\hat{y}$ plotted against two orthogonal predictor axes $X_1$ and $X_2$, with a best-fit plane drawn through the point cloud.
:::

If $X_1$ and $X_2$ are **orthogonal** (at right angles):
- they give independent information
- their inner product is 0
- if they span the space, they "form a basis"

If $X_1$ and $X_2$ are **not orthogonal**:
- they are not independent and share redundant information

The usefulness of independent predictors motivates EOF analysis (to be discussed later).

### Generalized normal equations

For multiple predictors $x_1, x_2, \dots, x_n$, minimize:

$$Q = \sum_{i=1}^{N}\left(a_0 + a_1 x_{1,i} + a_2 x_{2,i} + \dots + a_n x_{n,i} - y_i\right)^2$$

Setting $\partial Q/\partial a_i = 0$ for $i = 0, \dots, n$ gives $n+1$ equations. If the mean has been removed from all variables ($a_0 = 0$), the $j$th equation is:

$$\overline{x_j y} = \sum_{i=1}^{n} a_i \overline{x_j x_i}$$

In matrix form:

$$\underbrace{\begin{bmatrix}
\overline{x_1^2} & \overline{x_1 x_2} & \cdots \\
\overline{x_2 x_1} & \overline{x_2^2} & \cdots \\
\vdots & & \ddots
\end{bmatrix}}_{\mathbf{C}_{xx}}
\underbrace{\begin{bmatrix} a_1 \\ a_2 \\ \vdots \end{bmatrix}}_{\mathbf{a}}
=
\underbrace{\begin{bmatrix} \overline{x_1 y} \\ \overline{x_2 y} \\ \vdots \end{bmatrix}}_{\mathbf{C}_{xy}}$$

Or compactly: $C_{x_i x_j}\, a_j = C_{x_i y}$

Key observations:
1. The left-hand side $\mathbf{C}_{xx}$ is the **covariance matrix** of the predictors (diagonal = variances, off-diagonal = covariances)
2. The right-hand side $\mathbf{C}_{xy}$ is the **covariance vector** between predictors and predictand
3. If each variable is standardized, $\mathbf{C}_{xx}$ becomes the **correlation matrix** and $\mathbf{C}_{xy}$ the **correlation vector**
4. If predictors are linearly independent, off-diagonal elements are 0 and the $a_j$'s can be found algebraically
5. Otherwise, solve via matrix inversion:

$$\mathbf{a} = \mathbf{C}_{xx}^{-1}\,\mathbf{C}_{xy}$$

### Multiple regression — how many predictors should I use?

Assuming all variables are standardized, the normal equations become $r(x_i, x_j)\,a_i = r(x_j, y)$. For two predictors:

$$\begin{bmatrix} 1 & r_{1,2} \\ r_{1,2} & 1 \end{bmatrix}
\begin{bmatrix} a_1 \\ a_2 \end{bmatrix}
=
\begin{bmatrix} r_{1,y} \\ r_{2,y} \end{bmatrix}$$

Solving:

$$a_1 = \frac{r_{1,y} - r_{1,2}r_{2,y}}{1 - r_{1,2}^2}, \qquad a_2 = \frac{r_{2,y} - r_{1,2}r_{1,y}}{1 - r_{1,2}^2}$$

The total fraction of explained variance ($R^2$) with two predictors:

$$R^2 = \frac{r_{1,y}^2 + r_{2,y}^2 - 2r_{1,y}r_{2,y}r_{1,2}}{1 - r_{1,2}^2}$$

:::{admonition} Example: does adding a second predictor help?
:class: note
Say $r_{1,y} = r_{2,y} = r_{1,2} = 0.5$.

With only $x_1$: $R_1^2 = r_{1,y}^2 = 0.25$

Adding $x_2$:

$$R_{1,2}^2 = \frac{0.5^2 + 0.5^2 - 2\times 0.5\times 0.5\times 0.5}{1-0.5^2} = 0.33$$

Adding $x_2$ increases the explained variance from 25% to 33%. ✓

Now suppose $r_{2,y} = 0.25$ (everything else the same):

$$R_{1,2}^2 = \frac{0.5^2 + 0.25^2 - 2\times 0.5\times 0.25\times 0.5}{1-0.5^2} = 0.25$$

Adding $x_2$ adds **nothing**. The minimum useful correlation for $x_2$ is:

$$|r(x_2,y)|_{\text{min}} > |r(x_1,y)\cdot r(x_1,x_2)| = 0.5\times 0.5 = 0.25$$

So $r_{2,y}$ must exceed 0.25 to be worth adding.
:::

Key guidelines:
- Ideal case: $r_{1,2} = 0$ — two completely independent predictors
- Worst case: $r_{1,2} = 1$ — $x_2$ provides no new information
- Adding too many predictors can lead to **overfitting** — fitting the noise rather than the signal. Always use as few predictors as possible and test the fit on **independent data**.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: Adjusted $R^2$* — see <https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2>
:::

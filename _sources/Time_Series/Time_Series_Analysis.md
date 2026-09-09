(Time_Series)=
# Week 12-15: Time Series Analysis

In Weeks 8-11 we decomposed data into patterns that the data itself chose — the EOFs were *empirical*, derived from the covariance matrix, and their shapes were whatever maximized variance. This block does the same kind of decomposition, but onto a basis we *impose*: sines and cosines.

That sounds like a step backwards, and in one sense it is — a prescribed basis cannot adapt to your data. But it buys something EOFs cannot give you. Sines and cosines have a **frequency**, so the coefficients tell you not just "how much structure" but "on what timescale". Questions like *does this signal have a 40-day periodicity?*, *how much variance lives at interannual timescales?*, and *does the wind at 120°E lead the wind at 150°E by a quarter cycle?* are all frequency questions.

The thread from Weeks 4-6 continues to hold: **harmonic analysis is multiple regression** where the predictors happen to be sines and cosines. Because those predictors are mutually orthogonal, the messy normal equations from multiple regression collapse into one-line formulas. Almost everything in this block follows from that.

---

## Part 1 — Harmonic Analysis

### Brief review of sines and cosines

$$y(t) = A\sin(2\pi f t + \phi) = A\sin(\omega t + \phi)$$

where

- $A$ is the **amplitude**, the deviation of the curve from zero
- $f$ is the **frequency**, in oscillations per time step
- $\omega = 2\pi f$ is the **angular frequency**, in radians per time step
- $\phi$ is the **phase** in radians. When $\phi \ne 0$ the wave is shifted in time by $\phi/\omega$ time steps

Four related quantities that are easy to confuse:

| Quantity | Symbol | Definition |
|---|---|---|
| wavenumber | $k$ | an integer; number of oscillations that fit inside the domain |
| frequency | $f$ | oscillations per time step, $f = k/T$ |
| angular frequency | $\omega$ | radians per time step, $\omega = 2\pi f = 2\pi k/T$ |
| period | $T$ | time steps per oscillation, $1/f$ |

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure:* $\cos(x)$ and $\sin(x)$ plotted over $-2\pi \le x \le 2\pi$, to fix the relationship between amplitude, phase and period.
:::

### Fourier sums

**The basic idea.** Interpret a time (or space) series as a summation of contributions from harmonic functions, each with a unique temporal or spatial scale. These harmonic functions form a **basis** — exactly the sense of "basis" from the linear algebra review in Weeks 8-11.

This is analogous to multiple regression, where you describe one variable as a combination of others. Here the predictors are

$$\cos\left(\frac{2\pi k t}{T}\right), \qquad \sin\left(\frac{2\pi k t}{T}\right)$$

where $k$ is an integer from 1 to $N/2$ ($N$ being the number of values you have), often called the **wavenumber**. Note $k/T = f$, the number of oscillations per time step.

So $y(t)$ decomposes as

$$\boxed{y(t) = A_0 + \sum_{k=1}^{N/2}A_k\cos\left(\frac{2\pi kt}{T}\right) + \sum_{k=1}^{N/2}B_k\sin\left(\frac{2\pi kt}{T}\right)}$$

with the following properties:

- $y$ is resolved on the interval $0 \le t \le T$
- the length of the time series is $N+1$. ($N$ must be even so the highest wavenumber is an integer, but the record length must be odd so that all harmonic functions fit an integer number of times into the full record.)
- $A_k$ and $B_k$ are the **regression coefficients** for each predictor
- each predictor has frequency $2\pi k/T$ and so fits into $t_1 = 0$ to $t_{N+1} = T$ an integral $k$ number of times
- if $k = 1$, the function completes $2\pi$ radians in the full length of the data — the **lowest resolved frequency**
- if $k = N/2$, the function completes $2\pi$ radians in 2 time steps — the highest resolvable frequency, the **Nyquist frequency**
- if $k = N/2$ then $B_k = 0$. (Try drawing it: you can fit a cosine at two points per cycle, but not a sine.)
- for a given $k$ you fit one cosine and one sine of variable amplitude. This is equivalent to fitting one cosine with variable amplitude *and* variable phase. Either way, **two predictors per $k$**.

:::{admonition} Example / deeper dive
:class: note
**Why $2\pi$ radians in 2 time steps, and not 1?**

To define an oscillation you need to see it go up *and* come back down. With one sample per cycle every sample lands at the same phase, and the wave is indistinguishable from a constant. Two samples per cycle is the bare minimum that can register an alternation — and even then only if your samples don't land exactly on the zero crossings.

Counting the degrees of freedom gives the same answer. $N$ data points can support exactly $N$ independent predictors. Each $k$ costs two (a sine and a cosine), so $k$ can run to $N/2$ — no further. The Nyquist limit is not a property of waves; it is a bookkeeping constraint.
:::

### More on the Nyquist frequency and aliasing

The Nyquist frequency is the highest frequency resolvable in your data:

- one oscillation per 2 time steps, i.e. 1 oscillation per $2\Delta t$
- wavenumber $k = (\text{length of time series} - 1)/2$

What if higher frequencies are present? They get **aliased** onto lower frequencies — a real signal at high frequency masquerades as a spurious signal at low frequency. This is a serious problem when there is a large amount of variance above the Nyquist frequency.

:::{admonition} A point that is easy to get backwards
:class: warning
The Nyquist frequency depends **entirely on your data time step** $\Delta t$ — *not* on how much data you have. Collecting a longer record does not raise the Nyquist frequency; it only lowers the lowest resolved frequency. The smaller the time step, the higher the Nyquist frequency.

So the two ends of your spectrum are set by two different things: $\Delta t$ sets the high-frequency limit, record length $T$ sets the low-frequency limit.
:::

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `harmonic_analysis_nyquist.m`
:::

### Discrete Fourier Transform

For discrete, evenly spaced data on $0 \le t \le T$, where $0$ and $T$ coincide with data points 1 and $N+1$ and $N$ is even:

1. the functions are $\cos(2\pi kt/T)$ and $\sin(2\pi kt/T)$

2. the average of each sine/cosine on $t_1 = 0$ to $t_{N+1} = T$ is zero, so $A_0 = \overline{y}$

3. the harmonic functions are **mutually orthogonal** on the interval, e.g.

$$r(\cos x, \sin x) = 0, \qquad r(\cos x, \cos 2x) = 0$$

4. the regression coefficients of $y(t)$ onto harmonic function $x_k(t)$ follow from the normal equations of multiple regression:

$$A_k = \frac{\overline{x_k'y'}}{\overline{x_k'^2}}$$

$$\begin{bmatrix}
\overline{x_1^2} & \overline{x_1x_2} & \overline{x_1x_3} & \cdots \\
\overline{x_2x_1} & \overline{x_2^2} & \overline{x_2x_3} & \cdots \\
\overline{x_3x_1} & \overline{x_3x_2} & \overline{x_3^2} & \cdots \\
\vdots & \vdots & \vdots &
\end{bmatrix}
\begin{bmatrix}A_1\\A_2\\A_3\\\vdots\end{bmatrix}
=
\begin{bmatrix}\overline{x_1y}\\\overline{x_2y}\\\overline{x_3y}\\\vdots\end{bmatrix}$$

**This is the payoff of orthogonality.** In general multiple regression that matrix is dense and must be inverted. Because the predictors here are orthogonal, all cross-covariances vanish, the matrix is diagonal, and each $A_k$ is obtained independently — no inversion, no interaction between predictors.

5. each harmonic function has variance $\overline{x_k'^2} = \tfrac12$, except at $k = N/2$ where the cosine term has variance 1 and the sine term 0

6. the solutions therefore reduce to

$$A_k = 2\cdot\overline{\cos\left(\frac{2\pi kt}{T}\right)y'(t)}, \qquad B_k = 2\cdot\overline{\sin\left(\frac{2\pi kt}{T}\right)y'(t)}$$

$$A_{N/2} = \overline{\cos\left(\frac{2\pi kt}{T}\right)y'(t)}, \qquad B_{N/2} = 0$$

where the factor 2 comes from $\overline{x_k'^2} = \tfrac12$.

7. thus

$$y(t) = \overline{y} + \sum_{k=1}^{N/2-1}\left[A_k\cos\left(\frac{2\pi kt}{T}\right) + B_k\sin\left(\frac{2\pi kt}{T}\right)\right] + A_{N/2}\cos\left(\frac{\pi Nt}{T}\right)$$

### Amplitude-phase form

The same result in terms of cosine alone:

$$y(t) = \overline{y} + \sum_{k=1}^{N/2}C_k\cos\left(\frac{2\pi kt}{T} - \phi_k\right), \qquad \boxed{C_k^2 = A_k^2 + B_k^2}$$

which follows from the Pythagorean theorem with $A_k$ on the $x$-axis, $B_k$ on the $y$-axis, and $\phi$ the angle between.

The fraction of variance explained by a particular harmonic is

$$R^2(y, x_k) = \frac{\left(\overline{x_k'y'}\right)^2}{\overline{x_k'^2}\cdot\overline{y'^2}}$$

which for the cosine and sine components gives

$$R^2 = \frac{A_k^2}{2\overline{y'^2}}, \qquad R^2 = \frac{B_k^2}{2\overline{y'^2}}, \qquad R^2 = \frac{A_{N/2}^2}{\overline{y'^2}}\ \ (k = N/2)$$

and the **actual** (not fractional) variance explained by harmonic $k$ is

$$\boxed{\text{variance explained} = \frac{A_k^2 + B_k^2}{2} = \frac{C_k^2}{2}}$$

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `example_harmonic_coeff_calculator.m`
:::

---

## Part 2 — The Power Spectrum

The plot of $C_k^2/2$ versus $k$ is the **spectrum** of $y(t)$:

- if $t$ is time, it is the **frequency spectrum**
- if $t$ is space, it is the **wavenumber spectrum**

This is a **line spectrum** — defined only at integer $k$. That is fine for an infinite sample, but has three serious drawbacks for a finite record:

1. integer values of $k$ have no special relationship to the population being sampled; they depend entirely on the record length $T$, which is presumably arbitrary
2. each spectral line contains about **2 degrees of freedom** — $N/2$ estimates from $N$ data points
3. except for the diurnal and seasonal cycles, most geophysical phenomena are only *quasi*-periodic, so they do not conform to a single integer wavenumber and are better represented as **spectral bands of finite width**

These considerations motivate a **continuous power spectrum**, where the variance of $y(t)$ is given per unit frequency:

$$\overline{y'^2} = \int_0^{k^*}\Phi(k)\,dk$$

where $\Phi(k)$ is the continuous power spectrum and $k^*$ is the Nyquist frequency (one cycle per $2\Delta t$). $\Phi(k)$ is just a continuous version of $C_k^2/2$, estimated either by smoothing the line spectrum or by averaging many line spectra together.

### The resolution–reliability trade-off

Two ways to gain degrees of freedom:

1. **Average adjacent spectral estimates.** Lose resolution, gain dof. Averaging 10 adjacent estimates decreases frequency resolution by a factor of 10 and increases dof by a factor of 10.
2. **Average realizations of the spectra.** Lose information at low frequencies, gain dof. Subdividing into 10 series of equal length raises the lowest resolved frequency by a factor of 10 and increases dof by a factor of 10.

:::{admonition} The trade-off you cannot escape
:class: warning
There is a permanent tension between **resolution** and **reproducibility**:

> High resolution / high information / **low quality statistics**
> $\Downarrow$
> Lower resolution / lower information / **higher quality statistics**

It is somewhat akin to Heisenberg's Uncertainty Principle. Every choice below — chunk length, window, overlap, band averaging — is a position on this trade-off, not a solution to it.
:::

Degrees of freedom per spectral estimate:

$$\text{D.O.F.} = N/M^*$$

where $M^*$ is the number of independent spectral estimates and $N$ the number of actual data points. For line estimates $M^* = N/2$, giving $\text{D.O.F.} = N/(N/2) = 2$ — as claimed above.

Note: **as long as we use a red-noise fit as the null hypothesis, we do not need to reduce $N$ to account for autocorrelation.** The autocorrelation is already built into the null.

### Methods for computing the power spectrum

**Method 1: the direct method.** Calculate $C_k$ through harmonic analysis, i.e. via the analytic DFT solutions for $A_k$ and $B_k$.

- easy to code, but inefficient — there are numerous redundancies in the $A_k$, $B_k$ calculations
- the **Fast Fourier Transform (FFT)** is widely available and exploits those redundancies
- the transform assumes **cyclic continuity** (waves fit fully into the domain); it is typical to taper the edges of the series — see Part 5

**Method 2: the lag correlation method.** Less often used, but worth knowing. It rests on a theorem of Norbert Wiener:

> The autocorrelation function and the power spectrum are Fourier transforms of each other.

So the power spectrum can be found by performing harmonic analysis on the autocorrelation function within $-T_L \le \tau \le T_L$, where $T_L$ is the maximum lag. This connects directly to the autocorrelation machinery of Weeks 4-6.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure:* slides showing autocorrelation / power spectrum pairs.

*Worked example:* `autocorrelation_of_linear_trend.m`
:::

### Plotting the power spectrum

Two schemes, both preserving the property that **area under the curve is proportional to variance**.

**Linear scale.** Plot frequency against power spectral density $\Phi(\omega)$, so that

$$\int_{\omega_1}^{\omega_2}\Phi(\omega)\,d\omega$$

is the total variance between $\omega_1$ and $\omega_2$.

**Logarithmic scale.** When the frequency range spans several orders of magnitude, plot $\ln\omega$ on the $x$-axis. Then the $y$-axis must be $\omega\cdot\Phi(\omega)$ so that visual area remains proportional to variance:

$$\int_{\ln\omega_1}^{\ln\omega_2}\omega\,\Phi(\omega)\,d(\ln\omega) = \int_{\omega_1}^{\omega_2}\Phi(\omega)\,d\omega$$

since $dx = x\,d(\ln x)$. On the log scale the low-frequency end is stretched relative to the high-frequency end. It is often useful to normalize total area to 1, so area between two frequencies is directly the fraction of variance explained.

:::{admonition} How spectra mislead
:class: warning
Be careful with axis choices — each combination creates a specific visual illusion:

- **log power, linear frequency** → weak variances look more important than they are
- **log frequency, linear power** → high frequencies are much more important than they look
- **log–log** → both illusions at once

You will see log–log plots constantly. They are not wrong, but you must read them knowing that the area-equals-variance intuition no longer applies.
:::

:::{admonition} Figure / in-class demonstration
:class: tip
*Figures 6.11 and 6.12 (Hartmann Ch. 6b):* (left) schematic of continuous power versus frequency, where variance between two frequencies is the area under the curve; (right) schematic of frequency-times-power versus log-frequency, with the same area-variance property.
:::

---

## Part 3 — The Complex Fourier Transform

So far the Fourier transform has been written in real arithmetic. It is far more efficient in complex form:

$$f(t) = \frac{1}{2\pi}\int_{-\infty}^{\infty}F(\omega)e^{i\omega t}\,d\omega, \qquad F(\omega) = \int_{-\infty}^{\infty}f(t)e^{-i\omega t}\,dt$$

where $\omega$ is the radial frequency (radians per unit time) and $F(\omega)$ is generally complex. Different textbooks place the $1/(2\pi)$ in different places — check the convention before comparing formulas.

The bridge back to the real form is Euler's identity:

$$e^{i\theta} = \cos\theta + i\sin\theta, \qquad \cos\theta = \frac{e^{i\theta}+e^{-i\theta}}{2}, \qquad \sin\theta = \frac{e^{i\theta}-e^{-i\theta}}{2i}$$

The practical advantage: a complex coefficient carries **amplitude and phase in a single number**, so the $A_k$/$B_k$ pair becomes one object. This is what makes the cross-spectrum algebra of Part 7 tractable.

---

## Part 4 — Statistical Significance of Spectral Peaks

Spectral peaks are tested against the null hypothesis that the data is just **red noise** — that it has low-frequency variability but is not periodic. The amplitude of a peak is tested against the spectrum of the red-noise fit to the data.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure 6.18 (Hartmann Ch. 6b):* schematic of a spectrum, its null-hypothesis spectrum, and the ratio of variances at a frequency of interest. $\Phi_0$ is the "background" spectrum forming the null hypothesis.
:::

### WOSA and the FFT method

The most common approach uses the FFT, which exploits the fact that the series length can be chosen as an integer power of two, $M_{ch} = 2^n$ (mixed-radix FFTs allow $M_{ch} = 2^n3^m5^j$). The total series of length $N\Delta t$ is broken into smaller chunks of length $M_{ch}$, and the spectra of those chunks are averaged into a grand spectrum with some statistical reliability. This is **Welch's Overlapping Segment Analysis (WOSA)**.

Because a tapered window is normally applied, data near the ends of each chunk is weighted less heavily than data in the center. So that all data is weighted equally, the chunks should **overlap**. In Matlab this is the `noverlap` parameter. If you overlap by exactly half the chunk length and use a Hanning window, every data point gets exactly the same weight in the averaged spectrum — except the first and last $M_{ch}/2$ points of the record.

### The F-test

The significance of the ratio of variances $\Phi_1$ to $\Phi_0$ is assessed with the F-statistic:

$$\boxed{F = \frac{S_1^2}{S_0^2}}$$

where $S_1^2$ and $S_0^2$ are the variances of the original spectrum and the red-noise fit. **Remember — the power spectrum *is* this variance.**

The F-statistic requires degrees of freedom $\nu_1$ and $\nu_0$:

- for red noise ($\nu_0$), we typically assume a large number
- for the real series ($\nu_1$), use $N/(M_{chunk}/2) = N/(\text{number of spectral estimates})$

:::{admonition} A priori versus a posteriori — the multiple comparisons trap
:class: warning
You may only use *a priori* statistics if you **expected a peak at that frequency in advance**. Otherwise you must use *a posteriori* statistics, accounting for the probability that one frequency out of $M_{chunk}/2$ gave a significant peak by chance. Your true significance level is

$$(0.99)^{M_{chunk}/2}$$

if you test each peak at 99% confidence. With 64 spectral estimates that is $0.99^{64} \approx 0.53$ — barely better than a coin flip. Hunting for the largest peak in a spectrum and then testing it as though you had predicted it is one of the most common errors in spectral analysis.
:::

:::{admonition} Example / deeper dive
:class: note
**Example 4.2 — F-statistic and spectral peaks**

You have a daily time series 512 days long, split into 128-day chunks (4 realizations), and a peak in the spectrum 3 times higher than the red-noise value. Is it significant at 99%?

1. assume the peak is where we expected it *a priori*
2. we need $F_{crit} < 3$ to reject the null that the peak is a random fluctuation of red noise
3. determine the dof for the data and the red-noise series, $\nu_{data}$, $\nu_{red}$
4. $\nu_{data} = N/(M_{chunk}/2) = 512/(128/2) = 8$
5. $\nu_{red} = $ large (typically $> 100$)
6. for the F-statistic, the data is the numerator and red noise the denominator
7. from the table, the 99% level is $F_{crit} = 2.69$
8. since $F_{crit} = 2.69 < 3$, the peak passes and we reject the null hypothesis
:::

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `run_testing_spectral_significance.m`
:::

### Finding the power spectrum of red noise

The autocorrelation function for red noise is

$$\gamma(\tau) = \exp\left(-\frac{\tau}{T_e}\right)$$

with $T_e$ the e-folding timescale. By Wiener's theorem the power spectrum is its Fourier transform,

$$\Phi(\omega) = \int_{-\infty}^{\infty}e^{-\tau/T_e}e^{-i\omega\tau}\,d\tau$$

which integrates to

$$\boxed{\Phi(\omega) = \frac{2T_e}{1 + T_e^2\omega^2}}$$

This is largest at $\omega = 0$ and falls off more rapidly for large $T_e$.

Practical notes:

- it is typical to set the area under the spectrum to either 1 or the total variance of the data being fitted
- the amplitude looks different depending on whether you plot against $k$, period, $\omega$, etc.
- in the real world we only have an *estimate* of red noise from a limited sample, while the above is for an infinite sample. In practice use the **discrete red-noise spectrum** of Gilman et al. (1963), whose unnormalized form is

$$\Phi(\omega)_{\text{discrete}} = \frac{1-\rho^2}{1 - 2\rho\cos\left(\frac{h\pi}{N/2}\right) + \rho^2}$$

with $h = [0:1:N/2]$ and $\rho$ the lag-1 autocorrelation of the red-noise process
- the discrete version does not differ much from the theoretical one, and the theoretical version is more **conservative**
- **the FFT assumes stationarity — remove any trend before analysis.** A linear trend is not a periodic signal, but it will deposit large amounts of spurious power at the lowest frequencies.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure:* slide comparing the theoretical red-noise spectrum with Gilman et al. (1963).
:::

---

## Part 5 — Windows Due to Finite Data

Two key problems arise with spectral analysis of a finite discrete record:

1. **Aliasing** — we resolve only frequencies below the Nyquist frequency $1/(2\Delta t)$; higher frequencies get aliased down (Part 1)
2. **Leakage** — we assume all waveforms start and stop at $t = 0$ and $t = T$, but in reality many wavenumbers do not complete an integer number of cycles in the domain, so power leaks to neighbouring wavenumbers

Leakage is handled by applying a **window** to the data. To understand windows we first need convolution.

### Convolution

The convolution of two functions is

$$(w * g)(t) \overset{\text{def}}{=} \int_{-\infty}^{\infty}w(\tau)g(t-\tau)\,d\tau$$

:::{admonition} Example / deeper dive
:class: note
**Convolution is just a sliding weighted average**

Let $w(\tau) = [0,\ 1/3,\ 1/3,\ 1/3,\ 0]$ for $\tau = -2,-1,0,1,2$, and let

$$g(t) = [9,\ 7,\ 8,\ 6,\ 5,\ 4,\ 3,\ 8,\ 9,\ 10,\ 11] \quad \text{for } t = 1,2,3,4\ldots$$

Then

$$(w*g)(t=4) = 7\cdot0 + 8\cdot\tfrac13 + 6\cdot\tfrac13 + 5\cdot\tfrac13 + 4\cdot0 = \tfrac{19}{3}$$

That is all a convolution is — slide the weights along, multiply, sum. Every filter in Part 6 is this operation with different weights.
:::

### The Convolution Theorem

$$\boxed{\mathcal{F}\big(f(t)g(t)\big) = F(\omega) * G(\omega)}$$
$$\boxed{\mathcal{F}\big(f(t) * g(t)\big) = F(\omega)\cdot G(\omega)}$$

where $\mathcal{F}$ is the Fourier transform operator and capital letters denote transforms of the corresponding lower-case functions.

In English: **the Fourier transform of a product is the convolution of the individual transforms, and vice versa.** Multiplication in one domain is convolution in the other.

Now consider what this means for a finite record. We are in effect applying a window $g(t)$ to our data so that we only see $0 \le t \le T$ rather than $-\infty \le t \le \infty$:

- **infinite window:** $g(t) = 1$ for all $t$
- **finite window:** $g(t) = 1$ for $0 \le t \le T$, and $0$ elsewhere

So when we think we are computing the Fourier transform of $f(t)$, we are actually computing the transform of $f(t)\cdot g(t)$, which by the convolution theorem is

$$F(\omega) * G(\omega) = \int_{-\infty}^{\infty}F(\omega_0)G(\omega - \omega_0)\,d\omega_0$$

One way to think about this: **the true infinite-record spectrum $F(\omega)$ gets smeared by a moving weighted average whose shape is $G(\omega)$**, the transform of the window. You never see the true spectrum; you see it convolved with your window.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure 6.13 (Hartmann Ch. 6b):* schematic comparing the window function of an infinite time series (left, $w(t) = 1$ for all $t$) with that of a finite time series (right, $w(t) = 1$ for $0 \le t \le T$ and 0 elsewhere).
:::

### Response functions for various windows

$G(\omega - \omega_0)$ is called the **response function** of the window $g(t)$.

**The ideal case.** We would like

$$G(\omega - \omega_0) = 1 \ \text{ when } \omega = \omega_0, \qquad G(\omega-\omega_0) = 0 \ \text{ when } \omega \ne \omega_0$$

i.e. a delta function — each frequency estimate contaminated by no other frequency.

**The boxcar window.** This is the window you apply automatically, without choosing to, whenever you analyze a finite record: $g(t) = 1$ for $0 \le t \le T$ and 0 otherwise. Its transform is

$$G(\omega) = \frac{\sin(\omega T/2)}{\omega T/2} = 2\,\text{sinc}\!\left(\frac{\omega T}{2\pi}\right), \qquad \text{sinc}(x) = \frac{\sin(\pi x)}{\pi x}$$

Consequences:

- **smoothing**: the window smooths and distorts the true spectrum
- **side lobes**: these add frequencies to your spectrum that were not in the data. They exist because high-frequency harmonics are needed to fit the sharp corners of the boxcar
- **small $T$**: the response function is much wider, so there is more smoothing across frequencies
- **large $T$**: the response function approaches a delta function, closer to ideal

That last pair makes intuitive sense: for large $T$ the edges of your boxcar are a tiny part of the data; for small $T$ the edges are close to most of your data.

**The Hanning window** (also Hann, cosine bell, or Tukey window) tapers the ends rather than cutting sharply. It is applied automatically in Matlab routines like `spectrum.m` and `spectrum.welch`:

$$w(t) = \frac12\left(1 - \cos\left(\frac{2\pi t}{T}\right)\right) = \frac{2}{T}\cos^2\left(\frac{\pi t}{T}\right), \qquad 0 \le t \le T$$

with response function

$$G(\omega) = \text{sinc}\!\left(\frac{\omega T}{2\pi}\right) + \frac12\left[\text{sinc}\!\left(\frac{\omega T}{2\pi}+1\right) + \text{sinc}\!\left(\frac{\omega T}{2\pi}-1\right)\right]$$

The first term is identical to the boxcar; the last two **cancel the side lobes**.

Once again you cannot get something for nothing: the central lobe widens in order to remove the side lobes. Usually removing side lobes outweighs the extra smoothing.

That extra smoothing is not captured in the dof calculation, so a correction factor is added:

$$n = \frac{N}{M_{chunk}/2}\cdot f_w$$

where $f_w$ compensates for the window's smoothing, generally between 1 and 1.5. **For a Hanning window $f_w = 1.2$.**

**The Hamming window** is similar to Hanning with a slight modification, giving slightly better side-lobe reduction and slightly more smoothing of the central lobe.

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked examples:* `boxcar_hanning_response_functions.m`, `convolution_example.m`; slide of the boxcar response function.
:::

### Applying overlaps of the windows

If you have 400 data points split into 4 chunks of 100 to increase dof, the tapered window means you barely use the data at the chunk edges — but that data is no more or less important than data at the centers.

To weight each point more equally, apply the window to one chunk, then take the next chunk with a **50% overlap**. With 400 points and 100-point chunks you now have 7 chunks rather than 4. This is **WOSA** (Welch's Overlapping Segment Analysis).

:::{admonition} Overlapping does not multiply your dof
:class: warning
Do **not** conclude that 7 chunks means 7× the degrees of freedom. You are using most of the data twice, and you cannot double count. Typically dof is kept at $4\times2$ to be safe, though Welch's 1967 paper shows you actually gain $4\times2.8$ for 50% overlap.
:::

---

## Part 6 — Filtering

Windowing in the time domain has an exact counterpart in the frequency domain, usually called **frequency filtering** — though we can also filter in the time domain (e.g. a moving average).

- **low pass filter**: isolates low frequencies
- **high pass filter**: isolates high frequencies
- **bandpass filter**: removes both high and low frequencies, isolating a band in the middle

:::{admonition} The single most important warning in this block
:class: warning
**A bandpass filter will make even noise look periodic and interesting.**

If you band-pass white noise at 30–60 days and plot the result, you will see a beautiful oscillation with a ~40-day period. It is not there in the data; you put it there by discarding everything else. Always ask whether the periodicity survives a significance test against red noise (Part 4) *before* filtering, not after.
:::

### General idea of frequency filtering

Recall harmonic analysis in amplitude/phase form. The filtered series is

$$y(t) = \overline{y} + \sum_{k=1}^{N/2}R(k)\,C_k\cos\left(\frac{2\pi kt}{T} - \phi_k\right)$$

where $R(k)$ is the **response function** of the filter,

$$R(k) = \frac{\text{amplitude of output at wavenumber }k}{\text{amplitude of input at wavenumber }k} = \frac{C(k)_{\text{filtered}}}{C(k)_{\text{unfiltered}}}$$

Generally $R$ is given as a function of frequency $R(\omega)$. In general $R(\omega)$ may be **complex**, since a filter can change both amplitude and phase. Mostly we care about amplitude, so we consider the real case:

$$R(\omega) = \frac{C(\omega)_{\text{filtered}}}{C(\omega)_{\text{original}}}, \qquad |R(\omega)|^2 = \frac{\Phi(\omega)_{\text{filtered}}}{\Phi(\omega)_{\text{original}}}$$

### Filtering in the frequency domain

Transform to the frequency domain (e.g. FFT), multiply $C(\omega)_{\text{original}}$ by your desired $R(\omega)$, then reconstruct. This is very effective and guarantees you keep exactly the frequencies you want.

Possible problems:

- **cannot be done in real time**
- a response function with a **sharp cutoff introduces spurious behaviour in the time domain** (and vice versa), particularly at the ends. So a boxcar filter on the FFT coefficients is not a good idea unless variance at higher frequencies is known to be very small

Note that selecting a single spectral estimate, say $k=3$, is applying a boxcar in the frequency domain of width one wavenumber.

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `seasonal_cycle_filtering_FFT.m`
:::

### Filtering in the time domain

Apply a weighting function (e.g. a weighted running average):

$$y(t) = \sum_{\tau=-J}^{J}g(\tau)\cdot x(t+\tau)$$

where $g(\tau)$ is the weighting function, $x(t)$ the original series, $y(t)$ the filtered series, and the length of $g$ is $2J+1$, i.e. $J = (L-1)/2$. Note that **$J$ data points are lost at each end**.

A moving average of width $L=3$ (the 1-1-1 filter):

$$y(t) = \tfrac13x(t-\tau) + \tfrac13x(t) + \tfrac13x(t+\tau)$$

and of width $L=5$ (1-1-1-1-1):

$$y(t) = \tfrac15x(t-2\tau) + \tfrac15x(t-\tau) + \tfrac15x(t) + \tfrac15x(t+\tau) + \tfrac15x(t+2\tau)$$

A moving average can only smooth, so it is **always a low-pass filter**. To get the high-pass component, subtract:

$$R(\omega)_H = 1 - R(\omega)_L, \qquad y(t)_{H} = x(t) - y(t)_{L}$$

and a bandpass is the difference of two low-pass series with different cutoffs:

$$y(t)_{3-6\,\text{days}} = y(t)_{L>3} - y(t)_{L>6}$$

**What does $g(\tau)$ do to the frequencies?** The filter is just a convolution,

$$y(t) = \int_{-\infty}^{\infty}x(t+\tau)g(\tau)\,d\tau = x * g$$

so by the convolution theorem

$$\mathcal{F}(y) = \mathcal{F}(x)\cdot\mathcal{F}(g) = X\cdot G$$

and therefore

$$\boxed{R(\omega) = \frac{Y_{\text{filtered}}(\omega)}{X(\omega)} = G(\omega)}$$

**The response function of a time-domain filter is the Fourier transform of its weights.** The weights and the response function are a Fourier transform pair. This single fact lets you design filters in whichever domain is more convenient.

**Example: the running mean.** $g(\tau) = 1/L$ for $-\frac{L-1}{2} \le \tau \le \frac{L-1}{2}$, zero elsewhere. Its response function is

$$R(\omega) = G(\omega) = \frac{\sin(\omega L/2)}{\omega L/2} = \text{sinc}\!\left(\frac{\omega L}{2\pi}\right)$$

— the same sinc we met as the boxcar *window* in Part 5, now appearing as a moving *average*. Note that $R(\omega) = 0$ when $\omega L = 2\pi, 4\pi, 6\pi\ldots$, i.e. at 1 cycle per $L$, 2 cycles per $L$, and so on.

So a 12-month running mean applied to monthly data exactly removes the 12-month cycle, the 6-month cycle, the 4-month cycle, the 2-month cycle… **But the problem is the side lobes of the sinc function**, which phase-shift frequencies between those perfectly removed harmonics, leaking power elsewhere.

**Example: the "perfect" response function.** A rectangle around the frequencies you wish to keep. In the time domain that window is again a sinc — extending over a wide range of the series, with negative lobes that *subtract* time steps. This causes severe problems at the ends of the record. As with windowing, the best solution is something in between: a tapered filter.

**Example: the Gaussian bell.** $f(x) = e^{-x^2}$ forms a Fourier transform pair with itself, so a Gaussian moving average in time gives a Gaussian response in frequency — no side lobes at all. Frequently used for spatial averaging, though the cutoff is not sharp.

### Creating your own non-recursive filters

In general: (1) specify a desired $R(\omega)$, then (2) Fourier transform it to get the time-domain weights $g(\tau)$.

**Applying filters in a row.** If you apply $w_1(\tau)$ and then $w_2(\tau)$ to the result,

$$\mathcal{F}(y_2) = \mathcal{F}(y)\cdot\mathcal{F}(w_1)\cdot\mathcal{F}(w_2)$$

so **the total response function is simply the product of the individual response functions.**

**The 1-2-1 filter.** A symmetric but non-rectangular moving average:

$$y(t) = \tfrac14x(t-\tau) + \tfrac12x(t) + \tfrac14x(t+\tau)$$

How do we get its response function? Either run it through an FFT, or use the formula below.

:::{admonition} Example / deeper dive
:class: note
**Deriving the response function of a symmetric filter**

For a symmetric filter $y(t) = \sum_{\tau=-J}^{J}C_\tau x(t-\tau)$ with $C_{-\tau} = C_\tau$:

$$Y(\omega) = \int_{-\infty}^{\infty}\left(\sum_{\tau=-J}^{J}C_\tau x(t-\tau)\right)e^{-i\omega t}dt$$

Take the sum outside the integral, since it is over $\tau$ and not $t$:

$$Y(\omega) = \sum_{\tau=-J}^{J}\left(C_\tau\int_{-\infty}^{\infty}x(t-\tau)e^{-i\omega t}dt\right)$$

The integrand almost looks like the transform of $x$, except for that pesky $\tau$. This is where the **Time Shifting Theorem** rescues us:

$$\mathcal{F}\big(x(t-a)\big) = \mathcal{F}\big(x(t)\big)\cdot e^{-i\omega a}$$

Thus

$$Y(\omega) = \left(\sum_{\tau=-J}^{J}C_\tau e^{-i\omega\tau}\right)X(\omega)$$

and the response function is

$$R(\omega) = \sum_{\tau=-J}^{J}C_\tau e^{-i\omega\tau}$$

Using the symmetry of the $C_\tau$ and $\cos x = \tfrac12(e^{ix}+e^{-ix})$, this simplifies to the form you will actually use:

$$\boxed{R(\omega) = C_0 + 2\sum_{\tau=1}^{J}C_\tau\cos(\omega\tau)}$$
:::

:::{admonition} Example / deeper dive
:class: note
**Example 4.3 — response functions of simple filter weights**

**Running mean 3 (1-1-1):** $C_{-1} = C_0 = C_{+1} = 1/3$

$$R(\omega) = \tfrac13 + 2\cdot\tfrac13\cos(\omega) = \tfrac13 + \tfrac23\cos(\omega)$$

**Running mean 5 (1-1-1-1-1):** all $C_\tau = 1/5$

$$R(\omega) = \tfrac15 + \tfrac25\cos(\omega) + \tfrac25\cos(2\omega)$$

**1-2-1 filter:** $C_{-1} = 1/4$, $C_0 = 1/2$, $C_{+1} = 1/4$

$$R(\omega) = \tfrac12 + 2\cdot\tfrac14\cos(\omega) = \tfrac12 + \tfrac12\cos(\omega)$$

Comparing these on a plot of response versus cycles per time step shows why the 1-2-1 filter is often preferred to a 3-point running mean of the same width: it has no negative lobes, so it never inverts the sign of a frequency component.
:::

**Lanczos smoothing.** A low-pass filter designed for as sharp a cutoff as possible, with tweaks to remove the harmonic lobes of the sinc function (**Gibbs phenomenon**). Writing $J$ weights that give the most abrupt cutoff,

$$R(\omega) = C_0 + 2\sum_{\tau=1}^{J}C_\tau\cos(\omega\tau), \qquad C_k = \frac{1}{\pi}\int_0^{\pi}\cos(k\omega')R(\omega')\,d\omega'$$

For an abrupt cutoff at $\omega = \alpha\pi$,

$$C_k = \frac{1}{k\pi}\sin(\alpha k\pi)$$

The sharper the cutoff, the more weights required — and in all cases you get wiggles from the sharp edges (Gibbs). The **Lanczos weights** are sinc functions applied to the ideal weights:

$$\tilde{C}_k = \text{sinc}\!\left(\frac{\pi k}{J}\right)C_k$$

allowing a small number of weights (fast computation) while suppressing Gibbs oscillations.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figures:* `15_boxcar_runmean_responsefcn.pdf`, Lanczos figures from Hartmann's notes.

*Worked example:* `response_functions_moving_agv.m`
:::

### Recursive filters

Recursive filters have properties that often make them ideal:

- less computation time, since you use already-filtered data rather than raw data at each step
- fewer weights required
- **can be done in real time**, using only previous data

However, construction is not always easy and does not always give stable results.

Non-recursive:

$$y(t) = \sum_{\tau=-J}^{J}g(\tau)\cdot x(t+\tau)$$

Recursive:

$$y(t) = \sum_{\tau=-J}^{0}b(\tau)\cdot x(t+\tau) + \sum_{\tau=-J}^{-1}a(\tau)\cdot y(t+\tau)$$

:::{admonition} Spin-up contamination
:class: warning
A critical aspect of a recursive filter is **how long the influence of a given data point takes to die away** — calculable via the *impulse response*. Because the filter feeds on its own output, the beginning of your series (while the filter spins up) is corrupted. Always discard the first $M$ data points, where $M$ is the number of steps for the impulse response to become negligible.
:::

**Butterworth filters** are a well-known family of recursive filters, designed to be smooth, to have high tangency at the origin and infinity, and to have as flat a passband as possible:

$$|R(\omega)|^2 = R(\omega)R^*(\omega) = \frac{1}{1 + \left(\dfrac{\omega}{\omega_c}\right)^{2N}}$$

where $N+1$ is the number of weights ($a$'s and $b$'s) and $\omega_c$ the cutoff frequency. As $N \to \infty$ the response approaches a boxcar. In Matlab, `butter.m` outputs the weights $a$ and $b$ for use with `filter.m`.

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `butterworth_example.m`
:::

---

## Part 7 — Cross-Spectrum Analysis

Cross-spectrum analysis is used when you have two time series and want their relationship **as a function of frequency**. Sometimes both series have spectral peaks near the same frequency and you want the phase relationship between them. But even when *neither* series has any spectral peak, coherent modes can be lurking in the background — cross-spectrum analysis can draw them out.

The goal: given $x(t)$ and $y(t)$, find the frequency and phase of their relationship. For example, for a propagating wave such as the MJO, winds at one location and another further east may oscillate at the same frequency but be 90° out of phase.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure:* slides `17_intro_cross_spectral_analysis.pdf`.
:::

### Cross-spectrum, co-spectrum and quadrature spectrum

To shorten notation define $\alpha_k = 2\pi kt/T$. From the DFT, a pair of series decompose as

$$x(t) = \overline{x} + \left[\sum_{k=1}^{N/2-1}A_{xk}\cos(\alpha_k) + \sum_{k=1}^{N/2-1}B_{xk}\sin(\alpha_k)\right] + A_{xN/2}\cos(\alpha_{N/2})$$

and similarly for $y$. Using the orthogonality of sines and cosines, the covariance of $x$ and $y$ is

$$\overline{x'y'} = \frac12\sum_{k=1}^{N/2-1}(A_{xk}A_{yk} + B_{xk}B_{yk}) + A_{xN/2}A_{yN/2} = \sum_{k=1}^{N/2}CO(k)$$

where $CO(k)$ is the **co-spectrum** of $x$ and $y$. In the special case $x=y$ this reduces to the ordinary power spectrum, $\overline{x'^2} = \sum_k C_k^2/2$.

Now define the complex quantities

$$F_x(k) = (A_{xk} - iB_{xk})e^{i\alpha_k} = C_{xk}e^{-i\theta_{xk}}e^{i\alpha_k}$$

and similarly $F_y(k)$. Using $e^{i\theta} = \cos\theta + i\sin\theta$, the real and imaginary parts are

$$F_x(k) = \underbrace{A_{xk}\cos(\alpha_k) + B_{xk}\sin(\alpha_k)}_{\text{Re}} + i\underbrace{\big(A_{xk}\sin(\alpha_k) - B_{xk}\cos(\alpha_k)\big)}_{\text{Im}}$$

so that $x = \overline{x} + \sum_k\text{Re}(F_x(k))$.

Define three products, with $*$ denoting complex conjugate:

$$F_{xx} = F_xF_x^*, \qquad F_{yy} = F_yF_y^*, \qquad F_{xy} = F_xF_y^*$$

Since $ZZ^* = (\text{Re}\,Z)^2 + (\text{Im}\,Z)^2$, working through the algebra gives

$$F_{xx}(k) = A_{xk}^2 + B_{xk}^2 = C_{xk}^2, \qquad F_{yy}(k) = C_{yk}^2$$

which is just twice the power spectrum. And for the cross term,

$$\boxed{F_{xy}(k) = \underbrace{(A_{xk}A_{yk} + B_{xk}B_{yk})}_{CO(k)} + i\underbrace{(A_{xk}B_{yk} - A_{yk}B_{xk})}_{Q(k)}}$$

In words:

- $F_{xy}(k)$ is the **cross-spectrum** of $x$ and $y$
- the real part is the **co-spectrum** $CO(k)$ — the **in-phase** signal
- the imaginary part is the **quadrature spectrum** $Q(k)$ — the **out-of-phase** signal
- the amplitude of the cross-spectrum is $|F_{xy}| = C_{xk}C_{yk}$, and the phase is $\theta_{xyk} = \theta_{yk} - \theta_{xk}$

### Coherence squared

$$\boxed{\text{Coh}^2(k) = \frac{|F_{xy}|^2}{F_{xx}F_{yy}}}$$

:::{admonition} Coherence at a single k is always exactly 1
:class: warning
Evaluate $\text{Coh}^2$ at one wavenumber and the algebra collapses:

$$\text{Coh}^2(k=k_1) = \frac{CO(k)^2 + Q(k)^2}{C_x^2C_y^2} = \frac{(A_xA_y + B_xB_y)^2 + (A_xB_y - A_yB_x)^2}{(A_x^2+B_x^2)(A_y^2+B_y^2)} = 1$$

Every cross-term cancels. **An unaveraged coherence carries no information whatsoever** — any two time series are perfectly coherent at every individual frequency. Coherence only becomes meaningful once averaged over adjacent wavenumbers or over multiple chunks.
:::

Averaged over two spectral coefficients $k_1$ and $k_2$ it becomes

$$\text{Coh}^2 = \frac{C_{x1}^2C_{y1}^2 + 2C_{x1}C_{y1}C_{x2}C_{y2}\cos(\theta_1-\theta_2) + C_{x2}^2C_{y2}^2}{C_{x1}^2C_{y1}^2 + C_{x1}^2C_{y2}^2 + C_{x2}^2C_{y1}^2 + C_{x2}^2C_{y2}^2}$$

Things to notice:

- in practice, for two unrelated series, $\text{Coh}^2$ **drops off rapidly** as it is averaged over a range of $k$ or over more than one subset
- coherence is largest when the phase difference is smallest, $\theta_1 = \theta_2$. If the phase relationships at $k=1$ and $k=2$ differ greatly, the average coherence is small — which makes physical sense, since we expect similar phase relationships for closely spaced wavenumbers
- even with identical phase difference, coherence is not automatically 1. It is maximized when $C_{x1}/C_{y1} = C_{x2}/C_{y2}$, i.e. when the **amplitude ratios** match

> The coherence will be 1 when realizations are averaged only under the conditions that the two realizations show the **same phase difference** between the two variables and the **same amplitude ratio** between them. The coherence thus shows how well these two conditions are satisfied.

**In practice**, having divided your data into $M$ chunks:

1. calculate the FFT of $x$ and $y$ to get the $A$'s and $B$'s
2. calculate $CO(k)$, $Q(k)$, $C_x(k)^2$ and $C_y(k)^2$ for each subset
3. compute

$$\text{Coh}^2(k) = \frac{\left(\frac{1}{M}\sum_{i=1}^{M}CO(k)\right)^2 + \left(\frac{1}{M}\sum_{i=1}^{M}Q(k)\right)^2}{\left(\frac{1}{M}\sum_{i=1}^{M}C_x^2(k)\right)\left(\frac{1}{M}\sum_{i=1}^{M}C_y^2(k)\right)}$$

Averaged this way, $\text{Coh}^2(k)$ gives the **fraction of the variance of $y(t)$ explainable using $x(t)$ as a predictor over that frequency band** — note how closely it parallels $R^2$ from Weeks 4-6:

$$R^2 = \frac{\left(\overline{x'y'}\right)^2}{\overline{x'^2}\cdot\overline{y'^2}}$$

Coherence squared is $R^2$, computed frequency by frequency.

:::{admonition} Degrees of freedom for coherence
:class: warning
The significance of $\text{Coh}^2$ depends on the dof of the individual spectral estimates. Most tables are given in terms of $M$ — the number of $k$ averaged over, or the number of chunks. So although the dof is roughly $2M$ as usual, **you should use only $M$ in these lookup tables.**
:::

### Getting the phase information

The phase of the cross-spectrum reveals the phase lag between $x$ and $y$ at a given $k$ — but **is only meaningful if the $\text{Coh}^2$ is meaningful**. By the Pythagorean theorem,

$$\boxed{\tan(\theta_{xk} - \theta_{yk}) = \frac{Q(k)}{CO(k)} = \frac{\text{Im}(F_{xy})}{\text{Re}(F_{xy})}}$$

:::{admonition} Figure / in-class demonstration
:class: tip
*Worked example:* `cross_spectrum_example.m`
:::

---

## Part 8 — Mixed Space-Time Analysis

Spectral analysis forms a family of four methods:

| Method | Variables | FFT dimensions |
|---|---|---|
| spectral analysis | one | one |
| cross-spectral analysis | two | one |
| mixed space/time spectral analysis | one | two |
| mixed space/time cross-spectral analysis | two | two |

The algebra is developed in a series of papers by Hayashi, Y. — *J. Meteor. Soc. Japan* (1971), *J. Atmos. Sci.* (1977), *J. Meteor. Soc. Japan* (1982).

### Space-time spectral analysis

In ordinary spectral analysis we decompose $y(t)$ into temporal harmonics, which in amplitude/phase form is

$$y(t) = \overline{y} + \sum_{k=1}^{N/2}C_k\cos(2\pi kt/T - \phi_k)$$

or, converting to frequency,

$$y(t) = \overline{y} + \sum_{\omega=2\pi/T}^{\pi}W_\omega\cos(\omega t + \phi_\omega)$$

In mixed space/time analysis we start from $y(x,t)$ and decompose into **both** temporal and spatial harmonics. The most common application uses frequency for time and waves along latitude circles for space. (A nice bonus: if you use latitude circles you do **not** need to window in that dimension, because the data is genuinely periodic.)

$$y(t,\lambda) = \overline{y} + \sum_{k=1}^{n/2}\ \sum_{\omega=\pm2\pi/T}^{\pi}W_{k,\pm\omega}\cos(\pm\omega t + \lambda k + \phi_{k,\pm\omega})$$

where $T$ is the record length, $\lambda$ is longitude in radians, $k$ is the zonal wavenumber, and $n$ is the number of grid points along a latitude circle.

**$\omega$ can be positive or negative at each $k$**, because the phase speed of a wave can be positive (westward) or negative (eastward), with

$$c_{ph} = \omega/k$$

This is the essential gain over a one-dimensional spectrum: a plain frequency spectrum tells you a 40-day oscillation exists, but cannot tell you whether it propagates east or west. Splitting $\omega$ into $\pm$ at each $k$ separates the two.

The goal is to find the weights $W(\omega, k)$. In steps:

1. **Decompose $y$ into zonal wavenumbers** $k = 1$ to $n/2$, giving sine and cosine functions at each time step. With $Y(\lambda,t) = (144, 3650)$ you obtain two matrices, $X_C(k,t) = (72,3650)$ and $X_S(k,t) = (72,3650)$.

2. **Fourier transform $X_C$ and $X_S$ in time:**

$$X_C(k,t) = \sum_\omega A_{k,\omega}\cos(\omega t) + B_{k,\omega}\sin(\omega t)$$
$$X_S(k,t) = \sum_\omega a_{k,\omega}\cos(\omega t) + b_{k,\omega}\sin(\omega t)$$

3. **Combine the coefficients.** They relate to $W_{k,\pm\omega}$ and $\phi_{k,\pm\omega}$ through

$$4W_{k,\pm\omega}^2 = (A \mp b)^2 + (\mp B - a)^2, \qquad \phi_{k,\pm\omega} = \arctan\left(\frac{\mp B - a}{A \mp b}\right)$$

and the power spectrum is

$$P_{k,\pm\omega} = \tfrac12W_{k,\pm\omega}^2$$

4. **Plot** the resulting power spectrum, conventionally with frequency on the $x$-axis and spatial wavenumber on the $y$-axis.

### Space-time cross-spectral analysis

The methodology is the same, but now investigating the relationship between two variables, say $u'$ and $v'$.

1. Calculate $W_{k,\pm\omega}$ and $\phi_{k,\pm\omega}$ for $X$ and $Y$ separately, giving $W^X$, $W^Y$, $\phi^X$, $\phi^Y$.

2. For each $(k,\omega)$ pair, the co-spectrum and quadrature spectrum are

$$Co_{k,\pm\omega} = \tfrac12W^X_{k,\pm\omega}W^Y_{k,\pm\omega}\cos\left(\phi^X_{k,\pm\omega} - \phi^Y_{k,\pm\omega}\right)$$
$$Q_{k,\pm\omega} = \tfrac12W^X_{k,\pm\omega}W^Y_{k,\pm\omega}\sin\left(\phi^X_{k,\pm\omega} - \phi^Y_{k,\pm\omega}\right)$$

3. Average the co- and quadrature spectra over a series of spatial, or more commonly temporal, harmonics — over subsets of the data or adjacent frequencies.

4. The coherence squared is then

$$\text{Coh}^2(k,\pm\omega) = \frac{\overline{Co}^2_{k,\pm\omega} + \overline{Q}^2_{k,\pm\omega}}{\overline{P^X_{k,\pm\omega}}\ \overline{P^Y_{k,\pm\omega}}}$$

:::{admonition} Example / deeper dive
:class: note
**Worked application: $\overline{u'v'}$ by phase speed, after Randel & Held (1991)**

Say you want winter-time plots of $\overline{u'v'}$ as a function of phase speed and latitude.

1. For each winter (DJFM) you have 120 days. These are your chunks.
2. Focus on a single latitude. Using NCEP data, $u = u(144,120)$ and $v = v(144,120)$.
3. FFT the data **in space** to get four matrices (two for $u$, two for $v$), each $(72,120)$ since 72 wavenumbers are resolved.
4. FFT the sine and cosine matrices **in time**. Now you have four matrices for $u$ and four for $v$ — a temporal sine and cosine matrix for each spatial sine matrix, and so on.
5. From the temporal cosine and sine matrices, solve for $W$ and $\phi$ for $u$ and $v$ separately.
6. Solve for $CO$ and $Q$ as above.
7. You now have $CO$ and $Q$ at a single latitude. Repeat for all latitudes to get $CO(\text{lat},k,\omega)$ and $Q(\text{lat},k,\omega)$.
8. Repeat for all winters (all chunks) and average the $CO$'s and $Q$'s to get mean $CO$ and mean $Q$.
9. Calculate $\text{Coh}^2$ from the averaged $CO$ and $Q$.
10. Plot either $\text{Coh}^2$, or the in-phase co-spectral power density

$$\text{CPSD}_{k,\pm\omega} = 2\,\text{Re}\left(\overline{F_{X,k,\pm\omega}F^*_{Y,k,\pm\omega}}\right) = 2\cdot\overline{CO}_{\text{lat},\pm\omega}$$

where the overline denotes averages over chunks or bandwidths. Note the CPSD **can be negative**, implying for example a negative zonal momentum flux.
11. You now have $\text{Coh}^2$ or CPSD as a function of latitude, $k$ and $\omega$. Convert to (latitude, phase speed) using $c_{ph} = \omega/k$.
:::

---

## Summary

The whole block rests on a small number of transferable facts:

| Concept | Statement |
|---|---|
| **Harmonic analysis** | multiple regression onto an orthogonal sine/cosine basis, so each coefficient is independent |
| **Nyquist** | set by $\Delta t$ alone; record length sets only the *low*-frequency limit |
| **Variance** | harmonic $k$ explains $C_k^2/2$; total area under the spectrum is the total variance |
| **Convolution theorem** | multiplication in one domain is convolution in the other |
| **Window / filter duality** | a filter's response function is the Fourier transform of its weights |
| **Coherence²** | $R^2$ computed frequency by frequency — meaningless until averaged |

Three cautions to carry into your own analysis:

1. **You are always looking through a window.** A finite record silently applies a boxcar, whose sinc response smears the true spectrum and adds side lobes. Choosing a Hanning window is not adding a distortion — it is replacing an unavoidable bad one with a better one.

2. **Bandpass filtering manufactures periodicity.** Filtered white noise oscillates beautifully at whatever band you chose. Test significance against red noise *before* filtering, and remember the a-posteriori penalty $(0.99)^{M/2}$ if you went looking for the largest peak.

3. **Every knob is the same trade-off.** Chunk length, band averaging, window choice, overlap — all of them buy degrees of freedom with resolution, or resolution with degrees of freedom. There is no setting that gives you both.

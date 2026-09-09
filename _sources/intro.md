# Syllabus
This is the course handout for Objective Analysis — a graduate-level introduction to the statistical and mathematical tools used in the atmospheric and oceanic sciences.

## How to read these notes

Throughout the notes, two types of highlighted boxes are used:

:::{admonition} Figure / in-class demonstration
:class: tip
**Green boxes** mark figures, in-class demonstrations, or Python/Matlab examples to be shown during lecture. They are placeholders for visual material that accompanies the written derivations.
:::

:::{admonition} Example / deeper dive
:class: note
**Blue boxes** contain worked examples, mathematical proofs, or deeper conceptual discussions that go beyond the main narrative. These are important — they are where most of the practice problems and physical intuition live.
:::

## Course Outline

The course runs over 16 weeks, with the mid-term in Week 7 and the final in Week 16. There is no class in the final-exam week.

| Weeks | Topic |
| --- | --- |
| 1-3 | {ref}`Rule 101 — foundations of statistics <Origin_of_Statisics>` |
| 4-6 | {ref}`Regression & AR1 <Regression>` |
| **7** | **Mid-term exam** |
| 8-11 | {ref}`Seeking structure in data — EOFs & clustering <Seeking_Structure>` |
| 12-15 | {ref}`Time series analysis — spectra & filtering <Time_Series>` |
| **16** | **Final exam** |

### Part I: Foundations of Statistics

**{ref}`Origin_of_Statisics`** — Week 1-3: Rule 101

* **The origin of statistics: why the bell curve?**
* **Mean, variance and higher moments**
* **Basic probabilities; unions; intersections; conditional probabilities; Bayes theorem**
* **Tying it together: Moments, PDFs, and Bayes' Theorem**
* **Statistical significance testing**
* **Hypothesis testing**
* **Monte Carlo and resampling techniques**
* **Compositing**
* **Other common distributions**
* **Non-parametric tests**

### Part II: Regression & Autocorrelation

**{ref}`Regression`** — Week 4-6: Regression & AR1

* **Linear regression: least squares, slope & intercept**
* **Theory of correlation: Pearson's r, Fisher-Z, Spearman's rank**
* **Autocorrelation & effective sample size: AR1/red noise, Leith & Bretherton**
* **Multiple regression: generalized normal equations, overfitting**

---

**Week 7 — Mid-term exam**, covering Parts I and II.

---

### Part III: Seeking Structure in Data

**{ref}`Seeking_Structure`** — Week 8-11: EOFs & Clustering

* **Linear algebra review: inner products, covariance matrices, inverse, rank & null space**
* **Eigenvalues and eigenvectors: diagonalizing the covariance matrix**
* **EOFs via eigenanalysis: maximizing explained variance, orthogonality, principal components**
* **EOFs via Singular Value Decomposition, and its equivalence to eigenanalysis**
* **EOFs with real data: weighting, standardization, presentation in physical units**
* **How many EOFs to retain: North et al. (1982) and degeneracy**
* **Cluster analysis: k-means and Lloyd's algorithm**
* **Self-organizing maps (SOMs): training, mapping, quantization & topographic error**

### Part IV: Time Series Analysis

**{ref}`Time_Series`** — Week 12-15: Spectral Analysis & Filtering

* **Harmonic analysis: Fourier sums, the Nyquist frequency and aliasing**
* **The discrete Fourier transform as orthogonal multiple regression**
* **The power spectrum: line vs. continuous, and the resolution–reliability trade-off**
* **The complex Fourier transform**
* **Significance of spectral peaks: red-noise null, the F-test, a priori vs. a posteriori**
* **Windows and finite data: convolution, the convolution theorem, boxcar & Hanning, WOSA**
* **Filtering: response functions, non-recursive & Lanczos, recursive & Butterworth**
* **Cross-spectrum analysis: co-spectrum, quadrature spectrum, coherence² and phase**
* **Mixed space-time analysis: propagating waves and phase speed**

---

**Week 16 — Final exam.** No class this week.

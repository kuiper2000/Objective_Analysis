# Syllabus
This is the course handout for Objective Analysis — a graduate-level introduction to the statistical and mathematical tools used in the atmospheric and oceanic sciences.

## Course Outline

### Part I: Foundations of Statistics

**{ref}`Origin_of_Statisics`** — Week 1-4: Rule 101

* **The origin of statistics: why the bell curve?**
  * A practical problem: combining noisy measurements
  * Setting up the likelihood
  * The key constraint: forcing the mean to be optimal
  * From the constraint to an ODE — deriving the Gaussian

* **Mean, variance and higher moments**
  * Sample mean, median, and their properties
  * Sample variance and the N−1 correction
  * Higher moments: skewness and kurtosis

* **Basic probabilities; unions; intersections; conditional probabilities; Bayes theorem**
  * Unions and intersections; Venn diagrams
  * Conditional probability and the multiplicative law
  * Bayes' theorem: worked examples (chemical detection, cab accidents, Monty Hall problem)
  * Probability philosophy: frequentist vs. Bayesian approaches

* **Tying it together: Moments, PDFs, and Bayes' Theorem**
  * Moments as expectations over a PDF — the population limit
  * From the product rule of probability to Bayes' theorem in continuous form
  * Worked example: Gaussian prior × Gaussian likelihood → Gaussian posterior
  * Physical interpretation: precision-weighted averaging and the Kalman Filter

* **Statistical significance testing**
  * The normal distribution, PDF, and CDF
  * The z-statistic and z-score
  * Estimating the significance of the sample mean
  * Confidence intervals
  * The Central Limit Theorem
  * The t-statistic: when sample sizes are small
  * When to use (and not use) the t-test
  * A note on independence and effective sample size

* **Hypothesis testing**
  * Terminology: significance level, critical value, p-value
  * Five-step framework for hypothesis testing
  * Comparison of means
  * Type I and Type II errors
  * A priori vs. a posteriori statistics; the multiple-testing problem
  * Bayesian vs. frequentist approaches revisited

* **Monte Carlo and resampling techniques**
  * Why use resampling and Monte Carlo?
  * Bootstrap resampling
  * Jackknife resampling
  * Monte Carlo simulation

* **Compositing (superposed epoch analysis)**
  * Steps to compositing; advantages and disadvantages
  * Significance of composites

* **Other common distributions**
  * Chi-square distribution: tests of variance
  * F-statistic: comparing two sample variances
  * Binomial distribution and its normal approximation

* **Non-parametric tests**
  * Signs test (Wilcoxon test)
  * Runs test (Wald–Wolfowitz test)

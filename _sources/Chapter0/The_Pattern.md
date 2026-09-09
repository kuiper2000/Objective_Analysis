(The_Pattern)=
# Chapter 0: The Idea Behind the Course

:::{admonition} How to read this chapter
:class: warning
This chapter states, in one place, the single idea that every later part of the course is an instance of. It is deliberately the most abstract thing you will read all semester.

**You are not expected to absorb it on first reading.** Nothing in Weeks 1-3 depends on it. Skim it now so the shape is familiar, then come back at the end of each Part — there is a short *"The pattern"* box at the head of every Part that points back here and fills in that Part's row of the table below. The idea is meant to accumulate, not to land all at once.

**There are no prerequisites.** The one tool everything runs on — the **Lagrange multiplier** — is explained from scratch in the interlude below, in plain language, before it is used. Every derivation here is then done in ordinary partial derivatives. Where the continuum version of an argument would need calculus of variations, it is quarantined in a blue box you can skip.
:::

## Two questions, one answer

Almost everything we do this semester is one of two questions:

1. **Which probability distribution should I assume?** You want to test whether a trend is significant, or whether two composites differ. Before you can test anything you must commit to a PDF, because the null hypothesis *is* a PDF — it says "here is the distribution the data would follow if nothing interesting were happening," and the test asks which aspect of that distribution the data contradicts.

2. **Which model best represents the data?** You want a line through a scatterplot, a spatial pattern that explains the most variance, or a set of Fourier coefficients. In each case there are infinitely many candidate answers and you must pick one.

These look like different problems. They are the same problem, and it has the same shape in both cases:

> **Write down a quantity that measures how good a candidate answer is. Then find the candidate that optimizes it, subject to whatever you are unwilling to give up.**

The quantity you optimize is a **functional** — a function whose input is itself a whole function (a PDF, a fitted curve, a spatial pattern) and whose output is a single number. The things you are unwilling to give up are **constraints**. The machinery that handles "optimize this, subject to that" is the **method of Lagrange multipliers**.

That word *functional* sounds like it demands new mathematics. It does not, and you will never have to handle one abstractly here: chop the function into finitely many values and a functional becomes an ordinary function of finitely many variables, at which point the ordinary multiplier recipe applies. That is exactly how the derivation below proceeds.

That is the whole idea. The rest of this chapter shows it working twice, and the rest of the course is the same move applied to progressively richer objects.

Before either, though, we need the one tool both of them use.

---

## Interlude: what a Lagrange multiplier actually is

Every derivation in this chapter — and the EOF derivation in Part III — runs on Lagrange multipliers. If the phrase means nothing to you, or means only "a symbol you introduce and then cancel," this section is the one to read properly. Everything else follows easily once this is clear.

### The problem it solves

Maximizing something is easy when you are free to try anything: differentiate, set to zero, done. The hard case is when you are **not** free — when the answer must also satisfy some requirement. Find the highest point on a hill, *but you must stay on the marked path*. Fit the data as well as possible, *but the pattern must have unit length*. Spread the probability out as much as possible, *but the variance must come out to $\sigma^2$*.

Formally: maximize $f(x,y)$ subject to $g(x,y) = c$.

### The picture: contour lines and a fence

Think of $f$ as the height of a hill, drawn as contour lines — curves along which the height is constant. The constraint $g = c$ is a fence you must walk along.

Now walk. **As long as your fence *crosses* a contour line, you are still going uphill or downhill**, so you cannot be at the highest reachable point yet — take another step and you improve.

So at the best point, the fence cannot cross a contour. It must just **graze** it — the fence and the contour are *tangent*. And two curves are tangent exactly when their gradients point along the same line:

$$\nabla f = \lambda\,\nabla g$$

**That equation is the whole method, and $\lambda$ is nothing more exotic than the number that makes the two arrows match.** The gradients must be *parallel*; $\lambda$ says by what factor they differ in length.

### What you actually do

Rather than reason geometrically each time, bundle the condition into a single object — the **Lagrangian** — by subtracting the constraint with a multiplier attached:

$$\mathcal{L}(x,y,\lambda) = f(x,y) - \lambda\big(g(x,y) - c\big)$$

Then set *all* of its partial derivatives to zero, treating $\lambda$ as just another variable:

- $\partial\mathcal{L}/\partial x = 0$ and $\partial\mathcal{L}/\partial y = 0$ together give $\nabla f = \lambda\nabla g$, the tangency condition
- $\partial\mathcal{L}/\partial\lambda = 0$ gives back $g = c$, the constraint itself

So the trick is this: **you trade a constrained problem in $n$ variables for an unconstrained one in $n+1$.** Paying one extra unknown buys you the right to differentiate freely, as if the constraint were not there. That is why it is worth doing.

:::{admonition} Example / deeper dive
:class: note
**A two-line example**

Maximize $xy$ subject to $x + y = 10$. (Largest rectangle for a fixed perimeter.)

$$\mathcal{L} = xy - \lambda(x + y - 10)$$

$$\frac{\partial\mathcal{L}}{\partial x} = y - \lambda = 0, \qquad \frac{\partial\mathcal{L}}{\partial y} = x - \lambda = 0, \qquad \frac{\partial\mathcal{L}}{\partial\lambda} = -(x+y-10) = 0$$

The first two give $x = y = \lambda$; the third then gives $x = y = 5$, so $\lambda = 5$ and the maximum is $25$.

Notice we never solved the constraint for one variable and substituted. For this problem you easily could; for the ones later in this chapter — where the "variables" are the values of a probability density at every point — you cannot, and the multiplier is the only practical route.
:::

### What $\lambda$ means: an exchange rate

Here is the part usually left out, and it is the part that matters most for this course.

$\lambda$ is not scaffolding to be discarded once the algebra is done. It has a meaning: **$\lambda$ is the rate at which the best achievable value improves if you loosen the constraint.**

$$\boxed{\frac{df^*}{dc} = \lambda}$$

where $f^*$ is the maximum attainable given the constraint level $c$. Economists call this a **shadow price** — *"how much more profit would one extra dollar of budget buy me?"* If the constraint is a budget, $\lambda$ is what one more unit of budget is worth.

Check it on the example above. Maximizing $xy$ subject to $x+y=c$ gives $x=y=c/2$ and a maximum of $c^2/4$. Then

$$\frac{d}{dc}\left(\frac{c^2}{4}\right) = \frac{c}{2} = 5 \quad\text{at } c = 10$$

which is exactly the $\lambda = 5$ we found.

:::{admonition} Why the multiplier keeps turning out to be the thing you wanted
:class: seealso
Three times in this course a Lagrange multiplier will come back as a quantity with clear physical meaning:

| Constraint | Multiplier turns out to be |
| --- | --- |
| variance $= \sigma^2$ (Chapter 0) | $\lambda_2 = 1/2\sigma^2$, the inverse variance |
| $\mathbf{e}^T\mathbf{e} = 1$ (Part III, EOFs) | $\lambda$ = the eigenvalue = variance explained |

This is not luck. The exchange-rate property *guarantees* it: $\lambda$ is the derivative of the thing you maximized with respect to the thing you constrained, so it is automatically a meaningful rate.

You can verify the first row right now. The entropy of a Gaussian is $S = \tfrac12\ln(2\pi e\sigma^2)$, so the rate at which entropy grows as you permit more variance is

$$\frac{dS}{d(\sigma^2)} = \frac{1}{2\sigma^2}$$

which is precisely the $\lambda_2$ that will fall out of the derivation below. **So when a multiplier appears in this course, always ask what it turned out to be — it is usually the answer, not the bookkeeping.**
:::

---

## Question 1: which PDF? — maximize entropy

Suppose all you know about a quantity is its mean and its variance. Infinitely many distributions have that mean and that variance. Which should you assume?

The answer that has survived is: **assume the one that commits you to as little as possible beyond what you actually know.** Any other choice smuggles in structure you have no evidence for. The measure of "how little is committed" is the **entropy** of the distribution,

$$S[f] = -\int f(x)\ln f(x)\,dx$$

which is large when $f$ is spread out and noncommittal, and small when $f$ is sharply peaked and therefore highly specific. So the recipe is: **maximize $S[f]$ subject to the constraints that encode what you actually know.**

### Step 1: chop the axis into bins

"Maximize over all functions $f$" sounds like it needs machinery you may not have met. It does not. **Discretize first, and it becomes an ordinary Lagrange multiplier problem of the kind you already know.**

Chop the $x$-axis into $N$ narrow bins, and let $p_i$ be the probability of landing in bin $i$ (centred at $x_i$). The entropy is now an ordinary function of $N$ ordinary variables,

$$S(p_1,\dots,p_N) = -\sum_{i=1}^{N}p_i\ln p_i$$

and the three things we know become three ordinary constraints:

$$\sum_i p_i = 1, \qquad \sum_i x_i p_i = \mu, \qquad \sum_i (x_i-\mu)^2 p_i = \sigma^2$$

### Step 2: attach one multiplier per constraint

This is exactly the recipe from multivariable calculus — one multiplier for each thing you are unwilling to give up:

$$\mathcal{L}(p_1,\dots,p_N) = -\sum_i p_i\ln p_i - \lambda_0\Big(\sum_i p_i - 1\Big) - \lambda_1\Big(\sum_i x_ip_i - \mu\Big) - \lambda_2\Big(\sum_i (x_i-\mu)^2p_i - \sigma^2\Big)$$

### Step 3: differentiate with respect to one bin

Take an ordinary partial derivative with respect to a single $p_j$. Only the $i = j$ term of each sum survives, and $\frac{d}{dp}(p\ln p) = \ln p + 1$, so

$$\frac{\partial\mathcal{L}}{\partial p_j} = -\ln p_j - 1 - \lambda_0 - \lambda_1 x_j - \lambda_2(x_j-\mu)^2 = 0$$

Rearranging and exponentiating,

$$p_j = \exp\left(-1-\lambda_0-\lambda_1x_j-\lambda_2(x_j-\mu)^2\right)$$

**That is the whole derivation.** Nothing beyond a partial derivative and the multiplier recipe. Now shrink the bins: $p_j \to f(x)\,dx$, the index $j$ becomes the continuous label $x$, and

$$f(x) = \exp\left(-1-\lambda_0-\lambda_1x-\lambda_2(x-\mu)^2\right)$$

:::{admonition} Example / deeper dive
:class: note
**What the continuum version is really doing**

Taken directly in the continuum, the object being maximized is a **functional** — its input is a whole function and its output is one number — and the derivative you need is the **functional derivative** $\delta\mathcal{L}/\delta f$. The bin picture is the honest intuition for it: $\delta/\delta f$ is nothing more than $\partial/\partial p_j$ in the limit of infinitely many bins.

If you want the continuum argument on its own terms: write everything under one integral with integrand

$$G(f,x) = -f\ln f - \lambda_0 f - \lambda_1 xf - \lambda_2(x-\mu)^2f$$

perturb the maximizer by an arbitrary test function, $f \to f + \varepsilon\eta$, and demand stationarity:

$$\left.\frac{d}{d\varepsilon}\mathcal{L}[f+\varepsilon\eta]\right|_{\varepsilon=0} = \int\frac{\partial G}{\partial f}\,\eta(x)\,dx = 0$$

This has to hold for *every* $\eta$, and the **fundamental lemma of the calculus of variations** says that if $\int g(x)\eta(x)dx = 0$ for all smooth $\eta$, then $g\equiv0$. Hence $\partial G/\partial f = 0$ pointwise, which is the same condition Step 3 produced.

Note that it is the multipliers that let $\eta$ be *arbitrary*. Without them you could only perturb in directions that preserve normalization and both moments, and the argument would be far messier. That is what buying the multipliers actually buys you.

**One loose end in the bin argument.** The discrete and continuous entropies are not quite equal. Putting $p_i = f(x_i)\Delta x$,

$$-\sum_i p_i\ln p_i = -\int f\ln f\,dx - \ln\Delta x$$

and that $-\ln\Delta x$ blows up as the bins shrink. It is harmless here because it is a *constant* — it depends only on the bin width, not on the shape of $f$ — so it shifts the value of the entropy without moving the location of its maximum. The two routes therefore find the same $f$. (This offset is why the continuous quantity is properly called *differential* entropy and, unlike the discrete version, can be negative.)
:::

### Step 4: fix the multipliers using the constraints

We have the shape but not yet the numbers. Three multipliers, three constraints.

**The linear term vanishes.** Complete the square in $u = x-\mu$:

$$-\lambda_1 u - \lambda_2u^2 = -\lambda_2\left(u + \frac{\lambda_1}{2\lambda_2}\right)^2 + \frac{\lambda_1^2}{4\lambda_2}$$

so $f$ is a bell curve centred at $x = \mu - \lambda_1/(2\lambda_2)$. But the mean constraint says the centre is at $\mu$. That **forces** $\lambda_1 = 0$ — it is a consequence of the constraints, not an assumption. Absorbing the remaining constants into one prefactor $A$,

$$f(x) = A\exp\left(-\lambda_2(x-\mu)^2\right)$$

**The other two follow from two standard integrals:**

$$\int_{-\infty}^{\infty}e^{-au^2}du = \sqrt{\frac{\pi}{a}}, \qquad \frac{\int_{-\infty}^{\infty}u^2e^{-au^2}du}{\int_{-\infty}^{\infty}e^{-au^2}du} = \frac{1}{2a}$$

The second is just the variance of $e^{-au^2}$, so the variance constraint reads $\sigma^2 = 1/(2\lambda_2)$, giving

$$\lambda_2 = \frac{1}{2\sigma^2}$$

and normalization then gives $A = \sqrt{\lambda_2/\pi} = 1/(\sigma\sqrt{2\pi})$. Therefore

$$\boxed{f(x) = \frac{1}{\sigma\sqrt{2\pi}}\exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)}$$

:::{admonition} The prediction from the interlude, confirmed
:class: seealso
The exchange-rate argument said $\lambda_2$ *had* to equal $dS/d(\sigma^2) = 1/2\sigma^2$. The derivation just produced exactly that, without ever using the shadow-price property. The multiplier was never bookkeeping.
:::

:::{admonition} Two principles, one curve
:class: note
This is the **same Gaussian** that Week 1-3 derives from Gauss's 1809 argument — and the two derivations share no assumptions at all.

- **Gauss:** *postulate that the sample mean is the maximum-likelihood estimate of the true value.* That forces $f'/f$ to be linear, which is a differential equation, whose solution is the Gaussian.
- **Maximum entropy:** *postulate nothing beyond a known mean and variance, and assume as little as possible otherwise.* That is a constrained optimization, whose solution is the Gaussian.

One argument is about estimation, the other about ignorance. They converge on the same curve. This is the strongest reason to regard the normal distribution as genuinely special rather than merely convenient — and it is why so much of the course is built on it.
:::

### Other constraints, other distributions

The Gaussian is not privileged by the method — it is what you get from *those particular* constraints. Change what you know and the same procedure returns a different distribution:

| What you know | MaxEnt distribution |
|---|---|
| support is bounded, $x\in[a,b]$; nothing else | **uniform** |
| $x \ge 0$, mean $\mu$ fixed | **exponential**, $f = \frac{1}{\mu}e^{-x/\mu}$ |
| $x\in\mathbb{R}$, mean and **variance** $E[(x-\mu)^2]$ fixed | **Gaussian** |
| $x\in\mathbb{R}$, mean and **absolute deviation** $E\lvert x-\mu\rvert$ fixed | **Laplace**, $f \propto e^{-\lvert x-\mu\rvert/b}$ |
| $x \in \{0,1,2,\dots\}$, mean fixed | **geometric** |

Hold onto the Gaussian/Laplace pair — the difference between constraining the *squared* deviation and the *absolute* deviation is about to reappear as the difference between two loss functions.

---

## Question 2: which model? — minimize a loss

Now the second question. You have data and a family of candidate models, and you need to pick one. The engineering move is to define a **loss function** measuring how bad a candidate is, then choose the candidate that minimizes it.

This looks like a different kind of activity — pragmatic, not principled. It is not. **A loss function is a likelihood in disguise.**

Suppose your model makes predictions $\hat{y}_i$ and the errors $e_i = y_i - \hat{y}_i$ are independent draws from some density $f$. The likelihood of your data is $L = \prod_i f(e_i)$, and maximizing it is the same as minimizing

$$-\ln L = -\sum_i \ln f(e_i)$$

So the loss function **is** the negative log-likelihood, which means **choosing a loss is choosing an error distribution**. Work it through for the two densities highlighted above:

$$f(e) \propto e^{-e^2/2\sigma^2} \;\Longrightarrow\; -\ln L \propto \sum_i e_i^2$$

$$f(e) \propto e^{-\lvert e\rvert/b} \;\Longrightarrow\; -\ln L \propto \sum_i \lvert e_i\rvert$$

which gives the correspondence:

| Loss | Implied error PDF | Optimal estimator | Behaviour |
|---|---|---|---|
| $\sum e_i^2$ (L2) | Gaussian | **mean** / ordinary least squares | efficient, but outliers dominate |
| $\sum \lvert e_i\rvert$ (L1) | Laplace | **median** / least absolute deviations | robust to outliers, harder to solve |

:::{admonition} Why we square the errors — the honest answer
:class: note
Week 4-6 introduces least squares with the remark that squaring keeps errors positive and makes the minimization linear. Both are true, and both are conveniences rather than reasons.

The real reason is the table above. **Squaring the errors is equivalent to asserting that the errors are Gaussian.** That assertion is usually reasonable — by the Central Limit Theorem, errors that accumulate from many small independent causes tend to be Gaussian — but it is an assumption about your data, not a mathematical necessity, and it is the assumption that fails when you have outliers.

This also explains a fact that would otherwise look like coincidence: the mean is the estimator that minimizes squared error, and the Gaussian is the maximum-entropy distribution for fixed *squared* deviation. The median is the estimator that minimizes absolute error, and the Laplace is the maximum-entropy distribution for fixed *absolute* deviation. Mean-with-Gaussian and median-with-Laplace are two self-consistent packages, and the loss you pick selects which one you are in.
:::

---

## The same move, everywhere

Both questions turned out to be: *optimize a functional subject to constraints.* So does everything that follows.

| Part | Functional optimized | Subject to | Solution |
|---|---|---|---|
| **1-3** Statistics | entropy $-\int f\ln f$ | known moments | the PDF you assume |
| **1-3** Statistics | likelihood $\prod f(x_i-\mu)$ | mean is optimal | the Gaussian |
| **4-6** Regression | loss $=-\ln L$ | choice of model family | regression coefficients |
| **8-11** EOFs | variance explained $\mathbf{e}^T\mathbf{C}\mathbf{e}$ | $\mathbf{e}^T\mathbf{e} = 1$ | eigenvectors of $\mathbf{C}$ |
| **8-11** Clustering | within-cluster sum of squares | exactly $k$ clusters | cluster centroids |
| **12-15** Fourier | squared projection onto each harmonic | orthonormal basis | Fourier coefficients |

:::{admonition} Example / deeper dive
:class: note
**A preview: the Lagrange multiplier you will meet again**

The EOF row of that table is worth doing now, because it shows the multiplier is not just bookkeeping.

In Weeks 8-11 we look for the spatial pattern $\mathbf{e}$ explaining the most variance. "Most variance" means maximizing $\mathbf{e}^T\mathbf{C}\mathbf{e}$, where $\mathbf{C}$ is the covariance matrix. But that can be inflated without limit just by lengthening $\mathbf{e}$, so we constrain $\mathbf{e}$ to unit length. Same shape as before — attach a multiplier:

$$\mathcal{L}(\mathbf{e}) = \mathbf{e}^T\mathbf{C}\mathbf{e} - \lambda\left(\mathbf{e}^T\mathbf{e} - 1\right)$$

Differentiate and set to zero:

$$\frac{\partial\mathcal{L}}{\partial\mathbf{e}} = 2\mathbf{C}\mathbf{e} - 2\lambda\mathbf{e} = 0 \qquad\Longrightarrow\qquad \boxed{\mathbf{C}\mathbf{e} = \lambda\mathbf{e}}$$

That is the eigenvalue equation. **The Lagrange multiplier turns out to be the eigenvalue — which is the variance explained by that pattern.** The multiplier is not an algebraic device to be discarded at the end; it is the answer you were looking for.

You will meet exactly the same structure in the Fourier chapter, where the "constraint" is that the basis functions are orthonormal, which is what collapses the multiple-regression normal equations into one-line formulas for the coefficients.
:::

---

## Where the pattern stops

A framing is only useful if you know its edges. Two things this one does **not** cover.

**Derived distributions are not chosen.** Maximum entropy explains the distributions you *assume* — Gaussian, exponential, uniform, Laplace. It does not explain the $t$, $\chi^2$ and $F$ distributions you will meet in Weeks 1-3. Those are not separate maximum-entropy solutions and nobody chose them: they are **consequences** of assuming a Gaussian and then sampling it finitely. The $t$ distribution is what you get when you estimate $\sigma$ from the same small sample you are testing; $\chi^2$ is what a sum of squared Gaussians does; $F$ is a ratio of two such sums. Keep the categories separate:

- **distributions you choose** — fixed by what you know, via maximum entropy
- **distributions that follow** — forced on you by the sampling, given that choice

**Non-parametric methods refuse the question.** Bootstrap, jackknife and Monte Carlo tests are not an odds-and-ends appendix to Weeks 1-3. They are the deliberate opposite pole of this whole chapter: instead of choosing a functional and deriving a distribution, they **decline to assume any PDF at all** and build the null distribution by resampling the data itself. That is the right move when you have no basis for the constraints maximum entropy needs, and the price is that you need enough data for the resampling to be meaningful.

---

## What to carry forward

Three sentences:

1. **The null hypothesis is a PDF, and you can derive rather than guess it** — maximize entropy subject to what you actually know.
2. **A loss function is a negative log-likelihood, so choosing a loss is choosing an error distribution** — L2 asserts Gaussian errors, L1 asserts Laplace.
3. **Whenever you want "the best" of something subject to a restriction, attach a Lagrange multiplier** — and check what the multiplier turns out to mean, because it is often the quantity you actually wanted.

Everything else this semester is these three sentences applied to bigger objects.

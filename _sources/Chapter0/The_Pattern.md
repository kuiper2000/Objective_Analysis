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

That is the whole idea. The rest of this chapter shows it working twice, then shows the two answers merging into one — the merged version is **Bayes' theorem**. After that, the rest of the course is the same move applied to progressively richer objects.

Before either, though, we need the one tool both of them use.

---

## Interlude: what a Lagrange multiplier actually is

Every derivation in this chapter — and the EOF derivation in Weeks 8-11 — runs on Lagrange multipliers. If the phrase means nothing to you, or means only "a symbol you introduce and then cancel," this section is the one to read properly. Everything else follows easily once this is clear.

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
| $\sum_i q_i = 1$ when updating a prior (Chapter 0, Bayes) | $e^{1+\lambda_0} = p(\text{data})$, the **evidence** |
| $\mathbf{e}^T\mathbf{e} = 1$ (Weeks 8-11, EOFs) | $\lambda$ = the eigenvalue = variance explained |

This is not luck. The exchange-rate property *guarantees* it: $\lambda$ is the derivative of the thing you maximized with respect to the thing you constrained, so it is automatically a meaningful rate.

You can verify the first row right now. The entropy of a Gaussian is $S = \tfrac12\ln(2\pi e\sigma^2)$, so the rate at which entropy grows as you permit more variance is

$$\frac{dS}{d(\sigma^2)} = \frac{1}{2\sigma^2}$$

which is precisely the $\lambda_2$ that will fall out of the derivation below. **So when a multiplier appears in this course, always ask what it turned out to be — it is usually the answer, not the bookkeeping.**
:::

### But why does *defining* $\mathcal{L}$ get us what we want?

This is the step that most often feels like a sleight of hand, so it is worth saying plainly what is and is not happening.

**Nothing is proved by writing down $\mathcal{L}$.** It is a container, not an argument. You already believe two things about the answer — tangency, $\nabla f = \lambda\nabla g$, and feasibility, $g = c$. That is $n+1$ equations in $n+1$ unknowns ($x_1,\dots,x_n$ and $\lambda$), and you could solve them and never write a Lagrangian at all.

The object

$$\mathcal{L} = f - \lambda(g-c)$$

is built precisely so that differentiating it hands those same two facts back to you — the $x_j$ derivatives give the components of $\nabla f = \lambda\nabla g$, and the $\lambda$ derivative gives $g = c$. Note *why* that last one is so simple: $\mathcal{L}$ is **linear** in $\lambda$, so $\partial\mathcal{L}/\partial\lambda$ is nothing but $-(g-c)$. That is not deep; it is engineered. The payoff is only that you no longer have to remember two different *kinds* of condition, one geometric and one algebraic — you turn the same crank on every variable, $\lambda$ included.

**Where the tangency condition actually comes from.** The fence-and-contour picture above is correct but is the two-dimensional shadow of a cleaner statement, and it is the cleaner statement that survives into the $N$-bin problem below, where nothing can be drawn:

At the best point, $f$ must stop changing **in every direction you are allowed to move** — every direction tangent to the surface $g=c$. "No change along $\mathbf{v}$" means $\nabla f\cdot\mathbf{v}=0$. So $\nabla f$ is perpendicular to every tangent direction. But $\nabla g$ is perpendicular to that same surface too, since a gradient is always normal to its own level set. The directions perpendicular to the surface form a single line, and two vectors on one line are multiples of each other:

$$\nabla f = \lambda\,\nabla g$$

In words: **$\nabla f$ has no component along the permitted directions, so all of it lies along the forbidden one — which is where $\nabla g$ already points.** $\lambda$ is only the ratio of their lengths.

**What the multiplier actually buys.** Without it you are *not entitled* to set $\partial f/\partial x_j = 0$, because the $x_j$ are not free: nudging one forces another to move to stay on the constraint. Plain stationarity is simply false for a constrained problem. The term $-\lambda(g-c)$ is exactly the correction that removes the component of the gradient pointing off the constraint surface, so that after the subtraction the "as if unconstrained" condition $\partial\mathcal{L}/\partial x_j = 0$ becomes legitimate — it now encodes only the tangential part, which is the part that genuinely has to vanish. That is the trade stated above: one extra unknown in exchange for the right to differentiate every variable independently.

:::{admonition} Example / deeper dive
:class: note
**A reading that makes the definition feel inevitable**

If $\mathcal{L}$ still looks pulled from a hat, read $-\lambda(g-c)$ as a **charge for violating the constraint**, at price $\lambda$ per unit. Then $\mathcal{L}$ is simply *payoff minus bill*.

- Price too low → cheating is cheap, and the free optimizer of $\mathcal{L}$ overshoots, $g > c$.
- Price too high → cheating is over-penalized, and it undershoots, $g < c$.
- There is exactly one price $\lambda^*$ at which the **unconstrained** optimizer of $\mathcal{L}$ lands on $g=c$ of its own accord.

At that price the constraint stops being a wall and becomes a cost you have already internalized, so you may optimize freely and the constraint enforces itself. The equation $\partial\mathcal{L}/\partial\lambda = 0$ is precisely the statement "$\lambda$ has been tuned to that price." This is also why $\lambda$ *has* to be the shadow price $df^*/dc$ of the previous section — by construction it is a price. The two readings are one fact.

**One caveat.** You are finding a **stationary point** of $\mathcal{L}$, not a maximum of it. Because $\mathcal{L}$ is linear in $\lambda$ it is unbounded in that direction — the solution is a saddle (a maximum over $x$, a minimum over $\lambda$), never a maximum in all the variables at once. So "maximize the Lagrangian" is sloppy; the honest recipe is *set all partials to zero, then confirm separately that the stationary point is the extremum you wanted.* In the entropy derivation below, the concavity of $-\sum p\ln p$ is what supplies that confirmation.
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

**The linear term vanishes.** Substitute $u = x-\mu$. The exponent becomes $-1-\lambda_0-\lambda_1\mu-\lambda_1u-\lambda_2u^2$; the first three terms carry no $u$, so set them aside into an overall prefactor and complete the square in what remains:

$$-\lambda_1 u - \lambda_2u^2 = -\lambda_2\left(u + \frac{\lambda_1}{2\lambda_2}\right)^2 + \frac{\lambda_1^2}{4\lambda_2}$$

so $f$ is a bell curve centred at $x = \mu - \lambda_1/(2\lambda_2)$ (the leftover $\lambda_1^2/4\lambda_2$ is another $u$-free constant, and joins the prefactor too). But the mean constraint says the centre is at $\mu$. That **forces** $\lambda_1 = 0$ — it is a consequence of the constraints, not an assumption. Absorbing the remaining constants into one prefactor $A$,

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
This is the **same Gaussian** that Weeks 1-3 derive from Gauss's 1809 argument — and the two derivations share no assumptions at all.

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
Weeks 4-6 introduce least squares with the remark that squaring keeps errors positive and makes the minimization linear. Both are true, and both are conveniences rather than reasons.

The real reason is the table above. **Squaring the errors is equivalent to asserting that the errors are Gaussian.** That assertion is usually reasonable — by the Central Limit Theorem, errors that accumulate from many small independent causes tend to be Gaussian — but it is an assumption about your data, not a mathematical necessity, and it is the assumption that fails when you have outliers.

This also explains a fact that would otherwise look like coincidence: the mean is the estimator that minimizes squared error, and the Gaussian is the maximum-entropy distribution for fixed *squared* deviation. The median is the estimator that minimizes absolute error, and the Laplace is the maximum-entropy distribution for fixed *absolute* deviation. Mean-with-Gaussian and median-with-Laplace are two self-consistent packages, and the loss you pick selects which one you are in.
:::

---

## Both questions at once — Bayes' theorem

Look at what each question quietly assumed.

**Question 1 had no prior opinion.** Among all PDFs with the right mean and variance it chose the most noncommittal one, deliberately refusing to favour any particular shape.

**Question 2 had no prior opinion either.** Among all candidate models it chose whichever one the data liked best — however implausible that model looked *before* the data arrived.

Neither is the usual situation. Usually you have **both**: something you believed beforehand, and data that bears on it. Putting the two together takes the same move and exactly one new ingredient.

### The functional: the cost of changing your mind

Entropy measures how noncommittal a distribution is in absolute terms. That was the right thing to maximize when you started from nothing. But now you start from somewhere — a **prior** $\pi_i$, your probability for hypothesis $i$ before the data — and the principle becomes: *move as little as possible, but move enough to accommodate what you saw.*

So the quantity to measure is not how noncommittal a candidate belief $q_i$ is, but how far it has **moved** from $\pi_i$. That is the **relative entropy**, also called the Kullback-Leibler divergence:

$$D[q\|\pi] = \sum_i q_i\ln\frac{q_i}{\pi_i}$$

(Here $\pi_i$ is a prior probability — nothing to do with $3.14159$.) Read it as a **cost of revision**, in the same units as entropy: zero when $q = \pi$, since you changed nothing and so paid nothing; positive otherwise; and very large when $q$ puts weight where $\pi$ had put almost none, which is what a drastic change of mind looks like.

This does not replace Question 1 — it contains it. Set the prior flat, $\pi_i = 1/N$, and

$$-D[q\|\pi] = -\sum_i q_i\ln q_i - \ln N$$

which is the entropy of Question 1 up to an additive constant, the same kind of harmless offset as the $-\ln\Delta x$ in Step 1's loose end, and it moves no maximum. **Maximizing entropy is minimizing the departure from a flat prior.** Question 1 was the special case of believing nothing in particular to begin with.

### The constraint: the data

Let $L_i$ be the **likelihood** — the probability of the data you actually observed, if hypothesis $i$ were the true one. A belief that takes the data seriously must assign it a decent log-likelihood on average, so the constraint is

$$\sum_i q_i\ln L_i = c$$

with $c$ stating how far you insist on deferring to what you saw. And, as always, that the answer is a probability distribution: $\sum_i q_i = 1$.

### Attach one multiplier per constraint and differentiate

Two constraints, two multipliers. We are minimizing rather than maximizing now, so the signs read the other way round; nothing else changes:

$$\mathcal{L}(q_1,\dots,q_N) = \sum_i q_i\ln\frac{q_i}{\pi_i} + \lambda_0\Big(\sum_i q_i - 1\Big) - \beta\Big(\sum_i q_i\ln L_i - c\Big)$$

Differentiate with respect to a single hypothesis $q_j$. Only the $i=j$ term of each sum survives, and $\frac{d}{dq}\big(q\ln(q/\pi)\big) = \ln(q/\pi) + 1$, so

$$\frac{\partial\mathcal{L}}{\partial q_j} = \ln\frac{q_j}{\pi_j} + 1 + \lambda_0 - \beta\ln L_j = 0$$

Rearranging and exponentiating, exactly as in Step 3 above,

$$q_j = \frac{\pi_j\,L_j^{\beta}}{Z}, \qquad Z = \sum_i \pi_i L_i^{\beta}$$

Set $\beta = 1$ and read it out loud — **posterior $\propto$ prior $\times$ likelihood**:

$$\boxed{q_j = \frac{\pi_j L_j}{\sum_i \pi_i L_i}}$$

That is **Bayes' theorem**. It was not assumed anywhere. It is what this chapter's one move returns when the functional is the cost of revision and the constraint is the data.

### What the two multipliers turn out to be

Now do what the interlude said to always do.

**$\lambda_0$, on normalization, is the evidence.** Normalizing $q_j = \pi_j L_j^{\beta}e^{-(1+\lambda_0)}$ gives $e^{1+\lambda_0} = Z$, and at $\beta = 1$,

$$Z = \sum_i \pi_i L_i = p(\text{data})$$

the **marginal likelihood**, or *evidence* — the denominator of Bayes' theorem and the quantity all Bayesian model comparison is built on. The multiplier that did nothing but enforce "probabilities sum to one" turned out to be the probability of your data. (The stray $-1$ is the same one that appeared in Step 3, and comes from the same derivative.)

**$\beta$, on the data constraint, is an exchange rate between the data and your prior.** By the shadow-price property, $dD^*/dc = \beta$: it is how many nats of revision you must pay per extra nat of fit demanded. $\beta = 1$ is ordinary Bayes — each observation counted once, at face value. $\beta < 1$ deliberately discounts the data, which is what you want when the likelihood is misspecified, or when the observations are autocorrelated and therefore are not really $n$ independent facts. So **Bayes' rule is one member of a one-parameter family, and the parameter is a Lagrange multiplier.**

:::{admonition} Example / deeper dive
:class: note
**Bayes with no free parameter — and what a revision costs**

Introducing the likelihood as a *constraint* is what produced the adjustable $\beta$. There is a tighter version in which no such freedom exists.

Work in the joint space of hypothesis *and* data. Your prior there is $\pi(\theta,x) = \pi(\theta)\,p(x\mid\theta)$ — it already contains the likelihood, because to have a model at all is to know what data each hypothesis would produce. Then observing $x = x_0$ is not a soft constraint but a hard one: *all* the probability must sit on $x = x_0$. Minimize $D[q\|\pi]$ subject to that and to normalization, and the single remaining multiplier delivers

$$q(\theta) = \frac{\pi(\theta)\,p(x_0\mid\theta)}{p(x_0)} = \pi(\theta\mid x_0)$$

with no $\beta$ anywhere: the likelihood enters with power exactly one because it was never a constraint to be priced. **Conditioning on data is the projection of your prior onto the set of beliefs compatible with what you saw**, measured in relative entropy. The $\beta$ of the main derivation is the price of having smuggled the likelihood in as a constraint instead.

The value of the minimum is worth computing too. Substituting back,

$$D^* = \sum_\theta \pi(\theta\mid x_0)\ln\frac{\pi(\theta\mid x_0)}{\pi(\theta,x_0)} = -\ln p(x_0)$$

**The minimum cost of revising your beliefs is exactly the surprisal of the data.** Unsurprising data is cheap to absorb; data your prior thought nearly impossible forces an expensive revision. Note that the evidence has now appeared twice over — once as the multiplier, once as the optimal value.
:::

### A penalty is a negative log prior

Question 2 established that a loss function is a negative log-likelihood. Bayes finishes the sentence. Take logs of posterior $\propto$ prior $\times$ likelihood:

$$-\ln p(\theta\mid y) = \underbrace{-\ln L(y\mid\theta)}_{\text{a loss}}\;\underbrace{-\ln\pi(\theta)}_{\text{a penalty}} + \text{const}$$

**A penalty term is a negative log prior**, exactly as a loss term is a negative log-likelihood. And what converts a stated belief into a penalty is a Lagrange multiplier. Suppose you fit coefficients $\mathbf{w}$ by least squares but refuse to let them grow large — because you have more predictors than you trust, or they are collinear. That is a constrained problem of precisely this chapter's shape:

$$\text{minimize }\lVert y - X\mathbf{w}\rVert^2 \quad\text{subject to}\quad \lVert\mathbf{w}\rVert^2 \le t \qquad\Longrightarrow\qquad \mathcal{L} = \lVert y - X\mathbf{w}\rVert^2 + \lambda\lVert\mathbf{w}\rVert^2$$

Now compare that with the negative log posterior for Gaussian errors of variance $\sigma^2$ and a Gaussian prior of variance $\tau^2$ on each coefficient, which is $\lVert y - X\mathbf{w}\rVert^2 + (\sigma^2/\tau^2)\lVert\mathbf{w}\rVert^2$. The two are the same problem, so

$$\boxed{\lambda = \frac{\sigma^2}{\tau^2}}$$

**The multiplier is the ratio of how noisy the data is to how confident the prior is** — an exchange rate between data and belief, in the most literal sense available. $\lambda\to0$ is a prior so vague it may as well be absent, and ordinary least squares returns; $\lambda\to\infty$ is a prior so rigid that no data can move it.

And the Gaussian/Laplace pair comes back a third time:

| Prior on the coefficients | Penalty it becomes | Effect |
|---|---|---|
| Gaussian | $\lambda\lVert\mathbf{w}\rVert^2$ (L2) | shrinks every coefficient smoothly toward zero |
| Laplace | $\lambda\lVert\mathbf{w}\rVert_1$ (L1) | sets most coefficients exactly to zero |

Squared-versus-absolute chose the error distribution in Question 2; the identical choice chooses the prior here. Same two constraints from the maximum-entropy table, same two densities, same two norms — applied now to the parameters instead of the residuals.

### Where the prior comes from

One honest gap. Minimum relative entropy tells you how to **revise** a belief, not how to **start** one: it needs $\pi$ handed to it from outside.

The usual answer is Question 1. Choose the prior by maximum entropy from whatever you genuinely know, then let the data revise it by minimum relative entropy — the two halves of this chapter, in that order. What the pattern will not do is excuse you from saying what you believed beforehand. It only forces that statement into the open, where it can be argued with, which is more than a method that leaves it implicit can offer.

---

## The same move, everywhere

Both questions — and their combination — turned out to be: *optimize a functional subject to constraints.* So does everything that follows.

| Part | Functional optimized | Subject to | Solution |
|---|---|---|---|
| **0** Bayes | revision cost $\sum q\ln(q/\pi)$ | the observed likelihood | the posterior |
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

Four sentences:

1. **The null hypothesis is a PDF, and you can derive rather than guess it** — maximize entropy subject to what you actually know.
2. **A loss function is a negative log-likelihood, so choosing a loss is choosing an error distribution** — L2 asserts Gaussian errors, L1 asserts Laplace.
3. **Bayes' theorem is this same move with a prior attached** — minimize the cost of revising what you believed, subject to honouring what you saw; a penalty term is then a negative log prior, just as a loss is a negative log-likelihood.
4. **Whenever you want "the best" of something subject to a restriction, attach a Lagrange multiplier** — and check what the multiplier turns out to mean, because it is often the quantity you actually wanted.

Everything else this semester is these four sentences applied to bigger objects.

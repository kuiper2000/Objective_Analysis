(Seeking_Structure)=
# Week 9-12: Seeking Structure in Data

In Weeks 5-8 we asked a two-variable question: given $x(t)$ and $y(t)$, how strongly are they related? That question has an obvious limitation. A climate field is not two time series — it is thousands of grid points, each with its own time series, all correlated with one another in ways we would like to summarize.

This block is about finding **structure** in that kind of data: coherent spatial patterns that vary together in time, and groupings of observations that resemble one another. The two main tools are **Empirical Orthogonal Function (EOF) analysis** and **cluster analysis**.

The good news is that nothing here is conceptually new. EOF analysis is regression — the same projection of one vector onto another that gave us the correlation coefficient — applied simultaneously to every grid point, and organized so that the resulting patterns are uncorrelated with each other. What we need first is the language to write "every grid point at once" compactly, which is linear algebra.

---

## Part 1 — Linear Algebra Review

### Matrix notation

A matrix $\mathbf{A}$ with $M$ rows and $N$ columns (dimensions $M \times N$) is ordered as

$$\mathbf{A} = \begin{bmatrix}
a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\
a_{21} & a_{22} & a_{23} & \cdots & a_{2n} \\
a_{31} & a_{32} & a_{33} & \cdots & a_{3n} \\
\vdots & \vdots & \vdots & & \vdots \\
a_{m1} & a_{m2} & a_{m3} & \cdots & a_{mn}
\end{bmatrix}$$

The notation $a_{ij}$ denotes the element in the $i$th row and $j$th column, with $1 \le i \le M$ and $1 \le j \le N$.

:::{admonition} A convention you must fix before you start
:class: warning
Almost every error in EOF analysis traces back to losing track of which dimension is **space** and which is **time**. Throughout these notes we will use

$$\mathbf{X} = [M \times N] = [\text{sampling} \times \text{space}]$$

so $M$ is the number of time steps and $N$ is the number of grid points — *except* in the SVD section, where we follow Prof. Hartmann's convention of space first. Every recipe below states its convention explicitly in the first line. Read it before you code.
:::

### Transpose

The transpose of $\mathbf{X}$ is written $\mathbf{X}^T$ and is obtained by flipping rows and columns, so an $M \times N$ matrix becomes $N \times M$. In index notation, if $\mathbf{X} = \mathbf{Y}^T$ then

$$X_{ij} = Y_{ji}$$

### Inner product — where regression re-enters

All vectors are matrices whose second dimension is 1. If $\mathbf{X}$ and $\mathbf{Y}$ are two vectors of length $M \times 1$, their inner product is

$$\mathbf{X}^T\mathbf{Y} = \sum_{i=1}^{M} x_i y_i = \text{scalar}, \qquad [1 \times M]\cdot[M \times 1] = [1\times1]$$

This is worth pausing on, because it is the bridge from Weeks 5-8:

- if $\mathbf{X}$ and $\mathbf{Y}$ have **mean zero**, then $\frac{1}{M}\mathbf{X}^T\mathbf{Y}$ is the **covariance** of the two vectors
- if $\mathbf{X}$ and $\mathbf{Y}$ are additionally **standardized** (unit standard deviation), then $\frac{1}{M}\mathbf{X}^T\mathbf{Y}$ is the **correlation coefficient** — and also the regression coefficient, since for standardized data $a_1 = r$

So an inner product *is* a regression. Everything in the EOF sections is built from this one operation.

### Multiplying matrices

If $\mathbf{A}$ is $M \times N$ and $\mathbf{B}$ is $N \times P$, then

$$\mathbf{A}\mathbf{B} = \mathbf{C}, \qquad [M\times N]\cdot[N\times P] = [M \times P]$$

The number of **columns** of the first matrix must equal the number of **rows** of the second. Unlike ordinary multiplication, matrix multiplication does not commute — you cannot re-order it.

To multiply, take the row of the first matrix against the column of the second, multiply element-wise and sum:

$$\begin{bmatrix} 1 & 2 & 3 \\ -1 & -1 & -3 \end{bmatrix}
\times
\begin{bmatrix} 2 & 3 \\ 4 & 2 \\ -1 & 1 \end{bmatrix}
=
\begin{bmatrix} 7 & 10 \\ -3 & -8 \end{bmatrix}$$

$$[2\times3]\times[3\times2] = [2\times2]$$

### Covariance and correlation matrices

This is the object the entire EOF machinery acts on. Given $\mathbf{X}$ of dimension $M \times N$ whose **columns have mean zero**,

$$\boxed{\mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X}}$$

is the **covariance matrix**, of dimension $N \times N$. Element $C_{ij}$ is the covariance between grid point $i$ and grid point $j$ — i.e. it is the inner product of column $i$ with column $j$, which by the previous section is exactly a covariance.

If the columns of $\mathbf{X}$ are **standardized** rather than merely centered, the same expression

$$\mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X}$$

is the **correlation matrix**. The formula does not change; what changes is the preprocessing. This distinction returns in Part 3, because standardizing your data means you are eigenanalyzing the correlation matrix whether you intended to or not.

$\mathbf{C}$ is also sometimes called the **dispersion matrix**. It is always square and always symmetric ($\mathbf{C} = \mathbf{C}^T$) — a fact we will lean on heavily.

### Inverse

The inverse of a square matrix $\mathbf{X}$ is written $\mathbf{X}^{-1}$. This does **not** mean raising the matrix to the $-1$ power; it denotes the inverse specifically, defined by

$$\mathbf{X}\mathbf{X}^{-1} = \mathbf{X}^{-1}\mathbf{X} = \mathbf{I}$$

where $\mathbf{I}$ is the identity matrix — square, with ones on the diagonal and zeros elsewhere.

**Only square matrices have inverses.** (There are "left" and "right" inverses for non-square matrices, but we will not use them here.)

In ordinary algebra you divide to remove a multiplied variable. With matrices, you multiply by the inverse instead — that is what the inverse is *for*.

### Types of matrices

- $\mathbf{X}$ is **diagonal** if it is square with values along the diagonal and zeros everywhere else
- $\mathbf{X}$ is **symmetric** if $\mathbf{X}^T = \mathbf{X}$
- $\mathbf{X}$ is **orthonormal** if $\mathbf{X}^T = \mathbf{X}^{-1}$, equivalently $\mathbf{X}^T\mathbf{X} = \mathbf{X}\mathbf{X}^T = \mathbf{I}$

That last one is the useful one. If all off-diagonal elements of $\mathbf{X}^T\mathbf{X}$ are zero, then $\mathbf{X}^T = \mathbf{X}^{-1}$, and all columns of $\mathbf{X}$ are linearly independent — i.e. uncorrelated with one another. When we find that our EOFs form an orthonormal matrix, transposing will be all the "inverting" we ever need to do.

### Vector spaces

- $\mathbb{R}^0$: a single point (e.g. 0)
- $\mathbb{R}^1$: the real line
- $\mathbb{R}^2$: where all 2-element vectors live
- $\mathbb{R}^3$: 3-D space

Vector spaces are **spanned** by a set of **basis vectors**. A basis for $\mathbb{R}^N$ is a set of $N$ vectors that, through linear combination, can describe any vector in $\mathbb{R}^N$.

There are infinitely many possible bases for any $\mathbb{R}^N$, but the conventional choice is the unit vectors. For $\mathbb{R}^2$:

$$\mathbf{e}_1 = (1, 0), \qquad \mathbf{e}_2 = (0, 1)$$

Any point in $\mathbb{R}^2$ is a linear combination of these — for example $(3,2) = 3\mathbf{e}_1 + 2\mathbf{e}_2$. Try it; you can reach anywhere in the plane.

Now a question. Do these form a basis for $\mathbb{R}^2$?

$$\mathbf{e}_1 = (1,0), \qquad \mathbf{e}_2 = (2,0)$$

No — any point with a non-zero second component is unreachable. The direction $(0,1)$ and all its multiples form the **null space**: the part of the space these vectors cannot reach.

For $\mathbb{R}^3$ the unit basis is

$$\mathbf{e}_1 = (1,0,0), \quad \mathbf{e}_2 = (0,1,0), \quad \mathbf{e}_3 = (0,0,1)$$

By definition basis vectors must be **linearly independent** — each supplies information the others do not.

### Testing for a basis: determinants

To check whether a set of vectors forms a basis, assemble them as the columns of $\mathbf{A}$ and ask whether there is a non-trivial $\mathbf{x}$ with

$$\mathbf{A}\mathbf{x} = 0$$

If such an $\mathbf{x}$ exists, the vectors fail to span the space and/or are not all independent — there is a non-trivial null space.

Rather than solving for every $\mathbf{x}$, we only need to know whether a solution exists, and there is a shortcut. If

$$\det(\mathbf{A}) = |\mathbf{A}| = 0$$

then a non-trivial solution to $\mathbf{A}\mathbf{x} = 0$ exists. Determinants are painful for large matrices but trivial for $2\times2$:

$$\det(\mathbf{A}) = a_{11}a_{22} - a_{12}a_{21}$$

If $\mathbf{A}$ is square and $\mathbf{A}\mathbf{x} = 0$ has a non-trivial solution, then $\mathbf{A}$ is called **singular** or **degenerate**, and it has no inverse.

### Rank

The **rank** $r$ of a matrix is the number of linearly independent columns. If $\mathbf{A}$ holds the basis vectors of $\mathbb{R}^2$,

$$\mathbf{A} = \begin{bmatrix} 1 & 0 \\ 0 & 1\end{bmatrix}, \qquad \det\mathbf{A} = 1 \ne 0$$

so the rank is 2. But if

$$\mathbf{A} = \begin{bmatrix} 1 & 2 \\ 0 & 0\end{bmatrix}, \qquad \det\mathbf{A} = 0 - 0 = 0$$

the columns are not independent. Here you can see it directly — twice the first column gives the second — but for larger matrices the determinant is what you rely on. The rank of this $2\times2$ matrix is 1.

Rank is bounded by the smaller dimension:

$$r \le \min(M, N)$$

If $\mathbf{A}$ is $N \times N$ and $r < N$, there is a non-trivial null space; the directions in $\mathbb{R}^N$ not spanned by the columns of $\mathbf{A}$ *define* that null space.

- columns linearly independent $\Rightarrow$ $\mathbf{A}$ is **full rank**
- columns not linearly independent $\Rightarrow$ $\mathbf{A}$ is **rank deficient**

This matters practically: a rank-deficient covariance matrix has fewer non-zero eigenvalues than its dimension, which means fewer meaningful EOFs than grid points. If you have 100 time steps and 40000 grid points, your covariance matrix is $40000\times40000$ but its rank is at most 100 — so at most 100 EOFs carry any information at all.

---

## Part 2 — Eigenvalues and Eigenvectors

Recall the covariance (dispersion) matrix

$$\mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X} \qquad \text{or} \qquad \mathbf{C} = \frac{1}{N}\mathbf{X}\mathbf{X}^T$$

depending on which dimension is the sampling dimension. Either way $\mathbf{C}$ is square and symmetric.

Principal Component Analysis is an **eigenanalysis of $\mathbf{C}$**. Any symmetric matrix $\mathbf{C}$ can be decomposed as

$$\mathbf{C}\mathbf{e}_i = \lambda_i\mathbf{e}_i$$

or in full matrix form

$$\boxed{\mathbf{C}\mathbf{E} = \mathbf{E}\boldsymbol{\Lambda}}$$

where $\mathbf{E}$ has the eigenvectors $\mathbf{e}_i$ as its columns, and $\boldsymbol{\Lambda}$ is diagonal with the eigenvalues $\lambda_i$ along it.

The set of eigenvectors and eigenvalues represents a **coordinate transformation into a space where $\mathbf{C}$ becomes diagonal**. This is the whole idea, and it is worth stating carefully:

> Because the covariance matrix is diagonal in this new coordinate space, the variations in these new directions are uncorrelated with each other. The eigenvectors define directions in the initial coordinate space along which the maximum possible variance can be explained, and in which variance in one direction is orthogonal to the variance explained by other directions defined by the other eigenvectors. The eigenvalues indicate how much variance is explained by each eigenvector.
>
> — D. Hartmann notes, Chapter 4

:::{admonition} Example / deeper dive
:class: note
**A worked $2\times2$ eigenanalysis**

Let the data be two points with the means already removed:

$$\mathbf{A} = \begin{bmatrix} -1 & 1 \\ 1 & -1 \end{bmatrix}$$

The covariance matrix is

$$\mathbf{C} = \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

and one can show that

$$\mathbf{e}_1 = \begin{bmatrix} -\sqrt{2}/2 \\ \sqrt{2}/2 \end{bmatrix} = \frac{\sqrt2}{2}\begin{bmatrix} -1 \\ 1\end{bmatrix}$$

Note that the norm of $\mathbf{e}_1$ is 1. Check it against the eigenvalue equation:

$$\begin{bmatrix} 1 & -1 \\ -1 & 1\end{bmatrix}\begin{bmatrix} -\sqrt2/2 \\ \sqrt2/2\end{bmatrix} = 2\cdot\begin{bmatrix} -\sqrt2/2 \\ \sqrt2/2\end{bmatrix}$$

so $\lambda_1 = 2$. It turns out that $\sum_i \lambda_i$ equals the total variance of $\mathbf{A}$ — the sum of the diagonal of $\mathbf{C}$, which here is 2. So in this problem there is only one non-zero eigenvalue ($\lambda_2 = 0$), and the first eigenvector explains **all** the variance.

That makes sense geometrically: both data points lie on the line $y = -x$, so a single direction describes them completely.

Finally, note that

$$\mathbf{e}_1 = \begin{bmatrix} \sqrt2/2 \\ -\sqrt2/2\end{bmatrix}$$

is an equally valid eigenvector. **The overall sign of an eigenvector is arbitrary** — which is why you should never over-interpret the sign of an EOF map without also looking at its PC.
:::

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure Example: `eigenanalysis_2d.py`* — dots denote the points of $\mathbf{A}$ with means removed, $(-1,1)$ and $(1,-1)$. Black coordinates denote the original $x$–$y$ axes; purple coordinates denote positions along the new axis $\mathbf{e}_1$, i.e. $(-\sqrt2, 0)$ and $(\sqrt2, 0)$. The two points and the dotted purple line (the direction of $\mathbf{e}_1$) lie along $y = -x$.

A second panel adds a third point $(-1.1, -0.2)$ so the points no longer fall exactly on a line. Eigenanalysis now returns two eigenvectors $\mathbf{e}_1$ and $\mathbf{e}_2$, drawn as perpendicular dotted lines through the origin. Eigenvector 1 explains most of the variance; eigenvector 2 accounts for the deviations from it.
:::

### Solving for eigenvectors

$$\begin{aligned}
\mathbf{C}\mathbf{e}_i &= \lambda_i\mathbf{e}_i \\
\mathbf{C}\mathbf{e}_i - \lambda_i\mathbf{e}_i &= 0 \\
(\mathbf{C} - \lambda_i\mathbf{I})\mathbf{e}_i &= 0
\end{aligned}$$

So the eigenvector $\mathbf{e}_i$ lies in the **null space** of $(\mathbf{C} - \lambda_i\mathbf{I})$ — which is exactly the null-space machinery from Part 1. That matrix only has a null space if it is singular, i.e. if

$$\det(\mathbf{C} - \lambda_i\mathbf{I}) = 0$$

Solving this is the **eigenvalue problem**, and it appears throughout science and engineering.

For a real, symmetric $\mathbf{C}$ the eigenvalues $\lambda_i$ are positive. If $\mathbf{C}$ is $M\times M$ there will be $M$ eigenvalues and $M$ eigenvectors — unless $\mathbf{C}$ is degenerate, in which case there are only $r$ of them, where $r$ is the rank.

### Properties of eigenvectors

- the eigenvectors of a matrix form a linearly independent set
- if $\mathbf{C}$ ($M\times M$) has full rank, the eigenvectors form a basis for $\mathbb{R}^M$
- if $\mathbf{C}$ is singular with rank $r = k$, the first $k$ eigenvectors span the subspace of $\mathbf{C}$ ($\mathbb{R}^k$), and eigenvectors $k+1 : M$ span the null space of $\mathbf{C}$
- if $\mathbf{C}$ is real and symmetric, the eigenvectors are not merely linearly independent, they are **orthogonal**

That last property is what puts the "O" in EOF, and we prove it in Part 3.

### Terminology, precisely

- **linearly independent**: you cannot reach one vector by a linear combination of the others
- **orthogonal**: the angle between the vectors is $90^\circ$; they share no vector space
- **orthonormal**: orthogonal *and* of unit length ($\sqrt{x^2+y^2} = 1$)

---

## Part 3 — EOFs via Eigenanalysis

:::{admonition} Figure / in-class demonstration
:class: tip
*Show intro EOF slides:* `eofintro.pdf`
:::

A vocabulary note first, since the literature is loose about it:

- **EOF** stands for Empirical Orthogonal Functions
- **PCA** stands for Principal Component Analysis
- the two are practically the same thing
- **PCA** is the more general term, and just implies splitting a matrix into linearly uncorrelated variables called principal components
- **EOF analysis** normally refers to finding the principal components that maximize variance explained — so, a subset of all the PCs

We will cover two routes to EOFs: eigenanalysis of the covariance matrix (this Part) and singular value decomposition (Part 4). **Both give identical answers.**

### Setting up the maximization

Let $\mathbf{X}$ be $M \times N$, rows being times and columns being spatial locations. We want to decompose $\mathbf{X}$ into orthogonal eigenvectors that explain the most variance. Equivalently, we want to find the $\mathbf{e}_1$ (dimension $[N\times1]$, i.e. a spatial pattern) that most resembles the data.

Building the measure of "resemblance", one requirement at a time:

1. the projection of $\mathbf{e}_1$ onto the data is measured by the inner product of $\mathbf{e}_1$ with $\mathbf{X}$ — this is the resemblance
2. the projection is **squared**, so that positive and negative resemblance count the same (exactly as in linear regression, where we squared the errors)
3. to make the measure independent of the number of observations, divide by $M$
4. to make it independent of the *magnitude* of $\mathbf{e}_1$, normalize to unit length, $\mathbf{e}_1^T\mathbf{e}_1 = 1$ — otherwise the projection could be made arbitrarily large just by lengthening $\mathbf{e}_1$
5. following these rules, the measure of resemblance is

$$(\mathbf{X}\mathbf{e}_1)^2\frac{1}{M} = \mathbf{e}_1^T\mathbf{X}^T\mathbf{X}\mathbf{e}_1\frac{1}{M}$$

which is equivalent to maximizing

$$\mathbf{e}_1^T\mathbf{C}\mathbf{e}_1 \qquad \text{subject to } \mathbf{e}_1^T\mathbf{e}_1 = 1, \qquad \mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X}$$

Call the maximum value of this squared projection $\lambda_1$ — the variance explained by $\mathbf{e}_1$:

$$\mathbf{e}_1^T\mathbf{C}\mathbf{e}_1 = \lambda_1$$

This has exactly the form of the eigenvalue problem from Part 2. Since the norm of $\mathbf{e}_1$ is one, it becomes

$$\mathbf{C}\mathbf{e}_1 = \mathbf{e}_1\lambda_1$$

and in full matrix notation (with $\mathbf{C}$, $\boldsymbol{\Lambda}$, $\mathbf{E}$ all $N \times N$),

$$\mathbf{E}^T\mathbf{C}\mathbf{E} = \boldsymbol{\Lambda} \;\Rightarrow\; \mathbf{E}\mathbf{E}^T\mathbf{C}\mathbf{E} = \mathbf{E}\boldsymbol{\Lambda} \;\Rightarrow\; \mathbf{C}\mathbf{E} = \mathbf{E}\boldsymbol{\Lambda}$$

using $\mathbf{E}^T\mathbf{E} = \mathbf{E}\mathbf{E}^T = \mathbf{I}$.

So $\mathbf{e}_1$ **must be an eigenvector of $\mathbf{C}$**. We find it by eigenanalyzing $\mathbf{C}$ and taking the eigenvector belonging to the largest eigenvalue. The search for "the pattern that explains the most variance" and the eigenvalue problem are the same problem.

### What this buys you

- the 1st eigenvector explains the most variance in $\mathbf{X}$ and has the largest eigenvalue
- the 2nd eigenvector explains the second most, with the 2nd largest eigenvalue
- the eigenanalysis transforms $\mathbf{C}$ into a coordinate system where the new dispersion matrix is diagonal ($\boldsymbol{\Lambda}$)
- in that space all the variance is on the diagonal, since the vectors are orthogonal, so the **fraction of variance explained by the $j$th eigenvector** is

$$\boxed{\text{fraction of variance} = \frac{\lambda_j}{\sum_i \lambda_i}}$$

- EOFs are conventionally ordered by $\lambda$, so EOF 1 explains the largest variance and the last EOF the smallest
- **no other linear combination of $k$ predictors can explain a larger fraction of the variance than the first $k$ EOF/PC pairs**

That final point is the reason EOFs are used for compression and noise reduction: they are optimal in the variance sense, by construction.

:::{admonition} Example / deeper dive
:class: note
**Proof that eigenvectors are orthogonal**

Suppose we have two eigenvectors $\mathbf{e}_j$ and $\mathbf{e}_k$:

$$\mathbf{C}\mathbf{e}_j = \mathbf{e}_j\lambda_j, \qquad \mathbf{C}\mathbf{e}_k = \mathbf{e}_k\lambda_k$$

Multiply the $j$th equation by $\mathbf{e}_k^T$ and take the transpose; multiply the $k$th by $\mathbf{e}_j^T$; subtract. This gives

$$\mathbf{e}_j^T\mathbf{C}^T\mathbf{e}_k - \mathbf{e}_j^T\mathbf{C}\mathbf{e}_k = (\lambda_j - \lambda_k)\mathbf{e}_j^T\mathbf{e}_k$$

Since $\mathbf{C}$ is symmetric, $\mathbf{C} = \mathbf{C}^T$, so the left-hand side vanishes. Therefore

$$\mathbf{e}_j^T\mathbf{e}_k = 0 \qquad \text{unless} \qquad \lambda_j = \lambda_k$$

**The eigenvectors are orthogonal provided the eigenvalues are distinct.** Hold onto that caveat — it comes back in Part 5, where eigenvalues that are *not* well separated turn out to be the central practical problem in EOF analysis.
:::

### Getting your original data back

We want to keep all the original information in the new coordinate system, so we need the amplitudes that express each state in the EOF basis. Define

$$\boxed{\mathbf{Z} = \mathbf{X}\mathbf{E}}$$

Then, since $\mathbf{E}\mathbf{E}^T = \mathbf{I}$,

$$\mathbf{X} = \mathbf{Z}\mathbf{E}^T$$

The columns of $\mathbf{Z}$ are the **principal components (PCs)**, of length $M$. They tell you how much a given sample resembles a particular EOF structure.

To reconstruct a single observation at time $t$ and location $i$, combine the eigenvectors linearly with the PC amplitudes:

$$x_i(t) = \sum_{j=1}^{N} e_{ij}z_j(t)$$

Concretely, suppose the first eigenvector is

$$\mathbf{e}_1 = \begin{bmatrix} e_{11} \\ e_{21} \\ e_{31} \\ e_{41} \\ \vdots\end{bmatrix}$$

Project it onto the original data to get its amplitude at each time step:

$$\begin{bmatrix}
x_{11} & x_{12} & \cdots & x_{1n} \\
x_{21} & x_{22} & \cdots & x_{2n} \\
\vdots & \vdots & & \vdots \\
x_{m1} & x_{m2} & \cdots & x_{mn}
\end{bmatrix}
\begin{bmatrix} e_{11} \\ e_{21} \\ e_{31} \\ \vdots \end{bmatrix}
=
\begin{bmatrix} z_{11} \\ z_{21} \\ z_{31} \\ \vdots \\ z_{m1}\end{bmatrix}$$

that is, $\mathbf{X}\mathbf{E} = \mathbf{Z}$.

:::{admonition} A common mistake
:class: warning
Project $\mathbf{e}$ onto the **original data $\mathbf{X}$**, not onto the covariance matrix $\mathbf{C}$. The covariance matrix has already had the time dimension summed out of it — there is no time series left in it to recover.
:::

### Practical considerations: which covariance matrix?

The covariance matrix can be formed two ways:

$$\mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X} \to [N\times N] \qquad \text{or} \qquad \mathbf{C} = \frac{1}{N}\mathbf{X}\mathbf{X}^T \to [M\times M]$$

So far we have used the first. Which you *should* use is purely a question of computational cost:

- $M = 20000$ days, $N = 400$ grid points $\Rightarrow$ use $\frac{1}{M}\mathbf{X}^T\mathbf{X}$, a $400\times400$ problem. Easy.
- $M = 100$ days, $N = 40000$ grid points $\Rightarrow$ $\frac{1}{M}\mathbf{X}^T\mathbf{X}$ would be $40000\times40000$. Not easy.

In the second case you have two options:

1. interpolate your data onto a coarser mesh
2. perform the eigenanalysis on the **spatial dispersion matrix**, $\mathbf{C} = \frac{1}{N}\mathbf{X}\mathbf{X}^T$, which is only $100\times100$

Option 2 gives the same answer, with one twist you must remember: **the eigenvectors are now the PCs, and the PCs are now the EOFs.** The roles swap. Recall from Part 1 that the rank of that giant $40000\times40000$ matrix was at most 100 anyway — so nothing is lost by working in the smaller space.

### Summary: EOFs by eigenanalysis of the covariance matrix, $\mathbf{X} = [M\times N]$

1. subtract the mean along the sampling dimension of $\mathbf{X}$
2. unless data size forces otherwise, define $\mathbf{C}$ along the sampling dimension so you are left with the *space* dimension:

$$\mathbf{C} = \frac{1}{M}\mathbf{X}^T\mathbf{X} \to [N\times N] \quad \text{if } M \text{ is the sampling dimension}$$
$$\mathbf{C} = \frac{1}{N}\mathbf{X}\mathbf{X}^T \to [M\times M] \quad \text{if } N \text{ is the sampling dimension}$$

3. eigenanalyze $\mathbf{C}$ by diagonalizing it (`eig` in Matlab, `EIGENQL` in IDL, `numpy.linalg.eigh` in Python)
4. the first EOF is the eigenvector corresponding to the largest eigenvalue
5. to find the PCs, project the data onto the eigenvectors — PCs should have the length of your sampling dimension
6. the fraction of variance explained by the $i$th EOF pair is $\lambda_i/\sum_j\lambda_j$

:::{admonition} Figure / in-class demonstration
:class: tip
*Show slides:* `eof_physicalvars.pdf`

*Show webpage:* Eigenstyle — <http://blog.thehackerati.com/post/126701202241/eigenstyle> (PCA applied to a wardrobe of dresses; a memorable demonstration that "EOF 1" is just "the direction of greatest variance", whatever the data happen to be).
:::

---

## Part 4 — EOFs via Singular Value Decomposition (SVD)

There is a second route to the same answer, and in practice it is the one you will use. SVD is more general, more numerically stable, and skips the explicit formation of $\mathbf{C}$ entirely.

Any $M\times N$ matrix can be decomposed into a product of three matrices:

$$\boxed{\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T}$$

where

- $\mathbf{U}$ is $M\times M$
- $\boldsymbol{\Sigma}$ is $M\times N$
- $\mathbf{V}$ is $N\times N$
- $\mathbf{U}$ and $\mathbf{V}$ are orthogonal
- $\boldsymbol{\Sigma}$ is diagonal

If we let the $M$ dimension be **space** and the $N$ dimension be **time** (following Prof. Hartmann's notes — note this is the opposite of Part 3), then:

- the columns of $\mathbf{U}$ are the eigenvectors of $\mathbf{X}\mathbf{X}^T$ — the **EOFs**
- the columns of $\mathbf{V}$ (rows of $\mathbf{V}^T$) are the eigenvectors of $\mathbf{X}^T\mathbf{X}$ — the **PCs**
- the $r$ singular values along the diagonal of $\boldsymbol{\Sigma}$ are proportional to the square roots of the non-zero eigenvalues of both $\mathbf{X}\mathbf{X}^T$ and $\mathbf{X}^T\mathbf{X}$: $s_i = \sqrt{\lambda_i}$
- to get variance explained, square the elements of $\boldsymbol{\Sigma}$

:::{admonition} Watch your dimension order
:class: warning
If instead you have **time along $M$** and **space along $N$** (the Part 3 convention), then the columns of $\mathbf{V}$ are the EOFs and the columns of $\mathbf{U}$ are the PCs. Everything else is identical. The mathematics does not care; you do.
:::

The matrix $\boldsymbol{\Sigma}$ is odd in that it is not square, yet we care about its diagonal. It can be viewed two ways:

1. if $M > N$,

$$\boldsymbol{\Sigma} = \begin{bmatrix} \mathbf{D} \\ \mathbf{0}\end{bmatrix}$$

where $\mathbf{D}$ is $N\times N$ diagonal and $\mathbf{0}$ is $(M-N)\times N$ zeros. The first $r$ elements of $\mathbf{D}$ are the singular values.

2. if $\mathbf{X}$ is full rank and $M = N$, then $\boldsymbol{\Sigma}$ is genuinely square and diagonal.

If the singular values are ordered with $\sigma_1$ largest, then the first column of $\mathbf{U}$ is the leading EOF and the first column of $\mathbf{V}$ is the leading PC.

:::{admonition} Example / deeper dive
:class: note
**Why SVD and eigenanalysis give the same EOFs**

Start from the SVD of $\mathbf{X}$ and square both sides:

$$\begin{aligned}
\mathbf{X}\mathbf{X}^T &= (\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T)(\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T)^T \\
&= (\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T)(\mathbf{V}\boldsymbol{\Sigma}^T\mathbf{U}^T)
\end{aligned}$$

Since $\mathbf{V}^T\mathbf{V} = \mathbf{I}$ (the columns of $\mathbf{U}$ and $\mathbf{V}$ are made orthonormal),

$$\mathbf{X}\mathbf{X}^T = \mathbf{U}\boldsymbol{\Sigma}^2\mathbf{U}^T$$

Multiplying both sides by $\mathbf{U}$,

$$N\cdot\mathbf{C}\mathbf{U} = \mathbf{U}\boldsymbol{\Sigma}^2$$

Compare with the eigenanalysis of the covariance matrix, $\mathbf{C}\mathbf{U} = \mathbf{U}\boldsymbol{\Lambda}$. Therefore

$$\frac{1}{N}\boldsymbol{\Sigma}^2 = \boldsymbol{\Lambda}$$

and $\mathbf{U}$ holds the eigenvectors of $\mathbf{C}$. **The two methods are the same calculation.**

Similarly for the PCs: starting again from $\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T$ and writing $\mathbf{X} = \mathbf{U}\mathbf{Z}$, we identify

$$\mathbf{Z} = \boldsymbol{\Sigma}\mathbf{V}^T$$

which is consistent with $\mathbf{X} = \mathbf{E}\mathbf{Z}$ from the eigenanalysis, since $\mathbf{U} = \mathbf{E}$.
:::

### Summary: EOFs with SVD, $\mathbf{X} = [M\times N]$

1. subtract the mean along the sampling dimension of $\mathbf{X}$
2. decompose $\mathbf{X}$ with a built-in SVD routine
3. if $M$ is the sample dimension (time), the EOFs (spatial patterns) are the columns of $\mathbf{V}$ and the PCs (time series) the columns of $\mathbf{U}$
4. if $N$ is the sample dimension (time), the EOFs are $\mathbf{U}$ and the PCs are $\mathbf{V}$
5. the fraction of variance explained by the $i$th EOF pair is $\sigma_i^2/\sum_j\sigma_j^2$

---

## Part 5 — EOFs with Actual Data

The mathematics above is exact. Applying it to real geophysical data introduces a set of choices that the mathematics will not make for you, and which change your answer.

### Preparing the data

**Removing the sampling mean**

- it is standard to remove the mean along the sampling dimension
- this makes your PCs have zero mean
- you may also want to remove the seasonal cycle (or diurnal cycle) if it is not of interest — otherwise EOF 1 will simply *be* the seasonal cycle, which you already knew about

**Removing the spatial mean**

- less common, since you usually want to allow EOFs to have the same sign everywhere
- exceptions: you care about temperature *gradients* rather than the full field
- or you care about the eddy component, so you remove the zonal mean

**Standardizing the data**

Sometimes you want to remove amplitude information before the analysis:

- your data combine parameters with different units, and you don't want the one with the largest units to dominate the variance
- the variance varies greatly from place to place and that variation is not of interest

To do this, after subtracting the sampling-dimension mean, divide by the standard deviation — your data are now in units of $\sigma$. Note the consequence flagged back in Part 1: **you are now eigenanalyzing the correlation matrix, not the covariance matrix.**

**Weighting**

Gridded data must be weighted by grid area, or high-latitude grid points — which represent much less area — will be over-represented:

- **temporal covariance matrix**: weight by $\sqrt{\cos\theta}$ (when grid points are multiplied against each other in forming $\mathbf{C}$, you recover the full $\cos\theta$)
- **spatial dispersion matrix**: weight by $\cos\theta$, but carefully — $\mathbf{C} = \mathbf{X}_{\cos}^T\mathbf{X}_{\cos}$ where the elements of $\mathbf{X}_{\cos}$ are $x_{t,\theta}\times\cos(\theta)$
- **SVD**: weight your data by $\cos\theta$

The resulting EOF will look unphysical and carry unphysical values. That is expected, and is fixed at the presentation stage below.

**Missing data**

- most routines will not perform eigenanalysis or SVD with NaNs present
- the best option is to use only locations with at least $\alpha\%$ of data present, and to compute EOFs via the covariance matrix so that $\mathbf{C}$ contains no missing values when passed to the eigenanalysis routine

### Presentation of EOFs

Let $\mathbf{X}$ be $[M\times N]$ with $M$ samples and $N$ locations. From here on we write equations for EOF 1 and PC 1 only; repeat for the others.

Both EOFs and PCs get plotted in the literature depending on interest. Generally the spatial domain is the more interesting one, since the PC time series is usually noisy.

**Problem:** the EOFs have no physical units, so we need a meaningful way to plot the patterns.

The procedure:

1. calculate the leading PC $\mathbf{z}_1$ (directly from SVD, or by projection in the eigenanalysis method)
2. standardize $\mathbf{z}_1$ — the mean is already zero, so just divide by $\sigma_{z_1}$
3. project the standardized PC back onto your **original anomaly data** and divide by the length of $\mathbf{z}$, giving $\mathbf{d}_1$: EOF 1 in the units of the original, *unweighted* data

If the original data is $\mathbf{X}$ and the weighted data is $\mathbf{X}_w$, then following eigenanalysis:

$$\mathbf{z}_1 = \mathbf{X}_w\mathbf{e}_1$$

$$\tilde{\mathbf{z}}_1 = \frac{\mathbf{z}_1 - \mu_{z_1}}{\sigma_{z_1}}$$

$$\boxed{\mathbf{d}_1 = \frac{1}{M}\tilde{\mathbf{z}}_1^T\mathbf{X} = \frac{1}{\text{length}(\mathbf{z}_1)}\tilde{\mathbf{z}}_1^T\mathbf{X}}$$

Now $\mathbf{d}_1$ is in physical units, and means: **the anomaly associated with a one-standard-deviation variation of $\mathbf{z}_1$.** That is a statement you can put in a figure caption and defend.

Note this works for both the covariance method and SVD. For SVD there is one fewer step, since the raw PC 1 is handed to you.

Note also that this last step is a **regression** — projecting a standardized time series onto data, exactly as in Week 5-8. The EOF gave us the pattern; regression gives it units.

### How many EOFs should be retained?

This is the question that separates careful EOF analysis from decorative EOF analysis.

**EOFs will always give you an answer.** EOFs of red noise will produce sines and cosines that look thoroughly physical. The rate at which eigenvalues drop off is a function of the redness of the data, not of whether anything meaningful is present.

Worse, **orthogonality itself imposes a constraint on shape**. Higher-order EOFs are often trying to represent noise while remaining orthogonal to everything before them — so their structure is partly an artifact of the EOFs that preceded them, not of the physics.

North et al. (1982) argue that the significance of an EOF is a function of how well separated the eigenvalues are. The 95% confidence bounds on an eigenvalue are

$$\boxed{\Delta\lambda = \lambda\sqrt{\frac{2}{N^*}}}$$

Here $N^*$ is our friend the **effective degrees of freedom** from Week 5-8 — autocorrelated data has fewer independent samples than it has data points. This is not easy to compute, and note you need a single value representative of the entire data set.

If the confidence bounds of two eigenvalues overlap, those two EOFs are **effectively degenerate**: they cannot be distinguished, and their individual patterns should not be interpreted. This is the practical face of the caveat from the orthogonality proof in Part 3 — the eigenvectors are only orthogonal (and only uniquely determined) if the eigenvalues are distinct.

When in doubt, test robustness other ways:

- subdivide your sample and compare EOF structures — are they the same?
- are your results explainable through theory?
- are results sensitive to the size of the domain? (high sensitivity suggests you are fitting sines and cosines inside your boundaries)

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure: Eigenvalue spectrum with North et al. (1982) bounds* — variance explained (%) on the $y$-axis (0–20%) against eigenvalue index (1–25), with error bars showing the 95% confidence bounds. Two series: an orange "random sample" curve decaying gradually and smoothly from about 15% to under 5%, and a purple "sample with structure" curve starting near 20% at index 1, dropping sharply to about 7% by index 2, then decaying like the random sample — illustrating a single dominant, well-separated leading EOF.
:::

**Wave-like phenomena.** If a travelling wave is present in your data, it appears as **two eigenvalues close together and separated from the rest**. These are the paired sine and cosine components of the wave. Do not interpret them as two independent patterns; they are one propagating feature.

:::{admonition} Figure / in-class demonstration
:class: tip
*Matlab examples:* `EOF_from_scratch.m`, `EOF_smoothing.m`, `eof_noise_reduction.m`
:::

### Recipe 1: EOFs from the temporal covariance matrix

Assume $\mathbf{X} = [M\times N] = [\text{sampling}\times\text{space}]$; $\mathbf{e}_1, \mathbf{d}_1 = [N\times1]$ and $\mathbf{z}_1 = [M\times1]$.

**1. Prepare the data**
- calculate anomalies along the sampling dimension, call this $\mathbf{X}$ — may be just the sample mean, or may include removing the seasonal cycle
- decide whether to standardize (divide by the standard deviation along the sampling dimension) — call this $\mathbf{X}_w$
- **save your original anomaly matrix $\mathbf{X}$ — you will need it in step 8**

**2. Calculate the temporal covariance matrix**
- $\mathbf{C} = \frac{1}{M}\mathbf{X}_w^T\mathbf{X}_w$
- if area weighting is needed, multiply $\mathbf{X}$ by $\sqrt{\cos\theta}$ first (call it $\mathbf{X}_w$), then form $\mathbf{C}$

**3. Perform eigenanalysis**
- use a built-in routine; it returns eigenvectors $\mathbf{E}$ and eigenvalues $\boldsymbol{\Lambda}$
- `[E, Lambda] = eig(C);` in Matlab

**4. Calculate and plot variance explained**
- $\lambda_i/\sum_j\lambda_j$
- plot it and determine which of the first few are distinct — use the North et al. (1982) bounds $\Delta\lambda = \lambda\sqrt{2/N^*}$ and look for non-overlap

**5. Retain the EOFs you are interested in**
- keep those with interesting variance explained; this may be only the first column of $\mathbf{E}$, denoted $\mathbf{e}_1$

**6. Calculate the PCs**
- project the weighted/standardized anomaly matrix onto the eigenvectors: $\mathbf{z}_1 = \mathbf{X}_w\mathbf{e}_1$, $\mathbf{z}_2 = \mathbf{X}_w\mathbf{e}_2$, …

**7. Standardize the PCs**
- $\tilde{\mathbf{z}}_1 = (\mathbf{z}_1 - \mu_{z_1})/\sigma_{z_1}$
- $\tilde{\mathbf{z}}_1$ is now in units of $\sigma$ — this is the PC to use from here on

**8. Calculate EOFs in physical units**
- $\mathbf{d}_1 = \frac{1}{M}\tilde{\mathbf{z}}_1^T\mathbf{X}$, projecting onto the **original** anomaly data
- $\mathbf{d}_1$ is in the units of your data: the anomaly associated with a one-$\sigma$ variation of the PC
- if $\mathbf{X} \equiv \mathbf{X}_w$, then $\mathbf{d}_1$ has the same pattern as $\mathbf{e}_1$ with different values; but if you area-weighted or standardized, **the spatial patterns will genuinely differ**

### Recipe 2: EOFs from the spatial dispersion matrix

Same dimensions as Recipe 1. Use this when $N \gg M$.

**1. Prepare the data** — as in Recipe 1.

**2. Calculate the spatial dispersion matrix**
- $\mathbf{C} = \frac{1}{N}\mathbf{X}_w\mathbf{X}_w^T$
- for area weighting, let $\mathbf{X}_{w,\cos}$ be the elements of $\mathbf{X}_w$ multiplied by $\cos\theta$, then $\mathbf{C} = \frac{1}{N}\mathbf{X}_{w,\cos}\mathbf{X}_w^T$

**3. Perform eigenanalysis**
- the routine now returns your **PCs** in $\mathbf{E}$ and eigenvalues in $\boldsymbol{\Lambda}$ — the roles have swapped

**4. Calculate and plot variance explained** — as in Recipe 1.

**5. Retain the PCs you are interested in**
- keep the first column of $\mathbf{E}$, denoted $\mathbf{z}_1$ — this is actually the first PC, because you analyzed the dispersion matrix

**6. Standardize the PCs**
- $\tilde{\mathbf{z}}_1 = (\mathbf{z}_1 - \mu_{z_1})/\sigma_{z_1}$

**7. Calculate EOFs in physical units**
- $\mathbf{d}_1 = \frac{1}{M}\tilde{\mathbf{z}}_1^T\mathbf{X}$, $\mathbf{d}_2 = \frac{1}{M}\tilde{\mathbf{z}}_2^T\mathbf{X}$, …
- again in the units of the original data

Note that this recipe has one step fewer than Recipe 1: you never have to project to *get* the PCs, because the eigenanalysis handed them to you directly.

### Recipe 3: EOFs via SVD

Assume $\mathbf{X} = [N\times M] = [\text{space}\times\text{sampling}]$; $\mathbf{e}_1, \mathbf{d}_1 = [N\times1]$; $\mathbf{z}_1 = [M\times1]$. **Note the transposed convention.**

**1. Prepare the data**
- calculate anomalies along the sampling dimension
- decide whether to standardize — call it $\mathbf{X}_w$
- if area weighting is needed, weight $\mathbf{X}$ by $\cos(\text{latitude})$ and call it $\mathbf{X}_w$
- save the original anomaly matrix $\mathbf{X}$

**2. Perform SVD**
- `[U, Sigma, V] = svd(Xw);` in Matlab
- the eigenvectors are the columns of $\mathbf{U}$, the PCs the columns of $\mathbf{V}$, the singular values the diagonal of $\boldsymbol{\Sigma}$

**3. Calculate and plot variance explained**
- $s_i^2/\sum_j s_j^2$ — note the **square**, since singular values are square roots of eigenvalues
- apply the North et al. bounds as before

**4. Retain the EOFs you are interested in**
- the first column of $\mathbf{U}$, denoted $\mathbf{e}_1$

**5. Calculate PCs**
- already done for you: PC 1 is the first column of $\mathbf{V}$, denoted $\mathbf{z}_1$

**6. Standardize the PCs**
- $\tilde{\mathbf{z}}_1 = (\mathbf{z}_1 - \mu_{z_1})/\sigma_{z_1}$

**7. Calculate EOFs in physical units**
- $\mathbf{d}_1 = \frac{1}{M}\mathbf{X}\tilde{\mathbf{z}}_1$, $\mathbf{d}_2 = \frac{1}{M}\mathbf{X}\tilde{\mathbf{z}}_2$, …
- note the ordering differs from Recipes 1–2 because of the transposed convention
- as before, if you weighted or standardized, the pattern of $\mathbf{d}_1$ will genuinely differ from $\mathbf{e}_1$

---

## Part 6 — Cluster Analysis

Cluster analysis is a data-mining technique whose goal is to group observations into a number of distinct clusters. Each observation ends up in **exactly one** cluster, with the assignment determined by a cost function. Finding the optimal set of clusters is usually iterative, starting from an initial guess.

The contrast with EOF analysis is the important part. Recall that a principal component tells you how much a given observation resembles a particular EOF pattern. An observation may look predominantly like one EOF — a large PC for that mode — but it will still have non-zero PCs for all the others. Every observation is a *mixture* of all the patterns.

In cluster analysis this is not so: each observation belongs to one cluster and no other. EOFs decompose; clusters partition. Which you want depends on whether you believe your system moves continuously through a state space or occupies discrete preferred regimes.

### k-means clustering

k-means is the most popular clustering technique in Earth science, largely for its simplicity. The aim is to group $N$ observations into $k$ clusters, each observation belonging to the cluster with the nearest center. This partitions the data into **Voronoi cells**.

The optimal distribution of cluster centers minimizes the within-cluster sum of squares — the sum of distances from each point to its cluster center:

$$\boxed{\underset{S}{\arg\min}\ \sum_{i=1}^{k}\sum_{\mathbf{x}\in S_i}\|\mathbf{x} - \boldsymbol{\mu}_i\|^2}$$

where $S = \{S_1, S_2, \ldots, S_k\}$ is the set of $k$ clusters, $\mathbf{x}$ the observations, $\boldsymbol{\mu}_i$ the center of cluster $S_i$ (the mean of all points in it), and $\arg\min$ indicates we seek the minimum over possible partitions $S$.

Note that this is a sum of **squared** distances — the same least-squares logic as linear regression in Week 5-8 and the same squared projection that defined EOF 1 in Part 3. Squaring keeps the cost positive definite and makes the minimization tractable, at the price of weighting outliers heavily.

While the idea is straightforward, computing the true optimum is very difficult and time-consuming for large data sets. Most implementations therefore use **heuristic algorithms** that converge quickly to a *local* optimum that may or may not be the global one. In practice you run k-means many times from different initial guesses and keep the iteration with the smallest within-cluster sum of squares.

The standard algorithm is **Lloyd's algorithm**: begin with an initial guess of the centers $\boldsymbol{\mu}$, then iteratively (i) assign each point to its nearest center and (ii) recompute each center as the mean of its assigned points, until nothing changes.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure: Lloyd's algorithm* — a sequence of six scatter panels ("Data and Initial Guess", "Initial Clusters", "Iteration 1" through "Iteration 4") applying Lloyd's algorithm to a small 2-D data set with 3 clusters (red, green, blue). Cluster centers are stars; outline stars denote the previous center position. The algorithm has converged by the third iteration, after which assignments and centers no longer change.
:::

Two decisions must be made by **you**, not by the algorithm:

**1. How many clusters?**
This is a user choice, not an output. You will get different results for different $k$ — and you will always get *an* answer, exactly as EOFs always give you an answer. Methods such as the **gap statistic** (Tibshirani et al. 2001) can help.

**2. How do you seed the initial guess?**
Many options exist. Running k-means many times and choosing the best solution is the standard way to avoid over-dependence on any single initialization.

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure: sensitivity to $k$* — three scatter panels using k-means on the same data set. Left: the raw unclustered data. Middle: coloured by 7 clusters. Right: coloured by 4 clusters. Stars denote centroids. The algorithm was run 20 times and the iteration with the optimal cost function retained. The same data supports either partition; the choice of $k$ is yours to justify.
:::

### Self-organizing maps (SOMs)

A self-organizing map — also called a Kohonen map — is another form of cluster analysis, using unsupervised machine learning to train an artificial neural network. Under certain conditions it is identical to k-means. The distinguishing feature is that **when one cluster center is updated, its neighbouring centers are updated too**.

The result is a set of cluster centers (called **nodes**) organized so that similar nodes sit near one another and dissimilar nodes sit far apart when arranged on an $m\times p$ grid. You gain an interpretable *layout* — neighbouring nodes are physically similar states — which k-means does not give you.

Like most neural networks, the algorithm has two phases:

1. **Training** — where the heavy lifting happens. Nodes are updated by comparison against input observations, ingested one at a time; often the same training data is iterated through several times. The result is a set of trained nodes.
2. **Mapping** — classifying a new observation as belonging to one of the nodes. Sometimes the training data differ from the data you want to map; sometimes they are the same.

As with k-means, you must decide how many nodes to create — but here the choice is an $m\times p$ **grid**, so you must choose both $m$ and $p$. The resulting nodes depend strongly on the dimensions chosen.

Training begins by initializing the $m\cdot p$ nodes, either randomly or using the principal components of the observations (a nice link back to Part 3). The algorithm then ingests the training data either one observation at a time (**online training**) or as a single group (**batch training**). Each observation is compared to all $m\cdot p$ nodes and the **best match unit (BMU)** is identified as the node with the smallest Euclidean distance to that observation. Writing the BMU as $\mathbf{n}_i$, the BMU and its neighbours are updated as

$$\boxed{\mathbf{n}_j(t+1) = \mathbf{n}_j(t) + \alpha(t)\,h_{ci}(t)\,[\mathbf{x}(t) - \mathbf{n}_i(t)]}$$

where $1 \le j \le m\cdot p$ and

- **training time ($t$)**: how many iterations have been performed
- **learning rate ($\alpha$)**: how strongly to update the nodes given the observation; typically decreases with training time
- **neighbourhood function ($h$)**: describes how many surrounding nodes to update besides the BMU; typically decreases with training time

The learning rate and neighbourhood function must be chosen by the user. Commonly $\alpha$ starts large and decreases linearly with training time. Common neighbourhood functions include a 2-D Gaussian or the Epanechnikov function. It is good practice to set the neighbourhood parameters initially so that the **entire** SOM space is included, ensuring all nodes get updated; this can be tightened at later training times. SOM training is typically implemented in multiple stages with different parameter combinations (e.g. a "rough" stage followed by a "finetune" stage).

**Choosing among SOMs.** Because SOM analysis has so many free parameters, many different SOMs can be computed from one data set. Two metrics guide the choice:

1. **quantization error**: the average distance between each observation and its BMU
2. **topographic error**: the fraction of observations whose first and second BMUs are *not* adjacent units

You want to minimize both. Typically the topographic error is held near zero while the quantization error is minimized. Good practice is to perform SOM analysis several times with different parameters and choose based on these two errors.

:::{admonition} Example / deeper dive
:class: note
**A worked SOM: Christman Field, Fort Collins**

Hourly observations from Christman Field, Colorado, 1 January – 31 December 2016. The input is a 2-D matrix $\mathbf{X}$ of dimension $(8784, 6)$: 8784 hourly observation periods and 6 variables — temperature (°F), relative humidity (%), wind speed (mph), pressure (mb), solar radiation (W m$^{-2}$) and precipitation (mm).

Training was performed in two phases, rough and finetune, each with different training lengths and initial/final radii of influence for the neighbourhood function. Nodes were initialized with random data on a $20\times20$ grid (400 nodes).

Plotting the trained node weights variable by variable shows the physical structure the SOM has found: for example, times of very strong wind are often associated with low relative humidity and precipitation.

Plotting the **frequency of occurrence** of each node — how many of the 8784 observations had that node as their BMU — shows which weather states are common. Nodes on the left-hand side of the grid occur most frequently, and correspond to low wind, high relative humidity and generally low temperatures.

Note the interpretive payoff: because neighbouring nodes are similar by construction, the frequency map is readable as a *landscape* of preferred weather regimes, not just an unordered list of cluster counts.
:::

:::{admonition} Figure / in-class demonstration
:class: tip
*Figure: SOM node weights* — six panels, one per variable (temperature °F, RH %, wind mph, pressure mb, solar W m$^{-2}$, precipitation mm), each a $20\times20$ grid with its own colour scale, showing the pattern of node weights across the trained SOM.

*Figure: SOM frequency of occurrence* — a $20\times20$ grid of cells, each shaded by intensity and annotated with a count, showing how many of the 8784 hourly observations were assigned to that node as their best-match unit.
:::

---

## Summary

| | EOF / PCA | Cluster analysis |
|---|---|---|
| **What it produces** | orthogonal patterns, ordered by variance explained | a partition of observations into groups |
| **Membership** | every observation is a mixture of all patterns | every observation belongs to exactly one cluster |
| **What you choose** | how many EOFs to retain; weighting; standardization | $k$ (or $m\times p$); initialization; cost function |
| **Optimality** | first $k$ EOFs explain more variance than any other $k$ linear combinations | only a *local* optimum is guaranteed in practice |
| **Main hazard** | EOFs of noise look physical; orthogonality constrains higher modes | you always get $k$ clusters, whether or not they exist |

Three threads to carry forward:

1. **Everything here is projection.** The inner product that gave us the correlation coefficient in Week 5-8 is the same operation that defines EOF 1, and the same operation that returns an EOF to physical units. Linear algebra is the bookkeeping that lets us do it for thousands of grid points at once.

2. **Both methods will always return an answer.** EOFs of red noise give plausible-looking sines and cosines; k-means will happily partition structureless data into $k$ tidy groups. The burden of demonstrating that structure is real — via North et al. bounds, subsampling, domain-sensitivity tests, or physical theory — is entirely on you.

3. **Degeneracy is the practical limit.** The orthogonality of EOFs holds only when eigenvalues are distinct. When they are not — closely spaced eigenvalues, or the paired sine/cosine of a travelling wave — the individual patterns are not separately meaningful, however clean they look on a map.

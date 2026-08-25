
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hermiteStats: Distribution-Robust Statistics via Hermite Polynomial Quantile Modeling

<!-- badges: start -->

[![R-CMD-check](https://github.com/WLenhard/hermiteStats/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WLenhard/hermiteStats/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/hermiteStats)](https://CRAN.R-project.org/package=hermiteStats)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**`hermiteStats`** provides distribution-robust, regularized estimators
for two of the most fundamental statistics in quantitative research:

1.  **The Hermite-Mehler Correlation (r_HM):** a distribution-robust
    estimator of the raw-scale Pearson correlation coefficient.
2.  **The Distribution-Free Effect Size (d_reg):** a robust standardized
    mean difference for independent and paired samples that retains the
    metric and interpretation of Cohen’s *d* / Hedges’ *g*.

Both methods address a long-standing dilemma in robust statistics:
classical remedies for non-normal data — trimming, Winsorizing, or
converting to ranks — typically restore stability only by silently
changing the quantity being estimated, discarding genuine tail variance
in the process. `hermiteStats` instead models the empirical quantile
function via **regularized monotone polynomial smoothing** and derives
exact distributional moments in closed algebraic form using
**Probabilists’ Hermite polynomials** and **Mehler’s (1866) identity** —
recovering the original, interpretable estimand while gaining the
stability of a robust method.

------------------------------------------------------------------------

## Table of Contents

- [Key Features](#key-features)
- [Mathematical Foundations](#mathematical-foundations)
  - [1. Rank-Based Inverse-Normal Quantile
    Modeling](#1-rank-based-inverse-normal-quantile-modeling)
  - [2. The Orthogonal Hermite Basis and Univariate
    Moments](#2-the-orthogonal-hermite-basis-and-univariate-moments)
  - [3. Mehler’s Bilinear Expansion for
    Covariance](#3-mehlers-bilinear-expansion-for-covariance)
  - [4. Disentangling Construct Association from Scale Attenuation
    (A)](#4-disentangling-construct-association-from-scale-attenuation-a)
  - [5. Distribution-Free Effect Size
    (d_reg)](#5-distribution-free-effect-size-d_reg)
- [Installation](#installation)
- [Quick Start and Examples](#quick-start-and-examples)
  - [Estimating Distributional
    Moments](#estimating-distributional-moments)
  - [Hermite-Mehler Correlation
    (r_HM)](#hermite-mehler-correlation-r_hm)
  - [Correlation Matrices and Multivariate
    Analysis](#correlation-matrices-and-multivariate-analysis)
  - [Distribution-Free Effect Sizes
    (d_reg)](#distribution-free-effect-sizes-d_reg)
  - [Paired and Repeated-Measures
    Designs](#paired-and-repeated-measures-designs)
  - [Diagnostic Visualizations](#diagnostic-visualizations)
- [Citation](#citation)
- [References](#references)

------------------------------------------------------------------------

## Key Features

- **No Estimand Shift:** Retains the raw-scale Pearson correlation and
  standardized mean difference as the target parameters — unlike rank
  correlations or Winsorized effect sizes, which silently estimate a
  different population quantity.
- **Variance Reduction in Small Samples (n \< 50):** Regularized,
  low-degree polynomial quantile smoothing prevents extreme tail
  observations from destabilizing sample covariances and pooled standard
  deviations.
- **Negligible Efficiency Loss under Normality:** Closely tracks the
  classical Pearson’s *r* / Cohen’s *d* when parametric assumptions
  hold, so there is essentially nothing to lose by using the regularized
  estimator as a default choice.
- **Shape Attenuation Diagnostic (A):** Quantifies how much of the
  observed association is suppressed purely by marginal shape mismatch
  (Hoeffding-Fréchet bounds) — a diagnostic with no counterpart in
  classical correlation analysis.
- **Flexible Monotonicity Enforcement:** Supports both a relaxed,
  rank-concordance-based check and a strict, analytically verified
  monotonicity constraint (root-finding on the fitted derivative).
- **Fully Closed-Form Computation:** All moments and covariances are
  obtained algebraically via Hermite/Mehler identities — no numerical
  integration, simulation, or iterative optimization is required.

------------------------------------------------------------------------

## Mathematical Foundations

### 1. Rank-Based Inverse-Normal Quantile Modeling

Let $X$ and $Y$ be continuous random variables with unknown cumulative
distribution functions $F_X$ and $F_Y$. The observed data are mapped to
latent standard normal scores $Z \sim \mathcal{N}(0, 1)$ using a
rank-based inverse-normal transformation (Normal Quantile
Transformation; e.g. Hazen, 1914):

$$Z_{x, i} = \Phi^{-1}\left(\frac{\operatorname{rank}(X_i) - 0.5}{n}\right), \quad Z_{y, i} = \Phi^{-1}\left(\frac{\operatorname{rank}(Y_i) - 0.5}{n}\right)$$

Unlike a plain rank transform, this step is only a means to an end: the
inverse quantile functions $X = f(Z_x)$ and $Y = g(Z_y)$ are then
modeled by regressing the raw values on powers of $Z$,

$$X = \sum_{j=0}^{k_x} \beta_{x, j} Z_x^j, \quad Y = \sum_{j=0}^{k_y} \beta_{y, j} Z_y^j,$$

so that the original metric and tail shape of the data are preserved
rather than discarded. To guarantee a valid quantile map, monotonicity
of $f$ (and $g$) is enforced over the observed range of $Z$,
automatically stepping down the polynomial degree whenever a violation
is detected.

### 2. The Orthogonal Hermite Basis and Univariate Moments

Standard powers $Z^j$ are not orthogonal under the Gaussian measure,
which makes their coefficients awkward for moment extraction.
Re-expressing the polynomial in the basis of **monic Probabilists’
Hermite polynomials** $He_m(z)$ (Hermite, 1864),

$$He_0(z) = 1, \quad He_1(z) = z, \quad He_2(z) = z^2 - 1, \quad He_3(z) = z^3 - 3z, \quad \dots$$

solves this: because $He_m(z)$ is orthogonal with respect to the
standard normal density,
$\mathbb{E}[He_m(Z) He_n(Z)] = m! \, \delta_{mn}$, the monomial
coefficients $\boldsymbol{\beta}$ map to Hermite coefficients
$\mathbf{a}$ via a fixed linear operator,
$\mathbf{a} = \mathbf{H}\boldsymbol{\beta}$, with

$$H_{m+1, n+1} = \binom{n}{2k} (2k - 1)!! \quad \text{for } 2k = n - m \ge 0 \text{ (and } 0 \text{ otherwise)}.$$

The mean and variance then follow immediately from the Hermite
coefficients:

$$\mathbb{E}[X] = a_0, \qquad \operatorname{Var}(X) = \sum_{m=1}^{k_x} a_m^2 \, m!$$

Skewness $\gamma_1$ and excess kurtosis $\gamma_2$ follow analogously
from the third and fourth central moments of the fitted polynomial,
evaluated using the raw Gaussian moments $\mathbb{E}[Z^j]$ — a direct
consequence of Isserlis’ (1918) theorem for products of jointly normal
variables:

$$\gamma_1 = \frac{\mathbb{E}[(X-\mu)^3]}{\sigma^3}, \qquad \gamma_2 = \frac{\mathbb{E}[(X-\mu)^4]}{\sigma^4} - 3$$

All four moments are obtained this way in exact closed algebraic form,
with no numerical integration at any step.

### 3. Mehler’s Bilinear Expansion for Covariance

Under a latent Gaussian copula with latent correlation
$\rho_z = \operatorname{cor}(Z_x, Z_y)$, **Mehler’s (1866) identity**
states that cross-products of Hermite polynomials of different order
vanish in expectation, and matching orders collapse to a single
deterministic power of $\rho_z$:

$$\mathbb{E}\left[ He_m(Z_x) He_n(Z_y) \right] = \delta_{mn} \, m! \, \rho_z^m$$

The manifest covariance is therefore the sum over matching Hermite
orders of the two fitted quantile polynomials ($\mathbf{a}$ for $X$,
$\mathbf{b}$ for $Y$):

$$\operatorname{Cov}(X, Y) = \sum_{m=1}^{\min(k_x, k_y)} a_m b_m \, m! \, \rho_z^m$$

which, combined with the univariate variances from Section 2, yields the
**Hermite-Mehler correlation (r_HM)**:

$$r_{\text{HM}} = \frac{\sum_{m=1}^{\min(k_x, k_y)} a_m b_m \, m! \, \rho_z^m}{\sqrt{\left(\sum_{m=1}^{k_x} a_m^2 m!\right)\left(\sum_{m=1}^{k_y} b_m^2 m!\right)}}$$

### 4. Disentangling Construct Association from Scale Attenuation (A)

This framework cleanly separates two conceptually distinct quantities:
the pure latent dependency $\rho_z$ (the Gaussian copula parameter,
unaffected by either variable’s shape), and the manifest linear
correlation $r_{\text{HM}}$ (which *is* affected by shape, exactly as
classical Pearson’s $r$ would be). Their ratio defines the **Shape
Attenuation Factor (A)**:

$$A = \frac{r_{\text{HM}}}{\rho_z}$$

- If $A \approx 1.0$, the marginal distributions do not meaningfully
  constrain the correlation.
- If $A \ll 1.0$, part of what looks like a “weak” Pearson correlation
  is in fact a mathematical ceiling imposed by mismatched marginal
  shapes (Hoeffding-Fréchet bounds; Carroll, 1961) rather than weak
  underlying dependence.

### 5. Distribution-Free Effect Size (d_reg)

For two independent groups, marginal quantile polynomials
$X_1 = f_1(Z_1)$ and $X_2 = f_2(Z_2)$ are fitted separately to obtain
regularized population moments
$\hat{\mu}_1, \hat{\mu}_2, \hat{\sigma}_1^2, \hat{\sigma}_2^2$. The
standardized mean difference is then

$$d_{\text{reg}} = \frac{\hat{\mu}_2 - \hat{\mu}_1}{\hat{\sigma}_{\text{avg}}}, \qquad \hat{\sigma}_{\text{avg}} = \sqrt{\frac{\hat{\sigma}_1^2 + \hat{\sigma}_2^2}{2}}$$

For **paired / repeated-measures designs** (`paired = TRUE`), the
observed difference scores $D = X_2 - X_1$ are modeled directly with
their own regularized quantile fit, giving the standardized mean change
$d_z = \hat{\mu}_D / \hat{\sigma}_D$. This $\hat{\sigma}_D$ is
theoretically linked to the marginal moments and the paired
Hermite-Mehler correlation via the general variance-of-a-difference
identity

$$\hat{\sigma}_D = \sqrt{\hat{\sigma}_1^2 + \hat{\sigma}_2^2 - 2\, r_{\text{HM}} \, \hat{\sigma}_1 \hat{\sigma}_2}$$

(no equal-variance assumption required), which makes explicit why a
stronger paired correlation mechanically shrinks the variance of the
change scores, and hence inflates $d_z$ relative to $d_{\text{reg}}$.

------------------------------------------------------------------------

## Installation

`hermiteStats` is not yet on CRAN. Install the development version
directly from GitHub:

``` r
# install.packages(\"remotes\")
remotes::install_github(\"WLenhard/hermiteStats\")
```

The package depends only on base R (`graphics`, `grDevices`, `stats`)
and requires R \>= 3.5.0 — there are no heavyweight dependencies to
resolve.

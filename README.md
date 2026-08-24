---
title: "README"
output: html_document
---

# hermiteStats: Distribution-Robust Statistics via Hermite Polynomial Quantile Modeling

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/yourusername/hermiteStats/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/hermiteStats)](https://CRAN.R-project.org/package=hermiteStats)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**`hermiteStats`** provides distribution-robust, regularized estimators for two of the most fundamental statistics in quantitative research:
1. **The Hermite–Mehler Correlation ($r_{\text{HM}}$):** A distribution-robust estimator of the raw-scale Pearson correlation coefficient.
2. **The Distribution-Free Effect Size ($d_{\text{reg}}$):** A robust standardized mean difference for independent and paired samples that retains the metric and interpretation of Cohen’s $d$ / Hedges’ $g$.

Both methods overcome the classical **"Procrustean bed"** dilemma in robust statistics: rather than forcing non-normal data into normal-theory assumptions or amputating genuine tail variation through trimming/Winsorization (which silently changes the target estimand), `hermiteStats` models the empirical quantile function via **monotone polynomial smoothing** and evaluates exact distributional moments analytically in closed algebraic form using **Probabilists' Hermite polynomials** and **Mehler’s (1866) identity**.

---

## Table of Contents

- [Key Features](#key-features)
- [Mathematical Foundations](#mathematical-foundations)
  - [1. Rank-Based Inverse-Normal Quantile Modeling](#1-rank-based-inverse-normal-quantile-modeling)
  - [2. The Orthogonal Hermite Basis & Univariate Moments](#2-the-orthogonal-hermite-basis--univariate-moments)
  - [3. Mehler’s Bilinear Expansion for Covariance](#3-mehlers-bilinear-expansion-for-covariance)
  - [4. Disentangling Construct Association from Scale Attenuation ($A$)](#4-disentangling-construct-association-from-scale-attenuation-a)
  - [5. Distribution-Free Effect Size ($d_{\text{reg}}$)](#5-distribution-free-effect-size-d_textreg)
- [Installation](#installation)
- [Quick Start & Examples](#quick-start--examples)
  - [Hermite–Mehler Correlation ($r_{\text{HM}}$)](#hermite-mehler-correlation-r_texthm)
  - [Correlation Matrices & Multivariate Analysis](#correlation-matrices--multivariate-analysis)
  - [Distribution-Free Effect Sizes ($d_{\text{reg}}$)](#distribution-free-effect-sizes-d_textreg)
  - [Paired / Repeated-Measures Designs](#paired--repeated-measures-designs)
  - [Diagnostic Visualizations](#diagnostic-visualizations)
- [Central Simulation Findings](#central-simulation-findings)
  - [Hermite–Mehler Correlation Simulation ($N = 1,500,000$ runs)](#hermite-mehler-correlation-simulation-n--1500000-runs)
  - [Effect Size $d_{\text{reg}}$ Simulation ($N = 726,000$ runs)](#effect-size-d_textreg-simulation-n--726000-runs)
- [Citation](#citation)
- [References](#references)

---

## Key Features

* **No Estimand Shift:** Retains the exact raw-scale Pearson correlation $\rho_{\text{Pearson}}$ and standardized mean difference $\delta_{\text{SMD}}$ as target parameters.
* **Massive Variance Reduction in Small Samples ($n < 50$):** Regularization via low-degree polynomial quantile smoothing prevents extreme tail observations from destabilizing sample covariances and pooled standard deviations.
* **Zero Efficiency Loss under Normality:** Achieves $\approx 100.3\%$ relative efficiency compared to sample Pearson's $r$ when parametric assumptions are met.
* **Novel Scale Attenuation Diagnostic ($A$):** Quantifies how much empirical association is suppressed purely due to marginal shape mismatch (Hoeffding–Fréchet ceiling bounds).
* **Flexible Monotonicity Enforcement:** Includes both `"relaxed"` (rank-concordance) and `"strict"` (analytical derivative root checks via `polyroot()`) monotonicity constraints.
* **Fast Closed-Form Computation:** Zero numerical quadrature or iterative optimization required for moment calculations.

---

## Mathematical Foundations

### 1. Rank-Based Inverse-Normal Quantile Modeling

Let $X$ and $Y$ be continuous random variables with unknown cumulative distribution functions $F_X$ and $F_Y$. The observed data are mapped to latent standard normal scores $Z \sim \mathcal{N}(0, 1)$ using the rank-based inverse-normal transformation (Normal Quantile Transformation; Fisher & Yates, 1938):

$$Z_{x, i} = \Phi^{-1}\left(\frac{\operatorname{rank}(X_i) - 0.5}{n}\right), \quad Z_{y, i} = \Phi^{-1}\left(\frac{\operatorname{rank}(Y_i) - 0.5}{n}\right)$$

The inverse quantile transformation functions $X = f(Z_x)$ and $Y = g(Z_y)$ are modeled via polynomial regression:

$$X = \sum_{j=0}^{k_x} \beta_{x, j} Z_x^j, \quad Y = \sum_{j=0}^{k_y} \beta_{y, j} Z_y^j$$

To guarantee valid quantile maps, the algorithm enforces monotonicity over the support $[z_{\min}, z_{\max}]$, automatically stepping down the degree if a violation is detected.

### 2. The Orthogonal Hermite Basis & Univariate Moments

Standard powers $Z^j$ are not orthogonal under the Gaussian measure. We re-express the polynomial in the basis of **monic Probabilists' Hermite polynomials** $He_m(z)$ (Hermite, 1864):

$$He_0(z) = 1, \quad He_1(z) = z, \quad He_2(z) = z^2 - 1, \quad He_3(z) = z^3 - 3z, \quad \dots$$

Because $He_m(z)$ forms an orthogonal system with respect to the standard normal density $\varphi(z)$, with $\mathbb{E}[He_m(Z) He_n(Z)] = m! \, \delta_{mn}$, the monomial coefficients $\boldsymbol{\beta}$ map to Hermite coefficients $\mathbf{a}$ via the linear operator $\mathbf{a} = \mathbf{H} \boldsymbol{\beta}$, where:

$$H_{m+1, n+1} = \binom{n}{2k} (2k - 1)!! \quad \text{for } 2k = n - m \ge 0$$

All distributional moments then emerge in exact closed algebraic form:

$$\mathbb{E}[X] = a_0, \quad \operatorname{Var}(X) = \sum_{m=1}^{k_x} a_m^2 \, m!$$

### 3. Mehler’s Bilinear Expansion for Covariance

Under a latent Gaussian copula with latent correlation $\rho_z = \operatorname{cor}(Z_x, Z_y)$, **Mehler's (1866) identity** states that the cross-covariances of Hermite polynomials collapse to a single deterministic power:

$$\mathbb{E}\left[ He_m(Z_x) He_n(Z_y) \right] = \delta_{mn} \, m! \, \rho_z^m$$

All cross-terms of different orders vanish. The total covariance is the sum of matching Hermite orders:

$$\operatorname{Cov}(X, Y) = \sum_{m=1}^{\min(k_x, k_y)} a_m b_m \, m! \, \rho_z^m$$

Combining the regularized moments yields the **Hermite–Mehler correlation ($r_{\text{HM}}$)**:

$$r_{\text{HM}} = \frac{\sum_{m=1}^{\min(k_x, k_y)} a_m b_m \, m! \, \rho_z^m}{\sqrt{\left(\sum_{m=1}^{k_x} a_m^2 m!\right)\left(\sum_{m=1}^{k_y} b_m^2 m!\right)}}$$

### 4. Disentangling Construct Association from Scale Attenuation ($A$)

The framework separates the pure latent dependency ($\rho_z$, the Gaussian copula parameter) from the manifest linear correlation ($r_{\text{HM}}$). The ratio defines the **Shape Attenuation Factor ($A$)**:

$$A = \frac{r_{\text{HM}}}{\rho_z}$$

* If $A \approx 1.0$, marginal distributions do not constrain the correlation.
* If $A \ll 1.0$, the observed Pearson correlation is attenuated purely by scale and distributional incompatibility (Hoeffding–Fréchet bounds; Carroll, 1961).

### 5. Distribution-Free Effect Size ($d_{\text{reg}}$)

For two groups, marginal quantile polynomials $X_1 = f_1(Z_1)$ and $X_2 = f_2(Z_2)$ are fitted separately to extract regularized population moments $\hat{\mu}_1, \hat{\mu}_2, \hat{\sigma}_1^2, \hat{\sigma}_2^2$. The standardized mean difference is calculated as:

$$d_{\text{reg}} = \frac{\hat{\mu}_2 - \hat{\mu}_1}{\hat{\sigma}_{\text{avg}}}, \quad \text{where } \hat{\sigma}_{\text{avg}} = \sqrt{\frac{\hat{\sigma}_1^2 + \hat{\sigma}_2^2}{2}}$$

In **paired/repeated-measures designs** (`paired = TRUE`), difference scores $D = X_2 - X_1$ are modeled with a Hermite quantile function to yield the **standardized mean change** $d_z = \hat{\mu}_D / \hat{\sigma}_D$, which is linked to the raw-metric effect size via the Hermite correlation:

$$\hat{\sigma}_D = \sqrt{\hat{\sigma}_1^2 + \hat{\sigma}_2^2 - 2 \, r_{\text{HM}} \, \hat{\sigma}_1 \hat{\sigma}_2} \implies d_{\text{reg}} = d_z \sqrt{1 - r_{\text{HM}}}$$

---

## Installation

You can install the development version of `hermiteStats` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("yourusername/hermiteStats")
```

---

## Quick Start & Examples

```r
library(hermiteStats)
```

### Hermite–Mehler Correlation ($r_{\text{HM}}$)

Estimate the linear correlation between a normally distributed variable and a skewed log-normal latency:

```r
set.seed(42)
n <- 40
z_latent <- rnorm(n)
x <- rnorm(n, mean = 100, sd = 15)                      # Normal construct
y <- exp(0.5 * z_latent + rnorm(n, sd = 0.5))           # Skewed latency (Lognormal)

# Calculate Hermite-Mehler correlation with Fisher-z 95% CI
fit_cor <- cor_hermite(x, y, conf_level = 0.95, ci_method = "fisher")
print(fit_cor)
#>   Hermite-Mehler Pearson Correlation (r_HM)
#> ------------------------------------------------
#>   Hermite Correlation (r_HM) :  0.428
#>   Latent Copula (rho_z)      :  0.514
#>   Shape Attenuation (A)      :  0.833
#>   Polynomial Degrees Fitted  :  X = 3, Y = 3
#>   Monotonicity Check         :  relaxed
#>   95% CI (fisher): [0.138, 0.651]
```

### Correlation Matrices & Multivariate Analysis

Compute full regularized correlation matrices for use in SEM, factor analysis, or meta-analysis:

```r
data(iris)
R_hermite <- cor_hermite(iris[, 1:4])
print(R_hermite)
```

### Distribution-Free Effect Sizes ($d_{\text{reg}}$)

#### Formula Interface (Independent Groups)

```r
set.seed(123)
df_study <- data.frame(
  score = c(rlnorm(25, meanlog = 2.0, sdlog = 0.5),   # Control (Skewed)
            rlnorm(25, meanlog = 2.4, sdlog = 0.5)),  # Treatment
  group = factor(rep(c("Control", "Treatment"), each = 25))
)

# Compute d_reg with 95% Percentile Bootstrap CI
fit_d <- d_reg(score ~ group, data = df_study, conf_level = 0.95, ci_method = "bootstrap")
print(fit_d)
#>   Distribution-Free Effect Size Estimation (d_reg)
#> ----------------------------------------------------
#>   Effect Size (d_reg)            :  0.642
#>   Averaged Model SD (sigma_avg)  :  4.812
#>   Hedges' g (Benchmark)          :  0.698
#>   Sample Sizes                   :  n1 = 25, n2 = 25
#>   Polynomial Degrees             :  g1 = 3, g2 = 3
#>   Monotonicity Check             :  relaxed
#>   95% CI (bootstrap): [0.084, 1.218]
```

### Paired / Repeated-Measures Designs

```r
pre  <- rlnorm(30, meanlog = 3.0, sdlog = 0.4)
post <- pre + rnorm(30, mean = 5.0, sd = 2.0)

fit_paired <- d_reg(pre, post, paired = TRUE, conf_level = 0.95)
print(fit_paired)
#>   Distribution-Free Effect Size Estimation (d_reg)
#> ----------------------------------------------------
#>   Effect Size (d_reg, raw scale)    :  0.518
#>   Standardized Mean Change (d_z)    :  1.942
#>   Paired Hermite Correlation (r_HM) :  0.865
#>   Averaged Model SD (sigma_avg)     :  9.814
#>   Difference Model SD (sigma_diff)  :  2.621
#>   Hedges' g (Benchmark)             :  0.509
#>   Sample Sizes                      :  n1 = 30, n2 = 30
#>   Polynomial Degrees                :  g1 = 3, g2 = 3
#>   Monotonicity Check                :  relaxed
#>   95% CI (bootstrap): [1.412, 2.580]
```

### Diagnostic Visualizations

Inspect the empirical quantile mapping and copula projections using base graphics:

```r
# Diagnostic plot for Correlation
plot(fit_cor)

# Diagnostic plot for Group Quantile Fits
plot(fit_d)
```

---

## Central Simulation Findings

The methods in `hermiteStats` were evaluated across extensive Monte Carlo simulation studies crossing multiple distribution families, sample sizes, effect sizes, and copula dependencies.

### Hermite–Mehler Correlation Simulation ($N = 1,500,000$ runs)

Benchmarked against sample Pearson's $r$, Spearman's $\rho$, and Gaussian rank correlation across 10 distribution families (Normal, Lognormal, Mixed Normal, Normal vs. Lognormal mismatch, 1% Outlier contamination, Beta, 1PL IRT sum scores, Exponential, Negative Binomial, and Student's $t_3$) across sample sizes $n = 10$ to $1000$:

* **Efficiency under Normality:** $r_{\text{HM}}$ achieved an empirical relative efficiency of **100.3%** compared to Pearson's $r$ (MSE ratio = 1.003), confirming zero efficiency loss when assumptions hold.
* **Accuracy under Non-Normality:** $r_{\text{HM}}$ produced lower Mean Squared Error (MSE) than Pearson's $r$ in **83.4%** of Gaussian copula conditions. Efficiency gains were largest for skewed and heavy-tailed distributions (MSE ratios: **1.584** for Lognormal, **1.509** for Exponential, **1.410** for Mixed Normal, **1.311** for $t_3$).
* **Variance Reduction:** Variance ratios reached **1.603** for lognormal and **1.554** for exponential data.
* **Confidence Interval Coverage:** Fisher $z$-transformed intervals applied to $r_{\text{HM}}$ achieved **95.8%** coverage across all non-normal continuous distributions (holding nominal 95% coverage across all sample sizes), whereas rank-based intervals collapsed toward 70%–76% due to asymptotic estimand shift.
* **Copula Misspecification:** Under an asymmetric **Clayton copula** (strict lower-tail dependence), $r_{\text{HM}}$ retained lower MSE than Pearson's $r$ in **74.5%** of all conditions.

### Effect Size $d_{\text{reg}}$ Simulation ($N = 726,000$ runs)

Benchmarked against small-sample corrected Hedges' $g$ (with noncentral-$t$ lambda-prime CIs) and robust Winsorized $d_{\text{AKP}}$ (Algina et al., 2005) across 6 distributions, sample sizes $n = 10$ to $100$, and true $\delta = 0.0$ to $1.0$:

* **Overall Superiority:** $d_{\text{reg}}$ yielded lower MSE in **76.0%** and lower sampling variance in **83.4%** of conditions compared to Hedges' $g$, and outperformed $d_{\text{AKP}}$ in **100.0%** of conditions.
* **Sweet Spot ($n < 50, |\delta| \le 0.8$):** The average MSE reduction over Hedges' $g$ was **11.5% at $n = 10$**, **8.1% at $n = 20$**, and **4.8% at $n = 30$**.
* **Variance Reduction:** Average variance reduction of **8.3%** for $n < 50$ over Hedges' $g$ and **30.5%** over $d_{\text{AKP}}$.
* **Coverage Stability:** Nominal 95% bootstrap intervals for $d_{\text{reg}}$ maintained **95.3%** empirical coverage across all distributions.

---

## Citation

If you use `hermiteStats` in published research, please cite the foundational manuscripts:

```bibtex
@article{lenhard2026hermite,
  title   = {The {Hermite--Mehler} Correlation: A Distribution-Robust Estimator of the {Pearson} Correlation Coefficient},
  author  = {Lenhard, Wolfgang and Lenhard, Alexandra},
  journal = {Behavior Research Methods},
  year    = {2026},
  doi     = {10.3758/s13428-xxx-xxxxx-x}
}

@article{lenhard2026dreg,
  title   = {Distribution-Free Effect Size Estimation: A Robust Alternative to {Cohen's d} and Other Effect Size Estimators},
  author  = {Lenhard, Wolfgang and Lenhard, Alexandra},
  journal = {Behavior Research Methods},
  year    = {2026},
  doi     = {10.3758/s13428-xxx-xxxxx-y}
}
```

---

## References

* Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to Cohen’s standardized mean difference effect size. *Psychological Methods*, 10(3), 317–328.
* Carroll, J. B. (1961). The nature of the data, or how to choose a correlation coefficient. *Psychometrika*, 26(4), 347–372.
* Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Lawrence Erlbaum Associates.
* Fisher, R. A., & Yates, F. (1938). *Statistical Tables for Biological, Agricultural, and Medical Research*. Oliver & Boyd.
* Hedges, L. V. (1981). Distribution theory for Glass’s estimator of effect size and related estimators. *Journal of Educational Statistics*, 6(2), 107–128.
* Hermite, C. (1864). Sur un nouveau développement en série des fonctions. *Comptes Rendus de l'Académie des Sciences, Paris*, 58, 93–100.
* Isserlis, L. (1918). On a formula for the product-moment coefficient of any order of a normal frequency distribution. *Biometrika*, 12(1/2), 134–139.
* Mehler, F. G. (1866). Ueber die Entwicklung einer Function von beliebig vielen Variabeln nach Laplaceschen Functionen höherer Ordnung. *Journal für die reine und angewandte Mathematik*, 66, 161–176.
* Sklar, M. (1959). Fonctions de répartition à $n$ dimensions et leurs marges. *Publications de l'Institut de Statistique de l'Université de Paris*, 8, 229–231.


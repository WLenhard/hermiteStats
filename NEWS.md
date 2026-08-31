# hermiteStats 0.4.0
31.08.2026

## New Features

* **New unified test interface `hermite_test()`** with
  `test = c("complete", "mean", "median", "variance", "asymmetry",
  "tailweight", "scale")`. The default `"complete"` performs a joint
  five-dimensional distributional profile test (mean, median, scale,
  asymmetry, tail weight) against a single joint permutation null with
  common permutation indices, Westfall-Young step-down family-wise error
  control across all five contrasts, and a global omnibus test of
  distributional equality (`minP` by default; `maxT` and Hotelling-type
  `T2` also available). `test = "variance"` maps to the scale contrast
  and additionally reports the variance ratio.
* **New `median_hermite()`**: a distribution-robust test of population
  median differences based on the regularized median `f(0) = beta_0`.
  Solves the Type I error inflation of naive permutation tests on raw
  sample medians and outperforms Welch's t-test in power on skewed and
  heavy-tailed data.
* **New omnibus option `"T2"`** in `shape_hermite()`: a permutation
  Hotelling-type Mahalanobis statistic, powerful when shape differences
  are spread across several contrasts.
* **New `predict.hermite_fit()`**: evaluates the fitted monotone quantile
  polynomial at arbitrary standard normal scores, e.g.
  `predict(fit, z = qnorm(c(.25, .5, .75)))` for regularized quartiles.
* `shape_hermite()` now reports the **variance ratio**
  (`exp(2 * delta log sigma)`) whenever the scale contrast is tested.

## Breaking / Behavioral Changes

* **Permutation inference is now the default** in `t_hermite()` and
  `median_hermite()` (`method = "permutation"`); the analytical mode
  remains available but is documented as an asymptotic approximation.
* Permutation defaults unified at `nperm = 1000L` across all tests
  (`shape_hermite()` previously used 2000).
* In analytical mode, `n_perm_success` is now `0` (previously it
  misleadingly reported `nperm`).
* All permutation p-values consistently use the add-one convention
  `(1 + #{|T*| >= |T|}) / (B + 1)`.

## Internal Refactoring

* Single shared profile extractor (`.hermite_profile()`) now supplies all
  test statistics (mean, median, sigma, log sigma, c2, c3, and their
  standard errors) from one quantile fit; `t_hermite()`,
  `median_hermite()`, `shape_hermite()`, and `hermite_test()` are thin
  layers over this common engine.
* Removed a dead OLS-covariance code path in the location statistics
  (computed and then overwritten in v0.3.x); permutation loops are
  substantially faster as a result.
* The sorted-sample QR fast path (fixed design per sample size and degree)
  is now used by all tests; the basis cache is size-capped.
* Deduplicated Hermite basis evaluators: `.hermite_eval_basis()` is the
  single canonical evaluator (the duplicate in `plots.R` was removed);
  polynomial evaluation is centralized in the engine and exposed via
  `predict.hermite_fit()`.
* `t_hermite()` / `median_hermite()` refactored into shared internal
  location-test, print, and plot engines; identical behavior, roughly half
  the code.
* `cor_hermite()` matrix method: marginals are now fitted exactly once per
  column in both the complete-data and pairwise-deletion branches.

# hermiteStats 0.3.5

* Last release of the 0.3 series: `hermite_fit()`, `cor_hermite()` /
  `pcor_hermite()`, `d_reg()`, and the initial versions of `t_hermite()`
  and `shape_hermite()` with analytical defaults.

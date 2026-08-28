#' @title Distribution-Robust Permutation Test for Mean Differences (Hermite t-Test)
#'
#' @description
#' Computes a permutation test for the difference in means between two
#' independent groups, or between two paired (repeated) measurements, using
#' the distribution-robust, polynomial-quantile-regularized effect size
#' (\eqn{d_{\mathrm{reg}}} or \eqn{d_z}) from \code{\link{d_reg}} as the test
#' statistic, in place of the classical sample mean and standard deviation
#' used by Student's or Welch's \eqn{t}-test.
#'
#' @param x A numeric vector of observations for the first group (or the
#'   pre-measurement / baseline in a paired design).
#' @param y A numeric vector of observations for the second group (or the
#'   post-measurement in a paired design).
#' @param paired Logical; \code{TRUE} for a paired-samples (within-subjects,
#'   pre-post) design. Default \code{FALSE}.
#' @param degree Integer scalar; maximum polynomial degree used for the
#'   marginal quantile fits, passed to \code{\link{d_reg}} (default \code{3L}).
#' @param copula Character string; copula mode for the paired Hermite
#'   correlation reported by \code{\link{d_reg}} when \code{paired = TRUE}:
#'   \code{"none"} (default) or \code{"gaussian"}. Ignored for independent groups.
#' @param monotonicity Monotonicity constraint passed to \code{\link{hermite_fit}}:
#'   \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param type Effect-size variant used as the test statistic, passed to
#'   \code{\link{d_reg}}: \code{"regularized"} (default; \eqn{d_{\mathrm{reg}}}
#'   for independent groups, \eqn{d_z} for paired designs), \code{"hedges"},
#'   \code{"glass"}, or \code{"combined"}. See \code{\link{d_reg}} for details.
#' @param alternative Character string specifying the alternative hypothesis,
#'   in the direction of \eqn{\mathrm{mean}(y) - \mathrm{mean}(x)}:
#'   \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param nperm Integer; number of permutation resamples (default \code{2000L}).
#'   Resamples are drawn randomly rather than enumerated exhaustively, so the
#'   resulting \eqn{p}-value is a Monte Carlo approximation to the exact
#'   permutation \eqn{p}-value; increase \code{nperm} for greater precision.
#' @param min_success Minimum acceptable proportion of successfully fitted
#'   permutation replicates before a warning is issued (default \code{0.90}).
#'   Occasional failures (e.g. a degenerate resample) are excluded from the
#'   reference distribution rather than aborting the test.
#' @param ... Additional arguments passed to \code{\link{d_reg}} \strong{for the
#'   observed statistic only} (e.g. \code{conf_level} to additionally attach a
#'   confidence interval for the effect size). These arguments are deliberately
#'   \emph{not} forwarded to the internal permutation replicates, to avoid
#'   triggering redundant nested bootstrap/CI computations on every resample.
#'
#' @details
#' \subsection{Test statistic}{
#' The observed statistic converts the regularized effect size from
#' \code{\link{d_reg}} to a \eqn{t}-like quantity using the classical
#' relationship between Cohen's \eqn{d}/\eqn{d_z} and Student's \eqn{t}:
#' \deqn{t_{\mathrm{Hermite}} = \hat d \times \sqrt{\tilde n}, \qquad
#'       \tilde n = \begin{cases} n_1 & \text{(paired)} \\
#'                                 n_1 n_2 / (n_1 + n_2) & \text{(independent)}
#'                  \end{cases}}
#' Because \eqn{\hat d} is derived from an adaptively fitted, degree-selected
#' polynomial quantile model, its exact finite-sample null distribution is not
#' analytically tractable, so \eqn{p}-values are obtained by permutation
#' instead of the classical (asymptotic) \eqn{t} reference distribution. Note
#' that multiplying by the constant \eqn{\sqrt{\tilde n}} does not change the
#' resulting \eqn{p}-value (it is the same positive scalar for the observed
#' statistic and every permutation replicate) — it is included purely to
#' report a familiar \eqn{t}-like magnitude.
#' }
#' \subsection{Resampling scheme}{
#' On every replicate, the \emph{entire} \code{\link{d_reg}} pipeline (rank-based
#' NQT, adaptive polynomial degree selection, monotonicity checking) is refit
#' from scratch, so that the permutation null distribution reflects the same
#' estimation variability as the observed statistic (Edgington & Onghena, 2007).
#' \describe{
#'   \item{Independent groups}{The pooled sample \eqn{(x, y)} is randomly
#'     re-split into groups of size \eqn{n_1} and \eqn{n_2} without
#'     replacement, consistent with the null hypothesis of exchangeable group
#'     labels.}
#'   \item{Paired design}{For each pair, the observed difference
#'     \eqn{D_i = y_i - x_i} has its sign flipped with probability 0.5
#'     (Rademacher / sign-flip resampling), and a pseudo-pair
#'     \eqn{(x_i,\; x_i + s_i D_i)} is reconstructed. This preserves \eqn{x_i}'s
#'     actual marginal distribution and calls \code{\link{d_reg}} in exactly
#'     the same configuration as the observed analysis, while randomizing the
#'     direction of the paired difference under the null hypothesis of a
#'     symmetric, zero-centered change distribution (Anderson & ter Braak, 2003).}
#' }
#' }
#'
#' @return An S3 object of class \code{"t_hermite"} containing:
#' \describe{
#'   \item{\code{statistic}}{The observed Hermite \eqn{t}-like statistic.}
#'   \item{\code{p_value}}{The Monte Carlo permutation \eqn{p}-value.}
#'   \item{\code{estimate}}{The underlying \code{\link{d_reg}} point estimate
#'     (\eqn{d_{\mathrm{reg}}} or \eqn{d_z}, depending on \code{paired}/\code{type}).}
#'   \item{\code{d_reg}, \code{d_z}}{Raw-scale and (if paired) standardized
#'     mean-change effect sizes from \code{\link{d_reg}}.}
#'   \item{\code{mean_diff}}{Regularized mean difference (\eqn{\hat\mu_2 - \hat\mu_1}).}
#'   \item{\code{n1}, \code{n2}}{Group sample sizes.}
#'   \item{\code{paired}, \code{alternative}, \code{nperm}}{Settings used.}
#'   \item{\code{n_perm_success}, \code{success_rate}}{Number and proportion of
#'     permutation replicates that fitted successfully.}
#'   \item{\code{perm_distribution}}{The vector of successful permutation statistics.}
#'   \item{\code{d_reg_fit}}{The observed \code{\link{d_reg}} object underlying the test.}
#'   \item{\code{method}}{Character label describing the test.}
#' }
#'
#' @references
#' Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for multi-factorial analysis of variance. \emph{Journal of Statistical Computation and Simulation}, 73(2), 85-113. \doi{10.1080/00949650215733}
#'
#' Edgington, E. S., & Onghena, P. (2007). \emph{Randomization Tests} (4th ed.). Chapman & Hall/CRC.
#'
#' Lenhard, W., & Lenhard, A. (submitted). Distribution-Free Effect Size Estimation: A Robust Alternative to Cohen's d and other effect size estimators. \emph{Behavior Research Methods}.
#'
#' @seealso \code{\link{d_reg}}, \code{\link{hermite_fit}}, \code{\link{cor_hermite}}
#'
#' @examples
#' # 1. Independent groups, skewed data
#' set.seed(42)
#' ctrl <- rlnorm(25, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(25, meanlog = 2.3, sdlog = 0.5)
#' test1 <- t_hermite(ctrl, trt, nperm = 500)
#' print(test1)
#'
#' # 2. Paired / repeated-measures design
#' pre  <- rlnorm(20, meanlog = 3.0, sdlog = 0.4)
#' post <- pre + rnorm(20, mean = 3.0, sd = 1.5)
#' test2 <- t_hermite(pre, post, paired = TRUE, nperm = 500)
#' print(test2)
#'
#' # 3. One-sided alternative
#' test3 <- t_hermite(ctrl, trt, alternative = "less", nperm = 500)
#' print(test3)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
t_hermite <- function(x, y, paired = FALSE, degree = 3L,
                      copula = c("none", "gaussian"),
                      monotonicity = c("relaxed", "strict", "none"),
                      type = c("regularized", "hedges", "glass", "combined"),
                      alternative = c("two.sided", "greater", "less"),
                      nperm = 2000L, min_success = 0.90, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  type         <- match.arg(type)
  alternative  <- match.arg(alternative)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (paired && length(x) != length(y)) {
    stop("'x' and 'y' must have identical length for a paired test.")
  }
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L) {
    stop("'nperm' must be a positive integer.")
  }

  # Observed effect: full Hermite pipeline fit once on the real data.
  # '...' is honored here only (e.g. to additionally request a CI on the effect size).
  obs_fit <- d_reg(x, y, degree = degree, copula = copula, monotonicity = monotonicity,
                   paired = paired, type = type, ...)

  n1 <- obs_fit$n1
  n2 <- obs_fit$n2
  n_tilde <- if (paired) n1 else (n1 * n2) / (n1 + n2)
  t_obs <- obs_fit$estimate * sqrt(n_tilde)

  x1 <- obs_fit$x1
  x2 <- obs_fit$x2

  perm_t <- rep(NA_real_, nperm)

  if (paired) {
    diffs <- x2 - x1
    for (b in seq_len(nperm)) {
      flips   <- sample(c(-1, 1), n1, replace = TRUE)
      x2_perm <- x1 + flips * diffs
      perm_t[b] <- tryCatch(
        d_reg(x1, x2_perm, degree = degree, copula = copula, monotonicity = monotonicity,
              paired = TRUE, type = type)$estimate * sqrt(n_tilde),
        error = function(e) NA_real_
      )
    }
  } else {
    combined <- c(x1, x2)
    N <- n1 + n2
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)
      perm_t[b] <- tryCatch(
        d_reg(combined[idx], combined[-idx], degree = degree, copula = copula,
              monotonicity = monotonicity, paired = FALSE, type = type)$estimate * sqrt(n_tilde),
        error = function(e) NA_real_
      )
    }
  }

  success      <- is.finite(perm_t)
  n_success    <- sum(success)
  success_rate <- n_success / nperm

  if (n_success < 1L) {
    warning("All permutation replicates failed; cannot compute a p-value.", call. = FALSE)
    p_val <- NA_real_
  } else {
    if (success_rate < min_success) {
      warning(sprintf("Only %d of %d permutation replicates succeeded (%.1f%%).",
                      n_success, nperm, 100 * success_rate), call. = FALSE)
    }
    perm_valid <- perm_t[success]
    p_val <- switch(
      alternative,
      "two.sided" = (1 + sum(abs(perm_valid) >= abs(t_obs))) / (n_success + 1),
      "greater"   = (1 + sum(perm_valid >= t_obs))           / (n_success + 1),
      "less"      = (1 + sum(perm_valid <= t_obs))           / (n_success + 1)
    )
  }

  res <- list(
    statistic         = t_obs,
    p_value           = p_val,
    estimate          = obs_fit$estimate,
    d_reg             = obs_fit$d_reg,
    d_z               = obs_fit$d_z,
    mean_diff         = obs_fit$group2$mean - obs_fit$group1$mean,
    n1                = n1,
    n2                = n2,
    paired            = paired,
    alternative       = alternative,
    nperm             = nperm,
    n_perm_success    = n_success,
    success_rate      = success_rate,
    perm_distribution = if (n_success > 0L) perm_t[success] else numeric(0),
    d_reg_fit         = obs_fit,
    method            = if (paired) "Paired Hermite Permutation Test" else "Two-Sample Hermite Permutation Test"
  )
  class(res) <- "t_hermite"
  res
}

#' @export
print.t_hermite <- function(x, digits = 3L, ...) {
  eff_lab <- if (x$paired) "d_z" else "d_reg"
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("  t_Hermite Statistic        :  %.*f\n", digits, x$statistic))
  cat(sprintf("  Effect Size (%s)          :  %.*f\n", eff_lab, digits, x$estimate))
  cat(sprintf("  Mean Difference            :  %.*f\n", digits, x$mean_diff))
  cat(sprintf("  Sample Sizes               :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Alternative Hypothesis     :  %s\n", x$alternative))
  cat(sprintf("  Permutations (successful)  :  %d / %d (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))
  cat(sprintf("  p-value (Monte Carlo perm.):  %.*f\n", digits, x$p_value))
  cat("\n")
  invisible(x)
}

#' @title Distribution-Robust Permutation Test for Shape Differences (Variance, Skewness, Kurtosis)
#'
#' @description
#' Tests for differences in the regularized distributional \emph{shape}
#' (variance, skewness, excess kurtosis) between two independent groups or two
#' paired measurements, using the closed-form Hermite moments from
#' \code{\link{hermite_moments}} as the target quantities and a permutation
#' null distribution for inference. This is the shape-focused counterpart to
#' \code{\link{t_hermite}}, which tests location (mean) differences only.
#'
#' @param x A numeric vector of observations for the first group (or the
#'   pre-measurement / baseline in a paired design).
#' @param y A numeric vector of observations for the second group (or the
#'   post-measurement in a paired design).
#' @param paired Logical; \code{TRUE} for a paired-samples (within-subjects,
#'   pre-post) design. Default \code{FALSE}.
#' @param degree Integer scalar; maximum polynomial degree requested for the
#'   marginal quantile fits (default \code{3L}). Both groups are refit, if
#'   necessary, at a common \emph{realized} degree so that their moments are
#'   directly comparable; see Details.
#' @param monotonicity Monotonicity constraint passed to \code{\link{hermite_fit}}:
#'   \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Rank tie-handling method passed to \code{\link{hermite_fit}}:
#'   \code{"average"} (default) or \code{"random"}.
#' @param nperm Integer; number of permutation resamples (default \code{2000L}).
#'   Resamples are drawn randomly rather than enumerated exhaustively, so the
#'   resulting \eqn{p}-values are Monte Carlo approximations.
#' @param min_success Minimum acceptable proportion of successfully fitted
#'   permutation replicates before a warning is issued (default \code{0.90}).
#' @param omnibus Logical; if \code{TRUE} (default), additionally computes a
#'   permutation-based omnibus test combining all three shape moments into a
#'   single global test statistic (see Details).
#'
#' @details
#' \subsection{Matched polynomial degree}{
#' Because variance, skewness, and excess kurtosis are algebraic functions of
#' the fitted polynomial's coefficients (see \code{\link{hermite_fit}}), a
#' difference in \emph{realized} degree between the two groups' independently
#' adaptive fits would confound genuine shape differences with differences in
#' model flexibility. Before extracting moments, both groups are therefore
#' iteratively refit at a shared common degree: starting from \code{degree},
#' both groups are fit, and if either one required a lower degree to satisfy
#' the monotonicity constraint, the common cap is reduced to the lower of the
#' two realized degrees and both groups are refit; this repeats until the
#' realized degree is stable for both groups.
#' }
#' \subsection{Resampling scheme}{
#' On every replicate, the matched-degree fitting procedure above is rerun on
#' resampled data, so the permutation null reflects the same degree-selection
#' variability as the observed test.
#' \describe{
#'   \item{Independent groups}{The pooled sample is randomly re-split into
#'     groups of size \eqn{n_1} and \eqn{n_2} without replacement.}
#'   \item{Paired design}{For each pair \eqn{(x_i, y_i)}, the two values are
#'     swapped with probability 0.5. Unlike the sign-flip scheme used by
#'     \code{\link{t_hermite}} for testing the paired \emph{difference} series
#'     (which never needs to re-randomize either marginal's identity), testing
#'     genuine \emph{marginal} shape equality requires actually re-randomizing
#'     which member of each pair is labeled "group 1" vs. "group 2" under the
#'     null hypothesis of within-pair exchangeability.}
#' }
#' }
#' \subsection{Omnibus test}{
#' If \code{omnibus = TRUE}, the three observed moment differences are
#' standardized using the mean and standard deviation of their own permutation
#' distributions, and combined as a sum of squares
#' (\eqn{\sum_k z_k^2}, \eqn{k \in \{\text{variance, skewness, kurtosis}\}}),
#' analogous to O'Brien's (1984) global test for multiple endpoints. Its
#' permutation \eqn{p}-value is obtained by recomputing the same standardized
#' sum-of-squares statistic on every permutation replicate, which implicitly
#' captures the joint null covariance of the three moment estimators without
#' requiring it to be derived analytically.
#' }
#' \subsection{Scope}{
#' This function intentionally does not test the mean/location. Use
#' \code{\link{t_hermite}} for location inference; combining both functions
#' yields a full "mean and shape" comparison, each with a test scheme matched
#' to its own null hypothesis.
#' }
#'
#' @return An S3 object of class \code{"shape_hermite"} containing:
#' \describe{
#'   \item{\code{estimate}}{Named numeric vector of observed shape differences
#'     (Group 2 \eqn{-} Group 1): \code{variance}, \code{skewness}, \code{excess_kurtosis}.}
#'   \item{\code{p_value}}{Named vector of raw (unadjusted) permutation \eqn{p}-values per moment.}
#'   \item{\code{p_adjusted}}{Holm-adjusted \eqn{p}-values across the three moments.}
#'   \item{\code{omnibus}}{A list with \code{statistic} and \code{p_value} for the
#'     global shape test (\code{NULL} if \code{omnibus = FALSE}).}
#'   \item{\code{n1}, \code{n2}}{Group sample sizes.}
#'   \item{\code{degree_matched}, \code{degree_requested}}{The common realized
#'     polynomial degree used for both groups, and the originally requested maximum.}
#'   \item{\code{paired}, \code{nperm}}{Settings used.}
#'   \item{\code{n_perm_success}, \code{success_rate}}{Number/proportion of
#'     permutation replicates that fitted successfully.}
#'   \item{\code{perm_distribution}}{The \code{n_perm_success x 3} matrix of
#'     successful permutation moment differences.}
#'   \item{\code{fit1}, \code{fit2}}{The observed, degree-matched \code{\link{hermite_fit}} objects.}
#'   \item{\code{method}}{Character label describing the test.}
#' }
#'
#' @references
#' Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for multi-factorial analysis of variance. \emph{Journal of Statistical Computation and Simulation}, 73(2), 85-113. \doi{10.1080/00949650215733}
#'
#' Edgington, E. S., & Onghena, P. (2007). \emph{Randomization Tests} (4th ed.). Chapman & Hall/CRC.
#'
#' O'Brien, P. C. (1984). Procedures for comparing samples with multiple endpoints. \emph{Biometrics}, 40(4), 1079-1087. \doi{10.2307/2531158}
#'
#' Lenhard, W., & Lenhard, A. (submitted). Distribution-Free Effect Size Estimation: A Robust Alternative to Cohen's d and other effect size estimators. \emph{Behavior Research Methods}.
#'
#' @seealso \code{\link{t_hermite}}, \code{\link{d_reg}}, \code{\link{hermite_fit}}, \code{\link{hermite_moments}}
#'
#' @examples
#' # 1. Independent groups differing in skewness but not in mean/variance
#' set.seed(1)
#' g1 <- rnorm(60)
#' g2 <- as.numeric(scale(rgamma(60, shape = 2))) # right-skewed, same mean/sd after scaling
#' test1 <- shape_hermite(g1, g2, nperm = 500)
#' print(test1)
#'
#' # 2. Paired design
#' set.seed(2)
#' pre  <- rnorm(40)
#' post <- pre + rt(40, df = 3) * 0.3   # heavier-tailed change
#' test2 <- shape_hermite(pre, post, paired = TRUE, nperm = 500)
#' print(test2)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
shape_hermite <- function(x, y, paired = FALSE, degree = 3L,
                          monotonicity = c("relaxed", "strict", "none"),
                          ties_method = c("average", "random"),
                          nperm = 2000L, min_success = 0.90, omnibus = TRUE) {

  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (paired && length(x) != length(y)) {
    stop("'x' and 'y' must have identical length for a paired design.")
  }

  if (paired) {
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
  } else {
    x1 <- x[is.finite(x)]; x2 <- y[is.finite(y)]
  }
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 4L || n2 < 4L) stop("Each group must contain at least 4 valid observations.")

  moment_names <- c("variance", "skewness", "excess_kurtosis")

  # Fit both groups at a shared, realized common polynomial degree
  fit_matched <- function(a, b) {
    cap <- degree
    repeat {
      fa <- hermite_fit(a, degree = cap, monotonicity = monotonicity,
                        ties_method = ties_method, force_odd = TRUE)
      fb <- hermite_fit(b, degree = cap, monotonicity = monotonicity,
                        ties_method = ties_method, force_odd = TRUE)
      common <- min(fa$degree, fb$degree)
      if (common >= cap) break
      cap <- common
    }
    list(fit1 = fa, fit2 = fb, degree = cap)
  }

  extract_diff <- function(fm) {
    m1 <- hermite_moments(fm$fit1)
    m2 <- hermite_moments(fm$fit2)
    c(variance        = m2$variance - m1$variance,
      skewness        = m2$skewness - m1$skewness,
      excess_kurtosis = m2$excess_kurtosis - m1$excess_kurtosis)
  }

  obs_fm   <- fit_matched(x1, x2)
  obs_diff <- extract_diff(obs_fm)

  perm_diff <- matrix(NA_real_, nrow = nperm, ncol = 3L,
                      dimnames = list(NULL, moment_names))

  if (paired) {
    for (b in seq_len(nperm)) {
      swap <- sample(c(FALSE, TRUE), n1, replace = TRUE)
      a  <- ifelse(swap, x2, x1)
      bb <- ifelse(swap, x1, x2)
      perm_diff[b, ] <- tryCatch(extract_diff(fit_matched(a, bb)),
                                 error = function(e) rep(NA_real_, 3L))
    }
  } else {
    combined <- c(x1, x2); N <- n1 + n2
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)
      perm_diff[b, ] <- tryCatch(extract_diff(fit_matched(combined[idx], combined[-idx])),
                                 error = function(e) rep(NA_real_, 3L))
    }
  }

  success      <- stats::complete.cases(perm_diff)
  n_success    <- sum(success)
  success_rate <- n_success / nperm

  if (n_success < 1L) stop("All permutation replicates failed; cannot compute p-values.")
  if (success_rate < min_success) {
    warning(sprintf("Only %d of %d permutation replicates succeeded (%.1f%%).",
                    n_success, nperm, 100 * success_rate), call. = FALSE)
  }
  perm_valid <- perm_diff[success, , drop = FALSE]

  p_raw <- vapply(seq_along(moment_names), function(k) {
    (1 + sum(abs(perm_valid[, k]) >= abs(obs_diff[k]))) / (n_success + 1)
  }, numeric(1))
  names(p_raw) <- moment_names
  p_adj <- stats::p.adjust(p_raw, method = "holm")

  omnibus_res <- NULL
  if (omnibus) {
    ctr <- colMeans(perm_valid)
    scl <- apply(perm_valid, 2, stats::sd)
    ok_m <- is.finite(scl) & scl > 0
    if (!all(ok_m)) {
      warning("Degenerate permutation variance for: ",
              paste(moment_names[!ok_m], collapse = ", "),
              " -- excluded from the omnibus statistic.", call. = FALSE)
    }
    z_perm <- sweep(sweep(perm_valid[, ok_m, drop = FALSE], 2, ctr[ok_m], "-"), 2, scl[ok_m], "/")
    z_obs  <- (obs_diff[ok_m] - ctr[ok_m]) / scl[ok_m]
    stat_perm <- rowSums(z_perm^2)
    stat_obs  <- sum(z_obs^2)
    p_omni <- (1 + sum(stat_perm >= stat_obs)) / (length(stat_perm) + 1)
    omnibus_res <- list(statistic = stat_obs, p_value = p_omni)
  }

  res <- list(
    estimate          = obs_diff,
    p_value           = p_raw,
    p_adjusted        = p_adj,
    omnibus           = omnibus_res,
    n1                = n1,
    n2                = n2,
    degree_matched    = obs_fm$degree,
    degree_requested  = degree,
    paired            = paired,
    nperm             = nperm,
    n_perm_success    = n_success,
    success_rate      = success_rate,
    perm_distribution = perm_valid,
    fit1              = obs_fm$fit1,
    fit2              = obs_fm$fit2,
    method            = if (paired) "Paired Hermite Shape-Difference Permutation Test"
    else "Two-Sample Hermite Shape-Difference Permutation Test"
  )
  class(res) <- "shape_hermite"
  res
}

#' @export
print.shape_hermite <- function(x, digits = 3L, ...) {
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 58), "\n", sep = "")
  cat(sprintf("  Sample Sizes              :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Matched Polynomial Degree :  %d (requested: %d)\n",
              x$degree_matched, x$degree_requested))
  cat(sprintf("  Permutations (successful) :  %d / %d (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))

  cat("\n  Shape-Difference Profile (Group 2 - Group 1):\n")
  tab <- data.frame(
    Moment     = names(x$estimate),
    Difference = round(unname(x$estimate), digits),
    p_value    = round(unname(x$p_value), digits),
    p_holm     = round(unname(x$p_adjusted), digits)
  )
  print(tab, row.names = FALSE)

  if (!is.null(x$omnibus)) {
    cat("\n  Omnibus Shape Test (sum of squared standardized differences):\n")
    cat(sprintf("    Statistic = %.*f,  p-value = %.*f\n",
                digits, x$omnibus$statistic, digits, x$omnibus$p_value))
  }
  cat("\n")
  invisible(x)
}

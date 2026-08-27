#' @title Classical Benchmark Effect-Size Estimators and Confidence Intervals
#'
#' @description
#' Classical, non-regularized standardized mean difference estimators (Cohen's
#' d family) and the noncentral-t ("lambda-prime") confidence interval method
#' used to bracket them analytically. These functions serve as the benchmark
#' comparisons reported alongside \code{\link{d_reg}}, and are also fully
#' documented and exported for standalone use.
#'
#' @name benchmarks
#' @keywords internal
NULL

# -----------------------------------------------------------------------------
# Exported Function: hedges_correction
# -----------------------------------------------------------------------------

#' Exact Small-Sample Bias Correction Factor for Hedges' g
#'
#' Cohen's \eqn{d} is slightly biased upward in small samples. Hedges (1981)
#' showed that multiplying by the exact correction factor
#' \deqn{J(df) = \frac{\Gamma(df/2)}{\sqrt{df/2}\,\Gamma((df-1)/2)}}
#' removes this bias, yielding Hedges' \eqn{g}. \eqn{J(df)} lies in
#' \eqn{(0, 1]} and approaches 1 as \code{df} grows, reflecting that the bias
#' vanishes asymptotically. The factor is evaluated on the log scale via
#' \code{\link[base]{lgamma}} for numerical stability at both very small and
#' very large degrees of freedom.
#'
#' @param df Numeric; degrees of freedom of the standard deviation estimate
#'   used to standardize the mean difference (e.g. \eqn{n_1+n_2-2} for the
#'   pooled-SD case). Values below \code{1} are floored, since \eqn{J(df)}
#'   is mathematically undefined (and tends to 0) as \code{df} approaches 1.
#'
#' @return Numeric scalar; the correction factor \eqn{J(df)}.
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, 6(2), 107-128. \doi{10.3102/10769986006002107}
#'
#' @examples
#' hedges_correction(8)    # noticeable correction in small samples
#' hedges_correction(100)  # negligible correction in large samples
#'
#' @export
hedges_correction <- function(df) {
  df <- max(df, 1.0001)
  exp(lgamma(df / 2.0) - 0.5 * log(df / 2.0) - lgamma((df - 1.0) / 2.0))
}

# -----------------------------------------------------------------------------
# Exported Functions: d_cohen, hedges_g
# -----------------------------------------------------------------------------

#' Standardized Mean Differences (Cohen's d Family)
#'
#' Computes classical, sample-moment-based standardized mean differences,
#' primarily as benchmark comparisons for the regularized
#' \code{\link{d_reg}} estimator.
#'
#' @param x1 Numeric vector of Group 1 (baseline) scores.
#' @param x2 Numeric vector of Group 2 (comparison) scores.
#' @param type Character string selecting the standardizer:
#'   \describe{
#'     \item{\code{"pooled"}}{(Default) Pooled within-group standard
#'       deviation, \eqn{s_p = \sqrt{[(n_1-1)s_1^2 + (n_2-1)s_2^2]/(n_1+n_2-2)}}
#'       — the classical Cohen's \eqn{d} / Hedges' \eqn{g} standardizer.}
#'     \item{\code{"avg"}}{Unweighted root-mean-square of the two group
#'       variances, \eqn{\sqrt{(s_1^2+s_2^2)/2}}, with degrees of freedom
#'       approximated via the Welch-Satterthwaite equation. Recommended when
#'       group variances are expected to differ.}
#'     \item{\code{"glass"}}{Glass's \eqn{\Delta}: standardizes by the Group 1
#'       (baseline) standard deviation only, appropriate when Group 1 alone
#'       represents the reference population variance (e.g. an untreated
#'       control group).}
#'   }
#' @param correct_bias Logical; if \code{TRUE} (default), multiplies the raw
#'   ratio by \code{\link{hedges_correction}}, yielding an (approximately)
#'   unbiased estimator.
#'
#' @return Numeric scalar. Returns \code{NA_real_} if either group has fewer
#'   than 2 observations, and \code{NA_real_} with a warning if the
#'   standardizer is zero, negative, or non-finite (i.e. the ratio is
#'   mathematically undefined).
#'
#' @examples
#' set.seed(1)
#' g1 <- rnorm(20, mean = 0)
#' g2 <- rnorm(20, mean = 0.5)
#' d_cohen(g1, g2)                       # pooled-SD Hedges' g (bias-corrected)
#' d_cohen(g1, g2, type = "glass")       # Glass's Delta
#' d_cohen(g1, g2, correct_bias = FALSE) # uncorrected Cohen's d
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, 6(2), 107-128. \doi{10.3102/10769986006002107}
#'
#' @seealso \code{\link{hedges_g}}, \code{\link{d_reg}}, \code{\link{hedges_correction}}
#' @export
d_cohen <- function(x1, x2, type = c("pooled", "avg", "glass"), correct_bias = TRUE) {
  type <- match.arg(type)
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 2L || n2 < 2L) return(NA_real_)
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  md <- mean(x2) - mean(x1)

  if (type == "pooled") {
    sdv <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
    df  <- n1 + n2 - 2L
  } else if (type == "avg") {
    sdv <- sqrt((s1^2 + s2^2) / 2.0)
    df  <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1 - 1L) + (s2^2/n2)^2/(n2 - 1L))
  } else {
    sdv <- s1
    df  <- n1 - 1L
  }

  if (!is.finite(sdv) || sdv <= 0) {
    warning("Standardizer is zero, negative, or undefined; effect size is not defined for this sample.")
    return(NA_real_)
  }

  d <- md / sdv
  if (correct_bias) d <- d * hedges_correction(df)
  d
}

#' @rdname d_cohen
#' @export
hedges_g <- function(x1, x2) {
  d_cohen(x1, x2, type = "pooled", correct_bias = TRUE)
}

# -----------------------------------------------------------------------------
# Exported Function: ci_nct
# -----------------------------------------------------------------------------

#' Noncentral-t (Lambda-Prime) Confidence Interval for a Standardized Mean Difference
#'
#' Computes an analytical confidence interval for a standardized mean
#' difference by inverting the noncentral \eqn{t} distribution (the
#' "lambda-prime" method: Steiger & Fouladi, 1997; Lecoutre, 2007). Given an
#' observed point estimate, the method numerically solves for the two
#' noncentrality parameters whose sampling distributions place the observed
#' test statistic at the \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles.
#'
#' @param d_point Numeric; the observed standardized mean difference.
#' @param n1 Integer; sample size of Group 1 (or the number of pairs, for
#'   paired designs).
#' @param n2 Integer; sample size of Group 2.
#' @param conf Numeric; confidence level (default \code{0.95}).
#' @param df Numeric; degrees of freedom of the standardizer. Defaults to
#'   \eqn{n_1+n_2-2} (independent two-sample, pooled- or averaged-SD case).
#' @param n_tilde Numeric or \code{NULL}; the effective sample size entering
#'   \eqn{t_{\mathrm{obs}} = d\sqrt{\tilde n}}. If \code{NULL} (default),
#'   the two-sample harmonic form \eqn{\tilde n = n_1 n_2/(n_1+n_2)} is used,
#'   appropriate for an independent-groups difference standardized by a
#'   pooled or averaged SD. Supply \code{n_tilde} explicitly for other
#'   designs: \code{n_tilde = n1} together with \code{df = n1 - 1} for a
#'   paired-sample \eqn{d_z} (\eqn{n_1} pairs) or for Glass's \eqn{\Delta}
#'   (standardized by the Group 1 SD alone).
#'
#' @return A named numeric vector \code{c(lower, upper)}.
#'
#' @references
#' Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. \emph{Journal of Modern Applied Statistical Methods}, 6(1), 107-116. \doi{10.22237/jmasm/1177992600}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of statistical models. In L. L. Harlow, S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no significance tests?} (pp. 221-257). Lawrence Erlbaum Associates.
#'
#' @examples
#' # Independent two-sample Cohen's d (default two-sample n_tilde)
#' ci_nct(d_point = 0.5, n1 = 30, n2 = 30)
#'
#' # Paired-sample d_z: n_tilde and df both based on the number of pairs
#' ci_nct(d_point = 0.5, n1 = 30, n2 = 30, df = 29, n_tilde = 30)
#'
#' @export
ci_nct <- function(d_point, n1, n2, conf = 0.95, df = n1 + n2 - 2L, n_tilde = NULL) {
  if (is.null(n_tilde)) n_tilde <- (n1 * n2) / (n1 + n2)
  t_obs   <- d_point * sqrt(n_tilde)
  a       <- (1 - conf) / 2

  ncp_from_t <- function(t_val, df_val, p_val) {
    f <- function(ncp) suppressWarnings(stats::pt(t_val, df = df_val, ncp = ncp)) - p_val
    stats::uniroot(f, c(t_val - 6, t_val + 6), extendInt = "yes")$root
  }

  lo <- ncp_from_t(t_obs, df, 1 - a) / sqrt(n_tilde)
  hi <- ncp_from_t(t_obs, df, a) / sqrt(n_tilde)
  c(lower = lo, upper = hi)
}

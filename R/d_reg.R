# =============================================================================
# d_reg.R -- Distribution-Free Effect Size Estimation
# =============================================================================

#' @title Distribution-Free Effect Size Estimation (d_reg)
#'
#' @description
#' Computes the distribution-robust standardized mean difference (\eqn{d_{\mathrm{reg}}})
#' for independent or paired two-group designs. \eqn{d_{\mathrm{reg}}} retains the
#' classical metric and interpretation of Cohen's \eqn{d} and Hedges' \eqn{g}, but
#' substantially reduces sampling variance and Mean Squared Error (MSE) under
#' non-normal continuous distributions by replacing raw sample means/variances
#' with regularized polynomial quantile models and closed-form Hermite moments.
#'
#' \code{d_reg()} is the package's \emph{estimation} tool for mean differences;
#' its permutation-based \emph{hypothesis test} companion is
#' \code{\link{t_hermite}} (or \code{\link{hermite_test}} for a joint test of
#' all distributional aspects).
#'
#' @param x A numeric vector of observations for the first group (or baseline),
#'   or a two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param y A numeric vector of observations for the second group (or intervention).
#'   Ignored when \code{x} is a formula.
#' @param data An optional data frame, list, or environment containing the
#'   variables referenced in the formula.
#' @param degree Integer scalar; maximum polynomial degree for the marginal
#'   quantile functions (default \code{3}). Automatically reduced if
#'   monotonicity is violated or too few unique values are available; see
#'   \code{\link{hermite_fit}}.
#' @param copula Character string specifying the copula mode used for the paired
#'   correlation (\code{\link{cor_hermite}}) when \code{paired = TRUE}:
#'   \code{"none"} (default, copula-free cross-moments) or \code{"gaussian"}
#'   (Mehler identity).
#' @param monotonicity Monotonicity constraint passed to \code{\link{hermite_fit}}:
#'   \describe{
#'     \item{\code{"relaxed"}}{(Default) Empirical rank-concordance criterion.
#'       Recommended for general applied use.}
#'     \item{\code{"strict"}}{Analytical guarantee that \eqn{f'(z) \ge 0}
#'       everywhere on the standard normal support.}
#'     \item{\code{"none"}}{Unconstrained fit.}
#'   }
#' @param paired Logical; \code{TRUE} for a paired-samples (within-subjects /
#'   pre-post) design. Default \code{FALSE}.
#' @param type Character string selecting the effect size estimator:
#'   \describe{
#'     \item{\code{"regularized"}}{(Default) \eqn{d_{\mathrm{reg}}} (independent
#'       groups) or \eqn{d_z} (paired), from regularized polynomial quantile
#'       modeling. Recommended for \eqn{n < 50} and/or non-normal data.}
#'     \item{\code{"hedges"}}{Classical bias-corrected Hedges' \eqn{g}, using
#'       the sample pooled standard deviation and the exact \eqn{J(df)}
#'       correction.}
#'     \item{\code{"glass"}}{Glass's \eqn{\Delta}, standardized by the baseline
#'       (Group 1) sample standard deviation only.}
#'     \item{\code{"combined"}}{Hybrid rule: reports \eqn{d_{\mathrm{reg}}} by
#'       default, but switches to Hedges' \eqn{g} whenever both \eqn{n_1, n_2 > 50}
#'       and \eqn{|d_{\mathrm{reg}}| > 0.8}. Note that this switching rule always
#'       evaluates the independent-groups-style \eqn{d_{\mathrm{reg}}}, even in
#'       paired designs; it does not fall back to \eqn{d_z}.}
#'   }
#' @param conf_level Numeric value in \eqn{(0,1)}, e.g. \code{0.95}, or
#'   \code{NULL} (default) to skip confidence interval calculation.
#' @param ci_method Confidence interval method:
#'   \describe{
#'     \item{\code{"bootstrap"}}{(Default) Non-parametric percentile bootstrap,
#'       re-running the entire regularization pipeline on each resample.}
#'     \item{\code{"nct"}}{Analytical inversion of the noncentral \eqn{t}
#'       (lambda-prime) distribution, applied to the regularized point estimate
#'       as an approximation to its sampling distribution.}
#'   }
#' @param B Integer; number of bootstrap resamples when \code{ci_method = "bootstrap"}
#'   (default \code{1000L}).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \subsection{Rationale}{
#' Classical Cohen's \eqn{d} / Hedges' \eqn{g} rely on unregularized sample means
#' and variances. In small samples, or under skewness and heavy tails, sample
#' variances are volatile, inflating the MSE of the effect size. Purely
#' rank-based alternatives (Cliff's \eqn{\delta}, rank-biserial correlation) or
#' Winsorized estimators (\eqn{d_{\mathrm{AKP}}}; Algina et al., 2005) regain
#' stability but change the estimand away from the mean difference on the
#' original, interpretable metric.
#'
#' \eqn{d_{\mathrm{reg}}} avoids both extremes by (1) mapping each group to
#' standard normal scores via inverse-normal rank transformation, (2)
#' approximating each group's empirical quantile function with a low-degree
#' monotone polynomial, and (3) extracting population moments
#' (\eqn{\hat\mu_1, \hat\mu_2, \hat\sigma_1^2, \hat\sigma_2^2}) from these
#' polynomials in closed algebraic form via the orthogonal Hermite basis. For
#' independent groups the effect size is then standardized by the unweighted
#' root-mean-square of the two modeled variances (Cohen, 1988; Delacre et al., 2021):
#' \deqn{d_{\mathrm{reg}} = \frac{\hat\mu_2 - \hat\mu_1}{\hat\sigma_{\mathrm{avg}}},
#'       \qquad \hat\sigma_{\mathrm{avg}} = \sqrt{\frac{\hat\sigma_1^2 + \hat\sigma_2^2}{2}}.}
#' }
#'
#' \subsection{Paired / within-subject designs}{
#' When \code{paired = TRUE}, two complementary standardized metrics are reported:
#' \itemize{
#'   \item \strong{Raw-scale effect size} (\eqn{d_{\mathrm{reg}}}): standardized
#'         by \eqn{\hat\sigma_{\mathrm{avg}}}, directly comparable to
#'         independent-groups effect sizes (e.g. in meta-analysis).
#'   \item \strong{Standardized mean change} (\eqn{d_z}): standardized by the
#'         regularized SD of the observed change scores \eqn{D = X_2 - X_1},
#'         fitted directly as its own quantile model, so that
#'         \eqn{d_z = \hat\mu_D / \hat\sigma_D} with
#'         \eqn{\hat\sigma_D = \sqrt{\hat\sigma_1^2 + \hat\sigma_2^2 - 2 r_{\mathrm{Hermite}}\hat\sigma_1\hat\sigma_2}}.
#' }
#' }
#'
#' @return An S3 object of class \code{"d_reg"} containing:
#' \describe{
#'   \item{\code{estimate}}{The primary point estimate selected by \code{type}
#'     (\eqn{d_{\mathrm{reg}}} for independent groups, \eqn{d_z} for paired
#'     \code{type = "regularized"}).}
#'   \item{\code{d_reg}}{The regularized standardized mean difference on the
#'     raw scale (standardized by \eqn{\hat\sigma_{\mathrm{avg}}}).}
#'   \item{\code{d_z}}{Standardized mean change score (only if \code{paired = TRUE}).}
#'   \item{\code{hedges_g}}{Small-sample bias-corrected Hedges' \eqn{g} (benchmark).}
#'   \item{\code{glass_delta}}{Glass's \eqn{\Delta} using only the Group 1 SD (benchmark).}
#'   \item{\code{type}, \code{paired}, \code{copula}}{Settings used.}
#'   \item{\code{n1}, \code{n2}}{Group sample sizes.}
#'   \item{\code{group1}, \code{group2}}{Regularized and raw moments per group
#'     (\code{mean}, \code{sd}, \code{variance}, \code{raw_mean}, \code{raw_sd}).}
#'   \item{\code{diff_moments}}{Moments of the difference scores (paired only).}
#'   \item{\code{r_Hermite_paired}}{The \code{\link{cor_hermite}} object for the
#'     paired observations (paired only).}
#'   \item{\code{sd_standardizer}}{The standard deviation actually used to
#'     standardize \code{estimate}.}
#'   \item{\code{degrees}}{Realized polynomial degrees per group.}
#'   \item{\code{monotonicity}}{Monotonicity constraint applied.}
#'   \item{\code{fit1}, \code{fit2}}{Underlying \code{\link{hermite_fit}} objects.}
#'   \item{\code{ci}, \code{conf_level}, \code{ci_method}}{Confidence interval
#'     results, if requested.}
#' }
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to Cohen's standardized mean difference effect size: A robust parameter and confidence interval in the two independent groups case. \emph{Psychological Methods}, 10(3), 317-328. \doi{10.1037/1082-989X.10.3.317}
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral Sciences} (2nd ed.). Lawrence Erlbaum Associates.
#'
#' Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021). Why Hedges' gs based on the non-pooled standard deviation should be reported with Welch's t-test. \emph{Preprint}. \doi{10.31234/osf.io/tu6mp}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, 6(2), 107-128. \doi{10.3102/10769986006002107}
#'
#' Lenhard, W., & Lenhard, A. (submitted). Distribution-Free Effect Size Estimation: A Robust Alternative to Cohen's d and other effect size estimators. \emph{Behavior Research Methods}.
#'
#' @seealso \code{\link{t_hermite}} (matched hypothesis test),
#'   \code{\link{hermite_test}}, \code{\link{cor_hermite}},
#'   \code{\link{hermite_fit}}, \code{\link{d_cohen}}, \code{\link{hedges_g}}
#'
#' @examples
#' # 1. Independent groups with non-normal, skewed data
#' set.seed(42)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(30, meanlog = 2.3, sdlog = 0.5)
#' fit1 <- d_reg(ctrl, trt, conf_level = 0.95)
#' print(fit1)
#' summary(fit1)
#'
#' # 2. Formula interface
#' dat <- data.frame(
#'   rt   = c(ctrl, trt),
#'   cond = factor(rep(c("Control", "Treatment"), each = 30))
#' )
#' d_reg(rt ~ cond, data = dat, conf_level = 0.95)
#'
#' # 3. Paired / repeated-measures design
#' pre  <- rlnorm(25, meanlog = 3.0, sdlog = 0.4)
#' post <- pre + rnorm(25, mean = 4.0, sd = 1.5)
#' d_reg(pre, post, paired = TRUE, conf_level = 0.95)
#'
#' # 4. Matched hypothesis test for the same estimand
#' t_hermite(ctrl, trt, nperm = 500)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
d_reg <- function(x, ...) {
  UseMethod("d_reg")
}

#' @rdname d_reg
#' @export
d_reg.formula <- function(x, data = NULL, degree = 3L,
                          copula = c("none", "gaussian"),
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  if (missing(x) || length(x) != 3L) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  mf       <- stats::model.frame(formula = x, data = data)
  response <- mf[[1L]]
  group    <- as.factor(mf[[2L]])
  levels_g <- levels(group)
  if (length(levels_g) != 2L) stop("Grouping variable must have exactly two levels.")

  res <- d_reg.default(x = response[group == levels_g[1L]],
                       y = response[group == levels_g[2L]],
                       degree = degree, copula = copula,
                       monotonicity = monotonicity, paired = paired, type = type,
                       conf_level = conf_level, ci_method = ci_method, B = B, ...)
  res$group_labels <- levels_g
  res
}

#' @rdname d_reg
#' @export
d_reg.default <- function(x, y = NULL, degree = 3L,
                          copula = c("none", "gaussian"),
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  type         <- match.arg(type)
  ci_method    <- match.arg(ci_method)

  if (missing(y) || is.null(y)) stop("Vector 'y' must be supplied.")

  if (paired) {
    if (length(x) != length(y)) stop("Paired observations must have equal length.")
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
  } else {
    x1 <- x[is.finite(x)]
    x2 <- y[is.finite(y)]
  }

  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 3L || n2 < 3L) stop("Each group must contain at least 3 valid observations.")

  # Marginal quantile models and regularized moments
  fit1 <- hermite_fit(x1, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
  fit2 <- hermite_fit(x2, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
  m1 <- hermite_moments(fit1)
  m2 <- hermite_moments(fit2)

  # Averaged standardizer (Cohen, 1988; Delacre et al., 2021)
  sd_avg    <- sqrt((m1$variance + m2$variance) / 2.0)
  d_reg_val <- (m2$mean - m1$mean) / sd_avg

  # Paired: difference-score model and Hermite correlation
  if (paired) {
    fit_diff <- hermite_fit(x2 - x1, degree = degree,
                            monotonicity = monotonicity, force_odd = TRUE)
    m_diff <- hermite_moments(fit_diff)
    r_Hermite_fit <- .cor_hermite_assemble(fit1, fit2,
                                           poly_degree_requested = degree,
                                           copula = copula,
                                           monotonicity = monotonicity,
                                           ties_method = "average")
    d_z <- m_diff$mean / m_diff$sd
  } else {
    fit_diff <- NULL; m_diff <- NULL; r_Hermite_fit <- NULL; d_z <- NULL
  }

  # Classical benchmarks
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sd_pooled <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
  g_hedges    <- ((mean(x2) - mean(x1)) / sd_pooled) * hedges_correction(n1 + n2 - 2L)
  glass_delta <- (mean(x2) - mean(x1)) / s1

  combined_uses_hedges <- (n1 > 50L && n2 > 50L && abs(d_reg_val) > 0.8)

  primary_d <- switch(
    type,
    "regularized" = if (paired) d_z else d_reg_val,
    "hedges"      = g_hedges,
    "glass"       = glass_delta,
    "combined"    = if (combined_uses_hedges) g_hedges else d_reg_val
  )
  sd_standardizer <- switch(
    type,
    "glass"       = s1,
    "hedges"      = sd_pooled,
    "regularized" = sd_avg,
    "combined"    = if (combined_uses_hedges) sd_pooled else sd_avg
  )

  res <- list(
    estimate         = primary_d,
    d_reg            = d_reg_val,
    d_z              = d_z,
    hedges_g         = g_hedges,
    glass_delta      = glass_delta,
    type             = type,
    paired           = paired,
    copula           = copula,
    n1               = n1,
    n2               = n2,
    group1           = list(mean = m1$mean, sd = m1$sd, variance = m1$variance,
                            raw_mean = mean(x1), raw_sd = s1),
    group2           = list(mean = m2$mean, sd = m2$sd, variance = m2$variance,
                            raw_mean = mean(x2), raw_sd = s2),
    diff_moments     = m_diff,
    r_Hermite_paired = r_Hermite_fit,
    sd_standardizer  = sd_standardizer,
    degrees          = c(g1 = fit1$degree, g2 = fit2$degree),
    monotonicity     = monotonicity,
    fit1             = fit1,
    fit2             = fit2,
    x1               = x1,
    x2               = x2
  )
  class(res) <- "d_reg"

  if (!is.null(conf_level)) {
    res$ci <- stats::confint(res, level = conf_level, method = ci_method, B = B)
    res$conf_level <- conf_level
    res$ci_method  <- ci_method
  }
  res
}

# -----------------------------------------------------------------------------
# S3 Methods: confint, print, summary
# -----------------------------------------------------------------------------

#' Confidence Intervals for Distribution-Free Effect Sizes
#'
#' Computes non-parametric percentile bootstrap or analytical noncentral \eqn{t}
#' (lambda-prime) confidence intervals for a fitted \code{"d_reg"} object.
#'
#' @param object An object of class \code{"d_reg"} created by \code{\link{d_reg}}.
#' @param parm Ignored; included for S3 method consistency with \code{\link[stats]{confint}}.
#' @param level Numeric scalar in \eqn{(0, 1)}; the requested confidence level
#'   (default \code{0.95}).
#' @param method Character string specifying the confidence interval method:
#'   \describe{
#'     \item{\code{"bootstrap"}}{(Default) Non-parametric percentile bootstrap,
#'       refitting the entire regularized quantile pipeline on each resample
#'       (preserving the paired structure where applicable).}
#'     \item{\code{"nct"}}{Analytical inversion of the noncentral \eqn{t}
#'       (lambda-prime) distribution, treating the regularized point estimate
#'       as a pseudo-\eqn{t} pivot; see \code{\link{ci_nct}}.}
#'   }
#' @param B Integer scalar; number of bootstrap replications when
#'   \code{method = "bootstrap"} (default \code{1000L}).
#' @param ... Additional arguments passed to internal methods.
#'
#' @details
#' For \code{method = "nct"}, degrees of freedom and effective sample size
#' \eqn{\tilde{n}} are adapted to the design: independent groups use
#' \eqn{\tilde{n} = n_1 n_2 / (n_1 + n_2)} with \eqn{\mathrm{df} = n_1 + n_2 - 2};
#' paired designs and Glass's \eqn{\Delta} use \eqn{\tilde{n} = n_1} with
#' \eqn{\mathrm{df} = n_1 - 1}.
#'
#' @return A numeric matrix of dimension \code{1 x 2} containing the lower and
#'   upper confidence limits, with column names indicating the corresponding
#'   percentiles.
#'
#' @references
#' Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. \emph{Journal of Modern Applied Statistical Methods}, 6(1), 107-116. \doi{10.22237/jmasm/1177992600}
#'
#' @examples
#' set.seed(42)
#' ctrl <- rlnorm(25, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(25, meanlog = 2.3, sdlog = 0.5)
#' fit <- d_reg(ctrl, trt)
#'
#' confint(fit, level = 0.95, method = "bootstrap", B = 500)
#' confint(fit, level = 0.95, method = "nct")
#'
#' @seealso \code{\link{d_reg}}, \code{\link{ci_nct}}
#' @export
confint.d_reg <- function(object, parm, level = 0.95,
                          method = c("bootstrap", "nct"), B = 1000L, ...) {
  method <- match.arg(method)
  n1 <- object$n1; n2 <- object$n2

  if (method == "bootstrap") {
    x1 <- object$x1; x2 <- object$x2
    boot_d <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      fb <- tryCatch({
        if (object$paired) {
          idx <- sample.int(n1, n1, replace = TRUE)
          d_reg(x1[idx], x2[idx], degree = max(object$degrees),
                copula = object$copula, monotonicity = object$monotonicity,
                paired = TRUE, type = object$type)
        } else {
          d_reg(sample(x1, n1, replace = TRUE), sample(x2, n2, replace = TRUE),
                degree = max(object$degrees), copula = object$copula,
                monotonicity = object$monotonicity, paired = FALSE,
                type = object$type)
        }
      }, error = function(e) NULL)
      if (!is.null(fb)) boot_d[b] <- fb$estimate
    }
    alpha <- (1 - level) / 2
    ci <- stats::quantile(boot_d, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  } else {
    if (object$paired || object$type == "glass") {
      df <- n1 - 1L; n_tilde <- n1
    } else {
      df <- n1 + n2 - 2L; n_tilde <- (n1 * n2) / (n1 + n2)
    }
    ci <- ci_nct(object$estimate, n1 = n1, n2 = n2, conf = level,
                 df = df, n_tilde = n_tilde)
  }

  matrix(ci, nrow = 1L,
         dimnames = list("d_reg", c(paste0(100 * (1 - level) / 2, " %"),
                                    paste0(100 * (1 + level) / 2, " %"))))
}

#' @export
print.d_reg <- function(x, digits = 3L, ...) {
  cat("\n  Distribution-Free Effect Size Estimation (d_reg)\n")
  cat(strrep("-", 52), "\n", sep = "")

  if (x$paired) {
    cop_lbl <- if (!is.null(x$copula) && x$copula == "gaussian") {
      " [Gaussian]"
    } else " [Copula-Free]"
    cat(sprintf("  Effect Size (d_reg, raw scale)    :  %.*f\n", digits, x$d_reg))
    cat(sprintf("  Standardized Mean Change (d_z)    :  %.*f\n", digits, x$d_z))
    cat(sprintf("  Paired Hermite Correlation (r)    :  %.*f%s\n",
                digits, x$r_Hermite_paired$r_Hermite, cop_lbl))
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n",
                digits, sqrt((x$group1$variance + x$group2$variance) / 2)))
    cat(sprintf("  Difference Model SD (sigma_diff)  :  %.*f\n",
                digits, x$diff_moments$sd))
  } else {
    cat(sprintf("  Effect Size (d_reg)               :  %.*f\n", digits, x$estimate))
  }

  cat(sprintf("  Standardizer Used (denominator)   :  %.*f\n", digits, x$sd_standardizer))
  cat(sprintf("  Hedges' g (Benchmark)             :  %.*f\n", digits, x$hedges_g))
  cat(sprintf("  Sample Sizes                      :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Polynomial Degrees                :  g1 = %d, g2 = %d\n",
              x$degrees["g1"], x$degrees["g2"]))
  cat(sprintf("  Monotonicity Check                :  %s\n", x$monotonicity))

  if (!is.null(x$ci)) {
    cat(sprintf("  %s CI (%s): [%.*f, %.*f]\n",
                paste0(round(x$conf_level * 100), "%"), x$ci_method,
                digits, x$ci[1L], digits, x$ci[2L]))
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.d_reg <- function(object, digits = 3L, ...) {
  print(object, digits = digits, ...)
  cat("  Group Distributional Moments (Polynomial Modeled):\n")
  cat(sprintf("    Group 1: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$group1$mean, digits, object$group1$sd,
              digits, object$group1$variance))
  cat(sprintf("    Group 2: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$group2$mean, digits, object$group2$sd,
              digits, object$group2$variance))
  cat("\n")
  invisible(object)
}

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
#' @param x A numeric vector of observations for the first group (or baseline),
#'   or a two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param y A numeric vector of observations for the second group (or intervention).
#'   Ignored when \code{x} is a formula.
#' @param data An optional data frame, list, or environment containing the
#'   variables referenced in \code{formula}.
#' @param degree Integer scalar; maximum polynomial degree for the marginal
#'   quantile functions (default \code{3}). Automatically reduced if
#'   monotonicity is violated or too few unique values are available; see
#'   \code{\link{hermite_fit}}.
#' @param monotonicity Monotonicity constraint passed to
#'   \code{\link{hermite_fit}}:
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
#'         fitted directly as its own quantile model (rather than reconstructed
#'         from the marginal moments), so that
#'         \deqn{d_z = \frac{\hat\mu_D}{\hat\sigma_D}}
#'         with \eqn{\hat\sigma_D} empirically consistent with, though not
#'         algebraically forced to equal,
#'         \eqn{\sqrt{\hat\sigma_1^2 + \hat\sigma_2^2 - 2 r_{\mathrm{HM}}\hat\sigma_1\hat\sigma_2}}.
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
#'   \item{\code{hedges_g}}{Small-sample bias-corrected Hedges' \eqn{g} (always reported as a benchmark).}
#'   \item{\code{glass_delta}}{Glass's \eqn{\Delta} using only the Group 1 SD (always reported as a benchmark).}
#'   \item{\code{type}, \code{paired}}{Settings used.}
#'   \item{\code{n1, n2}}{Group sample sizes.}
#'   \item{\code{group1, group2}}{Regularized and raw moments per group (\code{mean}, \code{sd}, \code{variance}, \code{raw_mean}, \code{raw_sd}).}
#'   \item{\code{diff_moments}}{Moments of the difference scores (paired only).}
#'   \item{\code{r_hm_paired}}{The \code{\link{cor_hermite}} object for the paired observations (paired only).}
#'   \item{\code{sd_standardizer}}{The standard deviation actually used to standardize \code{estimate}.}
#'   \item{\code{degrees}}{Realized polynomial degrees per group.}
#'   \item{\code{monotonicity}}{Monotonicity constraint applied.}
#'   \item{\code{fit1, fit2}}{Underlying \code{\link{hermite_fit}} objects.}
#'   \item{\code{ci}, \code{conf_level}, \code{ci_method}}{Confidence interval results, if requested.}
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
#' Isserlis, L. (1918). On a formula for the product-moment coefficient of any order of a normal frequency distribution. \emph{Biometrika}, 12(1/2), 134-139. \doi{10.1093/biomet/12.1-2.134}
#'
#' Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. \emph{Journal of Modern Applied Statistical Methods}, 6(1), 107-116. \doi{10.22237/jmasm/1177992600}
#'
#' Lenhard, W., & Lenhard, A. (submitted). Distribution-Free Effect Size Estimation: A Robust Alternative to Cohen's d and other effect size estimators. \emph{Behavior Research Methods}.
#'
#' @seealso \code{\link{cor_hermite}}, \code{\link{hermite_fit}}, \code{\link{d_cohen}}, \code{\link{hedges_g}}
#'
#' @examples
#' # ---------------------------------------------------------
#' # 1. Independent groups with non-normal, skewed data
#' # ---------------------------------------------------------
#' set.seed(42)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(30, meanlog = 2.3, sdlog = 0.5)
#'
#' fit1 <- d_reg(ctrl, trt, conf_level = 0.95, ci_method = "bootstrap")
#' print(fit1)
#' summary(fit1)
#'
#' # ---------------------------------------------------------
#' # 2. Formula interface
#' # ---------------------------------------------------------
#' dat <- data.frame(
#'   rt = c(ctrl, trt),
#'   cond = factor(rep(c("Control", "Treatment"), each = 30))
#' )
#' fit_form <- d_reg(rt ~ cond, data = dat, conf_level = 0.95)
#' print(fit_form)
#'
#' # ---------------------------------------------------------
#' # 3. Paired / repeated-measures design
#' # ---------------------------------------------------------
#' pre  <- rlnorm(25, meanlog = 3.0, sdlog = 0.4)
#' post <- pre + rnorm(25, mean = 4.0, sd = 1.5)
#'
#' fit_paired <- d_reg(pre, post, paired = TRUE, conf_level = 0.95)
#' print(fit_paired)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
d_reg <- function(x, ...) {
  UseMethod("d_reg")
}

#' @rdname d_reg
#' @export
d_reg.formula <- function(formula, data = NULL, degree = 3L,
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  if (missing(formula) || (length(formula) != 3L)) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  mf <- stats::model.frame(formula = formula, data = data)
  response <- mf[[1L]]
  group <- as.factor(mf[[2L]])

  levels_g <- levels(group)
  if (length(levels_g) != 2L) {
    stop("Grouping variable must have exactly two levels.")
  }

  x1 <- response[group == levels_g[1L]]
  x2 <- response[group == levels_g[2L]]

  res <- d_reg.default(x = x1, y = x2, degree = degree, monotonicity = monotonicity,
                       paired = paired, type = type, conf_level = conf_level,
                       ci_method = ci_method, B = B, ...)
  res$group_labels <- levels_g
  res
}

#' @rdname d_reg
#' @export
d_reg.default <- function(x, y = NULL, degree = 3L,
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  monotonicity <- match.arg(monotonicity)
  type         <- match.arg(type)
  ci_method    <- match.arg(ci_method)

  if (missing(y) || is.null(y)) stop("Vector 'y' must be supplied.")

  if (paired) {
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
    if (length(x1) != length(x2)) stop("Paired observations must have equal length.")
  } else {
    x1 <- x[is.finite(x)]
    x2 <- y[is.finite(y)]
  }

  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 3L || n2 < 3L) stop("Each group must contain at least 3 valid observations.")

  # Fit marginal polynomial quantile distributions
  fit1 <- hermite_fit(x1, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
  fit2 <- hermite_fit(x2, degree = degree, monotonicity = monotonicity, force_odd = TRUE)

  m1 <- hermite_moments(fit1)
  m2 <- hermite_moments(fit2)

  # Averaged standardizer for independent groups (Cohen 1988, Delacre 2021)
  sd_avg <- sqrt((m1$variance + m2$variance) / 2.0)

  # If paired, compute the difference distribution and the latent correlation.
  # r_hm_paired reuses fit1/fit2 (already fitted above) rather than refitting
  # both marginals a second time inside cor_hermite().
  if (paired) {
    diff_scores <- x2 - x1
    fit_diff <- hermite_fit(diff_scores, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
    m_diff <- hermite_moments(fit_diff)
    r_hm_fit <- .cor_hermite_assemble(fit1, fit2, poly_degree_requested = degree,
                                      monotonicity = monotonicity, ties_method = "average")

    # d_z: standardized mean change; d_reg_val: raw-scale repeated-measures SMD
    d_z <- m_diff$mean / m_diff$sd
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  } else {
    fit_diff <- NULL; m_diff <- NULL; r_hm_fit <- NULL; d_z <- NULL
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  }

  # Classical benchmarks (always computed and reported alongside the primary estimate)
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sd_pooled <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
  j_corr <- hedges_correction(n1 + n2 - 2L)
  g_hedges <- ((mean(x2) - mean(x1)) / sd_pooled) * j_corr
  glass_delta <- (mean(x2) - mean(x1)) / s1

  # "combined" always evaluates the switching rule on the independent-groups
  # d_reg (raw scale), even for paired designs; it never falls back to d_z.
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
    estimate = primary_d,
    d_reg = d_reg_val,
    d_z = d_z,
    hedges_g = g_hedges,
    glass_delta = glass_delta,
    type = type,
    paired = paired,
    n1 = n1, n2 = n2,
    group1 = list(mean = m1$mean, sd = m1$sd, variance = m1$variance, raw_mean = mean(x1), raw_sd = s1),
    group2 = list(mean = m2$mean, sd = m2$sd, variance = m2$variance, raw_mean = mean(x2), raw_sd = s2),
    diff_moments = m_diff,
    r_hm_paired = r_hm_fit,
    sd_standardizer = sd_standardizer,
    degrees = c(g1 = fit1$degree, g2 = fit2$degree),
    monotonicity = monotonicity,
    fit1 = fit1, fit2 = fit2,
    x1 = x1, x2 = x2
  )
  class(res) <- "d_reg"

  if (!is.null(conf_level)) {
    res$ci <- stats::confint(res, level = conf_level, method = ci_method, B = B)
    res$conf_level <- conf_level
    res$ci_method <- ci_method
  }
  res
}

# -----------------------------------------------------------------------------
# S3 Methods for d_reg: confint, print, summary
# -----------------------------------------------------------------------------

#' Confidence Intervals for Distribution-Free Effect Sizes
#'
#' @param object An object of class \code{"d_reg"}.
#' @param parm Ignored (present for S3 consistency with \code{\link[stats]{confint}}).
#' @param level Numeric scalar; confidence level (default \code{0.95}).
#' @param method Character; \code{"bootstrap"} (default) or \code{"nct"}.
#' @param B Integer; number of bootstrap replications (default \code{1000L}).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' For \code{method = "nct"}, the analytical noncentral-t inversion
#' (\code{\link{ci_nct}}) requires an effective sample size \eqn{\tilde n}
#' matched to the design:
#' \itemize{
#'   \item \strong{Paired designs} (\code{object$paired = TRUE}): \eqn{\tilde n = n_1}
#'         (the number of pairs) and \code{df = n1 - 1}, treating \eqn{d_z}
#'         as a one-sample/paired-t statistic.
#'   \item \strong{Glass's \eqn{\Delta}} (\code{object$type == "glass"}):
#'         \eqn{\tilde n = n_1} and \code{df = n1 - 1}, consistent with
#'         standardizing by the Group 1 SD alone.
#'   \item \strong{All other independent-groups cases}: the classical
#'         two-sample harmonic form \eqn{\tilde n = n_1 n_2/(n_1+n_2)} and
#'         \code{df = n1 + n2 - 2}.
#' }
#' The \code{"bootstrap"} method makes no such distributional assumption and
#' is recommended as the default for general use, especially for skewed or
#' small-sample data.
#'
#' @return A \code{1 x 2} matrix of confidence limits.
#' @export
confint.d_reg <- function(object, parm, level = 0.95,
                          method = c("bootstrap", "nct"), B = 1000L, ...) {
  method <- match.arg(method)
  n1 <- object$n1; n2 <- object$n2

  if (method == "bootstrap") {
    x1 <- object$x1; x2 <- object$x2
    boot_d <- numeric(B)
    for (b in seq_len(B)) {
      if (object$paired) {
        idx <- sample.int(n1, n1, replace = TRUE)
        fb <- d_reg(x1[idx], x2[idx], degree = max(object$degrees),
                    monotonicity = object$monotonicity, paired = TRUE, type = object$type)
      } else {
        s1 <- sample(x1, n1, replace = TRUE)
        s2 <- sample(x2, n2, replace = TRUE)
        fb <- d_reg(s1, s2, degree = max(object$degrees),
                    monotonicity = object$monotonicity, paired = FALSE, type = object$type)
      }
      boot_d[b] <- fb$estimate
    }
    alpha <- (1 - level) / 2
    ci <- stats::quantile(boot_d, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  } else {
    if (object$paired) {
      df <- n1 - 1L
      n_tilde <- n1
    } else if (object$type == "glass") {
      df <- n1 - 1L
      n_tilde <- n1
    } else {
      df <- n1 + n2 - 2L
      n_tilde <- (n1 * n2) / (n1 + n2)
    }
    ci <- ci_nct(object$estimate, n1 = n1, n2 = n2, conf = level, df = df, n_tilde = n_tilde)
  }

  matrix(ci, nrow = 1L, dimnames = list("d_reg", c(paste0(100 * (1 - level)/2, " %"),
                                                   paste0(100 * (1 + level)/2, " %"))))
}

#' @export
print.d_reg <- function(x, digits = 3L, ...) {
  cat("\n  Distribution-Free Effect Size Estimation (d_reg)\n")
  cat(strrep("-", 52), "\n", sep = "")

  if (x$paired) {
    cat(sprintf("  Effect Size (d_reg, raw scale)    :  %.*f\n", digits, x$d_reg))
    cat(sprintf("  Standardized Mean Change (d_z)    :  %.*f\n", digits, x$d_z))
    cat(sprintf("  Paired Hermite Correlation (r_HM) :  %.*f\n", digits, x$r_hm_paired$r_hm))
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n", digits, sqrt((x$group1$variance + x$group2$variance) / 2)))
    cat(sprintf("  Difference Model SD (sigma_diff)  :  %.*f\n", digits, x$diff_moments$sd))
  } else {
    cat(sprintf("  Effect Size (d_reg)               :  %.*f\n", digits, x$estimate))
  }

  cat(sprintf("  Standardizer Used (denominator)   :  %.*f\n", digits, x$sd_standardizer))
  cat(sprintf("  Hedges' g (Benchmark)             :  %.*f\n", digits, x$hedges_g))
  cat(sprintf("  Sample Sizes                      :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Polynomial Degrees                :  g1 = %d, g2 = %d\n", x$degrees["g1"], x$degrees["g2"]))
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
              digits, object$group1$mean, digits, object$group1$sd, digits, object$group1$variance))
  cat(sprintf("    Group 2: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$group2$mean, digits, object$group2$sd, digits, object$group2$variance))
  cat("\n")
  invisible(object)
}

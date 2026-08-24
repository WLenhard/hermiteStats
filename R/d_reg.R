#' @title Distribution-Free Effect Size Estimation (d_reg)
#'
#' @description
#' Computes the distribution-robust standardized mean difference (\eqn{d_{\mathrm{reg}}})
#' for independent or paired two-group designs. The method retains the classical
#' metric and interpretation of Cohen's \eqn{d} and Hedges' \eqn{g} while
#' achieving substantial variance reduction and lower Mean Squared Error (MSE)
#' across non-normal continuous distributions through regularized polynomial
#' quantile modeling and analytical Hermite moment extraction.
#'
#' @param x A numeric vector of observations for the first group (or baseline),
#'   or a two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param y A numeric vector of observations for the second group (or intervention).
#'   Ignored when \code{x} is a formula.
#' @param data An optional data frame, list, or environment containing the
#'   variables in \code{formula}.
#' @param degree Integer scalar; maximum polynomial degree for the empirical
#'   quantile functions (default is \code{3}). Automatically steps down if
#'   monotonicity is violated or unique values are limited.
#' @param monotonicity Character string specifying the monotonicity constraint:
#'   \describe{
#'     \item{\code{"relaxed"}}{(Default) Empirical rank-concordance criterion
#'       (\eqn{\rho_{\text{Spearman}} \ge 0.95}). Recommended for general empirical applications.}
#'     \item{\code{"strict"}}{Analytical root check ensuring the first derivative
#'       \eqn{f'(z) \ge 0} across the entire standard normal support.}
#'     \item{\code{"none"}}{Unconstrained polynomial fitting.}
#'   }
#' @param paired Logical; if \code{TRUE}, evaluates a paired-samples (within-subjects /
#'   pre-post) design. Default is \code{FALSE}.
#' @param type Character string specifying the effect size estimator:
#'   \describe{
#'     \item{\code{"regularized"}}{(Default) Distribution-free effect size \eqn{d_{\mathrm{reg}}}
#'       using regularized polynomial quantile modeling. Optimal for \eqn{n < 50} and non-normal data.}
#'     \item{\code{"hedges"}}{Traditional small-sample bias-corrected Hedges' \eqn{g}
#'       using sample pooled standard deviation and exact \eqn{J(df)} correction.}
#'     \item{\code{"glass"}}{Glass's \eqn{\Delta} using only the baseline (Group 1) standard deviation.}
#'     \item{\code{"combined"}}{Hybrid estimator: utilizes \eqn{d_{\mathrm{reg}}} by default,
#'       but falls back to Hedges' \eqn{g} when both \eqn{n_1, n_2 > 50} and \eqn{|d_{\mathrm{reg}}| > 0.8}.}
#'   }
#' @param conf_level Numeric value in \eqn{(0, 1)} specifying the confidence level
#'   (e.g., \code{0.95} for a 95\% CI), or \code{NULL} (default) for point estimation only.
#' @param ci_method Character string specifying the confidence interval calculation method:
#'   \describe{
#'     \item{\code{"bootstrap"}}{(Default) Non-parametric percentile bootstrap refitting
#'       the full quantile regularization pipeline in each replication.}
#'     \item{\code{"nct"}}{Analytical inversion of the noncentral \eqn{t} (lambda-prime) distribution.}
#'   }
#' @param B Integer scalar; number of bootstrap resamples if \code{ci_method = "bootstrap"}.
#'   Default is \code{1000L}.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \subsection{Theoretical Foundation}{
#' Traditional standardized mean differences (Cohen's \eqn{d}, Hedges' \eqn{g}) rely on
#' unregularized sample means and variances. In small samples or in the presence of
#' skewness, heavy tails, and outliers, sample variances become volatile, severely
#' inflating the Mean Squared Error (MSE) of the effect size. Conversely, rank-based
#' alternatives (e.g., Cliff's \eqn{\delta}, rank-biserial correlation) or 20\% Winsorization
#' (\eqn{d_{\mathrm{AKP}}}; Algina et al., 2005) alter the target estimand, shifting away
#' from the population mean difference on the observed metric.
#'
#' The \eqn{d_{\mathrm{reg}}} estimator avoids both extremes:
#' \enumerate{
#'   \item \strong{Inverse Normal Quantile Transformation:} Observations within each group
#'         are mapped to standard normal scores \eqn{Z \sim \mathcal{N}(0, 1)}.
#'   \item \strong{Monotone Polynomial Smoothing:} The empirical quantile functions
#'         \eqn{X_1 = f_1(Z_1)} and \eqn{X_2 = f_2(Z_2)} are approximated with low-degree
#'         monotone polynomials via OLS regression.
#'   \item \strong{Analytical Moment Derivation:} Distributional population moments
#'         (\eqn{\hat{\mu}_1, \hat{\mu}_2, \hat{\sigma}_1^2, \hat{\sigma}_2^2}) are derived
#'         in exact closed form using orthogonal Probabilists' Hermite polynomials and
#'         Isserlis' (1918) theorem without numerical integration.
#'   \item \strong{Standardization:} For independent groups, the effect size is standardized
#'         by the unweighted root-mean-square of the population variances:
#'         \deqn{d_{\mathrm{reg}} = \frac{\hat{\mu}_2 - \hat{\mu}_1}{\hat{\sigma}_{\text{avg}}}, \quad \text{where } \hat{\sigma}_{\text{avg}} = \sqrt{\frac{\hat{\sigma}_1^2 + \hat{\sigma}_2^2}{2}}}
#' }
#' }
#'
#' \subsection{Paired / Within-Subject Designs}{
#' When \code{paired = TRUE}, \code{d_reg} provides two complementary standardized metrics:
#' \itemize{
#'   \item \strong{Raw-scale Effect Size (\eqn{d_{\mathrm{reg}}} or \eqn{d_{\mathrm{rm}}}):} Standardized
#'         by the average marginal standard deviation \eqn{\hat{\sigma}_{\text{avg}}}. Retains the
#'         original metric and is directly comparable to independent-group effect sizes in meta-analyses.
#'   \item \strong{Standardized Mean Change (\eqn{d_z}):} Standardized by the regularized standard
#'         deviation of the change scores \eqn{D = X_2 - X_1}:
#'         \deqn{d_z = \frac{\hat{\mu}_D}{\hat{\sigma}_D}}
#'         where \eqn{\hat{\sigma}_D = \sqrt{\hat{\sigma}_1^2 + \hat{\sigma}_2^2 - 2 r_{\mathrm{HM}} \hat{\sigma}_1 \hat{\sigma}_2}},
#'         linking the change score variance directly to the Hermite–Mehler correlation \eqn{r_{\mathrm{HM}}}.
#' }
#' }
#'
#' @return An S3 object of class \code{"d_reg"} containing:
#' \describe{
#'   \item{\code{estimate}}{The primary effect size point estimate (\eqn{d_{\mathrm{reg}}} for independent groups, \eqn{d_z} for paired change if \code{type = "regularized"}).}
#'   \item{\code{d_reg}}{The regularized standardized mean difference based on the average group SD (\eqn{\hat{\sigma}_{\text{avg}}}).}
#'   \item{\code{d_z}}{The standardized mean change score (only present if \code{paired = TRUE}).}
#'   \item{\code{hedges_g}}{The benchmark small-sample bias-corrected Hedges' \eqn{g}.}
#'   \item{\code{glass_delta}}{Glass's \eqn{\Delta} using only Group 1 standard deviation.}
#'   \item{\code{type}}{Character; the estimation method specified.}
#'   \item{\code{paired}}{Logical; indicates whether a paired analysis was conducted.}
#'   \item{\code{n1, n2}}{Sample sizes of Group 1 and Group 2.}
#'   \item{\code{group1, group2}}{Named lists containing regularized and raw moments (\code{mean}, \code{sd}, \code{variance}, \code{raw_mean}, \code{raw_sd}).}
#'   \item{\code{diff_moments}}{Named list containing moments of the difference scores (if \code{paired = TRUE}).}
#'   \item{\code{r_hm_paired}}{The fitted \code{\link{cor_hermite}} object for paired observations (if \code{paired = TRUE}).}
#'   \item{\code{sd_standardizer}}{The standardizer value used in the denominator.}
#'   \item{\code{degrees}}{Realized polynomial degrees for each group.}
#'   \item{\code{monotonicity}}{The monotonicity constraint applied.}
#'   \item{\code{fit1, fit2}}{Underlying \code{\link{hermite_fit}} objects for each group.}
#'   \item{\code{ci}}{Matrix containing confidence limits (if \code{conf_level} was specified).}
#'   \item{\code{conf_level}}{Requested confidence level.}
#'   \item{\code{ci_method}}{Method used for confidence interval estimation.}
#' }
#'
#' @references
#' Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to Cohen’s standardized mean difference effect size: A robust parameter and confidence interval in the two independent groups case. \emph{Psychological Methods}, 10(3), 317–328. \doi{10.1037/1082-989X.10.3.317}
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral Sciences} (2nd ed.). Lawrence Erlbaum Associates.
#'
#' Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021). Why Hedges' gs based on the non-pooled standard deviation should be reported with Welch's t-test. \emph{Preprint}. \doi{10.31234/osf.io/tu6mp}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass’s estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, 6(2), 107–128. \doi{10.3102/10769986006002107}
#'
#' Isserlis, L. (1918). On a formula for the product-moment coefficient of any order of a normal frequency distribution. \emph{Biometrika}, 12(1/2), 134–139. \doi{10.1093/biomet/12.1-2.134}
#'
#' Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. \emph{Journal of Modern Applied Statistical Methods}, 6(1), 107–116. \doi{10.22237/jmasm/1177992600}
#'
#' Lenhard, W., & Lenhard, A. (submitted). Distribution-Free Effect Size Estimation: A Robust Alternative to Cohen’s d and other effect size estimators. \emph{Behavior Research Methods}.
#'
#' @seealso \code{\link{cor_hermite}}, \code{\link{hermite_fit}}, \code{\link{d_cohen}}, \code{\link{hedges_g}}
#'
#' @examples
#' # ---------------------------------------------------------
#' # 1. Independent Groups with Non-Normal Skewed Data
#' # ---------------------------------------------------------
#' set.seed(42)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(30, meanlog = 2.3, sdlog = 0.5)
#'
#' # Default regularized effect size with 95% bootstrap CI
#' fit1 <- d_reg(ctrl, trt, conf_level = 0.95, ci_method = "bootstrap")
#' print(fit1)
#' summary(fit1)
#'
#' # ---------------------------------------------------------
#' # 2. Formula Interface
#' # ---------------------------------------------------------
#' dat <- data.frame(
#'   rt = c(ctrl, trt),
#'   cond = factor(rep(c("Control", "Treatment"), each = 30))
#' )
#' fit_form <- d_reg(rt ~ cond, data = dat, conf_level = 0.95)
#' print(fit_form)
#'
#' # ---------------------------------------------------------
#' # 3. Paired / Repeated-Measures Design
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

  # If paired, compute difference distribution and latent correlation
  if (paired) {
    diff_scores <- x2 - x1
    fit_diff <- hermite_fit(diff_scores, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
    m_diff <- hermite_moments(fit_diff)
    r_hm_fit <- cor_hermite(x1, x2, poly_degree = degree, monotonicity = monotonicity)

    # d_z: standardized mean change; d_reg_val: raw-scale repeated measures SMD
    d_z <- m_diff$mean / m_diff$sd
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  } else {
    fit_diff <- NULL; m_diff <- NULL; r_hm_fit <- NULL; d_z <- NULL
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  }

  # Classical benchmarks
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sd_pooled <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
  j_corr <- .J_correction(n1 + n2 - 2L)
  g_hedges <- ((mean(x2) - mean(x1)) / sd_pooled) * j_corr
  glass_delta <- (mean(x2) - mean(x1)) / s1

  # Type routing
  primary_d <- switch(
    type,
    "regularized" = if (paired) d_z else d_reg_val,
    "hedges"      = g_hedges,
    "glass"       = glass_delta,
    "combined"    = if (n1 > 50L && n2 > 50L && abs(d_reg_val) > 0.8) g_hedges else d_reg_val
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
    sd_standardizer = if (type == "glass") s1 else if (type == "hedges") sd_pooled else sd_avg,
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
#' Computes analytical noncentral \eqn{t} (lambda-prime) or non-parametric
#' percentile bootstrap confidence intervals for \code{d_reg} objects.
#'
#' @param object An object of class \code{"d_reg"}.
#' @param parm Ignored.
#' @param level Numeric scalar; confidence level (default is \code{0.95}).
#' @param method Character; \code{"bootstrap"} (default) or \code{"nct"}.
#' @param B Integer; number of bootstrap replications (default is \code{1000L}).
#' @param ... Additional arguments.
#'
#' @return A matrix of dimension \code{1 x 2} with confidence limits.
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
    df <- if (object$paired) n1 - 1L else (if (object$type == "glass") n1 - 1L else n1 + n2 - 2L)
    ci <- ci_nct(object$estimate, n1 = n1, n2 = n2, conf = level, df = df)
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
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n", digits, x$sd_standardizer))
    cat(sprintf("  Difference Model SD (sigma_diff)  :  %.*f\n", digits, x$diff_moments$sd))
  } else {
    cat(sprintf("  Effect Size (d_reg)               :  %.*f\n", digits, x$estimate))
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n", digits, x$sd_standardizer))
  }

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

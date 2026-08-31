# =============================================================================
# tests.R -- Distribution-Robust Hypothesis Testing via Hermite Quantile Models
## =============================================================================

#' @title Distribution-Robust Hypothesis Testing (Location, Scale, and Shape)
#'
#' @description
#' Permutation-based (default) and approximate analytical hypothesis tests for
#' two-sample and paired designs, all built on the same regularized polynomial
#' quantile engine used throughout \pkg{hermiteStats}:
#'
#' \describe{
#'   \item{\code{\link{t_hermite}}}{Mean difference test (regularized analogue
#'     of Welch's \eqn{t}-test).}
#'   \item{\code{\link{median_hermite}}}{Median difference test based on the
#'     regularized median \eqn{\widehat{\mathrm{Med}} = f(0) = \hat\beta_0}.}
#'   \item{\code{\link{shape_hermite}}}{Shape difference test for scale
#'     (\eqn{\Delta\log\sigma}), asymmetry (\eqn{\Delta c_2}), and tail weight
#'     (\eqn{\Delta c_3}), with Westfall-Young multiplicity control and an
#'     omnibus test.}
#'   \item{\code{\link{hermite_test}}}{Unified interface dispatching to the
#'     tests above, including a joint five-dimensional "complete" profile test.}
#' }
#'
#' All test statistics are functionals of a monotone polynomial quantile model
#' \eqn{X = f(Z) = \sum_j \beta_j Z^j} re-expressed in the orthogonal
#' Probabilists' Hermite basis \eqn{f(Z) = \sum_m a_m He_m(Z)}, from which
#' \deqn{\mu = a_0, \quad \mathrm{Med} = f(0) = \beta_0, \quad
#'       \sigma^2 = \sum_{m \ge 1} a_m^2\, m!, \quad
#'       c_m = a_m \sqrt{m!} / \sigma}
#' follow in closed form. Because the exact finite-sample distribution of these
#' regularized functionals has no closed form, permutation inference (refitting
#' the entire pipeline on every resample) is the recommended default.
#'
#' @name hermite_tests
#' @keywords internal
NULL

# =============================================================================
# 0. Internal Engine
# =============================================================================

# --- 0.1 Basis cache ---------------------------------------------------------

#' Environment caching the fixed QR decomposition and Hermite operator per
#' (sample size, degree) pair. For tie-free samples, the normal-score design
#' matrix depends only on n and degree, so each permutation replicate reduces
#' to a single qr.coef() call on the sorted values. The cache is size-capped.
#' @noRd
.hermite_test_cache <- new.env(parent = emptyenv())

#' Cached normal-score design (QR) and Hermite basis operator
#'
#' @param n Integer; sample size (tie-free case).
#' @param degree Integer; polynomial degree.
#' @return List with components \code{QR} (QR decomposition of the Vandermonde
#'   matrix on the Hazen normal-score grid), \code{Hmat} (monomial-to-Hermite
#'   change-of-basis operator), and \code{fac} (\code{factorial(0:degree)}).
#' @noRd
.hermite_test_basis <- function(n, degree) {
  key <- sprintf("n%d_d%d", n, degree)
  obj <- .hermite_test_cache[[key]]
  if (!is.null(obj)) return(obj)

  if (length(ls(envir = .hermite_test_cache)) >= 128L) {
    rm(list = ls(envir = .hermite_test_cache), envir = .hermite_test_cache)
  }

  z    <- stats::qnorm((seq_len(n) - 0.5) / n)
  Zmat <- outer(z, 0:degree, `^`)
  obj  <- list(
    QR   = qr(Zmat),
    Hmat = .hermite_basis_matrix(degree),
    fac  = factorial(0:degree)
  )
  assign(key, obj, envir = .hermite_test_cache)
  obj
}

# --- 0.2 Unified profile extractor -------------------------------------------

#' Full regularized distributional profile of a sample
#'
#' Fits the polynomial quantile model X = f(Z) by OLS on rank-based normal
#' scores and extracts every functional required by the hypothesis tests in a
#' single pass:
#' \itemize{
#'   \item \code{mean} = a_0, with asymptotic SE = sigma / sqrt(n);
#'   \item \code{median} = f(0) = beta_0, with asymptotic SE
#'         |beta_1| * sqrt(pi / (2n)) (density-based approximation via
#'         1/f(med) = sqrt(2*pi) * f'(0); approximate for the regularized
#'         median -- one reason permutation inference is the default);
#'   \item \code{variance} / \code{sd} / \code{log_sd} via Parseval's identity;
#'   \item \code{c2}, \code{c3}: Parseval-standardized Hermite shape weights
#'         (asymmetry and tail weight), bounded in [-1, 1] and location-scale
#'         invariant.
#' }
#' Tie-free samples use the cached QR fast path; tied samples fall back to an
#' explicit rank-based fit.
#'
#' @param v Numeric vector.
#' @param degree Integer; polynomial degree (odd).
#' @return A list of profile statistics, or \code{NULL} if the sample is too
#'   small or the fit is degenerate.
#' @noRd
.hermite_profile <- function(v, degree = 3L) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < degree + 2L) return(NULL)

  vs <- sort(v)
  if (anyDuplicated(vs) == 0L) {
    bas  <- .hermite_test_basis(n, degree)
    beta <- tryCatch(qr.coef(bas$QR, vs), error = function(e) NULL)
    Hmat <- bas$Hmat
    fac  <- bas$fac
  } else {
    z    <- stats::qnorm((rank(v, ties.method = "average") - 0.5) / n)
    Zmat <- outer(z, 0:degree, `^`)
    beta <- tryCatch(qr.solve(Zmat, v), error = function(e) NULL)
    Hmat <- .hermite_basis_matrix(degree)
    fac  <- factorial(0:degree)
  }
  if (is.null(beta) || anyNA(beta)) return(NULL)

  a <- as.vector(Hmat %*% beta)
  m <- seq_len(degree)

  var_val <- sum(a[m + 1L]^2 * fac[m + 1L])
  if (!is.finite(var_val) || var_val <= .Machine$double.eps) return(NULL)
  sdv <- sqrt(var_val)
  cm  <- a[m + 1L] * sqrt(fac[m + 1L]) / sdv

  list(
    n        = n,
    beta     = beta,
    hermite  = a,
    mean     = a[1L],
    se_mean  = sdv / sqrt(n),
    median   = beta[1L],
    se_med   = abs(beta[2L]) * sqrt(pi / (2.0 * n)),
    sd       = sdv,
    variance = var_val,
    log_sd   = log(sdv),
    c2       = if (degree >= 2L) cm[2L] else 0.0,
    c3       = if (degree >= 3L) cm[3L] else 0.0,
    df       = max(1L, n - degree - 1L)
  )
}

#' Shape contrast vector (Group 2 - Group 1) from two profiles
#' @noRd
.shape_contrast <- function(p1, p2) {
  c(scale      = p2$log_sd - p1$log_sd,
    asymmetry  = p2$c2     - p1$c2,
    tailweight = p2$c3     - p1$c3)
}

# --- 0.3 Multiplicity, omnibus, formatting ------------------------------------

#' Westfall-Young step-down min-P adjustment from a joint permutation matrix
#'
#' @param Tabs Numeric matrix of *absolute* standardized statistics with the
#'   observed statistics in row 1 and permutation replicates below; one column
#'   per contrast.
#' @return List with \code{p_raw} (marginal permutation p-values, observed row
#'   included in the reference set, i.e. the add-one convention), \code{p_adj}
#'   (step-down adjusted p-values, monotonicity enforced), and \code{P} (the
#'   full matrix of marginal p-values used for minP omnibus tests).
#' @references Westfall, P. H., & Young, S. S. (1993). \emph{Resampling-Based
#'   Multiple Testing}. Wiley.
#' @noRd
.wy_stepdown <- function(Tabs) {
  M <- nrow(Tabs); K <- ncol(Tabs)
  P <- matrix(NA_real_, M, K)
  for (k in seq_len(K)) {
    r <- rank(Tabs[, k], ties.method = "min")
    P[, k] <- (M - r + 1) / M
  }
  p_raw <- P[1L, ]
  if (K == 1L) return(list(p_raw = p_raw, p_adj = p_raw, P = P))

  o <- order(p_raw)
  Q <- P[, o, drop = FALSE]
  for (j in seq(K - 1L, 1L)) Q[, j] <- pmin(Q[, j], Q[, j + 1L])
  adj <- vapply(seq_len(K), function(j) mean(Q[, j] <= Q[1L, j] + 1e-12), numeric(1L))
  adj <- cummax(adj)
  out <- numeric(K); out[o] <- adj
  list(p_raw = p_raw, p_adj = out, P = P)
}

#' Omnibus test from a standardized joint permutation matrix
#'
#' @param Zall Standardized statistics ((B+1) x K), observed in row 1.
#' @param P Marginal p-value matrix from \code{.wy_stepdown}.
#' @param type One of "minP" (Tippett minimum-p), "maxT" (maximum absolute
#'   standardized statistic), or "T2" (Mahalanobis / permutation Hotelling
#'   type, using the permutation covariance of the standardized contrasts).
#' @return List with \code{type}, \code{statistic}, \code{p_value}, or
#'   \code{NULL} if \code{type == "none"} or K < 2.
#' @noRd
.omnibus_test <- function(Zall, P, type) {
  if (type == "none" || ncol(Zall) < 2L) return(NULL)
  if (type == "minP") {
    st <- apply(P, 1L, min)
    return(list(type = "minP", statistic = st[1L],
                p_value = mean(st <= st[1L] + 1e-12)))
  }
  if (type == "maxT") {
    st <- apply(abs(Zall), 1L, max)
    return(list(type = "maxT", statistic = st[1L],
                p_value = mean(st >= st[1L] - 1e-10)))
  }
  # T2: Mahalanobis distance under the permutation covariance
  S <- stats::cov(Zall[-1L, , drop = FALSE])
  Sinv <- tryCatch(solve(S), error = function(e) {
    solve(S + diag(1e-8, ncol(S)))
  })
  st <- rowSums((Zall %*% Sinv) * Zall)
  list(type = "T2", statistic = st[1L],
       p_value = mean(st >= st[1L] - 1e-10))
}

#' Format p-values for print methods
#' @noRd
.p_fmt <- function(p, digits = 3L) {
  vapply(p, function(v) {
    if (is.na(v)) "NA"
    else if (v < 10^(-digits)) sprintf("< %.*f", digits, 10^(-digits))
    else sprintf("%.*f", digits, v)
  }, character(1L))
}

#' Two-group input validation and cleaning shared by all tests
#' @noRd
.two_sample_data <- function(x, y, paired, degree, min_extra = 2L) {
  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (paired && length(x) != length(y)) stop("Paired vectors must have equal length.")

  if (paired) {
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
  } else {
    x1 <- x[is.finite(x)]
    x2 <- y[is.finite(y)]
  }
  if (min(length(x1), length(x2)) < degree + min_extra) {
    stop(sprintf("Each group needs at least %d valid observations for degree %d.",
                 degree + min_extra, degree))
  }
  list(x1 = x1, x2 = x2)
}

#' Extract the two groups from a response ~ group formula
#' @noRd
.formula_groups <- function(formula, data) {
  if (missing(formula) || length(formula) != 3L) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  mf <- stats::model.frame(formula = formula, data = data)
  response <- mf[[1L]]
  group    <- as.factor(mf[[2L]])
  levels_g <- levels(group)
  if (length(levels_g) != 2L) stop("Grouping variable must have exactly two levels.")
  list(x = response[group == levels_g[1L]],
       y = response[group == levels_g[2L]],
       labels = levels_g)
}

# --- 0.4 Shared location-test engine ------------------------------------------

#' Internal engine shared by t_hermite (statistic = "mean") and
#' median_hermite (statistic = "median")
#'
#' Both tests are studentized location tests: T = (theta_2 - theta_1) / SE,
#' where theta is the regularized mean (a_0) or regularized median (beta_0).
#'
#' Permutation null (default): each group is aligned by subtracting its own
#' regularized location estimate; the pooled residuals are exchangeable under
#' H0: theta_1 = theta_2, and the studentized statistic is recomputed on every
#' resample (refitting the full quantile pipeline). Paired designs use
#' sign-flipping of the difference scores.
#'
#' Analytical mode: Student-t reference with Welch-Satterthwaite df based on
#' the asymptotic SEs. Approximate; see the exported documentation.
#'
#' @noRd
.hermite_location_test <- function(x1, x2, paired, statistic, method, degree,
                                   alternative, nperm) {

  n1 <- length(x1); n2 <- length(x2)

  get_theta <- function(p) if (statistic == "mean") p$mean else p$median
  get_se    <- function(p) if (statistic == "mean") p$se_mean else p$se_med

  f1 <- .hermite_profile(x1, degree)
  f2 <- .hermite_profile(x2, degree)
  if (is.null(f1) || is.null(f2)) stop("Degenerate quantile fit in at least one group.")

  if (paired) {
    diffs  <- x2 - x1
    f_diff <- .hermite_profile(diffs, degree)
    if (is.null(f_diff)) stop("Degenerate quantile fit for the difference scores.")
    est     <- get_theta(f_diff)
    se_diff <- get_se(f_diff)
    df_val  <- f_diff$df
  } else {
    est     <- get_theta(f2) - get_theta(f1)
    se1 <- get_se(f1); se2 <- get_se(f2)
    se_diff <- sqrt(se1^2 + se2^2)
    df_val  <- (se1^2 + se2^2)^2 / (se1^4 / f1$df + se2^4 / f2$df)
  }

  sd_avg <- sqrt((f1$variance + f2$variance) / 2.0)
  t_stat <- est / se_diff

  if (method == "analytical") {
    p_val <- switch(
      alternative,
      "two.sided" = 2.0 * stats::pt(-abs(t_stat), df = df_val),
      "greater"   = stats::pt(t_stat, df = df_val, lower.tail = FALSE),
      "less"      = stats::pt(t_stat, df = df_val, lower.tail = TRUE)
    )
    perm_dist <- numeric(0)
    n_success <- 0L
    success_rate <- NA_real_

  } else {
    perm_t <- rep(NA_real_, nperm)

    if (paired) {
      for (b in seq_len(nperm)) {
        sgn <- ifelse(stats::runif(n1) < 0.5, -1.0, 1.0)
        fb  <- .hermite_profile(diffs * sgn, degree)
        if (!is.null(fb)) perm_t[b] <- get_theta(fb) / get_se(fb)
      }
    } else {
      pool <- c(x1 - get_theta(f1), x2 - get_theta(f2))
      N <- n1 + n2
      for (b in seq_len(nperm)) {
        idx <- sample.int(N, n1, replace = FALSE)
        pa  <- .hermite_profile(pool[idx],  degree)
        pb  <- .hermite_profile(pool[-idx], degree)
        if (!is.null(pa) && !is.null(pb)) {
          perm_t[b] <- (get_theta(pb) - get_theta(pa)) /
            sqrt(get_se(pa)^2 + get_se(pb)^2)
        }
      }
    }

    perm_dist    <- perm_t[is.finite(perm_t)]
    n_success    <- length(perm_dist)
    success_rate <- n_success / nperm
    if (n_success < 10L) stop("Too few usable permutation replicates.")

    tol <- 1e-10
    p_val <- switch(
      alternative,
      "two.sided" = (1.0 + sum(abs(perm_dist) >= abs(t_stat) - tol)) / (n_success + 1.0),
      "greater"   = (1.0 + sum(perm_dist >= t_stat - tol)) / (n_success + 1.0),
      "less"      = (1.0 + sum(perm_dist <= t_stat + tol)) / (n_success + 1.0)
    )
  }

  list(
    statistic         = t_stat,
    p_value           = p_val,
    estimate          = est,
    d_standardized    = est / sd_avg,
    se_diff           = se_diff,
    df                = df_val,
    loc1              = get_theta(f1),
    loc2              = get_theta(f2),
    mean1             = f1$mean,   mean2 = f2$mean,
    median1           = f1$median, median2 = f2$median,
    sd1               = f1$sd,     sd2 = f2$sd,
    c2_1              = f1$c2,     c2_2 = f2$c2,
    sd_avg            = sd_avg,
    n1                = n1, n2 = n2,
    degree            = degree,
    paired            = paired,
    statistic_type    = statistic,
    method_type       = method,
    alternative       = alternative,
    nperm             = if (method == "permutation") nperm else NA_integer_,
    n_perm_success    = n_success,
    success_rate      = success_rate,
    perm_distribution = perm_dist
  )
}

# =============================================================================
# 1. t_hermite: Mean Difference Test
# =============================================================================

#' Distribution-Robust Mean Difference Test (Hermite t-Test)
#'
#' Tests the null hypothesis of equal population means,
#' \eqn{H_0: \mu_2 = \mu_1}, for two independent groups or paired
#' measurements, using the regularized mean \eqn{\hat\mu = \hat a_0} of the
#' fitted polynomial quantile model in place of the raw sample mean.
#'
#' @details
#' \subsection{Statistic}{
#' The test statistic is the studentized regularized mean difference
#' \deqn{t_{\mathrm{Hermite}} = \frac{\hat\mu_2 - \hat\mu_1}
#'       {\sqrt{\hat\sigma_1^2/n_1 + \hat\sigma_2^2/n_2}},}
#' where \eqn{\hat\mu = \hat a_0} and
#' \eqn{\hat\sigma^2 = \sum_{m \ge 1} \hat a_m^2\, m!} are the closed-form
#' Hermite moments (Parseval's identity). For paired designs the statistic is
#' \eqn{\hat\mu_D / (\hat\sigma_D / \sqrt{n})}, where the difference scores
#' \eqn{D = X_2 - X_1} receive their own quantile fit.
#' }
#'
#' \subsection{Inference (\code{method})}{
#' \describe{
#'   \item{\code{"permutation"} (default)}{Each group is aligned by
#'     subtracting its regularized mean; the pooled residuals are permuted
#'     across groups (or difference-score signs are flipped, if paired), and
#'     the full quantile pipeline is refitted on every resample. This yields
#'     well-calibrated Type I error across normal, skewed, and heavy-tailed
#'     distributions, and is the recommended mode.}
#'   \item{\code{"analytical"}}{Student-\eqn{t} reference with
#'     Welch-Satterthwaite degrees of freedom based on the asymptotic standard
#'     errors \eqn{\hat\sigma/\sqrt{n}}. Fast, but only approximate in small
#'     or strongly non-normal samples; prefer permutation for reported
#'     results.}
#' }
#' }
#'
#' @param x A numeric vector of observations for Group 1 (or the baseline in
#'   paired designs), or a two-sided \code{\link[stats]{formula}} of the form
#'   \code{response ~ group}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"t_hermite"} containing, among others:
#' \describe{
#'   \item{\code{statistic}}{Observed studentized statistic \eqn{t_{\mathrm{Hermite}}}.}
#'   \item{\code{p_value}}{Permutation or analytical \eqn{p}-value.}
#'   \item{\code{mean_diff}}{Regularized mean difference (Group 2 \eqn{-} Group 1,
#'     or mean of the difference scores if paired).}
#'   \item{\code{d_reg}}{Standardized effect size (mean difference divided by
#'     \eqn{\hat\sigma_{\mathrm{avg}}}).}
#'   \item{\code{mean1}, \code{mean2}, \code{sd1}, \code{sd2}}{Regularized group moments.}
#'   \item{\code{se_diff}, \code{df}}{Standard error and (Welch) degrees of freedom.}
#'   \item{\code{perm_distribution}}{Permutation null statistics (permutation mode).}
#'   \item{\code{n_perm_success}, \code{success_rate}}{Permutation bookkeeping.}
#' }
#'
#' @seealso \code{\link{median_hermite}}, \code{\link{shape_hermite}},
#'   \code{\link{hermite_test}}, \code{\link{d_reg}}
#'
#' @examples
#' set.seed(42)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(30, meanlog = 2.4, sdlog = 0.5)
#'
#' # Permutation test (default)
#' res <- t_hermite(ctrl, trt, nperm = 500)
#' print(res)
#' plot(res)
#'
#' # Fast approximate analytical test
#' t_hermite(ctrl, trt, method = "analytical")
#'
#' # Formula interface
#' df <- data.frame(val = c(ctrl, trt), grp = factor(rep(c("A", "B"), each = 30)))
#' t_hermite(val ~ grp, data = df, nperm = 500)
#'
#' # Paired design
#' pre  <- rlnorm(25, meanlog = 3.0, sdlog = 0.4)
#' post <- pre + rnorm(25, mean = 4.0, sd = 1.5)
#' t_hermite(pre, post, paired = TRUE, nperm = 500)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
t_hermite <- function(x, ...) {
  UseMethod("t_hermite")
}

#' @rdname t_hermite
#' @param formula A two-sided formula of the form \code{response ~ group}.
#' @param data An optional data frame containing the variables in \code{formula}.
#' @export
t_hermite.formula <- function(formula, data = NULL, paired = FALSE,
                              method = c("permutation", "analytical"),
                              degree = 3L,
                              alternative = c("two.sided", "greater", "less"),
                              nperm = 1000L, ...) {
  fg  <- .formula_groups(formula, data)
  res <- t_hermite.default(x = fg$x, y = fg$y, paired = paired, method = method,
                           degree = degree, alternative = alternative,
                           nperm = nperm, ...)
  res$group_labels <- fg$labels
  res
}

#' @rdname t_hermite
#' @param y Numeric vector of observations for Group 2 (or the post-measurement
#'   in paired designs).
#' @param paired Logical; \code{TRUE} for a paired-samples design. Default \code{FALSE}.
#' @param method Character string; \code{"permutation"} (default) or
#'   \code{"analytical"} (approximate closed form). See Details.
#' @param degree Integer scalar; polynomial degree for the quantile models
#'   (default \code{3L}; must be an odd integer \eqn{\ge 1}).
#' @param alternative Alternative hypothesis: \code{"two.sided"} (default),
#'   \code{"greater"} (\eqn{\mu_2 > \mu_1}), or \code{"less"} (\eqn{\mu_2 < \mu_1}).
#' @param nperm Integer; number of permutation resamples (default \code{1000L}).
#' @export
t_hermite.default <- function(x, y, paired = FALSE,
                              method = c("permutation", "analytical"),
                              degree = 3L,
                              alternative = c("two.sided", "greater", "less"),
                              nperm = 1000L, ...) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)
  degree      <- as.integer(degree)
  if (degree < 1L || degree %% 2L == 0L) stop("'degree' must be an odd integer >= 1.")

  d   <- .two_sample_data(x, y, paired, degree)
  res <- .hermite_location_test(d$x1, d$x2, paired = paired, statistic = "mean",
                                method = method, degree = degree,
                                alternative = alternative, nperm = nperm)
  res$mean_diff <- res$estimate
  res$d_reg     <- res$d_standardized
  res$method    <- paste0(if (paired) "Paired" else "Two-Sample",
                          " Hermite Mean Difference Test")
  class(res) <- "t_hermite"
  res
}

# =============================================================================
# 2. median_hermite: Median Difference Test
# =============================================================================

#' Distribution-Robust Median Difference Test (Hermite Median Test)
#'
#' Tests the null hypothesis of equal population medians,
#' \eqn{H_0: \mathrm{Med}_2 = \mathrm{Med}_1}, for two independent groups or
#' paired measurements. The regularized median is obtained analytically as
#' \eqn{\widehat{\mathrm{Med}} = f(0) = \hat\beta_0}, i.e. the fitted quantile
#' polynomial evaluated at \eqn{Z = \Phi^{-1}(0.5) = 0}.
#'
#' @details
#' \subsection{Why the regularized median?}{
#' The empirical sample median is a step function of one or two middle order
#' statistics, and naive permutation tests on raw medians are severely
#' miscalibrated (Type I error of 9-12\% in simulations). The regularized
#' median \eqn{f(0) = \hat\beta_0} is a smooth linear functional of all order
#' statistics via OLS quantile smoothing, which restores exact calibration
#' under permutation while retaining the median's robustness against tail
#' variance. On skewed and heavy-tailed data, the permutation Hermite median
#' test achieves higher power than Welch's \eqn{t}-test while, unlike the
#' Wilcoxon test, remaining a genuine on-metric test of the population median.
#' }
#'
#' \subsection{Inference (\code{method})}{
#' \describe{
#'   \item{\code{"permutation"} (default)}{Groups are aligned by their
#'     regularized medians, pooled residuals are permuted (or difference-score
#'     signs flipped, if paired), and the studentized statistic is recomputed
#'     on every resample. Recommended: near-exact Type I error control across
#'     distribution families.}
#'   \item{\code{"analytical"}}{Student-\eqn{t} reference using the asymptotic
#'     density-based standard error
#'     \eqn{\widehat{\mathrm{SE}} = |\hat\beta_1| \sqrt{\pi / (2n)}}
#'     (from \eqn{1/f(\mathrm{med}) = \sqrt{2\pi}\,\hat\beta_1}). This is a
#'     large-sample approximation whose calibration degrades under heavy tails;
#'     prefer permutation for reported results.}
#' }
#' }
#'
#' @param x A numeric vector of observations for Group 1 (or the baseline in
#'   paired designs), or a two-sided \code{\link[stats]{formula}} of the form
#'   \code{response ~ group}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"median_hermite"} containing, among others:
#' \describe{
#'   \item{\code{statistic}}{Observed studentized statistic.}
#'   \item{\code{p_value}}{Permutation or analytical \eqn{p}-value.}
#'   \item{\code{median_diff}}{Regularized median difference (Group 2 \eqn{-}
#'     Group 1, or the median of the difference scores if paired).}
#'   \item{\code{d_med}}{Median difference standardized by \eqn{\hat\sigma_{\mathrm{avg}}}.}
#'   \item{\code{median1}, \code{median2}, \code{mean1}, \code{mean2},
#'     \code{c2_1}, \code{c2_2}}{Group profiles; the \eqn{c_2} asymmetry
#'     indices contextualize any mean-median divergence.}
#'   \item{\code{perm_distribution}}{Permutation null statistics (permutation mode).}
#' }
#'
#' @seealso \code{\link{t_hermite}}, \code{\link{shape_hermite}}, \code{\link{hermite_test}}
#'
#' @examples
#' set.seed(42)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.6)
#' trt  <- rlnorm(30, meanlog = 2.4, sdlog = 0.6)
#'
#' res <- median_hermite(ctrl, trt, nperm = 500)
#' print(res)
#' plot(res)
#'
#' # Formula interface
#' df <- data.frame(val = c(ctrl, trt), grp = factor(rep(c("Ctrl", "Trt"), each = 30)))
#' median_hermite(val ~ grp, data = df, nperm = 500)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
median_hermite <- function(x, ...) {
  UseMethod("median_hermite")
}

#' @rdname median_hermite
#' @param formula A two-sided formula of the form \code{response ~ group}.
#' @param data An optional data frame containing the variables in \code{formula}.
#' @export
median_hermite.formula <- function(formula, data = NULL, paired = FALSE,
                                   method = c("permutation", "analytical"),
                                   degree = 3L,
                                   alternative = c("two.sided", "greater", "less"),
                                   nperm = 1000L, ...) {
  fg  <- .formula_groups(formula, data)
  res <- median_hermite.default(x = fg$x, y = fg$y, paired = paired,
                                method = method, degree = degree,
                                alternative = alternative, nperm = nperm, ...)
  res$group_labels <- fg$labels
  res
}

#' @rdname median_hermite
#' @param y Numeric vector of observations for Group 2.
#' @param paired Logical; \code{TRUE} for paired samples. Default \code{FALSE}.
#' @param method Character string; \code{"permutation"} (default) or
#'   \code{"analytical"} (approximate closed form). See Details.
#' @param degree Integer scalar; polynomial degree (default \code{3L}; odd, \eqn{\ge 1}).
#' @param alternative Alternative hypothesis: \code{"two.sided"} (default),
#'   \code{"greater"} (\eqn{\mathrm{Med}_2 > \mathrm{Med}_1}), or \code{"less"}.
#' @param nperm Integer; number of permutation resamples (default \code{1000L}).
#' @export
median_hermite.default <- function(x, y, paired = FALSE,
                                   method = c("permutation", "analytical"),
                                   degree = 3L,
                                   alternative = c("two.sided", "greater", "less"),
                                   nperm = 1000L, ...) {
  method      <- match.arg(method)
  alternative <- match.arg(alternative)
  degree      <- as.integer(degree)
  if (degree < 1L || degree %% 2L == 0L) stop("'degree' must be an odd integer >= 1.")

  d   <- .two_sample_data(x, y, paired, degree)
  res <- .hermite_location_test(d$x1, d$x2, paired = paired, statistic = "median",
                                method = method, degree = degree,
                                alternative = alternative, nperm = nperm)
  res$median_diff <- res$estimate
  res$d_med       <- res$d_standardized
  res$method      <- paste0(if (paired) "Paired" else "Two-Sample",
                            " Hermite Median Difference Test")
  class(res) <- "median_hermite"
  res
}

# =============================================================================
# 3. shape_hermite: Scale, Asymmetry, and Tail-Weight Test
# =============================================================================

#' Distribution-Robust Shape Difference Test (Scale, Asymmetry, Tail Weight)
#'
#' Tests for differences in distributional shape between two independent
#' groups or paired measurements: \strong{scale} (\eqn{\Delta \log\sigma}),
#' \strong{asymmetry} (\eqn{\Delta c_2}, skewness direction), and \strong{tail
#' weight} (\eqn{\Delta c_3}, kurtosis direction), using Parseval-standardized
#' orthogonal Hermite weights and a location-aligned permutation null.
#'
#' @details
#' \subsection{The Hermite shape indices}{
#' For a quantile model \eqn{X = f(Z) = \sum_m a_m He_m(Z)} with regularized
#' variance \eqn{\sigma^2 = \sum_{m \ge 1} a_m^2\, m!} (Parseval), the
#' standardized weights
#' \deqn{c_m = \frac{a_m \sqrt{m!}}{\sigma}, \qquad \sum_{m \ge 1} c_m^2 = 1}
#' are direction cosines on the unit hypersphere: bounded in \eqn{[-1, 1]} and
#' location-scale invariant. Unlike classical moment ratios \eqn{g_1, g_2}
#' (volatile powers of order 3-4, Fleishman-boundary constrained), the
#' \eqn{c_m} are smooth linear functionals of the order statistics, yielding
#' well-behaved permutation distributions and substantially higher power (in
#' simulations, several times that of the Kolmogorov-Smirnov test for
#' asymmetry and tail departures).
#' }
#'
#' \subsection{Location alignment}{
#' With \code{align = "location"} (default), each group is centered at its
#' sample median before pooling, so that mean/median differences cannot leak
#' into the shape null (the tested null is equality of shape, not the omnibus
#' \eqn{F_1 = F_2}). Since \eqn{c_2}, \eqn{c_3}, and \eqn{\log\sigma} are
#' location-invariant, alignment does not alter the observed contrasts.
#' }
#'
#' \subsection{Multiplicity and omnibus testing}{
#' \describe{
#'   \item{Westfall-Young step-down}{Adjusted \eqn{p}-values controlling the
#'     family-wise error rate across the tested contrasts, exploiting the
#'     joint permutation distribution (uniformly more powerful than
#'     Bonferroni/Holm).}
#'   \item{\code{omnibus = "minP"} (default)}{Tippett minimum-\eqn{p} global
#'     test: powerful when the shape difference is concentrated in one
#'     contrast.}
#'   \item{\code{omnibus = "maxT"}}{Maximum standardized statistic.}
#'   \item{\code{omnibus = "T2"}}{Permutation Hotelling-type Mahalanobis
#'     statistic on the standardized contrast vector using the permutation
#'     covariance: powerful when differences are spread across several
#'     contrasts.}
#' }
#' }
#'
#' @param x A numeric vector of observations for Group 1 (or the baseline in
#'   paired designs), or a two-sided \code{\link[stats]{formula}} of the form
#'   \code{response ~ group}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"shape_hermite"} containing:
#' \describe{
#'   \item{\code{estimate}}{Named vector of observed contrasts (Group 2 \eqn{-}
#'     Group 1) among \code{scale} (\eqn{\Delta\log\sigma}), \code{asymmetry}
#'     (\eqn{\Delta c_2}), and \code{tailweight} (\eqn{\Delta c_3}).}
#'   \item{\code{variance_ratio}}{\eqn{\hat\sigma_2^2 / \hat\sigma_1^2 =
#'     \exp(2\,\Delta\log\sigma)}; reported when the scale contrast is tested.}
#'   \item{\code{p_value}, \code{p_adjusted}}{Raw and Westfall-Young adjusted
#'     permutation \eqn{p}-values.}
#'   \item{\code{omnibus}}{List with \code{type}, \code{statistic}, and
#'     \code{p_value} of the global shape test (or \code{NULL}).}
#'   \item{\code{indices}}{Matrix of group shape indices
#'     (\eqn{\log\sigma}, \eqn{c_2}, \eqn{c_3}).}
#'   \item{\code{perm_distribution}}{Matrix of permutation contrasts under \eqn{H_0}.}
#'   \item{\code{n1}, \code{n2}, \code{degree}, \code{align}, \code{contrasts},
#'     \code{paired}, \code{nperm}, \code{n_perm_success}, \code{success_rate},
#'     \code{method}}{Configuration and bookkeeping.}
#' }
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). \emph{Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment}. John Wiley & Sons.
#'
#' @seealso \code{\link{t_hermite}}, \code{\link{median_hermite}},
#'   \code{\link{hermite_test}}, \code{\link{hermite_fit}}
#'
#' @examples
#' # Asymmetry difference: normal vs. standardized gamma (equal mean/SD)
#' set.seed(42)
#' g1 <- rnorm(40)
#' g2 <- (rgamma(40, shape = 2) - 2) / sqrt(2)
#'
#' res <- shape_hermite(g1, g2, nperm = 500)
#' print(res)
#' plot(res)
#'
#' # Variance-only comparison, reported as a variance ratio
#' shape_hermite(g1, 1.5 * g1 + rnorm(40, sd = 0.1),
#'               contrasts = "scale", nperm = 500)
#'
#' # Formula interface
#' df <- data.frame(score = c(g1, g2),
#'                  group = factor(rep(c("Symmetric", "Skewed"), each = 40)))
#' shape_hermite(score ~ group, data = df, nperm = 500)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
shape_hermite <- function(x, ...) {
  UseMethod("shape_hermite")
}

#' @rdname shape_hermite
#' @param formula A two-sided formula of the form \code{response ~ group}.
#' @param data An optional data frame containing the variables in \code{formula}.
#' @export
shape_hermite.formula <- function(formula, data = NULL, paired = FALSE,
                                  degree = 3L,
                                  align = c("location", "none"),
                                  contrasts = c("asymmetry", "tailweight", "scale"),
                                  omnibus = c("minP", "maxT", "T2", "none"),
                                  nperm = 1000L, min_success = 0.90, ...) {
  fg  <- .formula_groups(formula, data)
  res <- shape_hermite.default(x = fg$x, y = fg$y, paired = paired,
                               degree = degree, align = align,
                               contrasts = contrasts, omnibus = omnibus,
                               nperm = nperm, min_success = min_success, ...)
  res$group_labels <- fg$labels
  res
}

#' @rdname shape_hermite
#' @param y Numeric vector of observations for Group 2.
#' @param paired Logical; \code{TRUE} for paired samples. Default \code{FALSE}.
#' @param degree Integer scalar; polynomial degree (default \code{3L}; odd, \eqn{\ge 3}).
#' @param align Alignment mode: \code{"location"} (default; median-centered
#'   pooling) or \code{"none"} (raw exchangeability; tests \eqn{F_1 = F_2}).
#' @param contrasts Character vector; any subset of \code{"asymmetry"}
#'   (\eqn{c_2}), \code{"tailweight"} (\eqn{c_3}), and \code{"scale"}
#'   (\eqn{\log\sigma}).
#' @param omnibus Omnibus method: \code{"minP"} (default), \code{"maxT"},
#'   \code{"T2"}, or \code{"none"}. See Details.
#' @param nperm Integer; number of permutation resamples (default \code{1000L}).
#' @param min_success Minimum proportion of usable permutation replicates
#'   before a warning is issued (default \code{0.90}).
#' @export
shape_hermite.default <- function(x, y, paired = FALSE, degree = 3L,
                                  align = c("location", "none"),
                                  contrasts = c("asymmetry", "tailweight", "scale"),
                                  omnibus = c("minP", "maxT", "T2", "none"),
                                  nperm = 1000L, min_success = 0.90, ...) {

  align     <- match.arg(align)
  omnibus   <- match.arg(omnibus)
  contrasts <- match.arg(contrasts, several.ok = TRUE)
  degree    <- as.integer(degree)
  if (degree < 3L || degree %% 2L == 0L) stop("'degree' must be an odd integer >= 3.")

  d  <- .two_sample_data(x, y, paired, degree)
  x1 <- d$x1; x2 <- d$x2
  n1 <- length(x1); n2 <- length(x2)

  # Observed profiles and contrasts
  p1 <- .hermite_profile(x1, degree)
  p2 <- .hermite_profile(x2, degree)
  if (is.null(p1) || is.null(p2)) stop("Degenerate quantile fit; shape contrasts undefined.")
  obs <- .shape_contrast(p1, p2)[contrasts]
  K   <- length(obs)

  # Location-aligned permutation pool
  if (align == "location") {
    r1 <- x1 - stats::median(x1)
    r2 <- x2 - stats::median(x2)
  } else {
    r1 <- x1; r2 <- x2
  }

  Tp <- matrix(NA_real_, nrow = nperm, ncol = K, dimnames = list(NULL, contrasts))

  if (paired) {
    for (b in seq_len(nperm)) {
      sw <- stats::runif(n1) < 0.5
      pa <- .hermite_profile(ifelse(sw, r2, r1), degree)
      pb <- .hermite_profile(ifelse(sw, r1, r2), degree)
      if (!is.null(pa) && !is.null(pb)) Tp[b, ] <- .shape_contrast(pa, pb)[contrasts]
    }
  } else {
    pool <- c(r1, r2); N <- n1 + n2
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)
      pa  <- .hermite_profile(pool[idx],  degree)
      pb  <- .hermite_profile(pool[-idx], degree)
      if (!is.null(pa) && !is.null(pb)) Tp[b, ] <- .shape_contrast(pa, pb)[contrasts]
    }
  }

  keep <- stats::complete.cases(Tp)
  n_success <- sum(keep); success_rate <- n_success / nperm
  if (n_success < 10L) stop("Too few usable permutation replicates.")
  if (success_rate < min_success) {
    warning(sprintf("Only %d of %d permutation replicates were usable (%.1f%%).",
                    n_success, nperm, 100 * success_rate), call. = FALSE)
  }
  Tp <- Tp[keep, , drop = FALSE]

  # Standardize by permutation moments; Westfall-Young on |Z|
  ctr <- colMeans(Tp)
  scl <- apply(Tp, 2L, stats::sd)
  scl[!is.finite(scl) | scl <= 0] <- 1

  Zall <- sweep(sweep(rbind(obs, Tp), 2L, ctr, "-"), 2L, scl, "/")
  wy   <- .wy_stepdown(abs(Zall))
  p_raw <- stats::setNames(wy$p_raw, contrasts)
  p_adj <- stats::setNames(wy$p_adj, contrasts)

  omni <- .omnibus_test(Zall, wy$P, omnibus)

  structure(
    list(
      estimate          = obs,
      variance_ratio    = if ("scale" %in% contrasts) exp(2 * obs[["scale"]]) else NULL,
      p_value           = p_raw,
      p_adjusted        = p_adj,
      omnibus           = omni,
      indices           = rbind(group1 = c(scale = p1$log_sd, asymmetry = p1$c2, tailweight = p1$c3),
                                group2 = c(scale = p2$log_sd, asymmetry = p2$c2, tailweight = p2$c3)),
      n1                = n1, n2 = n2,
      degree            = degree,
      align             = align,
      contrasts         = contrasts,
      paired            = paired,
      nperm             = nperm,
      n_perm_success    = n_success,
      success_rate      = success_rate,
      perm_distribution = Tp,
      method            = paste0(if (paired) "Paired" else "Two-Sample",
                                 " Hermite Shape Test (",
                                 if (align == "location") "location-aligned" else "exchangeable",
                                 " permutation)")
    ),
    class = "shape_hermite"
  )
}

# =============================================================================
# 4. hermite_test: Unified Interface and Complete Profile Test
# =============================================================================

#' Unified Distribution-Robust Hypothesis Test (Hermite Framework)
#'
#' Single entry point to the \pkg{hermiteStats} hypothesis testing suite.
#' Individual aspects (\code{"mean"}, \code{"median"}, \code{"variance"},
#' \code{"asymmetry"}, \code{"tailweight"}, \code{"scale"}) delegate to the
#' canonical test functions; \code{test = "complete"} (default) performs a
#' \strong{joint five-dimensional distributional profile test} with
#' Westfall-Young family-wise error control across all contrasts and a global
#' omnibus test of distributional equality.
#'
#' @details
#' \subsection{Test dispatch}{
#' \describe{
#'   \item{\code{"mean"}}{\code{\link{t_hermite}} (regularized mean \eqn{a_0}).}
#'   \item{\code{"median"}}{\code{\link{median_hermite}} (regularized median \eqn{f(0) = \beta_0}).}
#'   \item{\code{"variance"}, \code{"scale"}}{\code{\link{shape_hermite}} with
#'     the scale contrast \eqn{\Delta\log\sigma}. Both names address the same
#'     underlying test (a location-aligned permutation test of
#'     \eqn{H_0\!: \sigma_1 = \sigma_2}); the result additionally reports the
#'     variance ratio \eqn{\hat\sigma_2^2/\hat\sigma_1^2 = \exp(2\Delta\log\sigma)}.}
#'   \item{\code{"asymmetry"}, \code{"tailweight"}}{\code{\link{shape_hermite}}
#'     with the \eqn{c_2} or \eqn{c_3} contrast, respectively.}
#' }
#' For these single tests, the returned object is the native object of the
#' delegated function (\code{"t_hermite"}, \code{"median_hermite"}, or
#' \code{"shape_hermite"}), so all of its print/plot diagnostics remain
#' available.
#' }
#'
#' \subsection{The complete profile test (\code{test = "complete"})}{
#' All five contrasts -- \eqn{\Delta\mu} (studentized), \eqn{\Delta\mathrm{Med}}
#' (studentized), \eqn{\Delta\log\sigma}, \eqn{\Delta c_2}, \eqn{\Delta c_3} --
#' are evaluated against a \emph{single joint permutation null} built with
#' common permutation indices, so that the dependence structure among the
#' statistics is preserved:
#' \itemize{
#'   \item Location contrasts are computed on per-group \emph{aligned} pools
#'     (mean-centered for the mean contrast; median-centered for the median
#'     contrast), so each marginal null remains valid even when other aspects
#'     of the distributions differ (approximate subset pivotality).
#'   \item Shape contrasts are location-invariant and are computed on the
#'     median-centered pool at no additional cost (they share the median
#'     contrast's quantile fits).
#'   \item In paired designs, a single per-pair random swap vector drives the
#'     sign-flip of the difference scores (location contrasts) and the group
#'     swap (shape contrasts) simultaneously.
#' }
#' All contrasts are standardized by their permutation moments; Westfall-Young
#' step-down adjustment is applied across the full five-dimensional family,
#' and the omnibus test (default \code{"minP"}) provides a global test of
#' distributional equality across any aspect. The complete test is inherently
#' two-sided and always permutation-based; \code{method = "analytical"} is
#' only honored by the single mean/median tests.
#' }
#'
#' @param x A numeric vector of observations for Group 1 (or the baseline in
#'   paired designs), or a two-sided \code{\link[stats]{formula}} of the form
#'   \code{response ~ group}.
#' @param ... Additional arguments passed to methods (and, for single tests,
#'   forwarded to the delegated test function).
#'
#' @return
#' For \code{test = "complete"}, an S3 object of class \code{"hermite_test"} containing:
#' \describe{
#'   \item{\code{results}}{Data frame with one row per contrast: \code{test},
#'     \code{estimate}, \code{statistic} (studentized \eqn{t} for mean/median;
#'     permutation-standardized \eqn{z} for shape contrasts), \code{p_value}
#'     (raw permutation), and \code{p_adjusted} (Westfall-Young).}
#'   \item{\code{omnibus}}{Global test of distributional equality
#'     (\code{type}, \code{statistic}, \code{p_value}).}
#'   \item{\code{variance_ratio}}{\eqn{\hat\sigma_2^2/\hat\sigma_1^2}.}
#'   \item{\code{profiles}}{Matrix of regularized group profiles
#'     (mean, median, sd, \eqn{c_2}, \eqn{c_3}).}
#'   \item{\code{perm_distribution}}{Joint permutation contrast matrix.}
#'   \item{\code{n1}, \code{n2}, \code{degree}, \code{paired}, \code{nperm},
#'     \code{n_perm_success}, \code{success_rate}, \code{method}}{Configuration
#'     and bookkeeping.}
#' }
#' For single tests, the native object of the delegated function (see Details).
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). \emph{Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment}. John Wiley & Sons.
#'
#' @seealso \code{\link{t_hermite}}, \code{\link{median_hermite}},
#'   \code{\link{shape_hermite}}
#'
#' @examples
#' set.seed(42)
#' g1 <- rnorm(50, mean = 10, sd = 2)
#' g2 <- 10 + 2 * (rgamma(50, shape = 2) - 2) / sqrt(2)   # equal mean/SD, skewed
#'
#' # Complete distributional profile test (default)
#' res <- hermite_test(g1, g2, nperm = 500)
#' print(res)
#' plot(res)
#'
#' # Single aspects (delegated to the canonical functions)
#' hermite_test(g1, g2, test = "median", nperm = 500)
#' hermite_test(g1, g2, test = "variance", nperm = 500)
#'
#' # Formula interface
#' df <- data.frame(score = c(g1, g2),
#'                  group = factor(rep(c("A", "B"), each = 50)))
#' hermite_test(score ~ group, data = df, nperm = 500)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
hermite_test <- function(x, ...) {
  UseMethod("hermite_test")
}

#' @rdname hermite_test
#' @param formula A two-sided formula of the form \code{response ~ group}.
#' @param data An optional data frame containing the variables in \code{formula}.
#' @export
hermite_test.formula <- function(formula, data = NULL,
                                 test = c("complete", "mean", "median", "variance",
                                          "asymmetry", "tailweight", "scale"),
                                 paired = FALSE,
                                 method = c("permutation", "analytical"),
                                 degree = 3L,
                                 alternative = c("two.sided", "greater", "less"),
                                 nperm = 1000L,
                                 omnibus = c("minP", "maxT", "T2", "none"), ...) {
  fg  <- .formula_groups(formula, data)
  res <- hermite_test.default(x = fg$x, y = fg$y, test = test, paired = paired,
                              method = method, degree = degree,
                              alternative = alternative, nperm = nperm,
                              omnibus = omnibus, ...)
  res$group_labels <- fg$labels
  res
}

#' @rdname hermite_test
#' @param y Numeric vector of observations for Group 2.
#' @param test Character string selecting the aspect(s) to test; see Details.
#'   Default \code{"complete"}.
#' @param paired Logical; \code{TRUE} for paired samples. Default \code{FALSE}.
#' @param method Character string; \code{"permutation"} (default) or
#'   \code{"analytical"}. Only honored by the single mean/median tests; the
#'   complete and shape tests are always permutation-based.
#' @param degree Integer scalar; polynomial degree (default \code{3L}; odd, \eqn{\ge 3}).
#' @param alternative Alternative hypothesis for single mean/median tests;
#'   the complete and shape tests are inherently two-sided.
#' @param nperm Integer; number of permutation resamples (default \code{1000L}).
#' @param omnibus Omnibus method for the complete test: \code{"minP"}
#'   (default), \code{"maxT"}, \code{"T2"}, or \code{"none"}.
#' @export
hermite_test.default <- function(x, y,
                                 test = c("complete", "mean", "median", "variance",
                                          "asymmetry", "tailweight", "scale"),
                                 paired = FALSE,
                                 method = c("permutation", "analytical"),
                                 degree = 3L,
                                 alternative = c("two.sided", "greater", "less"),
                                 nperm = 1000L,
                                 omnibus = c("minP", "maxT", "T2", "none"), ...) {

  test        <- match.arg(test)
  method      <- match.arg(method)
  alternative <- match.arg(alternative)
  omnibus     <- match.arg(omnibus)
  degree      <- as.integer(degree)
  if (degree < 3L || degree %% 2L == 0L) stop("'degree' must be an odd integer >= 3.")

  # --- Single-test delegation ------------------------------------------------
  if (test == "mean") {
    return(t_hermite(x, y, paired = paired, method = method, degree = degree,
                     alternative = alternative, nperm = nperm, ...))
  }
  if (test == "median") {
    return(median_hermite(x, y, paired = paired, method = method, degree = degree,
                          alternative = alternative, nperm = nperm, ...))
  }
  if (test %in% c("variance", "scale", "asymmetry", "tailweight")) {
    contrast <- if (test == "variance") "scale" else test
    return(shape_hermite(x, y, paired = paired, degree = degree,
                         contrasts = contrast, omnibus = "none",
                         nperm = nperm, ...))
  }

  # --- Complete five-dimensional profile test ---------------------------------
  if (method == "analytical") {
    message("The complete profile test is permutation-based; 'method = \"analytical\"' is ignored.")
  }

  d  <- .two_sample_data(x, y, paired, degree)
  x1 <- d$x1; x2 <- d$x2
  n1 <- length(x1); n2 <- length(x2)

  p1 <- .hermite_profile(x1, degree)
  p2 <- .hermite_profile(x2, degree)
  if (is.null(p1) || is.null(p2)) stop("Degenerate quantile fit in at least one group.")

  contrast_names <- c("mean", "median", "scale", "asymmetry", "tailweight")

  if (paired) {
    d0 <- x2 - x1
    pd <- .hermite_profile(d0, degree)
    if (is.null(pd)) stop("Degenerate quantile fit for the difference scores.")

    obs <- c(mean   = pd$mean / pd$se_mean,
             median = pd$median / pd$se_med,
             .shape_contrast(p1, p2))
    est <- c(mean = pd$mean, median = pd$median, .shape_contrast(p1, p2))

    r1 <- x1 - stats::median(x1)
    r2 <- x2 - stats::median(x2)

    Tp <- matrix(NA_real_, nperm, 5L, dimnames = list(NULL, contrast_names))
    for (b in seq_len(nperm)) {
      sw  <- stats::runif(n1) < 0.5
      fb  <- .hermite_profile(ifelse(sw, -d0, d0), degree)
      pa  <- .hermite_profile(ifelse(sw, r2, r1), degree)
      pb  <- .hermite_profile(ifelse(sw, r1, r2), degree)
      if (!is.null(fb) && !is.null(pa) && !is.null(pb)) {
        Tp[b, ] <- c(fb$mean / fb$se_mean,
                     fb$median / fb$se_med,
                     .shape_contrast(pa, pb))
      }
    }
  } else {
    obs <- c(mean   = (p2$mean - p1$mean) / sqrt(p1$se_mean^2 + p2$se_mean^2),
             median = (p2$median - p1$median) / sqrt(p1$se_med^2 + p2$se_med^2),
             .shape_contrast(p1, p2))
    est <- c(mean   = p2$mean - p1$mean,
             median = p2$median - p1$median,
             .shape_contrast(p1, p2))

    pool_mean <- c(x1 - p1$mean,   x2 - p2$mean)     # mean-aligned pool
    pool_med  <- c(x1 - p1$median, x2 - p2$median)   # median-aligned pool
    N <- n1 + n2

    Tp <- matrix(NA_real_, nperm, 5L, dimnames = list(NULL, contrast_names))
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)      # common indices: joint null
      qa_m <- .hermite_profile(pool_mean[idx],  degree)
      qb_m <- .hermite_profile(pool_mean[-idx], degree)
      qa_d <- .hermite_profile(pool_med[idx],   degree)
      qb_d <- .hermite_profile(pool_med[-idx],  degree)
      if (!is.null(qa_m) && !is.null(qb_m) && !is.null(qa_d) && !is.null(qb_d)) {
        Tp[b, ] <- c((qb_m$mean - qa_m$mean) / sqrt(qa_m$se_mean^2 + qb_m$se_mean^2),
                     (qb_d$median - qa_d$median) / sqrt(qa_d$se_med^2 + qb_d$se_med^2),
                     .shape_contrast(qa_d, qb_d))
      }
    }
  }

  keep <- stats::complete.cases(Tp)
  n_success <- sum(keep); success_rate <- n_success / nperm
  if (n_success < 10L) stop("Too few usable permutation replicates.")
  if (success_rate < 0.90) {
    warning(sprintf("Only %d of %d permutation replicates were usable (%.1f%%).",
                    n_success, nperm, 100 * success_rate), call. = FALSE)
  }
  Tp <- Tp[keep, , drop = FALSE]

  ctr <- colMeans(Tp)
  scl <- apply(Tp, 2L, stats::sd)
  scl[!is.finite(scl) | scl <= 0] <- 1

  Zall <- sweep(sweep(rbind(obs, Tp), 2L, ctr, "-"), 2L, scl, "/")
  wy   <- .wy_stepdown(abs(Zall))
  omni <- .omnibus_test(Zall, wy$P, omnibus)

  results <- data.frame(
    test       = contrast_names,
    estimate   = unname(est),
    statistic  = unname(c(obs[1:2], Zall[1L, 3:5])),
    p_value    = unname(wy$p_raw),
    p_adjusted = unname(wy$p_adj),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      results           = results,
      omnibus           = omni,
      variance_ratio    = exp(2 * est[["scale"]]),
      profiles          = rbind(
        group1 = c(mean = p1$mean, median = p1$median, sd = p1$sd, c2 = p1$c2, c3 = p1$c3),
        group2 = c(mean = p2$mean, median = p2$median, sd = p2$sd, c2 = p2$c2, c3 = p2$c3)
      ),
      n1                = n1, n2 = n2,
      degree            = degree,
      paired            = paired,
      nperm             = nperm,
      n_perm_success    = n_success,
      success_rate      = success_rate,
      perm_distribution = Tp,
      method            = paste0(if (paired) "Paired" else "Two-Sample",
                                 " Hermite Distributional Profile Test (joint permutation)")
    ),
    class = "hermite_test"
  )
}

# =============================================================================
# 5. S3 Print and Plot Methods
# =============================================================================

#' Shared print engine for location tests
#' @noRd
.print_location_test <- function(x, digits, stat_label, est_label, extra_lines = NULL) {
  cat(sprintf("\n  %s (%s)\n", x$method, x$method_type))
  cat(strrep("-", 60), "\n", sep = "")
  cat(sprintf("  %s :  %.*f\n", format(stat_label, width = 26), digits, x$statistic))
  cat(sprintf("  %s :  %.*f\n", format(est_label, width = 26), digits, x$estimate))
  cat(sprintf("  %s :  %.*f\n", format("Standardized Effect (d)", width = 26), digits, x$d_standardized))
  if (!is.null(extra_lines)) extra_lines(x, digits)
  cat(sprintf("  %s :  %.*f (df = %.1f)\n", format("Standard Error (SE_diff)", width = 26),
              digits, x$se_diff, x$df))
  cat(sprintf("  %s :  n1 = %d, n2 = %d\n", format("Sample Sizes", width = 26), x$n1, x$n2))
  cat(sprintf("  %s :  %s\n", format("Alternative", width = 26), x$alternative))
  if (x$method_type == "permutation") {
    cat(sprintf("  %s :  %d / %d (%.1f%%)\n", format("Permutations (usable)", width = 26),
                x$n_perm_success, x$nperm, 100 * x$success_rate))
  }
  cat(sprintf("  %s :  %s\n", format("p-value", width = 26), .p_fmt(x$p_value, digits)))
  cat("\n")
  invisible(x)
}

#' Shared plot engine for location tests
#' @noRd
.plot_location_test <- function(x, stat_label, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  if (x$method_type == "permutation" && length(x$perm_distribution) > 0L) {
    d <- stats::density(x$perm_distribution)
    xlim_vals <- range(d$x, x$statistic, -x$statistic)
    xlim <- xlim_vals + c(-1, 1) * 0.1 * diff(xlim_vals)
    plot(d, main = paste0("Permutation Null (", x$method, ")"),
         xlab = paste0(stat_label, " (under H0)"), bty = "l",
         col = "darkblue", lwd = 2, xlim = xlim,
         panel.first = graphics::grid(col = "gray90", lty = 1), ...)
    graphics::polygon(d, col = grDevices::adjustcolor("steelblue", alpha.f = 0.2), border = NA)
    null_lab <- "Permutation Null"
  } else {
    x_seq <- seq(-max(4, abs(x$statistic) + 1), max(4, abs(x$statistic) + 1), length.out = 300L)
    y_seq <- stats::dt(x_seq, df = x$df)
    plot(x_seq, y_seq, type = "l", col = "darkblue", lwd = 2, bty = "l",
         main = sprintf("Theoretical Null: t(df = %.1f)", x$df),
         xlab = paste0(stat_label, " (under H0)"), ylab = "Density",
         panel.first = graphics::grid(col = "gray90", lty = 1), ...)
    graphics::polygon(c(x_seq, rev(x_seq)), c(y_seq, rep(0, length(y_seq))),
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.2), border = NA)
    null_lab <- "Theoretical t Null"
  }

  graphics::abline(v = x$statistic, col = "firebrick", lwd = 2.5)
  if (x$alternative == "two.sided") {
    graphics::abline(v = -x$statistic, col = "firebrick", lwd = 1.5, lty = 2)
  }
  graphics::mtext(sprintf("Observed t = %.3f  |  p = %.3f (%s)",
                          x$statistic, x$p_value, x$alternative),
                  side = 3, line = 0.3, cex = 0.85, col = "gray20")
  graphics::legend("topright",
                   legend = c(null_lab, sprintf("Observed t = %.3f", x$statistic)),
                   col = c("darkblue", "firebrick"), lwd = c(2, 2.5),
                   bty = "n", cex = 0.85)
  invisible(x)
}

#' @export
print.t_hermite <- function(x, digits = 3L, ...) {
  .print_location_test(
    x, digits,
    stat_label = "t_Hermite Statistic",
    est_label  = "Regularized Mean Diff",
    extra_lines = function(x, digits) {
      cat(sprintf("    Group 1 Mean             :  %.*f (SD = %.*f)\n",
                  digits, x$mean1, digits, x$sd1))
      cat(sprintf("    Group 2 Mean             :  %.*f (SD = %.*f)\n",
                  digits, x$mean2, digits, x$sd2))
    }
  )
}

#' @export
print.median_hermite <- function(x, digits = 3L, ...) {
  .print_location_test(
    x, digits,
    stat_label = "t_Median Statistic",
    est_label  = "Regularized Median Diff",
    extra_lines = function(x, digits) {
      cat(sprintf("    Group 1 Median           :  %.*f (Mean = %.*f, c2 = %.*f)\n",
                  digits, x$median1, digits, x$mean1, digits, x$c2_1))
      cat(sprintf("    Group 2 Median           :  %.*f (Mean = %.*f, c2 = %.*f)\n",
                  digits, x$median2, digits, x$mean2, digits, x$c2_2))
    }
  )
}

#' @export
plot.t_hermite <- function(x, ...) {
  .plot_location_test(x, stat_label = "t_Hermite", ...)
}

#' @export
plot.median_hermite <- function(x, ...) {
  .plot_location_test(x, stat_label = "t_Median", ...)
}

#' Print Method for Hermite Shape-Difference Test Objects
#'
#' Displays the group shape indices, contrast estimates, raw and
#' Westfall-Young adjusted permutation \eqn{p}-values, the variance ratio
#' (when the scale contrast is tested), and the omnibus result.
#'
#' @param x An object of class \code{"shape_hermite"}.
#' @param digits Integer; number of decimal places (default \code{3L}).
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.shape_hermite <- function(x, digits = 3L, ...) {
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 64), "\n", sep = "")
  cat(sprintf("  Sample Sizes        :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Polynomial Degree   :  %d\n", x$degree))
  cat(sprintf("  Permutations        :  %d / %d usable (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))

  lab <- c(asymmetry = "Asymmetry (c2)", tailweight = "Tail Weight (c3)",
           scale = "Scale (log sigma)")
  g1 <- x$indices["group1", ]; g2 <- x$indices["group2", ]

  cat("\n  Hermite Contrast Profile (Group 2 - Group 1):\n")
  tab <- data.frame(
    Contrast   = unname(lab[x$contrasts]),
    Group_1    = sprintf("%.*f", digits, unname(g1[x$contrasts])),
    Group_2    = sprintf("%.*f", digits, unname(g2[x$contrasts])),
    Difference = sprintf("%+.*f", digits, unname(x$estimate)),
    p_value    = .p_fmt(unname(x$p_value), digits),
    p_WY       = .p_fmt(unname(x$p_adjusted), digits),
    stringsAsFactors = FALSE
  )
  names(tab) <- c("Contrast", "Group 1", "Group 2", "Difference", "p-value", "p (WY)")
  print(tab, row.names = FALSE)

  if (!is.null(x$variance_ratio)) {
    cat(sprintf("\n  Variance Ratio (sigma2^2 / sigma1^2) :  %.*f\n",
                digits, x$variance_ratio))
  }
  if (!is.null(x$omnibus)) {
    cat(sprintf("\n  Omnibus Shape Test (%s):  Statistic = %.*f  |  p-value = %s\n",
                x$omnibus$type, digits, x$omnibus$statistic,
                .p_fmt(x$omnibus$p_value, digits)))
  }
  cat("\n")
  invisible(x)
}

#' Plot Method for Hermite Shape-Difference Test Objects
#'
#' Multi-panel visualization of the permutation null distributions per shape
#' contrast, with observed differences and (adjusted) \eqn{p}-values annotated.
#'
#' @param x An object of class \code{"shape_hermite"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#' @return Invisibly returns \code{x}.
#' @export
plot.shape_hermite <- function(x, ...) {
  if (is.null(x$perm_distribution) || nrow(x$perm_distribution) == 0L) {
    stop("No permutation distribution available to plot.")
  }
  K <- length(x$contrasts)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, K), mar = c(4.2, 4.0, 3.4, 1.0), font.main = 1)

  titles  <- c(asymmetry = "Asymmetry (c2)", tailweight = "Tail Weight (c3)",
               scale = "Scale (log sigma)")
  pal_col <- c(asymmetry = "#d95f02", tailweight = "#7570b3", scale = "#1b9e77")

  for (i in seq_along(x$contrasts)) {
    cn      <- x$contrasts[i]
    d       <- stats::density(x$perm_distribution[, cn])
    obs_val <- x$estimate[cn]
    xlim_v  <- range(d$x, obs_val, -obs_val)
    xlim    <- xlim_v + c(-1, 1) * 0.1 * diff(xlim_v)

    plot(d, main = titles[cn], xlab = paste0("Delta ", cn, " (under H0)"),
         ylab = if (i == 1L) "Density" else "", bty = "l",
         col = pal_col[cn], lwd = 2, xlim = xlim,
         panel.first = graphics::grid(col = "gray90", lty = 1), ...)
    graphics::polygon(d, col = grDevices::adjustcolor(pal_col[cn], alpha.f = 0.2), border = NA)
    graphics::abline(v = obs_val, col = "firebrick", lwd = 2.5)
    graphics::abline(v = -obs_val, col = "firebrick", lwd = 1.2, lty = 2)
    graphics::mtext(sprintf("Diff = %+.3f | p = %.3f (WY: %.3f)",
                            obs_val, x$p_value[cn], x$p_adjusted[cn]),
                    side = 3, line = 0.3, cex = 0.8, col = "gray20")
  }
  invisible(x)
}

#' Print Method for the Complete Hermite Profile Test
#'
#' Displays the regularized group profiles, the five-dimensional contrast
#' table with raw and Westfall-Young adjusted permutation \eqn{p}-values, the
#' variance ratio, and the global omnibus test.
#'
#' @param x An object of class \code{"hermite_test"}.
#' @param digits Integer; number of decimal places (default \code{3L}).
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.hermite_test <- function(x, digits = 3L, ...) {
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 68), "\n", sep = "")
  cat(sprintf("  Sample Sizes        :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Polynomial Degree   :  %d\n", x$degree))
  cat(sprintf("  Permutations        :  %d / %d usable (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))

  g_lab <- if (!is.null(x$group_labels)) x$group_labels else c("Group 1", "Group 2")
  cat("\n  Regularized Group Profiles:\n")
  prof <- as.data.frame(round(x$profiles, digits))
  rownames(prof) <- paste0("  ", g_lab)
  print(prof)

  lab <- c(mean = "Mean", median = "Median", scale = "Scale (log sigma)",
           asymmetry = "Asymmetry (c2)", tailweight = "Tail Weight (c3)")
  cat("\n  Contrasts (Group 2 - Group 1):\n")
  tab <- data.frame(
    Test      = unname(lab[x$results$test]),
    Estimate  = sprintf("%+.*f", digits, x$results$estimate),
    Statistic = sprintf("%.*f", digits, x$results$statistic),
    p_value   = .p_fmt(x$results$p_value, digits),
    p_WY      = .p_fmt(x$results$p_adjusted, digits),
    stringsAsFactors = FALSE
  )
  names(tab) <- c("Test", "Estimate", "Statistic", "p-value", "p (WY)")
  print(tab, row.names = FALSE)

  cat(sprintf("\n  Variance Ratio (sigma2^2 / sigma1^2) :  %.*f\n",
              digits, x$variance_ratio))
  if (!is.null(x$omnibus)) {
    cat(sprintf("  Omnibus Test (%s)   :  Statistic = %.*f  |  p-value = %s\n",
                x$omnibus$type, digits, x$omnibus$statistic,
                .p_fmt(x$omnibus$p_value, digits)))
  }
  cat("\n")
  invisible(x)
}

#' Plot Method for the Complete Hermite Profile Test
#'
#' Displays the joint permutation null distribution of every contrast
#' (mean, median, scale, asymmetry, tail weight) with the observed statistic
#' and adjusted \eqn{p}-values annotated.
#'
#' @param x An object of class \code{"hermite_test"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#' @return Invisibly returns \code{x}.
#' @export
plot.hermite_test <- function(x, ...) {
  Tp <- x$perm_distribution
  if (is.null(Tp) || nrow(Tp) == 0L) stop("No permutation distribution available to plot.")

  cn <- colnames(Tp)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(2, 3), mar = c(4.2, 4.0, 3.4, 1.0), font.main = 1)

  titles <- c(mean = "Mean (t)", median = "Median (t)",
              scale = "Scale (log sigma)", asymmetry = "Asymmetry (c2)",
              tailweight = "Tail Weight (c3)")
  pal <- c(mean = "#1f78b4", median = "#33a02c", scale = "#1b9e77",
           asymmetry = "#d95f02", tailweight = "#7570b3")

  # Observed statistics on the same scale as the permutation columns:
  # studentized t for mean/median, raw contrasts for shape.
  obs_map <- stats::setNames(x$results$statistic, x$results$test)
  obs_map[c("scale", "asymmetry", "tailweight")] <-
    stats::setNames(x$results$estimate, x$results$test)[c("scale", "asymmetry", "tailweight")]

  for (k in cn) {
    d       <- stats::density(Tp[, k])
    obs_val <- obs_map[[k]]
    xlim_v  <- range(d$x, obs_val, -obs_val)
    xlim    <- xlim_v + c(-1, 1) * 0.1 * diff(xlim_v)

    plot(d, main = titles[k], xlab = "Under H0", ylab = "Density", bty = "l",
         col = pal[k], lwd = 2, xlim = xlim,
         panel.first = graphics::grid(col = "gray90", lty = 1), ...)
    graphics::polygon(d, col = grDevices::adjustcolor(pal[k], alpha.f = 0.2), border = NA)
    graphics::abline(v = obs_val, col = "firebrick", lwd = 2.5)
    graphics::abline(v = -obs_val, col = "firebrick", lwd = 1.2, lty = 2)

    row_k <- x$results[x$results$test == k, ]
    graphics::mtext(sprintf("p = %.3f (WY: %.3f)", row_k$p_value, row_k$p_adjusted),
                    side = 3, line = 0.3, cex = 0.75, col = "gray20")
  }
  invisible(x)
}

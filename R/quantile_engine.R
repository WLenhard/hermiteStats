# =============================================================================
# quantile_engine.R -- Mathematical Core and Quantile Modeling Engine
# =============================================================================

#' @title Mathematical Core and Quantile Modeling Engine
#'
#' @description
#' Low-level routines that turn an empirical sample into a regularized,
#' continuous quantile function and extract its population moments in closed
#' algebraic form. These functions are the computational backbone shared by
#' \code{\link{hermite_fit}}, \code{\link{cor_hermite}}, \code{\link{d_reg}},
#' and the hypothesis-testing suite (\code{\link{hermite_test}},
#' \code{\link{t_hermite}}, \code{\link{median_hermite}},
#' \code{\link{shape_hermite}}).
#'
#' @details
#' Every distribution handled by \pkg{hermiteStats} is represented as a
#' polynomial quantile function \eqn{X = f(Z)}, where \eqn{Z \sim \mathcal{N}(0,1)}
#' is obtained from the data by rank-based inverse-normal transformation
#' (Normal Quantile Transformation, NQT). Because \eqn{Z} is standard normal by
#' construction, the coefficients of \eqn{f} can be re-expressed in the basis of
#' monic Probabilists' Hermite polynomials \eqn{He_m(z)}, which are orthogonal
#' under the standard normal measure:
#' \deqn{E[He_m(Z) He_n(Z)] = \delta_{mn}\, m!}
#' This orthogonality turns the extraction of the mean, variance, skewness, and
#' kurtosis of \eqn{f(Z)} into simple, closed-form algebraic sums over the
#' fitted coefficients, evaluated exactly via the raw moments \eqn{E[Z^j]} of
#' the standard normal distribution (a direct consequence of Isserlis' 1918
#' theorem for Gaussian variables). No numerical integration is required at any
#' stage.
#'
#' @name quantile_engine
#' @keywords internal
NULL

# -----------------------------------------------------------------------------
# Internal Helpers: Polynomial Evaluation, Raw Normal Moments, Hermite Operator
# -----------------------------------------------------------------------------

#' Evaluate a monomial polynomial at given points
#'
#' Canonical evaluator for \eqn{f(z) = \sum_{j=0}^k \beta_j z^j}, shared by the
#' fitting, checking, prediction, and plotting code so the power-basis
#' construction exists in exactly one place.
#'
#' @param beta Numeric vector of monomial coefficients \eqn{(\beta_0, \dots, \beta_k)}.
#' @param z Numeric vector of evaluation points.
#' @return Numeric vector of polynomial values, same length as \code{z}.
#' @noRd
.poly_eval <- function(beta, z) {
  deg <- length(beta) - 1L
  as.vector(outer(z, 0:deg, `^`) %*% beta)
}

#' Raw moments of the standard normal distribution, E[Z^j]
#'
#' Returns \eqn{E[Z^0], E[Z^1], \dots, E[Z^{\text{max\_j}}]}. Odd-order moments
#' are zero; even-order moments follow the double-factorial identity
#' \eqn{E[Z^{2k}] = (2k-1)!!}, a direct consequence of Isserlis' (1918) theorem.
#'
#' @param max_j Integer; highest moment order required.
#' @return Numeric vector of length \code{max_j + 1}, indexed from order 0.
#' @noRd
.z_raw_moments <- function(max_j) {
  moms <- numeric(max_j + 1L)
  moms[1L] <- 1.0  # E[Z^0]
  if (max_j >= 2L) {
    for (j in seq(2L, max_j, by = 2L)) {
      moms[j + 1L] <- prod(seq(j - 1L, 1L, by = -2L))
    }
  }
  moms
}

#' Change-of-basis matrix from monomials to monic Hermite polynomials
#'
#' Constructs the matrix \code{H} such that, for a polynomial given by monomial
#' coefficients \code{beta} (i.e. \eqn{f(z) = \sum_j \beta_j z^j}), the vector
#' \code{a = H \%*\% beta} holds the coefficients of the same polynomial
#' expressed in the orthogonal Probabilists' Hermite basis,
#' \eqn{f(z) = \sum_m a_m He_m(z)}. This re-expression is what enables the
#' closed-form moment formulas in \code{.extract_all_moments()}.
#'
#' @param degree Integer; polynomial degree.
#' @param z_moms Optional precomputed vector of standard normal raw moments
#'   (as returned by \code{.z_raw_moments()}), of length at least
#'   \code{2 * degree + 1}. Supplying this avoids redundant recomputation
#'   when the caller already has the moments available.
#' @return A \code{(degree + 1) x (degree + 1)} numeric matrix.
#' @noRd
.hermite_basis_matrix <- function(degree, z_moms = NULL) {
  if (is.null(z_moms)) z_moms <- .z_raw_moments(2L * degree)
  H_mat <- matrix(0.0, nrow = degree + 1L, ncol = degree + 1L)
  for (n in 0:degree) {
    for (m in 0:n) {
      if ((n - m) %% 2L == 0L) {
        k <- (n - m) / 2L
        H_mat[m + 1L, n + 1L] <- choose(n, 2L * k) * z_moms[2L * k + 1L]
      }
    }
  }
  H_mat
}

#' Evaluate the Probabilists' Hermite basis up to a given order
#'
#' Evaluates \eqn{He_0(z), He_1(z), \dots, He_{\text{degree}}(z)} at each
#' element of \code{z} via the stable three-term recurrence
#' \eqn{He_n(z) = z\,He_{n-1}(z) - (n-1)\,He_{n-2}(z)}. This is the single
#' canonical Hermite evaluator in the package (used by the correlation
#' diagnostics and the Mehler conditional-mean curve in the plot methods).
#'
#' @param z Numeric vector of evaluation points.
#' @param degree Integer; highest Hermite order required.
#' @return Numeric matrix with \code{length(z)} rows and \code{degree + 1}
#'   columns; column \code{m + 1} holds \eqn{He_m(z)}.
#' @noRd
.hermite_eval_basis <- function(z, degree) {
  n <- length(z)
  H <- matrix(0.0, nrow = n, ncol = degree + 1L)
  H[, 1L] <- 1.0
  if (degree >= 1L) H[, 2L] <- z
  if (degree >= 2L) {
    for (k in 2:degree) {
      H[, k + 1L] <- z * H[, k] - (k - 1L) * H[, k - 1L]
    }
  }
  H
}

#' Assemble a cor_hermite result from two already-fitted quantile models
#'
#' Used by \code{\link{d_reg}} in paired designs, where both marginals are
#' already fitted and refitting inside \code{cor_hermite()} would be wasteful.
#'
#' @noRd
.cor_hermite_assemble <- function(fit_x, fit_y, poly_degree_requested = 3L,
                                  copula = "none", monotonicity = "relaxed",
                                  ties_method = "average", trim = 0) {
  pair <- .cor_hermite_pair(fit_x, fit_y, copula = copula, trim = trim)
  res <- list(
    r_Hermite             = pair$r_Hermite,
    copula                = copula,
    rho_z                 = pair$rho_z,
    attenuation           = pair$attenuation,
    cov_xy                = pair$cov_xy,
    var_x                 = pair$var_x,
    var_y                 = pair$var_y,
    mean_x                = pair$mean_x,
    mean_y                = pair$mean_y,
    degrees               = pair$degrees,
    poly_degree_requested = poly_degree_requested,
    monotonicity          = monotonicity,
    ties_method           = ties_method,
    trim                  = trim,
    n                     = length(fit_x$x),
    fit_x                 = fit_x,
    fit_y                 = fit_y,
    cross_moments         = pair$cross_moments,
    rho_z_reference       = pair$rho_z_reference,
    x                     = fit_x$x,
    y                     = fit_y$x
  )
  class(res) <- "cor_hermite"
  res
}

#' Closed-form distributional moments of a fitted quantile polynomial
#'
#' Given the monomial coefficients of a fitted quantile function
#' \eqn{X = f(Z)}, returns its exact population mean, variance, skewness, and
#' kurtosis under \eqn{Z \sim \mathcal{N}(0,1)}. Mean and variance are obtained
#' directly from the orthogonal Hermite weights (Parseval's identity);
#' skewness and kurtosis follow from the third and fourth central moments of
#' the centered polynomial, evaluated via the raw Gaussian moments of \eqn{Z}.
#'
#' @param beta Numeric vector of monomial coefficients \eqn{(\beta_0, \dots, \beta_k)}.
#' @return A list with elements \code{mean}, \code{variance}, \code{sd},
#'   \code{skewness}, \code{excess_kurtosis}, \code{kurtosis}, and
#'   \code{hermite_coeffs} (the orthogonal Hermite weights \eqn{a_0, \dots, a_k}).
#' @noRd
.extract_all_moments <- function(beta) {
  deg <- length(beta) - 1L
  if (deg < 0L) {
    return(list(mean = NA_real_, variance = NA_real_, sd = NA_real_,
                skewness = NA_real_, excess_kurtosis = NA_real_,
                kurtosis = NA_real_, hermite_coeffs = numeric(0)))
  }

  # All Gaussian raw moments needed below (up to order 4*deg) are computed
  # once and shared between the basis change and the central-moment sums.
  z_moms <- .z_raw_moments(4L * deg)

  # Mean and variance via the orthogonal Hermite basis
  H <- .hermite_basis_matrix(deg, z_moms = z_moms)
  a <- as.vector(H %*% beta)
  mean_val <- a[1L]

  if (deg == 0L) {
    return(list(mean = mean_val, variance = 0.0, sd = 0.0,
                skewness = 0.0, excess_kurtosis = 0.0, kurtosis = 3.0,
                hermite_coeffs = a))
  }

  var_val <- max(0.0, sum((a[-1L]^2) * factorial(1L:deg)))
  sd_val  <- sqrt(var_val)

  if (var_val <= 1e-12) {
    return(list(mean = mean_val, variance = 0.0, sd = 0.0,
                skewness = 0.0, excess_kurtosis = 0.0, kurtosis = 3.0,
                hermite_coeffs = a))
  }

  # Centered monomial polynomial, used for the 3rd and 4th central moments
  beta_c <- beta
  beta_c[1L] <- beta_c[1L] - mean_val

  # 3rd central moment (mu3)
  idx3  <- outer(outer(0:deg, 0:deg, `+`), 0:deg, `+`)
  coef3 <- outer(outer(beta_c, beta_c), beta_c)
  mu3   <- sum(coef3 * z_moms[idx3 + 1L])
  skew_val <- mu3 / (sd_val^3)

  # 4th central moment (mu4)
  idx4  <- outer(outer(outer(0:deg, 0:deg, `+`), 0:deg, `+`), 0:deg, `+`)
  coef4 <- outer(outer(outer(beta_c, beta_c), beta_c), beta_c)
  mu4   <- sum(coef4 * z_moms[idx4 + 1L])
  kurt_raw    <- mu4 / (var_val^2)
  kurt_excess <- kurt_raw - 3.0

  list(
    mean            = mean_val,
    variance        = var_val,
    sd              = sd_val,
    skewness        = skew_val,
    excess_kurtosis = kurt_excess,
    kurtosis        = kurt_raw,
    hermite_coeffs  = a
  )
}

# -----------------------------------------------------------------------------
# Exported Function: check_monotonicity
# -----------------------------------------------------------------------------

#' Check Monotonicity of a Polynomial Quantile Function
#'
#' A fitted quantile function \eqn{X = f(Z)} must be non-decreasing to be a
#' valid quantile function. Because unconstrained polynomial regression does
#' not guarantee this, \code{check_monotonicity} verifies it after fitting, and
#' \code{\link{hermite_fit}} uses it internally to decide whether to accept a
#' given polynomial degree or step down to a lower one.
#'
#' @param coeffs Numeric vector of monomial coefficients \eqn{(\beta_0, \beta_1, \dots, \beta_k)}
#'   describing the polynomial \eqn{f(z) = \sum_{j=0}^k \beta_j z^j}.
#' @param z_range Numeric vector of length 2, \code{c(min_z, max_z)}, giving the
#'   domain over which monotonicity is evaluated. Default \code{c(-4, 4)}, which
#'   covers essentially the entire probability mass for common sample sizes.
#' @param method Character string specifying the verification criterion:
#'   \describe{
#'     \item{\code{"relaxed"}}{(Default) Checks that the fitted values
#'       \eqn{\hat{x} = f(z)} rank-correlate with \eqn{z} at
#'       \eqn{\rho_{\mathrm{Spearman}} \ge 0.95} and that the linear coefficient
#'       is positive (\eqn{\beta_1 > 0}). This tolerates tiny, practically
#'       irrelevant non-monotonicities (e.g. from heavy ties) while rejecting
#'       polynomials with genuine reversals.}
#'     \item{\code{"strict"}}{Requires the first derivative \eqn{f'(z)} to be
#'       non-negative everywhere on \code{z_range}. The global minimum of
#'       \eqn{f'(z)} is found analytically by locating the roots of
#'       \eqn{f''(z)} via \code{\link[base]{polyroot}} and evaluating
#'       \eqn{f'(z)} at these critical points and at the interval endpoints.}
#'     \item{\code{"none"}}{Skips the check and always returns \code{TRUE}.
#'       Useful for diagnostic purposes or when monotonicity is already
#'       guaranteed (e.g. linear fits).}
#'   }
#' @param z Optional numeric vector of standard normal scores at which to
#'   evaluate the \code{"relaxed"} criterion. If omitted, an evenly spaced
#'   grid of 100 points over \code{z_range} is used instead.
#'
#' @return A list with elements \code{is_monotonic} (logical),
#'   \code{min_derivative} (only computed for \code{method = "strict"}),
#'   \code{rank_concordance} (only computed for \code{method = "relaxed"}),
#'   and \code{method}.
#'
#' @examples
#' # A cubic polynomial that dips slightly but not enough to matter empirically
#' check_monotonicity(c(0, 1, 0, 0.02), method = "relaxed")
#'
#' # The same polynomial checked analytically over the full domain
#' check_monotonicity(c(0, 1, 0, 0.02), method = "strict")
#'
#' @export
check_monotonicity <- function(coeffs, z_range = c(-4, 4),
                               method = c("relaxed", "strict", "none"),
                               z = NULL) {
  method <- match.arg(method)
  coeffs[is.na(coeffs)] <- 0
  deg <- length(coeffs) - 1L

  if (method == "none" || deg <= 0L) {
    return(list(is_monotonic = TRUE, min_derivative = 0,
                rank_concordance = 1.0, method = method))
  }

  if (deg == 1L) {
    slope <- coeffs[2L]
    return(list(
      is_monotonic     = slope >= 0,
      min_derivative   = slope,
      rank_concordance = if (slope >= 0) 1.0 else -1.0,
      method           = method
    ))
  }

  if (method == "strict") {
    d1 <- coeffs[-1L] * seq_len(deg)          # coefficients of f'(z)
    pts <- z_range
    d2 <- d1[-1L] * seq_len(deg - 1L)         # coefficients of f''(z)
    if (any(d2 != 0)) {
      r <- tryCatch(polyroot(d2), error = function(e) complex(0))
      real_r <- Re(r)[abs(Im(r)) < 1e-8]
      pts <- c(pts, real_r[real_r >= z_range[1L] & real_r <= z_range[2L]])
    }
    slopes    <- vapply(unique(pts), function(v) .poly_eval(d1, v), numeric(1L))
    min_slope <- min(slopes)
    return(list(
      is_monotonic     = min_slope >= -1e-8,
      min_derivative   = min_slope,
      rank_concordance = NA_real_,
      method           = method
    ))
  }

  # method == "relaxed"
  if (is.null(z)) z <- seq(z_range[1L], z_range[2L], length.out = 100L)
  fitted_vals <- .poly_eval(coeffs, z)

  rank_cor <- tryCatch(stats::cor(z, fitted_vals, method = "spearman"),
                       error = function(e) 0)
  is_mono <- (!is.na(rank_cor)) && (rank_cor >= 0.95) && (coeffs[2L] > 0)

  list(
    is_monotonic     = is_mono,
    min_derivative   = NA_real_,
    rank_concordance = rank_cor,
    method           = method
  )
}

# -----------------------------------------------------------------------------
# Exported Function: hermite_fit
# -----------------------------------------------------------------------------

#' Fit a Regularized Monotone Quantile Polynomial
#'
#' Maps a numeric sample onto standard normal scores via rank-based
#' inverse-normal transformation (NQT) and fits a low-degree, monotone
#' polynomial quantile function \eqn{X = f(Z)} to these scores. The polynomial
#' is re-expressed in the orthogonal Probabilists' Hermite basis, from which
#' the population mean, variance, skewness, and excess kurtosis of the
#' modeled distribution are obtained in exact closed algebraic form (no
#' numerical integration, no distributional assumptions beyond the polynomial
#' quantile model itself).
#'
#' @param x Numeric vector of observations. Missing and infinite values are
#'   removed automatically.
#' @param degree Integer scalar; the maximum polynomial degree to attempt
#'   (default \code{3}). The realized degree may be lower: it is capped by the
#'   number of unique values in \code{x}, reduced further under heavy tying,
#'   and stepped down whenever the fit fails the monotonicity check.
#' @param monotonicity Monotonicity constraint passed to
#'   \code{\link{check_monotonicity}}: \code{"relaxed"} (default),
#'   \code{"strict"}, or \code{"none"}.
#' @param ties_method Method used to break rank ties before the inverse-normal
#'   transform: \code{"average"} (default; midranks) or \code{"random"}.
#' @param force_odd Logical; if \code{TRUE} (default), only odd polynomial
#'   degrees (1, 3, 5, ...) are considered. Odd-degree polynomials most
#'   naturally accommodate both skewed and symmetric quantile shapes without
#'   the boundary curvature artifacts even-degree terms tend to introduce.
#'
#' @details
#' The polynomial degree is chosen adaptively for each sample. Starting from
#' the highest admissible degree, coefficients are estimated by ordinary least
#' squares regression of \code{x} on powers of the normal scores \code{z}, and
#' the fit is accepted if it passes the requested monotonicity check;
#' otherwise the degree is reduced (by 2 if \code{force_odd = TRUE}, by 1
#' otherwise) and refit, down to a linear (degree-1) fallback if necessary.
#' This keeps the model as flexible as the data support while guaranteeing a
#' valid, monotone quantile function.
#'
#' @return An S3 object of class \code{"hermite_fit"} containing:
#' \describe{
#'   \item{\code{beta}}{Fitted monomial coefficients \eqn{(\beta_0, \dots, \beta_k)}.}
#'   \item{\code{hermite_coeffs}}{Orthogonal Hermite weights \eqn{(a_0, \dots, a_k)}.}
#'   \item{\code{degree}}{Realized polynomial degree.}
#'   \item{\code{degree_requested}}{Requested polynomial degree.}
#'   \item{\code{mean}}{Regularized population mean, \eqn{\mu = a_0}.}
#'   \item{\code{variance}}{Regularized population variance, \eqn{\sigma^2 = \sum_{m=1}^k a_m^2\, m!}.}
#'   \item{\code{sd}}{Regularized population standard deviation.}
#'   \item{\code{skewness}}{Regularized population skewness, \eqn{\gamma_1 = \mu_3 / \sigma^3}.}
#'   \item{\code{excess_kurtosis}}{Regularized population excess kurtosis, \eqn{\gamma_2 = \mu_4 / \sigma^4 - 3}.}
#'   \item{\code{kurtosis}}{Regularized raw (Pearson) kurtosis, \eqn{\beta_2 = \mu_4 / \sigma^4}.}
#'   \item{\code{n}}{Sample size.}
#'   \item{\code{n_unique}}{Number of unique observations.}
#'   \item{\code{tie_proportion}}{Proportion of tied values in the sample.}
#'   \item{\code{monotonicity}}{Monotonicity constraint applied.}
#'   \item{\code{ties_method}}{Tie-handling method applied.}
#'   \item{\code{x}}{Cleaned input observations.}
#'   \item{\code{z}}{Latent standard normal scores.}
#' }
#'
#' @references
#' Isserlis, L. (1918). On a formula for the product-moment coefficient of any
#' order of a normal frequency distribution. \emph{Biometrika}, 12(1/2),
#' 134-139. \doi{10.1093/biomet/12.1-2.134}
#'
#' @examples
#' set.seed(1)
#' x <- rlnorm(80, meanlog = 2, sdlog = 0.4)
#' fit <- hermite_fit(x)
#' print(fit)
#' plot(fit)
#'
#' # Predicted raw values at selected latent quantiles (median and quartiles)
#' predict(fit, z = qnorm(c(0.25, 0.50, 0.75)))
#'
#' @seealso \code{\link{hermite_moments}}, \code{\link{predict.hermite_fit}},
#'   \code{\link{check_monotonicity}}, \code{\link{cor_hermite}},
#'   \code{\link{d_reg}}, \code{\link{hermite_test}}
#' @export
hermite_fit <- function(x, degree = 3L,
                        monotonicity = c("relaxed", "strict", "none"),
                        ties_method = c("average", "random"),
                        force_odd = TRUE) {
  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  x_clean <- x[is.finite(x)]
  n <- length(x_clean)
  if (n < 3L) stop("Need at least 3 valid observations.")

  n_uniq   <- length(unique(x_clean))
  tie_prop <- 1 - (n_uniq / n)

  d_cap <- min(as.integer(degree), n_uniq - 1L)
  if (tie_prop > 0.30 && d_cap > 3L) {
    d_cap <- min(d_cap, max(3L, floor(n_uniq / 2L)))
  }
  if (force_odd && (d_cap %% 2L == 0L)) d_cap <- d_cap - 1L
  d_cap <- max(1L, d_cap)

  p <- (rank(x_clean, ties.method = ties_method) - 0.5) / n
  z <- stats::qnorm(p)
  z_range <- range(z)

  # Built once at the maximum candidate degree; lower-degree fits reuse the
  # relevant leading columns instead of recomputing the power basis.
  Zmat_full <- outer(z, 0:d_cap, `^`)

  cur_deg  <- d_cap
  step     <- if (force_odd) 2L else 1L
  beta_res <- NULL

  while (cur_deg >= 1L) {
    Zmat <- Zmat_full[, seq_len(cur_deg + 1L), drop = FALSE]
    beta <- tryCatch(qr.solve(Zmat, x_clean),
                     error = function(e) rep(NA_real_, cur_deg + 1L))

    if (!anyNA(beta)) {
      chk <- check_monotonicity(beta, z_range = z_range,
                                method = monotonicity, z = z)
      if (chk$is_monotonic) {
        beta_res <- beta
        break
      }
    }
    cur_deg <- cur_deg - step
  }

  if (is.null(beta_res)) {
    cur_deg  <- 1L
    Zmat     <- Zmat_full[, 1:2, drop = FALSE]
    beta_res <- tryCatch(qr.solve(Zmat, x_clean),
                         error = function(e) c(mean(x_clean), stats::sd(x_clean)))
  }

  moms <- .extract_all_moments(beta_res)

  res <- list(
    beta             = beta_res,
    hermite_coeffs   = moms$hermite_coeffs,
    degree           = cur_deg,
    degree_requested = degree,
    mean             = moms$mean,
    variance         = moms$variance,
    sd               = moms$sd,
    skewness         = moms$skewness,
    excess_kurtosis  = moms$excess_kurtosis,
    kurtosis         = moms$kurtosis,
    n                = n,
    n_unique         = n_uniq,
    tie_proportion   = tie_prop,
    monotonicity     = monotonicity,
    ties_method      = ties_method,
    x                = x_clean,
    z                = z
  )
  class(res) <- "hermite_fit"
  res
}

# -----------------------------------------------------------------------------
# Exported Functions: predict.hermite_fit, hermite_moments
# -----------------------------------------------------------------------------

#' Predict Raw Values from a Fitted Hermite Quantile Model
#'
#' Evaluates the fitted monotone quantile polynomial \eqn{X = f(Z)} at given
#' latent standard normal scores. Since \eqn{Z = \Phi^{-1}(p)}, this maps
#' probabilities to raw-scale quantiles of the modeled distribution; for
#' example, \code{predict(fit, z = 0)} returns the regularized median.
#'
#' @param object An object of class \code{"hermite_fit"}.
#' @param z Numeric vector of standard normal scores at which to evaluate the
#'   quantile function. If \code{NULL} (default), the fitted values at the
#'   model's own latent scores \code{object$z} are returned.
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric vector of predicted raw values, same length as \code{z}.
#'
#' @examples
#' fit <- hermite_fit(rlnorm(60, meanlog = 2, sdlog = 0.4))
#'
#' # Regularized median and quartiles
#' predict(fit, z = qnorm(c(0.25, 0.50, 0.75)))
#'
#' # Fitted values at the observed latent scores
#' head(predict(fit))
#'
#' @seealso \code{\link{hermite_fit}}, \code{\link{hermite_moments}}
#' @export
predict.hermite_fit <- function(object, z = NULL, ...) {
  if (is.null(z)) z <- object$z
  if (!is.numeric(z)) stop("'z' must be a numeric vector of standard normal scores.")
  .poly_eval(object$beta, z)
}

#' Extract Regularized Distributional Moments
#'
#' Convenience accessor that retrieves the regularized population moments
#' (\eqn{\mu, \sigma^2, \sigma, \gamma_1, \gamma_2}) from a fitted
#' \code{\link{hermite_fit}} object, without exposing the underlying
#' polynomial coefficients.
#'
#' @param object An object of class \code{"hermite_fit"} returned by
#'   \code{\link{hermite_fit}}.
#'
#' @return A named list with elements \code{mean}, \code{variance}, \code{sd},
#'   \code{skewness}, \code{excess_kurtosis}, and \code{kurtosis}.
#'
#' @examples
#' fit <- hermite_fit(rnorm(50))
#' hermite_moments(fit)
#'
#' @export
hermite_moments <- function(object) {
  if (!inherits(object, "hermite_fit")) {
    stop("Argument 'object' must be an object of class 'hermite_fit'.")
  }
  list(
    mean            = object$mean,
    variance        = object$variance,
    sd              = object$sd,
    skewness        = object$skewness,
    excess_kurtosis = object$excess_kurtosis,
    kurtosis        = object$kurtosis
  )
}

# -----------------------------------------------------------------------------
# S3 Methods: print and summary
# -----------------------------------------------------------------------------

#' Print a Fitted Hermite Quantile Model
#'
#' @param x An object of class \code{"hermite_fit"}.
#' @param digits Integer; number of decimal places to display.
#' @param ... Additional arguments (currently unused).
#' @return The object \code{x}, invisibly.
#' @export
print.hermite_fit <- function(x, digits = 3L, ...) {
  cat("\n  Regularized Quantile Model (Hermite Basis)\n")
  cat(strrep("-", 45), "\n", sep = "")
  cat(sprintf("  Modeled Mean (mu)      :  %.*f\n", digits, x$mean))
  cat(sprintf("  Modeled SD (sigma)     :  %.*f\n", digits, x$sd))
  cat(sprintf("  Modeled Variance       :  %.*f\n", digits, x$variance))
  cat(sprintf("  Modeled Skewness (g1)  :  %.*f\n", digits, x$skewness))
  cat(sprintf("  Excess Kurtosis (g2)   :  %.*f\n", digits, x$excess_kurtosis))
  cat(sprintf("  Fitted Polynomial Deg  :  %d (requested: %d)\n",
              x$degree, x$degree_requested))
  cat(sprintf("  Sample Size (n)        :  %d (unique: %d, ties: %.1f%%)\n",
              x$n, x$n_unique, x$tie_proportion * 100))
  cat(sprintf("  Monotonicity Check     :  %s\n", x$monotonicity))
  cat("\n")
  invisible(x)
}

#' Summarize a Fitted Hermite Quantile Model
#'
#' @param object An object of class \code{"hermite_fit"}.
#' @param digits Integer; number of decimal places to display.
#' @param ... Additional arguments (currently unused).
#' @return The object \code{object}, invisibly.
#' @export
summary.hermite_fit <- function(object, digits = 3L, ...) {
  print(object, digits = digits, ...)

  cat("  Polynomial Coefficients (Monomial raw beta):\n")
  beta_out <- object$beta
  names(beta_out) <- paste0("z^", 0:object$degree)
  print(round(beta_out, digits = digits))

  cat("\n  Orthogonal Hermite Weights (a_m):\n")
  herm_out <- object$hermite_coeffs
  names(herm_out) <- paste0("He_", 0:object$degree)
  print(round(herm_out, digits = digits))
  cat("\n")
  invisible(object)
}

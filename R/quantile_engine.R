#' @title Mathematical Core and Quantile Modeling Engine
#'
#' @description
#' Low-level mathematical routines for rank-based inverse-normal quantile
#' modeling, orthogonal Probabilists' Hermite polynomial basis transformations,
#' analytical moment extraction via Isserlis' theorem, and monotonicity verification.
#'
#' @details
#' The functions in this module form the foundational mathematical backbone
#' of \code{hermiteStats}. They provide the machinery to transform arbitrary
#' continuous or finely graded empirical distributions into regularized,
#' continuous quantile functions \eqn{X = f(Z)} defined over the standard normal
#' domain \eqn{Z \sim \mathcal{N}(0, 1)}, and to derive their population moments
#' (\eqn{\mu, \sigma^2, \sigma}) in exact closed algebraic form.
#'
#' @name quantile_engine
#' @keywords internal
NULL


# -----------------------------------------------------------------------------
# Internal Helpers: Raw Normal Moments & Hermite Basis Operator
# -----------------------------------------------------------------------------

#' Compute Standard Normal Raw Moments E[Z^j]
#' @noRd
.z_raw_moments <- function(max_j) {
  moms <- numeric(max_j + 1L)
  moms[1L] <- 1.0  # j = 0
  if (max_j >= 2L) {
    for (j in seq(2L, max_j, by = 2L)) {
      moms[j + 1L] <- prod(seq(j - 1L, 1L, by = -2L))
    }
  }
  moms
}

#' Matrix Operator for Monomial-to-Hermite Basis Conversion
#' @noRd
.hermite_basis_matrix <- function(degree) {
  H_mat <- matrix(0.0, nrow = degree + 1L, ncol = degree + 1L)
  z_moms <- .z_raw_moments(2L * degree)
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

#' Compute Distributional Moments from Monomial Coefficients
#' @noRd
.extract_all_moments <- function(beta) {
  deg <- length(beta) - 1L
  if (deg < 0L) {
    return(list(mean = NA_real_, variance = NA_real_, sd = NA_real_,
                skewness = NA_real_, excess_kurtosis = NA_real_, kurtosis = NA_real_,
                hermite_coeffs = numeric(0)))
  }

  # Mean & Variance via Orthogonal Hermite Basis
  H <- .hermite_basis_matrix(deg)
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

  # Centered Polynomial for 3rd and 4th Central Moments
  beta_c <- beta
  beta_c[1L] <- beta_c[1L] - mean_val

  # Precompute moments up to power 4 * deg
  z_moms <- .z_raw_moments(4L * deg)

  # 3rd Central Moment (mu3)
  idx3 <- outer(outer(0:deg, 0:deg, `+`), 0:deg, `+`)
  coef3 <- outer(outer(beta_c, beta_c), beta_c)
  mu3 <- sum(coef3 * z_moms[idx3 + 1L])
  skew_val <- mu3 / (sd_val^3)

  # 4th Central Moment (mu4)
  idx4 <- outer(outer(outer(0:deg, 0:deg, `+`), 0:deg, `+`), 0:deg, `+`)
  coef4 <- outer(outer(outer(beta_c, beta_c), beta_c), beta_c)
  mu4 <- sum(coef4 * z_moms[idx4 + 1L])
  kurt_raw <- mu4 / (var_val^2)
  kurt_excess <- kurt_raw - 3.0

  list(
    mean = mean_val,
    variance = var_val,
    sd = sd_val,
    skewness = skew_val,
    excess_kurtosis = kurt_excess,
    kurtosis = kurt_raw,
    hermite_coeffs = a
  )
}

# -----------------------------------------------------------------------------
# Exported Function: check_monotonicity
# -----------------------------------------------------------------------------

#' Check Monotonicity of a Polynomial Quantile Function
#'
#' Evaluates whether a polynomial quantile function \eqn{X = f(Z)} is strictly or
#' empirically monotonically increasing over a specified standard normal domain
#' \eqn{[z_{\min}, z_{\max}]}.
#'
#' @param coeffs Numeric vector of monomial coefficients \eqn{(\beta_0, \beta_1, \dots, \beta_k)}
#'   representing the polynomial \eqn{f(z) = \sum_{j=0}^k \beta_j z^j}.
#' @param z_range Numeric vector of length 2 specifying the interval \code{c(min_z, max_z)}
#'   over which monotonicity is evaluated. Default is \code{c(-4, 4)}.
#' @param method Character string specifying the verification criterion:
#'   \describe{
#'     \item{\code{"relaxed"}}{(Default) Evaluates empirical rank concordance
#'       (\eqn{\rho_{\text{Spearman}} \ge 0.95}) between latent normal scores \eqn{z} and fitted
#'       values \eqn{\hat{x} = f(z)}, combined with a positive linear slope (\eqn{\beta_1 > 0}).}
#'     \item{\code{"strict"}}{Analytically computes the global minimum of the first derivative
#'       \eqn{f'(z)} over \code{z_range} by extracting the roots of \eqn{f''(z)} via
#'       \code{\link[stats]{polyroot}}. Requires \eqn{\min f'(z) \ge 0}.}
#'     \item{\code{"none"}}{Bypasses monotonicity testing; always returns \code{TRUE}.}
#'   }
#' @param z Optional numeric vector of observed or evaluated standard normal scores.
#'
#' @return A list with logical element \code{is_monotonic} and diagnostic parameters.
#' @export
check_monotonicity <- function(coeffs, z_range = c(-4, 4),
                               method = c("relaxed", "strict", "none"),
                               z = NULL) {
  method <- match.arg(method)
  coeffs[is.na(coeffs)] <- 0
  deg <- length(coeffs) - 1L

  if (method == "none" || deg <= 0L) {
    return(list(is_monotonic = TRUE, min_derivative = 0, rank_concordance = 1.0, method = method))
  }

  if (deg == 1L) {
    slope <- coeffs[2L]
    return(list(
      is_monotonic = slope >= 0,
      min_derivative = slope,
      rank_concordance = if (slope >= 0) 1.0 else -1.0,
      method = method
    ))
  }

  if (method == "strict") {
    d1 <- coeffs[-1L] * seq_len(deg)
    eval_d1 <- function(val) sum(d1 * (val^(0:(deg - 1L))))
    pts <- z_range
    if (deg >= 2L) {
      d2 <- d1[-1L] * seq_len(deg - 1L)
      if (any(d2 != 0)) {
        r <- tryCatch(stats::polyroot(d2), error = function(e) complex(0))
        real_r <- Re(r)[abs(Im(r)) < 1e-8]
        pts <- c(pts, real_r[real_r >= z_range[1L] & real_r <= z_range[2L]])
      }
    }
    slopes <- vapply(unique(pts), eval_d1, numeric(1L))
    min_slope <- min(slopes)
    return(list(
      is_monotonic = min_slope >= -1e-8,
      min_derivative = min_slope,
      rank_concordance = NA_real_,
      method = method
    ))
  }

  if (is.null(z)) {
    z <- seq(z_range[1L], z_range[2L], length.out = 100L)
  }
  Zmat <- outer(z, 0:deg, `^`)
  fitted_vals <- as.vector(Zmat %*% coeffs)

  rank_cor <- tryCatch(stats::cor(z, fitted_vals, method = "spearman"), error = function(e) 0)
  is_mono <- (!is.na(rank_cor)) && (rank_cor >= 0.95) && (coeffs[2L] > 0)

  list(
    is_monotonic = is_mono,
    min_derivative = NA_real_,
    rank_concordance = rank_cor,
    method = method
  )
}

# -----------------------------------------------------------------------------
# Exported Function: hermite_fit
# -----------------------------------------------------------------------------

#' Fit a Regularized Monotone Quantile Polynomial
#'
#' Maps empirical observations onto standard normal scores via rank-based
#' inverse-normal transformation (Normal Quantile Transformation; NQT) and fits
#' a regularized monotone polynomial quantile function \eqn{X = f(Z)}.
#' Orthogonal Probabilists' Hermite polynomial weights and closed-form population
#' moments (\eqn{\mu, \sigma^2, \sigma, \text{skewness } \gamma_1, \text{excess kurtosis } \gamma_2})
#' are automatically extracted.
#'
#' @param x Numeric vector of observations. Missing and infinite values are
#'   automatically removed.
#' @param degree Integer scalar; maximum polynomial degree (default is \code{3}).
#' @param monotonicity Monotonicity constraint: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Method for handling rank ties: \code{"average"} (default) or \code{"random"}.
#' @param force_odd Logical; if \code{TRUE} (default), restricts polynomial degree to odd values.
#'
#' @return An S3 object of class \code{"hermite_fit"} containing:
#' \describe{
#'   \item{\code{beta}}{Fitted monomial coefficients \eqn{(\beta_0, \dots, \beta_k)}.}
#'   \item{\code{hermite_coeffs}}{Orthogonal Probabilists' Hermite weights \eqn{(a_0, \dots, a_k)}.}
#'   \item{\code{degree}}{Realized polynomial degree.}
#'   \item{\code{degree_requested}}{Requested polynomial degree.}
#'   \item{\code{mean}}{Regularized population mean \eqn{\mu = a_0}.}
#'   \item{\code{variance}}{Regularized population variance \eqn{\sigma^2 = \sum_{m=1}^k a_m^2 m!}.}
#'   \item{\code{sd}}{Regularized population standard deviation \eqn{\sigma = \sqrt{\sigma^2}}.}
#'   \item{\code{skewness}}{Regularized population skewness \eqn{\gamma_1 = \mu_3 / \sigma^3}.}
#'   \item{\code{excess_kurtosis}}{Regularized population excess kurtosis \eqn{\gamma_2 = \mu_4 / \sigma^4 - 3}.}
#'   \item{\code{kurtosis}}{Regularized raw Pearson kurtosis \eqn{\beta_2 = \mu_4 / \sigma^4}.}
#'   \item{\code{n}}{Sample size.}
#'   \item{\code{n_unique}}{Number of unique observations.}
#'   \item{\code{tie_proportion}}{Proportion of tied values.}
#'   \item{\code{monotonicity}}{Monotonicity check applied.}
#'   \item{\code{ties_method}}{Tie-handling method applied.}
#'   \item{\code{x}}{Cleaned input observations.}
#'   \item{\code{z}}{Latent standard normal scores.}
#' }
#'
#' @references
#' Isserlis, L. (1918). On a formula for the product-moment coefficient of any order of a normal frequency distribution. \emph{Biometrika}, 12(1/2), 134–139. \doi{10.1093/biomet/12.1-2.134}
#'
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

  n_uniq <- length(unique(x_clean))
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

  cur_deg <- d_cap
  step <- if (force_odd) 2L else 1L
  beta_res <- NULL

  while (cur_deg >= 1L) {
    Zmat <- outer(z, 0:cur_deg, `^`)
    beta <- tryCatch(qr.solve(Zmat, x_clean), error = function(e) rep(NA_real_, cur_deg + 1L))

    if (!any(is.na(beta))) {
      chk <- check_monotonicity(beta, z_range = z_range, method = monotonicity, z = z)
      if (chk$is_monotonic) {
        beta_res <- beta
        break
      }
    }
    cur_deg <- cur_deg - step
  }

  if (is.null(beta_res)) {
    cur_deg <- 1L
    Zmat <- outer(z, 0:1, `^`)
    beta_res <- tryCatch(qr.solve(Zmat, x_clean), error = function(e) c(mean(x_clean), stats::sd(x_clean)))
  }

  # Extract all regularized moments
  moms <- .extract_all_moments(beta_res)

  res <- list(
    beta = beta_res,
    hermite_coeffs = moms$hermite_coeffs,
    degree = cur_deg,
    degree_requested = degree,
    mean = moms$mean,
    variance = moms$variance,
    sd = moms$sd,
    skewness = moms$skewness,
    excess_kurtosis = moms$excess_kurtosis,
    kurtosis = moms$kurtosis,
    n = n,
    n_unique = n_uniq,
    tie_proportion = tie_prop,
    monotonicity = monotonicity,
    ties_method = ties_method,
    x = x_clean,
    z = z
  )
  class(res) <- "hermite_fit"
  res
}

# -----------------------------------------------------------------------------
# Exported Function: hermite_moments
# -----------------------------------------------------------------------------

#' Extract Regularized Distributional Moments
#'
#' Retrieves the regularized population moments (\eqn{\mu, \sigma^2, \sigma, \text{skewness } \gamma_1, \text{excess kurtosis } \gamma_2})
#' derived from a fitted \code{hermite_fit} object.
#'
#' @param object An object of class \code{"hermite_fit"} returned by \code{\link{hermite_fit}}.
#'
#' @return A named list containing \code{mean}, \code{variance}, \code{sd},
#'   \code{skewness}, \code{excess_kurtosis}, and \code{kurtosis}.
#' @export
hermite_moments <- function(object) {
  if (!inherits(object, "hermite_fit")) {
    stop("Argument 'object' must be an object of class 'hermite_fit'.")
  }
  list(
    mean = object$mean,
    variance = object$variance,
    sd = object$sd,
    skewness = object$skewness,
    excess_kurtosis = object$excess_kurtosis,
    kurtosis = object$kurtosis
  )
}

# -----------------------------------------------------------------------------
# S3 Methods: print and summary
# -----------------------------------------------------------------------------

#' @export
print.hermite_fit <- function(x, digits = 3L, ...) {
  cat("\n  Regularized Quantile Model (Hermite Basis)\n")
  cat(strrep("-", 45), "\n", sep = "")
  cat(sprintf("  Modeled Mean (mu)      :  %.*f\n", digits, x$mean))
  cat(sprintf("  Modeled SD (sigma)     :  %.*f\n", digits, x$sd))
  cat(sprintf("  Modeled Variance       :  %.*f\n", digits, x$variance))
  cat(sprintf("  Modeled Skewness (g1)  :  %.*f\n", digits, x$skewness))
  cat(sprintf("  Excess Kurtosis (g2)   :  %.*f\n", digits, x$excess_kurtosis))
  cat(sprintf("  Fitted Polynomial Deg  :  %d (requested: %d)\n", x$degree, x$degree_requested))
  cat(sprintf("  Sample Size (n)        :  %d (unique: %d, ties: %.1f%%)\n",
              x$n, x$n_unique, x$tie_proportion * 100))
  cat(sprintf("  Monotonicity Check     :  %s\n", x$monotonicity))
  cat("\n")
  invisible(x)
}

#' @export
summary.hermite_fit <- function(object, digits = 3L, ...) {
  print(object, digits = digits, ...)
  cat("  Polynomial Coefficients (Monomial raw beta):\n")
  beta_names <- paste0("z^", 0:object$degree)
  names(object$beta) <- beta_names
  print(round(object$beta, digits = digits))

  cat("\n  Orthogonal Hermite Weights (a_m):\n")
  herm_names <- paste0("He_", 0:object$degree)
  names(object$hermite_coeffs) <- herm_names
  print(round(object$hermite_coeffs, digits = digits))
  cat("\n")
  invisible(object)
}

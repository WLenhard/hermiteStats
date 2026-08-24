#' @title Internal Core Engine for Hermite-Mehler Modeling
#' @description Low-level mathematical functions for moment generation, basis conversion,
#'   monotonicity checking, and single-distribution quantile fitting.
#' @keywords internal
NULL

# Standard normal raw moments E[Z^j]
.z_raw_moments <- function(max_j) {
  moms <- numeric(max_j + 1L)
  moms[1L] <- 1.0
  if (max_j >= 2L) {
    for (j in seq(2L, max_j, by = 2L)) {
      moms[j + 1L] <- prod(seq(j - 1L, 1L, by = -2L))
    }
  }
  moms
}

# Basis transformation matrix: Monomial coefficients beta -> Hermite coefficients a
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

#' Check Monotonicity of a Polynomial Function
#'
#' Checks whether a fitted polynomial quantile function is strictly or empirically
#' monotonic over the observed standard normal domain.
#'
#' @param coeffs Numeric vector of polynomial coefficients (from degree 0 to k).
#' @param z_range Numeric vector of length 2 defining \code{c(min_z, max_z)}.
#' @param method Character; \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param z Optional numeric vector of sample z-scores for empirical checks.
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
    return(list(is_monotonic = TRUE, min_derivative = 0, method = method))
  }

  if (deg == 1L) {
    slope <- coeffs[2L]
    return(list(is_monotonic = slope >= 0, min_derivative = slope, method = method))
  }

  # Method 1: Strict Analytical Check (Root finding on derivative)
  if (method == "strict") {
    d1 <- coeffs[-1L] * seq_len(deg)
    eval_d1 <- function(val) sum(d1 * (val^(0:(deg - 1L))))
    pts <- z_range
    if (deg >= 2L) {
      d2 <- d1[-1L] * seq_len(deg - 1L)
      if (any(d2 != 0)) {
        r <- tryCatch(polyroot(d2), error = function(e) complex(0))
        real_r <- Re(r)[abs(Im(r)) < 1e-8]
        pts <- c(pts, real_r[real_r >= z_range[1L] & real_r <= z_range[2L]])
      }
    }
    slopes <- vapply(unique(pts), eval_d1, numeric(1L))
    min_slope <- min(slopes)
    return(list(
      is_monotonic = min_slope >= -1e-8,
      min_derivative = min_slope,
      method = method
    ))
  }

  # Method 2: Relaxed Rank-Concordance Check
  if (is.null(z)) {
    z <- seq(z_range[1L], z_range[2L], length.out = 100L)
  }
  Zmat <- outer(z, 0:deg, `^`)
  fitted_vals <- as.vector(Zmat %*% coeffs)

  rank_cor <- tryCatch(stats::cor(z, fitted_vals, method = "spearman"), error = function(e) 0)
  is_mono <- (!is.na(rank_cor)) && (rank_cor >= 0.95) && (coeffs[2L] > 0)

  list(
    is_monotonic = is_mono,
    rank_concordance = rank_cor,
    min_derivative = NA_real_,
    method = method
  )
}

#' Fit Monotone Quantile Polynomial to a Single Variable
#'
#' Maps observations to standard normal scores via rank-based inverse-normal
#' transformation and models the quantile function \eqn{X = f(Z)} using
#' a regularized monotone polynomial.
#'
#' @param x Numeric vector of observations.
#' @param degree Integer; requested polynomial degree (default = 3).
#' @param monotonicity Character; \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Method for rank ties: \code{"average"} (default) or \code{"random"}.
#' @param force_odd Logical; if \code{TRUE} (default), restricts stepping down to odd degrees.
#'
#' @return An S3 object of class \code{"hermite_fit"}.
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

  # Max allowable degree given discrete sample properties
  d_cap <- min(as.integer(degree), n_uniq - 1L)
  if (tie_prop > 0.30 && d_cap > 3L) {
    d_cap <- min(d_cap, max(3L, floor(n_uniq / 2L)))
  }
  if (force_odd && (d_cap %% 2L == 0L)) d_cap <- d_cap - 1L
  d_cap <- max(1L, d_cap)

  # Rank-based standard normal scores
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

  # Extract Hermite moments
  H <- .hermite_basis_matrix(cur_deg)
  a <- as.vector(H %*% beta_res)

  mean_val <- a[1L]
  var_val  <- if (cur_deg == 0L) 0.0 else sum((a[-1L]^2) * factorial(1L:cur_deg))

  res <- list(
    beta = beta_res,
    hermite_coeffs = a,
    degree = cur_deg,
    degree_requested = degree,
    mean = mean_val,
    variance = max(0.0, var_val),
    sd = sqrt(max(0.0, var_val)),
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

#' Extract Moments from a hermite_fit Object
#'
#' @param object An object of class \code{"hermite_fit"}.
#' @return A list with elements \code{mean}, \code{variance}, and \code{sd}.
#' @export
hermite_moments <- function(object) {
  if (!inherits(object, "hermite_fit")) stop("'object' must be of class 'hermite_fit'.")
  list(mean = object$mean, variance = object$variance, sd = object$sd)
}

#' Print Method for hermite_fit Objects
#'
#' @param x An object of class \code{"hermite_fit"}.
#' @param digits Integer; number of decimal places to display (default = 3).
#' @param ... Additional arguments.
#' @export
print.hermite_fit <- function(x, digits = 3L, ...) {
  cat("\n  Regularized Quantile Model (Hermite Basis)\n")
  cat(strrep("-", 45), "\n", sep = "")
  cat(sprintf("  Modeled Mean (mu)      :  %.*f\n", digits, x$mean))
  cat(sprintf("  Modeled SD (sigma)     :  %.*f\n", digits, x$sd))
  cat(sprintf("  Modeled Variance       :  %.*f\n", digits, x$variance))
  cat(sprintf("  Fitted Polynomial Deg  :  %d (requested: %d)\n", x$degree, x$degree_requested))
  cat(sprintf("  Sample Size (n)        :  %d (unique: %d, ties: %.1f%%)\n",
              x$n, x$n_unique, x$tie_proportion * 100))
  cat(sprintf("  Monotonicity Check     :  %s\n", x$monotonicity))
  cat("\n")
  invisible(x)
}

#' Summary Method for hermite_fit Objects
#'
#' @param object An object of class \code{"hermite_fit"}.
#' @param digits Integer; number of decimal places to display (default = 3).
#' @param ... Additional arguments.
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

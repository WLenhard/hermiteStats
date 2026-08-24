#' @title Hermite-Mehler Distribution-Robust Pearson Correlation
#'
#' @description Computes the distribution-robust Hermite-Mehler Pearson correlation (\eqn{r_{\mathrm{HM}}}),
#' the latent Gaussian copula correlation (\eqn{\rho_z}), and the shape attenuation factor (\eqn{A}).
#'
#' @param x A numeric vector, matrix, or data.frame.
#' @param y A numeric vector (if \code{x} is a vector).
#' @param poly_degree Integer; maximum polynomial degree (default = 3, odd values only).
#' @param monotonicity Monotonicity criterion: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Tie handling: \code{"average"} (default) or \code{"random"}.
#' @param conf_level Numeric in (0, 1) specifying the confidence level, or \code{NULL} for no CI.
#' @param ci_method Confidence interval method: \code{"fisher"} (analytical Fisher z) or \code{"bootstrap"}.
#' @param B Number of bootstrap resamples (default = 1000).
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"cor_hermite"} (or \code{"cor_hermite_matrix"}).
#' @export
cor_hermite <- function(x, ...) {
  UseMethod("cor_hermite")
}

#' @rdname cor_hermite
#' @export
cor_hermite.default <- function(x, y, poly_degree = 3L,
                                monotonicity = c("relaxed", "strict", "none"),
                                ties_method = c("average", "random"),
                                conf_level = NULL,
                                ci_method = c("fisher", "bootstrap"),
                                B = 1000L, ...) {

  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)
  ci_method    <- match.arg(ci_method)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied for two-vector correlation.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (length(x) != length(y)) stop("'x' and 'y' must have identical length.")

  ok <- is.finite(x) & is.finite(y)
  x_c <- x[ok]; y_c <- y[ok]
  n <- length(x_c)

  if (n < 4L) {
    warning("Insufficient complete observations.")
    return(structure(list(r_hm = NA_real_, rho_z = NA_real_, attenuation = NA_real_), class = "cor_hermite"))
  }

  # Fit marginal polynomial quantile functions
  fit_x <- hermite_fit(x_c, degree = poly_degree, monotonicity = monotonicity, ties_method = ties_method, force_odd = TRUE)
  fit_y <- hermite_fit(y_c, degree = poly_degree, monotonicity = monotonicity, ties_method = ties_method, force_odd = TRUE)

  # Latent Copula Correlation
  rho_z <- stats::cor(fit_x$z, fit_y$z)
  rho_z <- max(-1.0, min(1.0, rho_z))

  # Mehler Bilinear Covariance Integration
  dx <- fit_x$degree; dy <- fit_y$degree
  max_d <- max(dx, dy)

  H <- .hermite_basis_matrix(max_d)
  bx_pad <- rep(0.0, max_d + 1L); bx_pad[1:(dx + 1L)] <- fit_x$beta
  by_pad <- rep(0.0, max_d + 1L); by_pad[1:(dy + 1L)] <- fit_y$beta

  a <- as.vector(H %*% bx_pad)
  b <- as.vector(H %*% by_pad)
  fact_v <- factorial(1L:max_d)

  var_x <- sum((a[-1L]^2) * fact_v)
  var_y <- sum((b[-1L]^2) * fact_v)
  cov_xy <- sum(a[-1L] * b[-1L] * fact_v * (rho_z^(1L:max_d)))

  r_hm <- if (var_x <= 0 || var_y <= 0) NA_real_ else max(-1.0, min(1.0, cov_xy / sqrt(var_x * var_y)))
  attenuation <- if (abs(rho_z) > 1e-5) r_hm / rho_z else 1.0

  res <- list(
    r_hm = r_hm,
    rho_z = rho_z,
    attenuation = attenuation,
    cov_xy = cov_xy,
    var_x = var_x, var_y = var_y,
    mean_x = a[1L], mean_y = b[1L],
    degrees = c(x = dx, y = dy),
    poly_degree_requested = poly_degree,
    monotonicity = monotonicity,
    ties_method = ties_method,
    n = n,
    fit_x = fit_x,
    fit_y = fit_y,
    x = x_c, y = y_c
  )
  class(res) <- "cor_hermite"

  if (!is.null(conf_level)) {
    res$ci <- stats::confint(res, level = conf_level, method = ci_method, B = B)
    res$conf_level <- conf_level
    res$ci_method <- ci_method
  }
  res
}

#' @rdname cor_hermite
#' @export
cor_hermite.matrix <- function(x, poly_degree = 3L,
                               monotonicity = c("relaxed", "strict", "none"),
                               ties_method = c("average", "random"), ...) {
  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  p <- ncol(x)
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(p))

  R_hm <- matrix(1.0, nrow = p, ncol = p, dimnames = list(cnames, cnames))
  R_copula <- matrix(1.0, nrow = p, ncol = p, dimnames = list(cnames, cnames))
  Att_mat <- matrix(1.0, nrow = p, ncol = p, dimnames = list(cnames, cnames))

  for (i in 1:(p - 1L)) {
    for (j in (i + 1L):p) {
      fit <- cor_hermite.default(x[, i], x[, j], poly_degree = poly_degree,
                                 monotonicity = monotonicity, ties_method = ties_method)
      R_hm[i, j] <- R_hm[j, i] <- fit$r_hm
      R_copula[i, j] <- R_copula[j, i] <- fit$rho_z
      Att_mat[i, j] <- Att_mat[j, i] <- fit$attenuation
    }
  }

  structure(
    list(r_hm = R_hm, rho_z = R_copula, attenuation = Att_mat, n = nrow(x)),
    class = "cor_hermite_matrix"
  )
}

#' @rdname cor_hermite
#' @export
cor_hermite.data.frame <- function(x, poly_degree = 3L,
                                   monotonicity = c("relaxed", "strict", "none"),
                                   ties_method = c("average", "random"), ...) {
  cor_hermite.matrix(as.matrix(x), poly_degree = poly_degree,
                     monotonicity = monotonicity, ties_method = ties_method, ...)
}

#' Confidence Intervals for Hermite Correlation
#' @export
confint.cor_hermite <- function(object, parm, level = 0.95,
                                method = c("fisher", "bootstrap"), B = 1000L, ...) {
  method <- match.arg(method)
  r <- object$r_hm
  n <- object$n

  if (method == "fisher") {
    z_r <- atanh(r)
    se <- 1 / sqrt(n - 3L)
    crit <- stats::qnorm((1 + level) / 2)
    ci <- tanh(z_r + c(-1, 1) * crit * se)
  } else {
    boot_vals <- numeric(B)
    x <- object$x; y <- object$y
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      fb <- cor_hermite(x[idx], y[idx], poly_degree = object$poly_degree_requested,
                        monotonicity = object$monotonicity, ties_method = object$ties_method)
      boot_vals[b] <- fb$r_hm
    }
    alpha <- (1 - level) / 2
    ci <- stats::quantile(boot_vals, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  }
  matrix(ci, nrow = 1L, dimnames = list("r_HM", c(paste0(100 * (1 - level)/2, " %"),
                                                  paste0(100 * (1 + level)/2, " %"))))
}

#' @export
print.cor_hermite <- function(x, digits = 3L, ...) {
  cat("\n  Hermite-Mehler Pearson Correlation (r_HM)\n")
  cat(strrep("-", 48), "\n", sep = "")
  cat(sprintf("  Hermite Correlation (r_HM) :  %.*f\n", digits, x$r_hm))
  cat(sprintf("  Latent Copula (rho_z)      :  %.*f\n", digits, x$rho_z))
  cat(sprintf("  Shape Attenuation (A)      :  %.*f\n", digits, x$attenuation))
  cat(sprintf("  Polynomial Degrees Fitted  :  X = %d, Y = %d\n", x$degrees["x"], x$degrees["y"]))
  cat(sprintf("  Monotonicity Check         :  %s\n", x$monotonicity))
  if (!is.null(x$ci)) {
    cat(sprintf("  %s CI (%s): [%.*f, %.*f]\n",
                paste0(round(x$conf_level * 100), "%"), x$ci_method,
                digits, x$ci[1L], digits, x$ci[2L]))
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.cor_hermite <- function(object, digits = 3L, ...) {
  print(object, digits = digits, ...)
  cat("  Distributional Moments (Implied):\n")
  cat(sprintf("    X: Mean = %.*f, SD = %.*f\n", digits, object$mean_x, digits, sqrt(object$var_x)))
  cat(sprintf("    Y: Mean = %.*f, SD = %.*f\n", digits, object$mean_y, digits, sqrt(object$var_y)))
  cat(sprintf("    Covariance(X, Y) = %.*f\n\n", digits, object$cov_xy))
  invisible(object)
}

#' @export
print.cor_hermite_matrix <- function(x, digits = 3L, ...) {
  cat("\n  Hermite-Mehler Correlation Matrix (r_HM):\n")
  print(round(x$r_hm, digits = digits))
  cat("\n  Latent Copula Correlation Matrix (rho_z):\n")
  print(round(x$rho_z, digits = digits))
  cat("\n  Shape Attenuation Factor Matrix (A = r_HM / rho_z):\n")
  print(round(x$attenuation, digits = digits))
  cat("\n")
  invisible(x)
}

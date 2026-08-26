#' @title Hermite-Mehler Distribution-Robust Pearson Correlation
#'
#' @description
#' Computes the distribution-robust Hermite-Mehler Pearson correlation (\eqn{r_{\mathrm{HM}}}),
#' the latent Gaussian copula correlation (\eqn{\rho_z}), and the shape-attenuation
#' factor (\eqn{A}) for a pair of variables, or for every pair of columns in a
#' matrix or data frame.
#'
#' Unlike Spearman's \eqn{\rho}, Kendall's \eqn{\tau}, or Gaussian rank correlation,
#' \eqn{r_{\mathrm{HM}}} targets the same estimand as the classical Pearson
#' correlation of the raw, manifest variables. Unlike trimming or Winsorizing,
#' it does not discard authentic tail variance to achieve stability.
#'
#' @param x A numeric vector, numeric matrix, or data frame.
#' @param y A numeric vector, required when \code{x} is a vector; ignored
#'   when \code{x} is a matrix or data frame.
#' @param poly_degree Integer scalar; maximum polynomial degree used for the
#'   marginal quantile fits (default \code{3L}). Odd degrees (1, 3, 5) are
#'   recommended; see \code{\link{hermite_fit}}.
#' @param monotonicity Monotonicity constraint passed to
#'   \code{\link{hermite_fit}}: \code{"relaxed"} (default), \code{"strict"},
#'   or \code{"none"}.
#' @param ties_method Rank tie-handling method passed to
#'   \code{\link{hermite_fit}}: \code{"average"} (default) or \code{"random"}.
#' @param conf_level Numeric value in \eqn{(0, 1)}, e.g. \code{0.95}, or
#'   \code{NULL} (default) to skip confidence interval calculation. Only used
#'   for the two-vector method.
#' @param ci_method Confidence interval method for the two-vector case:
#'   \describe{
#'     \item{\code{"fisher"}}{(Default) Analytical Fisher \eqn{z}-transformation,
#'       treating \eqn{r_{\mathrm{HM}}} like a Pearson correlation with standard
#'       error \eqn{1/\sqrt{n-3}}.}
#'     \item{\code{"bootstrap"}}{Non-parametric percentile bootstrap that
#'       reruns the entire Hermite-Mehler pipeline (NQT, polynomial refit,
#'       Mehler covariance) on each resample.}
#'   }
#' @param B Integer; number of bootstrap resamples when \code{ci_method = "bootstrap"}
#'   (default \code{1000L}).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \subsection{The Hermite-Mehler pipeline}{
#' \enumerate{
#'   \item Both variables are mapped to standard normal scores \eqn{Z_x, Z_y}
#'         via rank-based inverse-normal transformation. Their sample correlation
#'         is the latent Gaussian copula parameter \eqn{\rho_z}.
#'   \item The inverse quantile functions \eqn{X = f(Z_x)} and \eqn{Y = g(Z_y)}
#'         are each approximated by a low-degree monotone polynomial
#'         (\code{\link{hermite_fit}}).
#'   \item The polynomial coefficients are re-expressed in the orthogonal
#'         Probabilists' Hermite basis. Assuming \eqn{(Z_x, Z_y)} is bivariate
#'         normal with correlation \eqn{\rho_z}, Mehler's (1866) bilinear
#'         expansion gives
#'         \deqn{E\left[He_m(Z_x) He_n(Z_y)\right] = \delta_{mn}\, m!\, \rho_z^m,}
#'         so all cross-order terms vanish and the manifest covariance follows
#'         in exact closed form:
#'         \deqn{\mathrm{Cov}(X,Y) = \sum_{m=1}^{\min(k_x,k_y)} a_m b_m\, m!\, \rho_z^m,}
#'         \deqn{\mathrm{Var}(X) = \sum_{m=1}^{k_x} a_m^2\, m!, \qquad
#'               \mathrm{Var}(Y) = \sum_{m=1}^{k_y} b_m^2\, m!,}
#'         \deqn{r_{\mathrm{HM}} = \frac{\mathrm{Cov}(X,Y)}{\sqrt{\mathrm{Var}(X)\,\mathrm{Var}(Y)}}.}
#' }
#' }
#'
#' \subsection{Shape attenuation factor A}{
#' When the two marginal shapes differ substantially (e.g. a roughly symmetric
#' variable paired with a strongly skewed one), the Hoeffding-Frechet bounds
#' imply that the achievable Pearson correlation is strictly below 1 in
#' absolute value. The ratio
#' \deqn{A = \frac{r_{\mathrm{HM}}}{\rho_z}}
#' quantifies this ceiling.
#' }
#'
#' \subsection{Matrices and missing data}{
#' For matrix or data-frame input, each column's marginal quantile model is
#' fitted once and cached for every pair it appears in, provided the input
#' contains no missing values. If missing values are present, \code{cor_hermite}
#' falls back to pairwise-complete-case fitting for every pair (matching
#' \code{stats::cor(..., use = "pairwise.complete.obs")}).
#' }
#'
#' @return
#' For two vectors, an S3 object of class \code{"cor_hermite"} containing:
#' \describe{
#'   \item{\code{r_hm}}{The Hermite-Mehler Pearson correlation.}
#'   \item{\code{rho_z}}{The latent Gaussian copula correlation.}
#'   \item{\code{attenuation}}{The shape attenuation factor \eqn{A = r_{\mathrm{HM}} / \rho_z}.}
#'   \item{\code{cov_xy}}{The regularized manifest covariance.}
#'   \item{\code{var_x, var_y}}{Regularized manifest variances of \eqn{X} and \eqn{Y}.}
#'   \item{\code{mean_x, mean_y}}{Regularized manifest means of \eqn{X} and \eqn{Y}.}
#'   \item{\code{degrees}}{Named integer vector of realized polynomial degrees.}
#'   \item{\code{poly_degree_requested}}{Requested maximum polynomial degree.}
#'   \item{\code{monotonicity}, \code{ties_method}}{Methods applied.}
#'   \item{\code{n}}{Number of complete paired observations.}
#'   \item{\code{fit_x}, \code{fit_y}}{The underlying \code{\link{hermite_fit}} objects.}
#'   \item{\code{ci}, \code{conf_level}, \code{ci_method}}{Confidence interval results, if requested.}
#' }
#'
#' For a matrix or data frame, an S3 object of class \code{"cor_hermite_matrix"} containing:
#' \describe{
#'   \item{\code{r_hm}}{The \eqn{p \times p} regularized Hermite-Mehler correlation matrix.}
#'   \item{\code{rho_z}}{The \eqn{p \times p} latent Gaussian copula correlation matrix.}
#'   \item{\code{attenuation}}{The \eqn{p \times p} shape-attenuation factor matrix.}
#'   \item{\code{cov}}{The \eqn{p \times p} regularized covariance matrix (with variances on the diagonal).}
#'   \item{\code{marginals}}{A data frame containing regularized univariate summary statistics (Mean, Variance, SD, Degree) for each column.}
#'   \item{\code{n}}{Sample size.}
#' }
#'
#' @references
#' Carroll, J. B. (1961). The nature of the data, or how to choose a correlation coefficient. \emph{Psychometrika}, 26(4), 347-372. \doi{10.1007/BF02289768}
#'
#' Hermite, C. (1864). Sur un nouveau développement en série des fonctions. \emph{Comptes Rendus de l'Académie des Sciences, Paris}, 58, 93-100.
#'
#' Lenhard, W., & Lenhard, A. (in preparation). The Hermite-Mehler Correlation: A Distribution-Robust Estimator of the Pearson Correlation Coefficient.
#'
#' Mehler, F. G. (1866). Ueber die Entwicklung einer Function von beliebig vielen Variabeln nach Laplaceschen Functionen höherer Ordnung. \emph{Journal für die reine und angewandte Mathematik}, 66, 161-176. \doi{10.1515/crll.1866.66.161}
#'
#' @seealso \code{\link{d_reg}}, \code{\link{hermite_fit}}, \code{\link{check_monotonicity}}
#'
#' @examples
#' # 1. Bivariate correlation
#' set.seed(123)
#' x <- rnorm(50, mean = 100, sd = 15)
#' y <- exp(0.5 * scale(x) + rnorm(50, sd = 0.5))
#' r_fit <- cor_hermite(x, y, conf_level = 0.95)
#' print(r_fit)
#'
#' # 2. Matrix input
#' data(iris)
#' r_mat <- cor_hermite(iris[, 1:4])
#' print(r_mat)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
cor_hermite <- function(x, ...) {
  UseMethod("cor_hermite")
}

# -----------------------------------------------------------------------------
# Internal Helper: Mehler bilinear covariance for two already-fitted quantile models
# -----------------------------------------------------------------------------

#' Mehler bilinear covariance for two already-fitted quantile models
#' @noRd
.cor_hermite_pair <- function(fit_x, fit_y) {
  rho_z <- stats::cor(fit_x$z, fit_y$z)
  if (!is.finite(rho_z)) return(list(
    r_hm = NA_real_, rho_z = NA_real_, attenuation = NA_real_,
    cov_xy = NA_real_, var_x = fit_x$variance, var_y = fit_y$variance,
    mean_x = fit_x$mean, mean_y = fit_y$mean,
    degrees = c(x = fit_x$degree, y = fit_y$degree)
  ))

  rho_z <- max(-1.0, min(1.0, rho_z))

  dx <- fit_x$degree; dy <- fit_y$degree
  max_d <- max(dx, dy)

  H <- .hermite_basis_matrix(max_d)
  bx_pad <- rep(0.0, max_d + 1L); bx_pad[1:(dx + 1L)] <- fit_x$beta
  by_pad <- rep(0.0, max_d + 1L); by_pad[1:(dy + 1L)] <- fit_y$beta

  a <- as.vector(H %*% bx_pad)
  b <- as.vector(H %*% by_pad)
  fact_v <- factorial(1L:max_d)

  cov_xy <- sum(a[-1L] * b[-1L] * fact_v * (rho_z^(1L:max_d)))
  var_x <- fit_x$variance
  var_y <- fit_y$variance

  r_hm <- if (var_x <= 0 || var_y <= 0) NA_real_ else max(-1.0, min(1.0, cov_xy / sqrt(var_x * var_y)))
  attenuation <- if (abs(rho_z) > 1e-5) r_hm / rho_z else 1.0

  list(
    r_hm = r_hm,
    rho_z = rho_z,
    attenuation = attenuation,
    cov_xy = cov_xy,
    var_x = var_x,
    var_y = var_y,
    mean_x = fit_x$mean,
    mean_y = fit_y$mean,
    degrees = c(x = dx, y = dy)
  )
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

  fit_x <- hermite_fit(x_c, degree = poly_degree, monotonicity = monotonicity, ties_method = ties_method, force_odd = TRUE)
  fit_y <- hermite_fit(y_c, degree = poly_degree, monotonicity = monotonicity, ties_method = ties_method, force_odd = TRUE)

  pair <- .cor_hermite_pair(fit_x, fit_y)

  res <- list(
    r_hm = pair$r_hm,
    rho_z = pair$rho_z,
    attenuation = pair$attenuation,
    cov_xy = pair$cov_xy,
    var_x = pair$var_x, var_y = pair$var_y,
    mean_x = pair$mean_x, mean_y = pair$mean_y,
    degrees = pair$degrees,
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

#' Assemble a cor_hermite result from two already-fitted quantile models
#' @noRd
.cor_hermite_assemble <- function(fit_x, fit_y, poly_degree_requested,
                                  monotonicity, ties_method) {
  pair <- .cor_hermite_pair(fit_x, fit_y)
  res <- list(
    r_hm = pair$r_hm,
    rho_z = pair$rho_z,
    attenuation = pair$attenuation,
    cov_xy = pair$cov_xy,
    var_x = pair$var_x, var_y = pair$var_y,
    mean_x = pair$mean_x, mean_y = pair$mean_y,
    degrees = pair$degrees,
    poly_degree_requested = poly_degree_requested,
    monotonicity = monotonicity,
    ties_method = ties_method,
    n = length(fit_x$x),
    fit_x = fit_x,
    fit_y = fit_y,
    x = fit_x$x, y = fit_y$x
  )
  class(res) <- "cor_hermite"
  res
}

#' @rdname cor_hermite
#' @export
cor_hermite.matrix <- function(x, poly_degree = 3L,
                               monotonicity = c("relaxed", "strict", "none"),
                               ties_method = c("average", "random"), ...) {
  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  if (!is.numeric(x)) stop("'x' must be a numeric matrix.")
  p <- ncol(x)
  if (is.null(p) || p < 2L) stop("'x' must contain at least two numeric columns.")

  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(p))

  # Matrices initialized with proper 1.0 diagonals for correlations and attenuation
  R_hm     <- diag(1.0, nrow = p, ncol = p)
  R_copula <- diag(1.0, nrow = p, ncol = p)
  Att_mat  <- diag(1.0, nrow = p, ncol = p)
  cov_hm   <- matrix(NA_real_, nrow = p, ncol = p)

  dimnames(R_hm)     <- list(cnames, cnames)
  dimnames(R_copula) <- list(cnames, cnames)
  dimnames(Att_mat)  <- list(cnames, cnames)
  dimnames(cov_hm)   <- list(cnames, cnames)

  fully_complete <- all(is.finite(x))

  if (fully_complete) {
    # Fit each column once and cache
    fits <- vector("list", p)
    means_v <- numeric(p)
    vars_v  <- numeric(p)
    sds_v   <- numeric(p)
    skew_v   <- numeric(p)
    degs_v  <- integer(p)

    for (i in seq_len(p)) {
      fits[[i]] <- hermite_fit(x[, i], degree = poly_degree, monotonicity = monotonicity,
                               ties_method = ties_method, force_odd = TRUE)
      means_v[i] <- fits[[i]]$mean
      vars_v[i]  <- fits[[i]]$variance
      sds_v[i]   <- fits[[i]]$sd
      skew_v[i]   <- fits[[i]]$skewness
      degs_v[i]  <- fits[[i]]$degree
    }

    # Diagonal of covariance matrix = variances
    diag(cov_hm) <- vars_v

    # Off-diagonals
    for (i in 1:(p - 1L)) {
      for (j in (i + 1L):p) {
        pair <- .cor_hermite_pair(fits[[i]], fits[[j]])
        R_hm[i, j]     <- R_hm[j, i]     <- pair$r_hm
        R_copula[i, j] <- R_copula[j, i] <- pair$rho_z
        Att_mat[i, j]  <- Att_mat[j, i]  <- pair$attenuation
        cov_hm[i, j]   <- cov_hm[j, i]   <- pair$cov_xy
      }
    }
  } else {
    # Pairwise complete fallback
    means_v <- numeric(p)
    vars_v  <- numeric(p)
    sds_v   <- numeric(p)
    skew_v   <- numeric(p)
    degs_v  <- integer(p)

    for (i in seq_len(p)) {
      xi <- x[is.finite(x[, i]), i]
      fit_i <- hermite_fit(xi, degree = poly_degree, monotonicity = monotonicity,
                           ties_method = ties_method, force_odd = TRUE)
      means_v[i] <- fit_i$mean
      vars_v[i]  <- fit_i$variance
      sds_v[i]   <- fit_i$sd
      skew_v[i]   <- fit_i$skewness
      degs_v[i]  <- fit_i$degree
    }

    diag(cov_hm) <- vars_v

    for (i in 1:(p - 1L)) {
      for (j in (i + 1L):p) {
        fit <- cor_hermite.default(x[, i], x[, j], poly_degree = poly_degree,
                                   monotonicity = monotonicity, ties_method = ties_method)
        R_hm[i, j]     <- R_hm[j, i]     <- fit$r_hm
        R_copula[i, j] <- R_copula[j, i] <- fit$rho_z
        Att_mat[i, j]  <- Att_mat[j, i]  <- fit$attenuation
        cov_hm[i, j]   <- cov_hm[j, i]   <- fit$cov_xy
      }
    }
  }

  marginals_df <- data.frame(
    Variable = cnames,
    Mean     = means_v,
    Variance = vars_v,
    SD       = sds_v,
    Skewness = skew_v,
    Degree   = degs_v,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      r_hm        = R_hm,
      rho_z       = R_copula,
      attenuation = Att_mat,
      cov         = cov_hm,
      marginals   = marginals_df,
      n           = nrow(x)
    ),
    class = "cor_hermite_matrix"
  )
}

#' @rdname cor_hermite
#' @export
cor_hermite.data.frame <- function(x, poly_degree = 3L,
                                   monotonicity = c("relaxed", "strict", "none"),
                                   ties_method = c("average", "random"), ...) {
  is_num <- vapply(x, is.numeric, logical(1L))
  if (!all(is_num)) {
    stop("All columns of 'x' must be numeric. Non-numeric columns found: ",
         paste(names(x)[!is_num], collapse = ", "))
  }
  cor_hermite.matrix(as.matrix(x), poly_degree = poly_degree,
                     monotonicity = monotonicity, ties_method = ties_method, ...)
}

# -----------------------------------------------------------------------------
# S3 Methods for cor_hermite: confint, print, summary
# -----------------------------------------------------------------------------

#' Confidence Intervals for the Hermite-Mehler Correlation
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
  cat(sprintf("    X: Mean = %.*f, SD = %.*f, Var = %.*f\n", digits, object$mean_x, digits, sqrt(object$var_x), digits, object$var_x))
  cat(sprintf("    Y: Mean = %.*f, SD = %.*f, Var = %.*f\n", digits, object$mean_y, digits, sqrt(object$var_y), digits, object$var_y))
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

  cat("\n  Hermite-Mehler Covariance Matrix (cov_HM; Diagonale: Variances):\n")
  print(round(x$cov, digits = digits))

  cat("\n  Regularized Marginal Moments:\n")
  marg <- x$marginals
  marg$Mean     <- round(marg$Mean, digits = digits)
  marg$Variance <- round(marg$Variance, digits = digits)
  marg$SD       <- round(marg$SD, digits = digits)
  marg$Skewness <- round(marg$Skewness, digits = digits)
  print(marg, row.names = FALSE)
  cat("\n")

  invisible(x)
}

#' @title Distribution-Robust Pearson Correlation based on Probabilists' Hermite Polynomials
#'
#' @description
#' Computes the distribution-robust Hermite correlation (\eqn{r_{\mathrm{Hermite}}})
#' for a pair of variables or for every pair of columns in a matrix or data frame.
#' Supports both a **copula-free** empirical cross-moment formulation (\code{copula = "none"}, default)
#' and a **parametric Gaussian-copula** formulation via Mehler's bilinear expansion
#' (\code{copula = "gaussian"}).
#'
#' Unlike rank correlations (Spearman's \eqn{\rho}, Kendall's \eqn{\tau}) or Gaussian rank
#' correlation, \eqn{r_{\mathrm{Hermite}}} preserves the raw manifest Pearson correlation estimand
#' without altering the metric or discarding authentic tail variance.
#'
#' @param x A numeric vector, numeric \code{\link[base]{matrix}}, or \code{\link[base]{data.frame}}.
#' @param y A numeric vector, required when \code{x} is a vector; ignored when \code{x} is a matrix or data frame.
#' @param poly_degree Integer scalar; maximum polynomial degree used for marginal quantile
#'   fits (default \code{3L}). Must be an odd integer (\code{1}, \code{3}, or \code{5});
#'   see \code{\link{hermite_fit}}.
#' @param copula Character string specifying the bivariate dependence structure:
#'   \describe{
#'     \item{\code{"none"}}{(Default) Copula-free estimator. Evaluates the manifest covariance
#'       directly from empirical cross-moments of the fitted polynomial quantile reconstructions.
#'       Makes no bivariate normality or copula assumptions.}
#'     \item{\code{"gaussian"}}{Parametric Gaussian-copula estimator. Evaluates covariance
#'       analytically in closed form via Mehler's (1866) bilinear expansion under a latent
#'       bivariate normal copula. Also yields the latent copula correlation \eqn{\rho_z} and
#'       the shape-attenuation factor \eqn{A}.}
#'   }
#' @param monotonicity Monotonicity constraint passed to \code{\link{hermite_fit}}:
#'   \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Rank tie-handling method passed to \code{\link{hermite_fit}}:
#'   \code{"average"} (default) or \code{"random"}.
#' @param trim Numeric in \eqn{[0, 0.5)}. Proportion of symmetric trimming applied to the
#'   pointwise product of centered fitted reconstructions when \code{copula = "none"}.
#'   Default is \code{0} (exact closed-form empirical mean).
#' @param conf_level Numeric value in \eqn{(0, 1)}, e.g. \code{0.95}, or \code{NULL} (default)
#'   to skip confidence interval calculation (two-vector method only).
#' @param ci_method Confidence interval method for two-vector input:
#'   \describe{
#'     \item{\code{"fisher"}}{(Default) Analytical Fisher \eqn{z}-transformation with asymptotic
#'       standard error \eqn{\mathrm{SE} = 1/\sqrt{n - 3}}.}
#'     \item{\code{"bootstrap"}}{Non-parametric percentile bootstrap that refits the entire
#'       marginal quantile and covariance pipeline on each resample.}
#'   }
#' @param B Integer; number of bootstrap resamples when \code{ci_method = "bootstrap"}
#'   (default \code{1000L}).
#' @param diagnostics Logical; if \code{TRUE} (and \code{copula = "none"}), attaches the
#'   empirical Hermite cross-moment matrix \eqn{\hat{\mathbf{M}}_{kl}} and reference latent
#'   correlation \eqn{\rho_z} as attributes. Default is \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \subsection{The Mathematical Framework}{
#' Writing the fitted marginal polynomials in the orthogonal Probabilists' Hermite basis:
#' \deqn{f(Z_x) - \mu_x = \sum_{k=1}^{d_x} a_k He_k(Z_x), \qquad g(Z_y) - \mu_y = \sum_{l=1}^{d_y} b_l He_l(Z_y)}
#' the manifest covariance of the fitted quantile reconstructions decomposes exactly as:
#' \deqn{\operatorname{Cov}(X, Y) = \sum_{k=1}^{d_x} \sum_{l=1}^{d_y} a_k b_l \, \mathbb{E}\left[He_k(Z_x) He_l(Z_y)\right]}
#'
#' \itemize{
#'   \item \strong{Copula-Free (\code{copula = "none"}):} Evaluates the expectation directly from
#'         the sample cross-moments, which is algebraically identical to the sample covariance of the
#'         monotonically smoothed reconstructions \eqn{\hat{w} = (f(Z_x) - \mu_x)(g(Z_y) - \mu_y)}.
#'   \item \strong{Gaussian Copula (\code{copula = "gaussian"}):} Evaluates the expectation analytically
#'         via Mehler's (1866) expansion \eqn{\mathbb{E}[He_k(Z_x)He_l(Z_y)] = \delta_{kl} k! \rho_z^k},
#'         yielding \eqn{\operatorname{Cov}(X, Y) = \sum_{k=1}^{\min(d_x, d_y)} a_k b_k k! \rho_z^k}.
#' }
#' }
#'
#' \subsection{Shape Attenuation Factor (\eqn{A})}{
#' When \code{copula = "gaussian"}, the framework yields the latent copula parameter \eqn{\rho_z}
#' and the shape attenuation factor:
#' \deqn{A = \frac{r_{\mathrm{Hermite}}}{\rho_z}}
#' which quantifies how much the manifest linear correlation is constrained by marginal shape
#' asymmetry (Hoeffding-Fréchet bounds). For \code{copula = "none"}, \eqn{\rho_z} and \eqn{A}
#' are not defined.
#' }
#'
#' @return
#' For two vectors, an S3 object of class \code{"cor_hermite"} containing:
#' \describe{
#'   \item{\code{r_Hermite}}{The regularized Hermite Pearson correlation.}
#'   \item{\code{copula}}{The copula specification used (\code{"none"} or \code{"gaussian"}).}
#'   \item{\code{rho_z}}{The latent Gaussian copula correlation (only if \code{copula = "gaussian"}).}
#'   \item{\code{attenuation}}{The shape attenuation factor \eqn{A} (only if \code{copula = "gaussian"}).}
#'   \item{\code{cov_xy}}{The regularized manifest covariance.}
#'   \item{\code{var_x, var_y}}{Regularized manifest variances of \eqn{X} and \eqn{Y}.}
#'   \item{\code{mean_x, mean_y}}{Regularized manifest means of \eqn{X} and \eqn{Y}.}
#'   \item{\code{degrees}}{Named integer vector of realized polynomial degrees.}
#'   \item{\code{poly_degree_requested}}{Requested maximum polynomial degree.}
#'   \item{\code{monotonicity}, \code{ties_method}, \code{trim}}{Configuration settings applied.}
#'   \item{\code{n}}{Number of complete paired observations.}
#'   \item{\code{fit_x}, \code{fit_y}}{The underlying \code{\link{hermite_fit}} objects.}
#'   \item{\code{cross_moments}, \code{rho_z_reference}}{Diagnostic cross-moments (if requested).}
#'   \item{\code{ci}, \code{conf_level}, \code{ci_method}}{Confidence interval results, if requested.}
#' }
#'
#' For a matrix or data frame, an S3 object of class \code{"cor_hermite_matrix"} containing:
#' \describe{
#'   \item{\code{r_Hermite}}{The \eqn{p \times p} regularized Hermite correlation matrix.}
#'   \item{\code{copula}}{The copula specification used.}
#'   \item{\code{rho_z}}{The \eqn{p \times p} latent copula matrix (only if \code{copula = "gaussian"}).}
#'   \item{\code{attenuation}}{The \eqn{p \times p} shape-attenuation matrix (only if \code{copula = "gaussian"}).}
#'   \item{\code{cov}}{The \eqn{p \times p} regularized covariance matrix (with variances on the diagonal).}
#'   \item{\code{marginals}}{A data frame of regularized univariate summary statistics for each column.}
#'   \item{\code{n}}{Sample size.}
#' }
#'
#' @references
#' Carroll, J. B. (1961). The nature of the data, or how to choose a correlation coefficient. \emph{Psychometrika}, 26(4), 347–372. \doi{10.1007/BF02289768}
#'
#' Hermite, C. (1864). Sur un nouveau développement en série des fonctions. \emph{Comptes Rendus de l'Académie des Sciences, Paris}, 58, 93–100.
#'
#' Lenhard, W., & Lenhard, A. (2026). The Hermite-Mehler Correlation: A Distribution-Robust Estimator of the Pearson Correlation Coefficient. \emph{Behavior Research Methods}.
#'
#' Mehler, F. G. (1866). Ueber die Entwicklung einer Function von beliebig vielen Variabeln nach Laplaceschen Functionen höherer Ordnung. \emph{Journal für die reine und angewandte Mathematik}, 66, 161–176. \doi{10.1515/crll.1866.66.161}
#'
#' @seealso \code{\link{d_reg}}, \code{\link{hermite_fit}}, \code{\link{cor_hermite_boot_ci}}
#'
#' @examples
#' # 1. Copula-free Hermite correlation (default)
#' set.seed(123)
#' x <- rnorm(50, mean = 100, sd = 15)
#' y <- exp(0.5 * scale(x) + rnorm(50, sd = 0.5))
#' r_free <- cor_hermite(x, y, conf_level = 0.95)
#' print(r_free)
#'
#' # 2. Gaussian-copula Hermite-Mehler correlation
#' r_gauss <- cor_hermite(x, y, copula = "gaussian", conf_level = 0.95)
#' print(r_gauss)
#' summary(r_gauss)
#'
#' # 3. Multivariate Correlation Matrix
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
# Internal Helpers: Evaluated Hermite Basis & Pairwise Kernel
# -----------------------------------------------------------------------------



#' Pairwise Hermite Covariance & Correlation Kernel
#' @noRd
.cor_hermite_pair <- function(fit_x, fit_y, copula = "none", trim = 0, diagnostics = FALSE) {
  var_x  <- fit_x$variance
  var_y  <- fit_y$variance
  mean_x <- fit_x$mean
  mean_y <- fit_y$mean
  dx     <- fit_x$degree
  dy     <- fit_y$degree

  if (is.na(var_x) || is.na(var_y) || var_x <= 0 || var_y <= 0) {
    return(list(
      r_Hermite = NA_real_, rho_z = NA_real_, attenuation = NA_real_,
      cov_xy = NA_real_, var_x = var_x, var_y = var_y,
      mean_x = mean_x, mean_y = mean_y,
      degrees = c(x = dx, y = dy),
      cross_moments = NULL, rho_z_reference = NULL
    ))
  }

  denom <- sqrt(var_x * var_y)

  if (copula == "none") {
    # Copula-Free: Empirical products of centered fitted reconstructions
    Zmat_x <- outer(fit_x$z, 0:dx, `^`)
    Zmat_y <- outer(fit_y$z, 0:dy, `^`)
    fitted_x <- as.vector(Zmat_x %*% fit_x$beta)
    fitted_y <- as.vector(Zmat_y %*% fit_y$beta)

    w <- (fitted_x - mean_x) * (fitted_y - mean_y)
    cov_xy <- mean(w, trim = trim)

    if (!is.finite(cov_xy)) {
      r_Hermite <- NA_real_
    } else {
      cov_xy <- max(-denom, min(denom, cov_xy))
      r_Hermite <- max(-1.0, min(1.0, cov_xy / denom))
    }

    rho_z <- NA_real_
    attenuation <- NA_real_
    cross_moments <- NULL
    rho_z_ref <- NULL

    if (diagnostics) {
      n <- length(fit_x$z)
      Hx <- .hermite_eval_basis(fit_x$z, dx)[, -1L, drop = FALSE]
      Hy <- .hermite_eval_basis(fit_y$z, dy)[, -1L, drop = FALSE]
      cross_moments <- crossprod(Hx, Hy) / n
      dimnames(cross_moments) <- list(paste0("He", seq_len(dx)), paste0("He", seq_len(dy)))
      rho_z_ref <- stats::cor(fit_x$z, fit_y$z)
    }

  } else {
    # Gaussian Copula: Analytical Mehler bilinear expansion
    rho_z <- stats::cor(fit_x$z, fit_y$z)
    if (!is.finite(rho_z)) {
      return(list(
        r_Hermite = NA_real_, rho_z = NA_real_, attenuation = NA_real_,
        cov_xy = NA_real_, var_x = var_x, var_y = var_y,
        mean_x = mean_x, mean_y = mean_y,
        degrees = c(x = dx, y = dy),
        cross_moments = NULL, rho_z_reference = NULL
      ))
    }
    rho_z <- max(-1.0, min(1.0, rho_z))

    max_d <- max(dx, dy)
    H <- .hermite_basis_matrix(max_d)
    bx_pad <- rep(0.0, max_d + 1L); bx_pad[1:(dx + 1L)] <- fit_x$beta
    by_pad <- rep(0.0, max_d + 1L); by_pad[1:(dy + 1L)] <- fit_y$beta

    a <- as.vector(H %*% bx_pad)
    b <- as.vector(H %*% by_pad)
    fact_v <- factorial(1L:max_d)

    cov_xy <- sum(a[-1L] * b[-1L] * fact_v * (rho_z^(1L:max_d)))
    r_Hermite <- max(-1.0, min(1.0, cov_xy / denom))
    attenuation <- if (abs(rho_z) > 1e-8) r_Hermite / rho_z else 1.0
    cross_moments <- NULL
    rho_z_ref <- rho_z
  }

  list(
    r_Hermite            = r_Hermite,
    rho_z           = rho_z,
    attenuation     = attenuation,
    cov_xy          = cov_xy,
    var_x           = var_x,
    var_y           = var_y,
    mean_x          = mean_x,
    mean_y          = mean_y,
    degrees         = c(x = dx, y = dy),
    cross_moments   = cross_moments,
    rho_z_reference = rho_z_ref
  )
}

# -----------------------------------------------------------------------------
# Methods: default, matrix, data.frame
# -----------------------------------------------------------------------------

#' @rdname cor_hermite
#' @export
cor_hermite.default <- function(x, y, poly_degree = 3L,
                                copula = c("none", "gaussian"),
                                monotonicity = c("relaxed", "strict", "none"),
                                ties_method = c("average", "random"),
                                trim = 0,
                                conf_level = NULL,
                                ci_method = c("fisher", "bootstrap"),
                                B = 1000L,
                                diagnostics = FALSE, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)
  ci_method    <- match.arg(ci_method)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied for two-vector correlation.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (length(x) != length(y)) stop("'x' and 'y' must have identical length.")

  if (length(trim) != 1L || is.na(trim) || trim < 0 || trim >= 0.5) {
    stop("'trim' must be a single number in [0, 0.5).")
  }

  ok <- is.finite(x) & is.finite(y)
  x_c <- x[ok]; y_c <- y[ok]
  n <- length(x_c)

  if (n < 4L || diff(range(x_c)) == 0 || diff(range(y_c)) == 0) {
    warning("Insufficient or constant observations.")
    return(structure(
      list(r_Hermite = NA_real_, copula = copula, rho_z = NA_real_, attenuation = NA_real_),
      class = "cor_hermite"
    ))
  }

  fit_x <- hermite_fit(x_c, degree = poly_degree, monotonicity = monotonicity,
                       ties_method = ties_method, force_odd = TRUE)
  fit_y <- hermite_fit(y_c, degree = poly_degree, monotonicity = monotonicity,
                       ties_method = ties_method, force_odd = TRUE)

  pair <- .cor_hermite_pair(fit_x, fit_y, copula = copula, trim = trim, diagnostics = diagnostics)

  res <- list(
    r_Hermite                  = pair$r_Hermite,
    copula                = copula,
    rho_z                 = pair$rho_z,
    attenuation           = pair$attenuation,
    cov_xy                = pair$cov_xy,
    var_x                 = pair$var_x,
    var_y                 = pair$var_y,
    mean_x                = pair$mean_x,
    mean_y                = pair$mean_y,
    degrees               = pair$degrees,
    poly_degree_requested = poly_degree,
    monotonicity          = monotonicity,
    ties_method           = ties_method,
    trim                  = trim,
    n                     = n,
    fit_x                 = fit_x,
    fit_y                 = fit_y,
    cross_moments         = pair$cross_moments,
    rho_z_reference       = pair$rho_z_reference,
    x                     = x_c,
    y                     = y_c
  )
  class(res) <- "cor_hermite"

  if (!is.null(conf_level)) {
    res$ci <- stats::confint(res, level = conf_level, method = ci_method, B = B)
    res$conf_level <- conf_level
    res$ci_method  <- ci_method
  }
  res
}

#' @rdname cor_hermite
#' @export
cor_hermite.matrix <- function(x, poly_degree = 3L,
                               copula = c("none", "gaussian"),
                               monotonicity = c("relaxed", "strict", "none"),
                               ties_method = c("average", "random"),
                               trim = 0, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  if (!is.numeric(x)) stop("'x' must be a numeric matrix.")
  p <- ncol(x)
  if (is.null(p) || p < 2L) stop("'x' must contain at least two numeric columns.")

  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(p))

  R_Hermite     <- diag(1.0, nrow = p, ncol = p)
  cov_Hermite   <- matrix(NA_real_, nrow = p, ncol = p)
  dimnames(R_Hermite)   <- list(cnames, cnames)
  dimnames(cov_Hermite) <- list(cnames, cnames)

  R_copula <- if (copula == "gaussian") diag(1.0, nrow = p, ncol = p) else NULL
  Att_mat  <- if (copula == "gaussian") diag(1.0, nrow = p, ncol = p) else NULL
  if (copula == "gaussian") {
    dimnames(R_copula) <- list(cnames, cnames)
    dimnames(Att_mat)  <- list(cnames, cnames)
  }

  fully_complete <- all(is.finite(x))

  if (fully_complete) {
    fits <- vector("list", p)
    means_v <- numeric(p)
    vars_v  <- numeric(p)
    sds_v   <- numeric(p)
    skew_v  <- numeric(p)
    degs_v  <- integer(p)

    for (i in seq_len(p)) {
      fits[[i]] <- hermite_fit(x[, i], degree = poly_degree, monotonicity = monotonicity,
                               ties_method = ties_method, force_odd = TRUE)
      means_v[i] <- fits[[i]]$mean
      vars_v[i]  <- fits[[i]]$variance
      sds_v[i]   <- fits[[i]]$sd
      skew_v[i]  <- fits[[i]]$skewness
      degs_v[i]  <- fits[[i]]$degree
    }

    diag(cov_Hermite) <- vars_v

    for (i in 1:(p - 1L)) {
      for (j in (i + 1L):p) {
        pair <- .cor_hermite_pair(fits[[i]], fits[[j]], copula = copula, trim = trim)
        R_Hermite[i, j]   <- R_Hermite[j, i]   <- pair$r_Hermite
        cov_Hermite[i, j] <- cov_Hermite[j, i] <- pair$cov_xy
        if (copula == "gaussian") {
          R_copula[i, j] <- R_copula[j, i] <- pair$rho_z
          Att_mat[i, j]  <- Att_mat[j, i]  <- pair$attenuation
        }
      }
    }
  } else {
    # Pairwise complete fallback
    means_v <- numeric(p)
    vars_v  <- numeric(p)
    sds_v   <- numeric(p)
    skew_v  <- numeric(p)
    degs_v  <- integer(p)

    for (i in seq_len(p)) {
      xi <- x[is.finite(x[, i]), i]
      fit_i <- hermite_fit(xi, degree = poly_degree, monotonicity = monotonicity,
                           ties_method = ties_method, force_odd = TRUE)
      means_v[i] <- fit_i$mean
      vars_v[i]  <- fit_i$variance
      sds_v[i]   <- fit_i$sd
      skew_v[i]  <- fit_i$skewness
      degs_v[i]  <- fit_i$degree
    }

    diag(cov_Hermite) <- vars_v

    for (i in 1:(p - 1L)) {
      for (j in (i + 1L):p) {
        fit <- cor_hermite.default(x[, i], x[, j], poly_degree = poly_degree,
                                   copula = copula, monotonicity = monotonicity,
                                   ties_method = ties_method, trim = trim)
        R_Hermite[i, j]   <- R_Hermite[j, i]   <- fit$r_Hermite
        cov_Hermite[i, j] <- cov_Hermite[j, i] <- fit$cov_xy
        if (copula == "gaussian") {
          R_copula[i, j] <- R_copula[j, i] <- fit$rho_z
          Att_mat[i, j]  <- Att_mat[j, i]  <- fit$attenuation
        }
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
      r_Hermite        = R_Hermite,
      copula      = copula,
      rho_z       = R_copula,
      attenuation = Att_mat,
      cov         = cov_Hermite,
      marginals   = marginals_df,
      n           = nrow(x)
    ),
    class = "cor_hermite_matrix"
  )
}

#' @rdname cor_hermite
#' @export
cor_hermite.data.frame <- function(x, poly_degree = 3L,
                                   copula = c("none", "gaussian"),
                                   monotonicity = c("relaxed", "strict", "none"),
                                   ties_method = c("average", "random"),
                                   trim = 0, ...) {
  is_num <- vapply(x, is.numeric, logical(1L))
  if (!all(is_num)) {
    stop("All columns of 'x' must be numeric. Non-numeric columns found: ",
         paste(names(x)[!is_num], collapse = ", "))
  }
  cor_hermite.matrix(as.matrix(x), poly_degree = poly_degree, copula = copula,
                     monotonicity = monotonicity, ties_method = ties_method,
                     trim = trim, ...)
}


#' Distribution-Robust Partial and Semipartial Hermite Correlation
#'
#' Computes the distribution-robust partial or semipartial correlation between
#' two variables controlling for one or more covariates, or computes the full
#' partial correlation matrix for a multivariate dataset via precision matrix inversion.
#'
#' @param x A numeric vector (variable X), a numeric matrix, or a data frame.
#' @param y A numeric vector (variable Y; required if \code{x} is a vector).
#' @param z A numeric vector, matrix, or data frame of controlling covariates \eqn{\mathbf{Z}}.
#' @param semi Logical; if \code{TRUE}, computes the semipartial (part) correlation
#'   where \code{z} is partialled out of \code{y} only. Default is \code{FALSE} (partial correlation).
#' @param copula Character string; \code{"none"} (default, copula-free cross-moments)
#'   or \code{"gaussian"} (Mehler bilinear identity).
#' @param poly_degree Integer; maximum polynomial degree (default \code{3L}).
#' @param monotonicity Monotonicity constraint: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param conf_level Numeric value in \eqn{(0, 1)} for confidence intervals, or \code{NULL} (default).
#' @param ci_method Method for CI: \code{"fisher"} (analytical Fisher z) or \code{"bootstrap"}.
#' @param B Integer; number of bootstrap replications (default \code{1000L}).
#' @param ... Additional arguments passed to \code{\link{cor_hermite}}.
#'
#' @return
#' If \code{x} is a matrix or data frame (and \code{y} is \code{NULL}), returns an
#' object of class \code{"pcor_hermite_matrix"} containing:
#' \describe{
#'   \item{\code{pcor}}{The \eqn{p \times p} partial correlation matrix.}
#'   \item{\code{r_Hermite}}{The \eqn{p \times p} pairwise Hermite correlation matrix.}
#'   \item{\code{cov}}{The \eqn{p \times p} regularized covariance matrix.}
#'   \item{\code{copula}}{The copula mode applied.}
#'   \item{\code{n}}{Sample size.}
#' }
#'
#' If \code{x}, \code{y}, and \code{z} are supplied, returns an object of class \code{"pcor_hermite"} containing:
#' \describe{
#'   \item{\code{estimate}}{The partial or semipartial correlation point estimate.}
#'   \item{\code{rho_z_pcor}}{The latent copula partial correlation (if \code{copula = "gaussian"}).}
#'   \item{\code{attenuation}}{The partial shape-attenuation factor (if \code{copula = "gaussian"}).}
#'   \item{\code{semi}}{Logical; whether semipartial correlation was computed.}
#'   \item{\code{copula}}{The copula mode applied.}
#'   \item{\code{k_controls}}{Number of controlled covariates.}
#'   \item{\code{n}}{Number of complete cases.}
#'   \item{\code{ci}}{Confidence limits (if \code{conf_level} was specified).}
#' }
#'
#' @examples
#' # 1. Full Partial Correlation Matrix (iris dataset)
#' data(iris)
#' pcor_mat <- pcor_hermite(iris[, 1:4])
#' print(pcor_mat)
#'
#' # 2. Specific Pair controlling for a Confounder
#' set.seed(42)
#' z <- rlnorm(60, meanlog = 1, sdlog = 0.5)
#' x <- 0.6 * z + rnorm(60, sd = 2)
#' y <- 0.7 * z + rnorm(60, sd = 2)
#'
#' fit_pcor <- pcor_hermite(x, y, z = z, conf_level = 0.95)
#' print(fit_pcor)
#'
#' @export
pcor_hermite <- function(x, y = NULL, z = NULL,
                         semi = FALSE,
                         copula = c("none", "gaussian"),
                         poly_degree = 3L,
                         monotonicity = c("relaxed", "strict", "none"),
                         conf_level = NULL,
                         ci_method = c("fisher", "bootstrap"),
                         B = 1000L, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  ci_method    <- match.arg(ci_method)

  # =========================================================================
  # Case 1: Full Matrix / Data Frame Precision Inversion
  # =========================================================================
  if (is.matrix(x) || is.data.frame(x)) {
    if (is.data.frame(x)) {
      is_num <- vapply(x, is.numeric, logical(1L))
      if (!all(is_num)) {
        stop("All columns must be numeric. Non-numeric columns found: ",
             paste(names(x)[!is_num], collapse = ", "))
      }
      mat <- as.matrix(x)
    } else {
      mat <- x
    }

    p <- ncol(mat)
    if (is.null(p) || p < 3L) {
      stop("Matrix or data frame must have at least 3 numeric columns for partial correlations.")
    }

    cnames <- colnames(mat)
    if (is.null(cnames)) cnames <- paste0("V", seq_len(p))

    # Compute full Hermite correlation matrix
    fit_mat <- cor_hermite(mat, poly_degree = poly_degree, copula = copula,
                           monotonicity = monotonicity, ...)
    R <- fit_mat$r_Hermite

    # Invert correlation matrix (with automatic ridge stabilization if ill-conditioned)
    inv_R <- tryCatch(
      solve(R),
      error = function(e) {
        solve(R + diag(1e-4, nrow = p, ncol = p))
      }
    )

    # Standardize precision matrix to partial correlation matrix
    diag_inv <- diag(inv_R)
    P_mat <- -inv_R / sqrt(outer(diag_inv, diag_inv))
    diag(P_mat) <- 1.0
    dimnames(P_mat) <- list(cnames, cnames)

    res <- list(
      pcor   = P_mat,
      r_Hermite   = R,
      cov    = fit_mat$cov,
      copula = copula,
      n      = fit_mat$n
    )
    class(res) <- "pcor_hermite_matrix"
    return(res)
  }

  # =========================================================================
  # Case 2: Specific Pair (X, Y) Controlling for Z
  # =========================================================================
  if (is.null(y) || is.null(z)) {
    stop("Arguments 'y' and 'z' must be supplied when 'x' is a numeric vector.")
  }

  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric vectors.")
  }

  Z_mat <- as.matrix(z)
  if (!is.numeric(Z_mat)) stop("'z' must be numeric.")

  if (length(x) != length(y) || length(x) != nrow(Z_mat)) {
    stop("'x', 'y', and 'z' must have identical row dimensions.")
  }

  df_all <- data.frame(X = x, Y = y, Z_mat)
  ok <- stats::complete.cases(df_all)
  df_all <- df_all[ok, , drop = FALSE]
  n <- nrow(df_all)
  k <- ncol(Z_mat)

  if (n < (k + 4L)) {
    warning("Insufficient complete observations for partial correlation.")
    return(structure(list(estimate = NA_real_, semi = semi, copula = copula,
                          k_controls = k, n = n), class = "pcor_hermite"))
  }

  # Compute regularized covariance matrix for the combined system
  fit_mat <- cor_hermite(df_all, poly_degree = poly_degree, copula = copula,
                         monotonicity = monotonicity, ...)
  Sigma <- fit_mat$cov

  # Partition Covariance Matrix: Block 1 = (X, Y), Block 2 = Z
  S11 <- Sigma[1:2, 1:2, drop = FALSE]
  S12 <- Sigma[1:2, 3:(2 + k), drop = FALSE]
  S21 <- Sigma[3:(2 + k), 1:2, drop = FALSE]
  S22 <- Sigma[3:(2 + k), 3:(2 + k), drop = FALSE]

  inv_S22 <- tryCatch(
    solve(S22),
    error = function(e) solve(S22 + diag(1e-4, nrow = k, ncol = k))
  )
  S_cond <- S11 - S12 %*% inv_S22 %*% S21

  # Compute partial or semipartial correlation
  if (semi) {
    denom <- sqrt(S11[1, 1] * S_cond[2, 2])
  } else {
    denom <- sqrt(S_cond[1, 1] * S_cond[2, 2])
  }

  pcor_val <- if (is.na(denom) || denom <= 0) NA_real_ else max(-1.0, min(1.0, S_cond[1, 2] / denom))

  # Latent Copula Partial Correlation (if copula = "gaussian")
  rho_z_pcor <- NA_real_
  att_val    <- NA_real_

  if (copula == "gaussian" && !is.null(fit_mat$rho_z)) {
    R_z <- fit_mat$rho_z
    R11 <- R_z[1:2, 1:2, drop = FALSE]
    R12 <- R_z[1:2, 3:(2 + k), drop = FALSE]
    R21 <- R_z[3:(2 + k), 1:2, drop = FALSE]
    R22 <- R_z[3:(2 + k), 3:(2 + k), drop = FALSE]

    inv_R22 <- tryCatch(
      solve(R22),
      error = function(e) solve(R22 + diag(1e-4, nrow = k, ncol = k))
    )
    R_cond <- R11 - R12 %*% inv_R22 %*% R21
    denom_z <- sqrt(R_cond[1, 1] * R_cond[2, 2])
    rho_z_pcor <- if (is.na(denom_z) || denom_z <= 0) NA_real_ else max(-1.0, min(1.0, R_cond[1, 2] / denom_z))
    att_val <- if (is.finite(rho_z_pcor) && abs(rho_z_pcor) > 1e-6) pcor_val / rho_z_pcor else 1.0
  }

  res <- list(
    estimate     = pcor_val,
    rho_z_pcor   = rho_z_pcor,
    attenuation  = att_val,
    semi         = semi,
    copula       = copula,
    k_controls   = k,
    n            = n,
    poly_degree  = poly_degree,
    monotonicity = monotonicity,
    data         = df_all
  )
  class(res) <- "pcor_hermite"

  # Confidence Intervals
  if (!is.null(conf_level) && is.finite(pcor_val)) {
    if (ci_method == "fisher") {
      z_val <- atanh(pcor_val)
      se    <- 1 / sqrt(max(1, n - k - 3L))
      crit  <- stats::qnorm((1 + conf_level) / 2)
      res$ci <- tanh(z_val + c(-1, 1) * crit * se)
    } else {
      boot_vals <- numeric(B)
      for (b in seq_len(B)) {
        idx <- sample.int(n, n, replace = TRUE)
        boot_vals[b] <- pcor_hermite(
          x = df_all$X[idx], y = df_all$Y[idx], z = df_all[idx, 3:(2 + k), drop = FALSE],
          semi = semi, copula = copula, poly_degree = poly_degree, monotonicity = monotonicity
        )$estimate
      }
      alpha <- (1 - conf_level) / 2
      res$ci <- stats::quantile(boot_vals, probs = c(alpha, 1 - alpha), na.rm = TRUE)
    }
    res$conf_level <- conf_level
    res$ci_method  <- ci_method
  }

  res
}

# -----------------------------------------------------------------------------
# S3 Print Methods for Partial Correlation
# -----------------------------------------------------------------------------

#' @export
print.pcor_hermite <- function(x, digits = 3L, ...) {
  type_lab <- if (x$semi) "Semipartial (Part)" else "Partial"
  cop_lab  <- if (x$copula == "gaussian") "Gaussian Copula" else "Copula-Free"

  cat(sprintf("\n  Distribution-Robust %s Correlation (r_Hermite)\n", type_lab))
  cat(strrep("-", 54), "\n", sep = "")
  cat(sprintf("  Estimate (r_Hermite.z)          :  %.*f\n", digits, x$estimate))
  cat(sprintf("  Copula Mode                :  %s\n", cop_lab))
  cat(sprintf("  Controlled Covariates (k)  :  %d\n", x$k_controls))
  cat(sprintf("  Sample Size (n)            :  %d\n", x$n))

  if (x$copula == "gaussian" && is.finite(x$rho_z_pcor)) {
    cat(sprintf("  Latent Copula Partial r    :  %.*f\n", digits, x$rho_z_pcor))
    cat(sprintf("  Shape Attenuation Factor   :  %.*f\n", digits, x$attenuation))
  }

  if (!is.null(x$ci)) {
    cat(sprintf("  %s CI (%s): [%.*f, %.*f]\n",
                paste0(round(x$conf_level * 100), "%"), x$ci_method,
                digits, x$ci[1L], digits, x$ci[2L]))
  }
  cat("\n")
  invisible(x)
}

#' @export
print.pcor_hermite_matrix <- function(x, digits = 3L, ...) {
  cop_lab <- if (x$copula == "gaussian") "Gaussian Copula" else "Copula-Free"

  cat(sprintf("\n  Hermite Partial Correlation Matrix (r_Hermite.z; %s):\n", cop_lab))
  cat(strrep("-", 56), "\n", sep = "")
  print(round(x$pcor, digits = digits))
  cat("\n")
  invisible(x)
}



#' Convenience Alias for Gaussian-Copula Hermite-Mehler Correlation
#' @rdname cor_hermite
#' @export
cor_hermitegauss <- function(x, y, poly_degree = 3L,
                             monotonicity = c("relaxed", "strict", "none"),
                             ties_method = c("average", "random"), ...) {
  cor_hermite(x = x, y = y, poly_degree = poly_degree, copula = "gaussian",
              monotonicity = monotonicity, ties_method = ties_method, ...)
}

# -----------------------------------------------------------------------------
# Inference & S3 Methods
# -----------------------------------------------------------------------------

#' Percentile Bootstrap Confidence Interval for Hermite Correlation
#'
#' @param object A fitted \code{"cor_hermite"} object.
#' @param x Numeric vector of original complete-case x data (optional if stored in \code{object}).
#' @param y Numeric vector of original complete-case y data (optional if stored in \code{object}).
#' @param B Integer; number of bootstrap replicates (default \code{1000L}).
#' @param conf Numeric scalar; confidence level (default \code{0.95}).
#' @param min_success Minimum proportion of successful replicates (default \code{0.90}).
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"cor_hermite_boot_ci"}.
#' @export
cor_hermite_boot_ci <- function(object, x = NULL, y = NULL,
                                B = 1000L, conf = 0.95, min_success = 0.90, ...) {
  if (!inherits(object, "cor_hermite")) {
    stop("'object' must be a fitted 'cor_hermite' object.")
  }

  poly_degree  <- if (!is.null(object$poly_degree_requested)) object$poly_degree_requested else 3L
  copula       <- if (!is.null(object$copula)) object$copula else "none"
  monotonicity <- if (!is.null(object$monotonicity)) object$monotonicity else "relaxed"
  ties_method  <- if (!is.null(object$ties_method)) object$ties_method else "average"
  trim         <- if (!is.null(object$trim)) object$trim else 0

  dots <- list(...)
  if (is.null(x)) x <- if (!is.null(object$x)) object$x else dots$x
  if (is.null(y)) y <- if (!is.null(object$y)) object$y else dots$y

  if (is.null(x) || is.null(y)) {
    stop("The original data are required for bootstrap resampling. Supply 'x' and 'y'.")
  }

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)

  if (n < poly_degree + 2L || diff(range(x)) == 0 || diff(range(y)) == 0) {
    return(structure(
      list(estimate = object$r_Hermite, lower = NA_real_, upper = NA_real_,
           conf = conf, B = B, n_success = 0L, n_fail = B, success_rate = 0,
           poly_degree = poly_degree, copula = copula, ties_method = ties_method,
           trim = trim, boot_estimates = numeric(0)),
      class = "cor_hermite_boot_ci"
    ))
  }

  boot_estimates <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_estimates[b] <- tryCatch(
      cor_hermite(x[idx], y[idx], poly_degree = poly_degree, copula = copula,
                  monotonicity = monotonicity, ties_method = ties_method,
                  trim = trim)$r_Hermite,
      error = function(e) NA_real_
    )
  }

  success <- is.finite(boot_estimates)
  n_success <- sum(success)
  success_rate <- n_success / B

  if (success_rate < min_success) {
    warning(sprintf("Only %d of %d bootstrap replicates succeeded (%.1f%%).",
                    n_success, B, 100 * success_rate), call. = FALSE)
    return(structure(
      list(estimate = object$r_Hermite, lower = NA_real_, upper = NA_real_,
           conf = conf, B = B, n_success = n_success, n_fail = B - n_success,
           success_rate = success_rate, poly_degree = poly_degree, copula = copula,
           ties_method = ties_method, trim = trim, boot_estimates = boot_estimates[success]),
      class = "cor_hermite_boot_ci"
    ))
  }

  alpha <- (1 - conf) / 2
  limits <- stats::quantile(boot_estimates[success], probs = c(alpha, 1 - alpha), names = FALSE, type = 7)

  structure(
    list(estimate = object$r_Hermite, lower = unname(limits[1L]), upper = unname(limits[2L]),
         conf = conf, B = B, n_success = n_success, n_fail = B - n_success,
         success_rate = success_rate, poly_degree = poly_degree, copula = copula,
         ties_method = ties_method, trim = trim, boot_estimates = boot_estimates[success]),
    class = "cor_hermite_boot_ci"
  )
}

#' Confidence Intervals for Hermite Correlation Objects
#'
#' Computes analytical Fisher \eqn{z} or non-parametric percentile bootstrap
#' confidence intervals for a fitted \code{"cor_hermite"} object.
#'
#' @param object An object of class \code{"cor_hermite"} created by \code{\link{cor_hermite}}.
#' @param parm Ignored; included for S3 method consistency with \code{\link[stats]{confint}}.
#' @param level Numeric scalar in \eqn{(0, 1)}; the requested confidence level (default is \code{0.95}).
#' @param method Character string specifying the confidence interval calculation method:
#'   \describe{
#'     \item{\code{"fisher"}}{(Default) Fast analytical Fisher \eqn{z}-transformation with
#'       asymptotic standard error \eqn{\mathrm{SE} = 1/\sqrt{n - 3}}.}
#'     \item{\code{"bootstrap"}}{Non-parametric percentile bootstrap refitting the full
#'       Hermite pipeline (NQT, polynomial quantile model, covariance estimation) across resamples.}
#'   }
#' @param B Integer scalar; number of bootstrap replications when \code{method = "bootstrap"}
#'   (default is \code{1000L}).
#' @param ... Additional arguments passed to internal methods.
#'
#' @details
#' When \code{method = "fisher"}, the interval is obtained by inverting the Fisher
#' \eqn{z}-transformation on \eqn{r_{\mathrm{Hermite}}}. When \code{method = "bootstrap"},
#' the function calls \code{\link{cor_hermite_boot_ci}}, refitting the model from scratch
#' in each bootstrap replicate while inheriting the original polynomial degree and
#' copula settings.
#'
#' @return A numeric matrix of dimension \code{1 x 2} containing the lower and upper
#'   confidence limits, with column names indicating the corresponding percentiles.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(40)
#' y <- 0.5 * x + rnorm(40)
#' fit <- cor_hermite(x, y)
#'
#' # Analytical Fisher-z 95% CI
#' confint(fit, level = 0.95, method = "fisher")
#'
#' # Percentile Bootstrap 95% CI
#' confint(fit, level = 0.95, method = "bootstrap", B = 500)
#'
#' @seealso \code{\link{cor_hermite}}, \code{\link{cor_hermite_boot_ci}}
#' @export
confint.cor_hermite <- function(object, parm, level = 0.95,
                                method = c("fisher", "bootstrap"), B = 1000L, ...) {
  method <- match.arg(method)
  r <- object$r_Hermite
  n <- object$n

  if (method == "fisher") {
    z_r <- atanh(r)
    se <- 1 / sqrt(n - 3L)
    crit <- stats::qnorm((1 + level) / 2)
    ci <- tanh(z_r + c(-1, 1) * crit * se)
  } else {
    boot_res <- cor_hermite_boot_ci(object, B = B, conf = level)
    ci <- c(boot_res$lower, boot_res$upper)
  }
  matrix(ci, nrow = 1L, dimnames = list("r_Hermite", c(paste0(100 * (1 - level)/2, " %"),
                                                  paste0(100 * (1 + level)/2, " %"))))
}

#' @export
print.cor_hermite <- function(x, digits = 3L, ...) {
  copula_label <- if (x$copula == "gaussian") "Gaussian Copula (Mehler's identity)" else "Copula-Free (empirical cross-moments)"

  cat("\n  Distribution-Robust Hermite Correlation\n")
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("  Hermite Correlation (r_Hermite) :  %.*f\n", digits, x$r_Hermite))
  cat(sprintf("  Copula Model               :  %s\n", copula_label))
  cat(sprintf("  Polynomial Degrees Fitted  :  X = %d, Y = %d\n", x$degrees["x"], x$degrees["y"]))
  cat(sprintf("  Monotonicity Check         :  %s\n", x$monotonicity))

  if (x$copula == "none" && !is.null(x$trim) && x$trim > 0) {
    cat(sprintf("  Trim Proportion            :  %.*f\n", digits, x$trim))
  }
  if (x$copula == "gaussian") {
    cat(sprintf("  Latent Copula (rho_z)      :  %.*f\n", digits, x$rho_z))
    cat(sprintf("  Shape Attenuation (A)      :  %.*f\n", digits, x$attenuation))
  }
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
  cat(sprintf("    X: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$mean_x, digits, sqrt(object$var_x), digits, object$var_x))
  cat(sprintf("    Y: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$mean_y, digits, sqrt(object$var_y), digits, object$var_y))
  cat(sprintf("    Covariance(X, Y) = %.*f\n\n", digits, object$cov_xy))

  if (!is.null(object$cross_moments)) {
    cat("  Empirical Hermite Cross-Moments (M_hat):\n")
    print(round(object$cross_moments, digits = digits))
    if (!is.null(object$rho_z_reference)) {
      cat(sprintf("  Reference Latent Correlation (rho_z): %.*f\n", digits, object$rho_z_reference))
    }
    cat("\n")
  }
  invisible(object)
}

#' @export
print.cor_hermite_matrix <- function(x, digits = 3L, ...) {
  copula_label <- if (x$copula == "gaussian") "Gaussian Copula" else "Copula-Free"

  cat(sprintf("\n  Hermite Correlation Matrix (r_Hermite; %s):\n", copula_label))
  cat(strrep("-", 52), "\n", sep = "")
  print(round(x$r_Hermite, digits = digits))

  if (x$copula == "gaussian") {
    cat("\n  Latent Copula Correlation Matrix (rho_z):\n")
    print(round(x$rho_z, digits = digits))

    cat("\n  Shape Attenuation Factor Matrix (A = r_Hermite / rho_z):\n")
    print(round(x$attenuation, digits = digits))
  }

  cat("\n  Hermite Covariance Matrix (cov_Hermite; Diagonal: Variances):\n")
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

#' @export
print.cor_hermite_boot_ci <- function(x, digits = 3L, ...) {
  cat("\n  Bootstrap Confidence Interval for Hermite Correlation\n")
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("  Estimate (r_Hermite)           :  %.*f\n", digits, x$estimate))
  cat(sprintf("  Copula Model              :  %s\n", x$copula))
  cat(sprintf("  %d%% Confidence Interval  : [%.*f, %.*f]\n",
              as.integer(x$conf * 100), digits, x$lower, digits, x$upper))
  cat(sprintf("  Bootstrap Replicates (B)  :  %d (%d successful, %.1f%%)\n",
              x$B, x$n_success, 100 * x$success_rate))
  cat("\n")
  invisible(x)
}

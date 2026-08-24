#' @title Hermite-Mehler Distribution-Robust Pearson Correlation
#'
#' @description
#' Computes the distribution-robust Hermite-Mehler Pearson correlation (\eqn{r_{\mathrm{HM}}}),
#' the latent Gaussian copula correlation (\eqn{\rho_z}), and the shape attenuation
#' factor (\eqn{A}) for bivariate pairs or multivariate datasets.
#'
#' The method preserves the raw manifest Pearson correlation estimand on the original
#' metric without changing the underlying parameter (unlike Spearman's \eqn{\rho},
#' Kendall's \eqn{\tau}, or Gaussian rank correlations) and without discarding authentic
#' tail variation (unlike trimming or Winsorization).
#'
#' @param x A numeric vector, numeric \code{\link[base]{matrix}}, or \code{\link[base]{data.frame}}.
#' @param y A numeric vector (required if \code{x} is a vector; ignored if \code{x} is a matrix or data frame).
#' @param poly_degree Integer scalar; maximum polynomial degree for quantile function
#'   smoothing (default is \code{3L}). Must be an odd integer (\code{1}, \code{3}, or \code{5}).
#' @param monotonicity Character string specifying the monotonicity constraint:
#'   \describe{
#'     \item{\code{"relaxed"}}{(Default) Empirical rank concordance (\eqn{\rho_{\text{Spearman}} \ge 0.95}).
#'       Preserves degrees well in discrete, psychometric, or count data.}
#'     \item{\code{"strict"}}{Analytical root check of the first derivative (\eqn{\min f'(z) \ge 0}).}
#'     \item{\code{"none"}}{Unconstrained OLS polynomial fit.}
#'   }
#' @param ties_method Character string specifying the rank tie-handling method:
#'   \code{"average"} (default, midranks) or \code{"random"}.
#' @param conf_level Numeric value in \eqn{(0, 1)} specifying the confidence level
#'   (e.g., \code{0.95} for a 95\% CI), or \code{NULL} (default) for point estimation only.
#' @param ci_method Character string specifying the confidence interval calculation method:
#'   \describe{
#'     \item{\code{"fisher"}}{(Default) Fast analytical Fisher \eqn{z}-transformation with
#'       asymptotic standard error \eqn{\mathrm{SE} = 1/\sqrt{n - 3}}.}
#'     \item{\code{"bootstrap"}}{Non-parametric percentile bootstrap refitting the full
#'       Hermite-Mehler pipeline across resamples.}
#'   }
#' @param B Integer scalar; number of bootstrap resamples if \code{ci_method = "bootstrap"}.
#'   Default is \code{1000L}.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' \subsection{The Estimand Dilemma in Robust Correlation}{
#' Classical Pearson's \eqn{r} suffers from severe sampling variance inflation when data
#' depart from normality (e.g., skewness, heavy tails, outliers). Traditional robust
#' alternatives solve this stability problem by changing the estimand:
#' \itemize{
#'   \item \strong{Spearman's \eqn{\rho} and Kendall's \eqn{\tau}} evaluate rank concordance
#'         rather than linear covariation, systematically underestimating Pearson's \eqn{r}
#'         under bivariate normality (\eqn{\mathbb{E}[\rho_s] = \frac{6}{\pi}\arcsin(\rho/2)}).
#'   \item \strong{Trimming / Winsorization} (e.g., percentage-bend, biweight midcorrelation)
#'         amputates structural tail variance, altering the population covariance structure.
#'   \item \strong{Gaussian Rank Correlation} evaluates the correlation on normalized latent
#'         scores, yielding the copula correlation \eqn{\rho_z} rather than the raw manifest \eqn{r}.
#' }
#' }
#'
#' \subsection{The Hermite-Mehler Pipeline}{
#' The \code{cor_hermite} estimator resolves this dilemma in three steps:
#' \enumerate{
#'   \item \strong{Rank-Based Inverse-Normal Transformation:} Both variables are mapped
#'         to standard normal scores \eqn{Z_x, Z_y \sim \mathcal{N}(0, 1)}. The sample
#'         correlation of these scores yields the latent Gaussian copula parameter \eqn{\rho_z}.
#'   \item \strong{Monotone Quantile Polynomial Fitting:} The inverse quantile functions
#'         \eqn{X = f(Z_x)} and \eqn{Y = g(Z_y)} are approximated via low-degree monotone
#'         polynomials via OLS regression.
#'   \item \strong{Closed-Form Mehler Expansion:} Polynomial coefficients are re-expressed
#'         in the basis of orthogonal Probabilists' Hermite polynomials \eqn{He_m(z)}. By
#'         Mehler's (1866) bilinear expansion:
#'         \deqn{\mathbb{E}\left[ He_m(Z_x) He_n(Z_y) \right] = \delta_{mn} \, m! \, \rho_z^m}
#'         All cross-terms of differing orders vanish. The covariance and variances are
#'         obtained in exact algebraic form:
#'         \deqn{\operatorname{Cov}(X, Y) = \sum_{m=1}^{\min(k_x, k_y)} a_m b_m \, m! \, \rho_z^m}
#'         \deqn{\operatorname{Var}(X) = \sum_{m=1}^{k_x} a_m^2 \, m!, \quad \operatorname{Var}(Y) = \sum_{m=1}^{k_y} b_m^2 \, m!}
#'         \deqn{r_{\mathrm{HM}} = \frac{\operatorname{Cov}(X, Y)}{\sqrt{\operatorname{Var}(X)\operatorname{Var}(Y)}}}
#' }
#' }
#'
#' \subsection{Shape Attenuation Factor (\eqn{A})}{
#' When marginal distribution shapes differ (e.g., a normal variable paired with a lognormal
#' latency), the mathematically achievable Pearson correlation is strictly bounded below \eqn{1.0}
#' (Hoeffding-Fréchet bounds; Carroll, 1961). The ratio:
#' \deqn{A = \frac{r_{\mathrm{HM}}}{\rho_z}}
#' provides an interpretable diagnostic index. If \eqn{A \approx 1.0}, marginal shapes do not
#' restrict the association. If \eqn{A \ll 1.0}, reporting a raw Pearson correlation severely
#' understates the substantive construct association purely due to scale mismatch.
#' }
#'
#' @return
#' For two vectors, an S3 object of class \code{"cor_hermite"} containing:
#' \describe{
#'   \item{\code{r_hm}}{Numeric; the regularized Hermite-Mehler Pearson correlation.}
#'   \item{\code{rho_z}}{Numeric; the latent Gaussian copula correlation (Gaussian rank correlation).}
#'   \item{\code{attenuation}}{Numeric; the shape attenuation factor \eqn{A = r_{\mathrm{HM}} / \rho_z}.}
#'   \item{\code{cov_xy}}{Numeric; the regularized manifest covariance.}
#'   \item{\code{var_x, var_y}}{Numeric; regularized manifest variances for \eqn{X} and \eqn{Y}.}
#'   \item{\code{mean_x, mean_y}}{Numeric; regularized manifest means for \eqn{X} and \eqn{Y}.}
#'   \item{\code{degrees}}{Named integer vector of realized polynomial degrees for \eqn{X} and \eqn{Y}.}
#'   \item{\code{poly_degree_requested}}{Integer; requested maximum polynomial degree.}
#'   \item{\code{monotonicity}}{Character; monotonicity method applied.}
#'   \item{\code{ties_method}}{Character; tie-handling method applied.}
#'   \item{\code{n}}{Integer; number of complete paired observations.}
#'   \item{\code{fit_x, fit_y}}{The underlying \code{\link{hermite_fit}} objects.}
#'   \item{\code{ci}}{Matrix containing confidence limits (if \code{conf_level} was specified).}
#'   \item{\code{conf_level}}{Requested confidence level.}
#'   \item{\code{ci_method}}{Method used for confidence interval estimation.}
#' }
#'
#' For a matrix or data frame, an S3 object of class \code{"cor_hermite_matrix"} containing:
#' \describe{
#'   \item{\code{r_hm}}{Matrix of pairwise Hermite-Mehler correlations.}
#'   \item{\code{rho_z}}{Matrix of pairwise latent copula correlations.}
#'   \item{\code{attenuation}}{Matrix of pairwise shape attenuation factors.}
#'   \item{\code{n}}{Integer; number of rows in the input dataset.}
#' }
#'
#' @references
#' Carroll, J. B. (1961). The nature of the data, or how to choose a correlation coefficient. \emph{Psychometrika}, 26(4), 347–372. \doi{10.1007/BF02289768}
#'
#' Fréchet, M. (1951). Sur les tableaux de corrélation dont les marges sont données. \emph{Annales de l'Université de Lyon}, 14, 53–77.
#'
#' Hermite, C. (1864). Sur un nouveau développement en série des fonctions. \emph{Comptes Rendus de l'Académie des Sciences, Paris}, 58, 93–100.
#'
#' Hoeffding, W. (1940). Masstabinvariante Korrelationstheorie. \emph{Schriften des Mathematischen Instituts und des Instituts für Angewandte Mathematik der Universität Berlin}, 5, 181–233.
#'
#' Lenhard, W., & Lenhard, A. (2026). The Hermite-Mehler Correlation: A Distribution-Robust Estimator of the Pearson Correlation Coefficient. \emph{Behavior Research Methods}. \doi{10.3758/s13428-xxx-xxxxx-x}
#'
#' Mehler, F. G. (1866). Ueber die Entwicklung einer Function von beliebig vielen Variabeln nach Laplaceschen Functionen höherer Ordnung. \emph{Journal für die reine und angewandte Mathematik}, 66, 161–176. \doi{10.1515/crll.1866.66.161}
#'
#' @seealso \code{\link{d_reg}}, \code{\link{hermite_fit}}, \code{\link{check_monotonicity}}
#'
#' @examples
#' # ---------------------------------------------------------
#' # 1. Bivariate Correlation with Marginal Shape Mismatch
#' # ---------------------------------------------------------
#' set.seed(123)
#' n <- 50
#' z <- rnorm(n)
#' x <- rnorm(n, mean = 100, sd = 15)           # Normally distributed
#' y <- exp(0.5 * z + rnorm(n, sd = 0.5))        # Skewed lognormal latency
#'
#' # Compute r_HM with 95% Fisher-z CI
#' r_fit <- cor_hermite(x, y, conf_level = 0.95, ci_method = "fisher")
#' print(r_fit)
#' summary(r_fit)
#'
#' # ---------------------------------------------------------
#' # 2. Multivariate Correlation Matrix
#' # ---------------------------------------------------------
#' data(iris)
#' r_mat <- cor_hermite(iris[, 1:4])
#' print(r_mat)
#'
#' # ---------------------------------------------------------
#' # 3. Percentile Bootstrap Confidence Interval
#' # ---------------------------------------------------------
#' r_boot <- cor_hermite(x, y, conf_level = 0.95, ci_method = "bootstrap", B = 500)
#' confint(r_boot)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
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

# -----------------------------------------------------------------------------
# S3 Methods for cor_hermite: confint, print, summary
# -----------------------------------------------------------------------------

#' Confidence Intervals for Hermite Correlation
#'
#' Computes analytical Fisher \eqn{z} or non-parametric percentile bootstrap
#' confidence intervals for \code{cor_hermite} objects.
#'
#' @param object An object of class \code{"cor_hermite"}.
#' @param parm Ignored.
#' @param level Numeric scalar; confidence level (default is \code{0.95}).
#' @param method Character; \code{"fisher"} (default) or \code{"bootstrap"}.
#' @param B Integer; number of bootstrap replications (default is \code{1000L}).
#' @param ... Additional arguments.
#'
#' @return A matrix of dimension \code{1 x 2} with confidence limits.
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

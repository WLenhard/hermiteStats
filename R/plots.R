#' @title Diagnostic and Interpretive Visualizations
#'
#' @description
#' Base-graphics plotting methods for the three primary model objects in
#' \code{hermiteStats}: \code{\link{hermite_fit}}, \code{\link{cor_hermite}},
#' and \code{\link{d_reg}}. Every plot is designed to make the regularization
#' step transparent by displaying the raw empirical data alongside the fitted or
#' implied regularized curves, allowing the plausibility of the underlying
#' polynomial quantile models to be assessed visually.
#'
#' @name hermiteStats_plots
#' @keywords internal
NULL

# -----------------------------------------------------------------------------
# Internal Helpers: shared color palette and quantile-curve prediction
# -----------------------------------------------------------------------------

#' Shared color palette for hermiteStats plots
#' @noRd
.hermite_pal <- list(
  g1      = "#1b9e77",  # teal   (Group 1 / X)
  g2      = "#d95f02",  # orange (Group 2 / Y)
  fit     = "firebrick",
  ref     = "gray70",
  grid    = "gray90",
  annot   = "gray20"
)

#' Evaluate a fitted hermite_fit polynomial on a grid of Z values
#'
#' @param fit An object of class \code{"hermite_fit"}.
#' @param z_grid Optional numeric vector of Z values. Defaults to an evenly
#'   spaced grid spanning the observed range of \code{fit$z}.
#' @param n_grid Integer; number of grid points if \code{z_grid} is not supplied.
#' @return A list with components \code{z} and \code{x} (predicted raw values).
#' @noRd
.hermite_fit_predict <- function(fit, z_grid = NULL, n_grid = 200L) {
  if (is.null(z_grid)) {
    z_grid <- seq(min(fit$z), max(fit$z), length.out = n_grid)
  }
  Zmat <- outer(z_grid, 0:fit$degree, `^`)
  list(z = z_grid, x = as.vector(Zmat %*% fit$beta))
}

#' Evaluate Probabilists' Hermite Polynomials on a Grid
#'
#' Evaluates \eqn{He_0(z), He_1(z), \dots, He_{\text{degree}}(z)} at each
#' element of \code{z} via the three-term recurrence
#' \eqn{He_n(z) = z\,He_{n-1}(z) - (n-1)\,He_{n-2}(z)}.
#'
#' @param z Numeric vector of evaluation points.
#' @param degree Integer; highest Hermite polynomial order required.
#' @return A numeric matrix with \code{length(z)} rows and \code{degree + 1}
#'   columns; column \code{m + 1} holds \eqn{He_m(z)}.
#' @noRd
.hermite_polynomial_matrix <- function(z, degree) {
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

#' Model-Implied Conditional Mean Curve Under the Fitted Gaussian Copula
#'
#' Computes the exact, closed-form conditional mean \eqn{\mathbb{E}[Y \mid X = x]}
#' implied by the fitted marginal quantile functions and latent copula correlation
#' \eqn{\rho_z} via the Hermite reproducing property of the Mehler kernel:
#' \deqn{\mathbb{E}\left[He_m(Z_y) \mid Z_x = z_x\right] = \rho_z^m\, He_m(z_x)}
#'
#' @param fit_x,fit_y Objects of class \code{"hermite_fit"}.
#' @param rho_z Numeric; the latent Gaussian copula correlation.
#' @param n_grid Integer; number of evaluation points along the latent range of \code{fit_x$z}.
#' @return A list with components \code{x} and \code{y}.
#' @noRd
.mehler_conditional_mean <- function(fit_x, fit_y, rho_z, n_grid = 200L) {
  z_grid <- seq(min(fit_x$z), max(fit_x$z), length.out = n_grid)
  x_vals <- .hermite_fit_predict(fit_x, z_grid)$x

  dy   <- fit_y$degree
  He_y <- .hermite_polynomial_matrix(z_grid, dy)
  y_vals <- as.vector(He_y %*% (fit_y$hermite_coeffs * (rho_z^(0:dy))))

  list(x = x_vals, y = y_vals)
}

# -----------------------------------------------------------------------------
# plot.hermite_fit
# -----------------------------------------------------------------------------

#' Plot a Fitted Hermite Quantile Model
#'
#' Displays the empirical (rank-based) standard normal scores against the raw
#' observations, overlaid with the fitted regularized polynomial quantile
#' function \eqn{X = f(Z)}, and annotates the regularized moments.
#'
#' @param x An object of class \code{"hermite_fit"}.
#' @param main Character; plot title (realized degree is appended automatically).
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' fit <- hermite_fit(rlnorm(60, meanlog = 2, sdlog = 0.4))
#' plot(fit)
#'
#' @export
plot.hermite_fit <- function(x, main = "Regularized Quantile Map", ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  pred <- .hermite_fit_predict(x)

  plot(x$z, x$x, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor("gray25", alpha.f = 0.45),
       xlab = "Standard Normal Score (Z)", ylab = "Raw Value (X)",
       main = sprintf("%s (degree = %d)", main, x$degree),
       bty = "l", panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1),
       ...)
  graphics::lines(pred$z, pred$x, col = .hermite_pal$fit, lwd = 2.5)

  graphics::mtext(
    sprintf("mu = %.2f  |  sigma = %.2f  |  skew = %.2f  |  exc. kurt. = %.2f",
            x$mean, x$sd, x$skewness, x$excess_kurtosis),
    side = 3, line = 0.3, cex = 0.8, col = .hermite_pal$annot
  )
  invisible(x)
}

# -----------------------------------------------------------------------------
# plot.cor_hermite
# -----------------------------------------------------------------------------

#' Plot a Hermite Correlation Fit
#'
#' Displays two diagnostic panels:
#' \enumerate{
#'   \item \strong{Latent Normal Score Association:} Plots \eqn{Z_x} against \eqn{Z_y}.
#'         For \code{copula = "gaussian"}, overlays the latent linear copula correlation \eqn{\rho_z};
#'         for \code{copula = "none"}, displays the empirical rank-based dependence.
#'   \item \strong{Manifest Association:} Displays the raw scatterplot annotated with
#'         \eqn{r_{\mathrm{Hermite}}}. For \code{copula = "gaussian"}, overlays both the
#'         model-implied conditional mean curve \eqn{\mathbb{E}[Y \mid X=x]} (exact via Mehler's identity)
#'         and a non-parametric lowess smooth, displaying the shape-attenuation factor \eqn{A}.
#'         For \code{copula = "none"}, displays the empirical lowess smooth without parametric assumptions.
#' }
#'
#' @param x An object of class \code{"cor_hermite"}.
#' @param ... Additional graphical parameters passed to the latent scatter plot.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' # 1. Copula-free fit
#' x <- rnorm(50, 100, 15)
#' y <- exp(0.5 * scale(x) + rnorm(50, sd = 0.5))
#' fit_free <- cor_hermite(x, y, copula = "none")
#' plot(fit_free)
#'
#' # 2. Gaussian copula fit with Mehler reference curve
#' fit_gauss <- cor_hermite(x, y, copula = "gaussian")
#' plot(fit_gauss)
#'
#' @export
plot.cor_hermite <- function(x, ...) {
  if (is.null(x$fit_x) || is.null(x$fit_y)) {
    stop("Cannot plot: this 'cor_hermite' object has no fitted quantile models ",
         "(likely created with fewer than 4 complete observations).")
  }

  copula_mode <- if (!is.null(x$copula)) x$copula else (if (!is.null(x$rho_z) && !is.na(x$rho_z)) "gaussian" else "none")

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  # --- Panel 1: Latent Normal Scores -----------------------------------------
  zx <- x$fit_x$z; zy <- x$fit_y$z
  lim_z <- range(zx, zy)

  main_p1 <- if (copula_mode == "gaussian") "Latent Gaussian Copula" else "Latent Normal Scores (Copula-Free)"

  plot(zx, zy, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor(.hermite_pal$g1, alpha.f = 0.45),
       xlim = lim_z, ylim = lim_z,
       xlab = "Normal Score Z(X)", ylab = "Normal Score Z(Y)",
       main = main_p1, bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1), ...)
  graphics::abline(a = 0, b = 1, col = .hermite_pal$ref, lty = 3)

  if (copula_mode == "gaussian" && is.finite(x$rho_z)) {
    graphics::abline(a = 0, b = x$rho_z, col = "darkblue", lwd = 2)
    graphics::mtext(sprintf("rho_z = %.3f", x$rho_z),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  } else {
    emp_rz <- stats::cor(zx, zy)
    graphics::abline(a = 0, b = emp_rz, col = "darkblue", lwd = 2, lty = 2)
    graphics::mtext(sprintf("Empirical r_z = %.3f", emp_rz),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  }

  # --- Panel 2: Manifest (Raw-Scale) Association ------------------------------
  plot(x$x, x$y, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor(.hermite_pal$g2, alpha.f = 0.45),
       xlab = "Raw X", ylab = "Raw Y",
       main = "Manifest Association", bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1))

  # 1. Non-parametric empirical trend
  sm <- stats::lowess(x$x, x$y)

  if (copula_mode == "gaussian" && is.finite(x$rho_z)) {
    # Gaussian Copula: Model-implied conditional curve (Mehler) vs. Lowess
    mh <- .mehler_conditional_mean(x$fit_x, x$fit_y, x$rho_z)
    ord <- order(mh$x)
    graphics::lines(mh$x[ord], mh$y[ord], col = "darkblue", lwd = 2.5, lty = 1)
    graphics::lines(sm, col = "gray30", lwd = 2, lty = 2)

    graphics::legend("topleft",
                     legend = c("Model-implied (Mehler)", "Empirical (lowess)"),
                     col = c("darkblue", "gray30"), lwd = c(2.5, 2), lty = c(1, 2),
                     bty = "n", cex = 0.8)

    graphics::mtext(sprintf("r_Hermite = %.3f  |  A = %.2f", x$r_Hermite, x$attenuation),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)

  } else {
    # Copula-Free: Regularized Linear Regression Line vs. Lowess
    beta_Hermite  <- x$cov_xy / x$var_x
    alpha_Hermite <- x$mean_y - beta_Hermite * x$mean_x

    # Draw regularized linear line across the observed range of X
    x_grid <- seq(min(x$x), max(x$x), length.out = 100L)
    y_pred <- alpha_Hermite + beta_Hermite * x_grid
    graphics::lines(x_grid, y_pred, col = "darkblue", lwd = 2.5, lty = 1)
    graphics::lines(sm, col = "gray30", lwd = 2, lty = 2)

    graphics::legend("topleft",
                     legend = c("Regularized Line (r_Hermite)", "Empirical (lowess)"),
                     col = c("darkblue", "gray30"), lwd = c(2.5, 2), lty = c(1, 2),
                     bty = "n", cex = 0.8)

    trim_txt <- if (!is.null(x$trim) && x$trim > 0) sprintf(" (trim = %.2f)", x$trim) else ""
    graphics::mtext(sprintf("r_Hermite = %.3f%s", x$r_Hermite, trim_txt),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  }
  invisible(x)
}

# -----------------------------------------------------------------------------
# plot.d_reg
# -----------------------------------------------------------------------------

#' Plot a Distribution-Free Effect Size Fit
#'
#' Displays two diagnostic panels:
#' \enumerate{
#'   \item \strong{Group Distributions:} Overlays the empirical kernel densities of
#'         both groups with regularized group means marked by dashed vertical lines,
#'         annotating the effect size estimate and confidence interval.
#'   \item \strong{Regularized Quantile Maps:} Displays both groups' fitted monotone
#'         polynomial quantile curves \eqn{X = f(Z)} on a shared latent normal axis,
#'         providing a direct visual inspection of differences in location, variance,
#'         and shape asymmetry.
#' }
#'
#' @param x An object of class \code{"d_reg"}.
#' @param ... Currently unused (present for S3 method consistency).
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' ctrl <- rlnorm(30, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(30, meanlog = 2.3, sdlog = 0.5)
#' fit <- d_reg(ctrl, trt, conf_level = 0.95)
#' plot(fit)
#'
#' @export
plot.d_reg <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  g1_name <- if (!is.null(x$group_labels)) x$group_labels[1L] else "Group 1"
  g2_name <- if (!is.null(x$group_labels)) x$group_labels[2L] else "Group 2"
  col1 <- .hermite_pal$g1; col2 <- .hermite_pal$g2

  # --- Panel 1: Empirical densities with regularized means -------------------
  d1 <- stats::density(x$x1); d2 <- stats::density(x$x2)
  xlim <- range(d1$x, d2$x)
  ylim <- range(0, d1$y, d2$y)

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Raw Value", ylab = "Density", main = "Group Distributions",
       bty = "l", panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1))
  graphics::polygon(d1, col = grDevices::adjustcolor(col1, alpha.f = 0.15), border = NA)
  graphics::polygon(d2, col = grDevices::adjustcolor(col2, alpha.f = 0.15), border = NA)
  graphics::lines(d1, col = col1, lwd = 2.5)
  graphics::lines(d2, col = col2, lwd = 2.5)
  graphics::rug(x$x1, side = 1, col = grDevices::adjustcolor(col1, alpha.f = 0.6))
  graphics::rug(x$x2, side = 3, col = grDevices::adjustcolor(col2, alpha.f = 0.6))
  graphics::abline(v = x$group1$mean, col = col1, lty = 2, lwd = 1.5)
  graphics::abline(v = x$group2$mean, col = col2, lty = 2, lwd = 1.5)
  graphics::legend("topright", legend = c(g1_name, g2_name), col = c(col1, col2),
                   lwd = 2.5, bty = "n", cex = 0.85)

  est_lab <- if (x$paired) "d_z" else "d_reg"
  est_txt <- sprintf("%s = %.3f", est_lab, x$estimate)
  if (!is.null(x$ci)) {
    est_txt <- sprintf("%s   [%.3f, %.3f]", est_txt, x$ci[1L], x$ci[2L])
  }
  graphics::mtext(est_txt, side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)

  # --- Panel 2: Combined regularized quantile maps ----------------------------
  zlim  <- range(x$fit1$z, x$fit2$z)
  ylim2 <- range(x$fit1$x, x$fit2$x)

  plot(NA, xlim = zlim, ylim = ylim2,
       xlab = "Standard Normal Score (Z)", ylab = "Raw Value",
       main = "Regularized Quantile Maps", bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1))
  graphics::points(x$fit1$z, x$fit1$x, pch = 19, cex = 0.7,
                   col = grDevices::adjustcolor(col1, alpha.f = 0.35))
  graphics::points(x$fit2$z, x$fit2$x, pch = 19, cex = 0.7,
                   col = grDevices::adjustcolor(col2, alpha.f = 0.35))
  p1 <- .hermite_fit_predict(x$fit1)
  p2 <- .hermite_fit_predict(x$fit2)
  graphics::lines(p1$z, p1$x, col = col1, lwd = 2.5)
  graphics::lines(p2$z, p2$x, col = col2, lwd = 2.5)
  graphics::legend("topleft", legend = c(g1_name, g2_name), col = c(col1, col2),
                   lwd = 2.5, bty = "n", cex = 0.85)

  invisible(x)
}

# -----------------------------------------------------------------------------
# plot.t_hermite
# -----------------------------------------------------------------------------

#' Plot a Hermite t-Test Object
#'
#' Displays the permutation null distribution of the Hermite \eqn{t}-statistic,
#' marking the observed test statistic with a solid vertical line, two-sided
#' critical boundaries with dashed lines, and annotating the \eqn{p}-value.
#'
#' @param x An object of class \code{"t_hermite"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' set.seed(42)
#' ctrl <- rlnorm(25, 2.0, 0.5)
#' trt  <- rlnorm(25, 2.3, 0.5)
#' res <- t_hermite(ctrl, trt, nperm = 500)
#' plot(res)
#'
#' @export
plot.t_hermite <- function(x, ...) {
  if (is.null(x$perm_distribution) || length(x$perm_distribution) == 0L) {
    stop("No permutation distribution available to plot.")
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  d <- stats::density(x$perm_distribution)
  xlim_vals <- range(d$x, x$statistic, -x$statistic)
  xlim <- xlim_vals + c(-1, 1) * 0.1 * diff(xlim_vals)

  plot(d, main = paste0("Permutation Null Distribution (", x$method, ")"),
       xlab = "t_Hermite (under H0)", bty = "l", col = "darkblue", lwd = 2,
       xlim = xlim, panel.first = graphics::grid(col = "gray90", lty = 1), ...)

  graphics::polygon(d, col = grDevices::adjustcolor("steelblue", alpha.f = 0.2), border = NA)
  graphics::abline(v = x$statistic, col = "firebrick", lwd = 2.5, lty = 1)

  if (x$alternative == "two.sided") {
    graphics::abline(v = -x$statistic, col = "firebrick", lwd = 1.5, lty = 2)
  }

  graphics::mtext(sprintf("Observed t = %.3f  |  p = %.3f (%s)",
                          x$statistic, x$p_value, x$alternative),
                  side = 3, line = 0.3, cex = 0.85, col = "gray20")

  graphics::legend("topright",
                   legend = c("Permutation Null", sprintf("Observed t = %.3f", x$statistic)),
                   col = c("darkblue", "firebrick"), lwd = c(2, 2.5), lty = c(1, 1),
                   bty = "n", cex = 0.85)

  invisible(x)
}

# -----------------------------------------------------------------------------
# plot.shape_hermite
# -----------------------------------------------------------------------------

#' Plot a Hermite Shape-Difference Test Object
#'
#' Displays a 3-panel visualization showing the permutation null distributions
#' for \strong{Variance}, \strong{Skewness}, and \strong{Excess Kurtosis} differences,
#' marking observed differences and annotating raw and Holm-adjusted \eqn{p}-values.
#'
#' @param x An object of class \code{"shape_hermite"}.
#' @param ... Additional graphical parameters.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' set.seed(42)
#' g1 <- rnorm(40)
#' g2 <- rgamma(40, shape = 2)
#' sres <- shape_hermite(g1, g2, nperm = 500)
#' plot(sres)
#'
#' @export
plot.shape_hermite <- function(x, ...) {
  if (is.null(x$perm_distribution) || nrow(x$perm_distribution) == 0L) {
    stop("No permutation distribution available to plot.")
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.0, 3.4, 1.0), font.main = 1)

  moment_titles <- c(variance        = "Variance Difference",
                     skewness        = "Skewness Difference",
                     excess_kurtosis = "Excess Kurtosis Difference")

  pal_col <- c("#1b9e77", "#d95f02", "#7570b3")

  for (i in seq_along(moment_titles)) {
    m_name   <- names(moment_titles)[i]
    perm_vals <- x$perm_distribution[, m_name]
    obs_val  <- x$estimate[m_name]
    p_val    <- x$p_value[m_name]
    p_adj    <- x$p_adjusted[m_name]

    d <- stats::density(perm_vals)
    xlim_vals <- range(d$x, obs_val, -obs_val)
    xlim <- xlim_vals + c(-1, 1) * 0.1 * diff(xlim_vals)

    plot(d, main = moment_titles[m_name],
         xlab = paste0("Delta ", m_name, " (under H0)"),
         ylab = if (i == 1L) "Density" else "",
         bty = "l", col = pal_col[i], lwd = 2, xlim = xlim,
         panel.first = graphics::grid(col = "gray90", lty = 1), ...)

    graphics::polygon(d, col = grDevices::adjustcolor(pal_col[i], alpha.f = 0.2), border = NA)
    graphics::abline(v = obs_val, col = "firebrick", lwd = 2.5, lty = 1)
    graphics::abline(v = -obs_val, col = "firebrick", lwd = 1.2, lty = 2)

    graphics::mtext(sprintf("Diff = %.3f | p = %.3f (adj: %.3f)", obs_val, p_val, p_adj),
                    side = 3, line = 0.3, cex = 0.8, col = "gray20")
  }

  invisible(x)
}

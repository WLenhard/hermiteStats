# =============================================================================
# plots.R -- Diagnostic and Interpretive Visualizations
# =============================================================================

#' @title Diagnostic and Interpretive Visualizations
#'
#' @description
#' Base-graphics plotting methods for the primary estimation objects in
#' \pkg{hermiteStats}: \code{\link{hermite_fit}}, \code{\link{cor_hermite}},
#' and \code{\link{d_reg}}. Every plot is designed to make the regularization
#' step transparent by displaying the raw empirical data alongside the fitted
#' or implied regularized curves, so the plausibility of the underlying
#' polynomial quantile models can be assessed visually.
#'
#' Plot methods for the hypothesis-test objects (\code{\link{t_hermite}},
#' \code{\link{median_hermite}}, \code{\link{shape_hermite}},
#' \code{\link{hermite_test}}) are documented alongside those functions.
#'
#' All numerical work (polynomial evaluation, Hermite basis construction) is
#' delegated to the quantile engine; this file contains visualization code
#' only.
#'
#' @name hermiteStats_plots
#' @keywords internal
NULL

# -----------------------------------------------------------------------------
# Internal Helpers: palette, prediction grid, Mehler conditional mean
# -----------------------------------------------------------------------------

#' Shared color palette for hermiteStats plots
#' @noRd
.hermite_pal <- list(
  g1    = "#1b9e77",   # teal   (Group 1 / X)
  g2    = "#d95f02",   # orange (Group 2 / Y)
  fit   = "firebrick",
  model = "darkblue",
  ref   = "gray70",
  grid  = "gray90",
  annot = "gray20"
)

#' Fitted quantile curve on an even grid spanning the observed latent range
#'
#' Thin wrapper around \code{\link{predict.hermite_fit}} used by all plot
#' methods that draw the regularized quantile function.
#'
#' @param fit An object of class \code{"hermite_fit"}.
#' @param n_grid Integer; number of grid points (default \code{200L}).
#' @return A list with components \code{z} (grid) and \code{x} (predictions).
#' @noRd
.hermite_fit_grid <- function(fit, n_grid = 200L) {
  z_grid <- seq(min(fit$z), max(fit$z), length.out = n_grid)
  list(z = z_grid, x = predict(fit, z = z_grid))
}

#' Model-implied conditional mean curve under the fitted Gaussian copula
#'
#' Computes the exact, closed-form conditional mean \eqn{E[Y | X = x]}
#' implied by the fitted marginal quantile functions and the latent copula
#' correlation \eqn{\rho_z}, via the Hermite reproducing property of the
#' Mehler kernel:
#' \deqn{E[He_m(Z_y) | Z_x = z_x] = \rho_z^m \, He_m(z_x)}
#'
#' @param fit_x,fit_y Objects of class \code{"hermite_fit"}.
#' @param rho_z Numeric; latent Gaussian copula correlation.
#' @param n_grid Integer; number of evaluation points.
#' @return A list with components \code{x} and \code{y}.
#' @noRd
.mehler_conditional_mean <- function(fit_x, fit_y, rho_z, n_grid = 200L) {
  z_grid <- seq(min(fit_x$z), max(fit_x$z), length.out = n_grid)
  x_vals <- predict(fit_x, z = z_grid)

  dy     <- fit_y$degree
  He_y   <- .hermite_eval_basis(z_grid, dy)
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
#' @param main Character; plot title (the realized degree is appended
#'   automatically).
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

  pred <- .hermite_fit_grid(x)

  plot(x$z, x$x, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor("gray25", alpha.f = 0.45),
       xlab = "Standard Normal Score (Z)", ylab = "Raw Value (X)",
       main = sprintf("%s (degree = %d)", main, x$degree),
       bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1),
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
#'   \item \strong{Latent Normal Score Association:} \eqn{Z_x} against
#'         \eqn{Z_y}. For \code{copula = "gaussian"}, overlays the latent
#'         linear copula correlation \eqn{\rho_z}; for \code{copula = "none"},
#'         the empirical rank-based dependence.
#'   \item \strong{Manifest Association:} the raw scatterplot annotated with
#'         \eqn{r_{\mathrm{Hermite}}}. For \code{copula = "gaussian"}, overlays
#'         both the model-implied conditional mean curve \eqn{E[Y | X = x]}
#'         (exact via Mehler's identity) and a non-parametric lowess smooth,
#'         and displays the shape-attenuation factor \eqn{A}. For
#'         \code{copula = "none"}, overlays the regularized linear regression
#'         line and the empirical lowess smooth.
#' }
#'
#' @param x An object of class \code{"cor_hermite"}.
#' @param ... Additional graphical parameters passed to the latent scatter plot.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(50, 100, 15)
#' y <- exp(0.5 * scale(x) + rnorm(50, sd = 0.5))
#'
#' # Copula-free fit
#' plot(cor_hermite(x, y, copula = "none"))
#'
#' # Gaussian copula fit with Mehler reference curve
#' plot(cor_hermite(x, y, copula = "gaussian"))
#'
#' @export
plot.cor_hermite <- function(x, ...) {
  if (is.null(x$fit_x) || is.null(x$fit_y)) {
    stop("Cannot plot: this 'cor_hermite' object has no fitted quantile models ",
         "(likely created with fewer than 4 complete observations).")
  }

  copula_mode <- if (!is.null(x$copula)) {
    x$copula
  } else if (!is.null(x$rho_z) && !is.na(x$rho_z)) "gaussian" else "none"

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.4, 1.2), font.main = 1)

  # --- Panel 1: Latent normal scores ----------------------------------------
  zx <- x$fit_x$z; zy <- x$fit_y$z
  lim_z <- range(zx, zy)

  main_p1 <- if (copula_mode == "gaussian") {
    "Latent Gaussian Copula"
  } else {
    "Latent Normal Scores (Copula-Free)"
  }

  plot(zx, zy, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor(.hermite_pal$g1, alpha.f = 0.45),
       xlim = lim_z, ylim = lim_z,
       xlab = "Normal Score Z(X)", ylab = "Normal Score Z(Y)",
       main = main_p1, bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1), ...)
  graphics::abline(a = 0, b = 1, col = .hermite_pal$ref, lty = 3)

  if (copula_mode == "gaussian" && is.finite(x$rho_z)) {
    graphics::abline(a = 0, b = x$rho_z, col = .hermite_pal$model, lwd = 2)
    graphics::mtext(sprintf("rho_z = %.3f", x$rho_z),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  } else {
    emp_rz <- stats::cor(zx, zy)
    graphics::abline(a = 0, b = emp_rz, col = .hermite_pal$model, lwd = 2, lty = 2)
    graphics::mtext(sprintf("Empirical r_z = %.3f", emp_rz),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  }

  # --- Panel 2: Manifest (raw-scale) association ------------------------------
  plot(x$x, x$y, pch = 19, cex = 0.75,
       col = grDevices::adjustcolor(.hermite_pal$g2, alpha.f = 0.45),
       xlab = "Raw X", ylab = "Raw Y",
       main = "Manifest Association", bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1))

  sm <- stats::lowess(x$x, x$y)   # non-parametric empirical trend

  if (copula_mode == "gaussian" && is.finite(x$rho_z)) {
    # Model-implied conditional mean curve (Mehler) vs. lowess
    mh  <- .mehler_conditional_mean(x$fit_x, x$fit_y, x$rho_z)
    ord <- order(mh$x)
    graphics::lines(mh$x[ord], mh$y[ord], col = .hermite_pal$model, lwd = 2.5)
    graphics::lines(sm, col = "gray30", lwd = 2, lty = 2)

    graphics::legend("topleft",
                     legend = c("Model-implied (Mehler)", "Empirical (lowess)"),
                     col = c(.hermite_pal$model, "gray30"),
                     lwd = c(2.5, 2), lty = c(1, 2), bty = "n", cex = 0.8)
    graphics::mtext(sprintf("r_Hermite = %.3f  |  A = %.2f",
                            x$r_Hermite, x$attenuation),
                    side = 3, line = 0.3, cex = 0.85, col = .hermite_pal$annot)
  } else {
    # Regularized linear regression line vs. lowess
    beta_h  <- x$cov_xy / x$var_x
    alpha_h <- x$mean_y - beta_h * x$mean_x
    x_grid  <- seq(min(x$x), max(x$x), length.out = 100L)
    graphics::lines(x_grid, alpha_h + beta_h * x_grid,
                    col = .hermite_pal$model, lwd = 2.5)
    graphics::lines(sm, col = "gray30", lwd = 2, lty = 2)

    graphics::legend("topleft",
                     legend = c("Regularized Line (r_Hermite)", "Empirical (lowess)"),
                     col = c(.hermite_pal$model, "gray30"),
                     lwd = c(2.5, 2), lty = c(1, 2), bty = "n", cex = 0.8)

    trim_txt <- if (!is.null(x$trim) && x$trim > 0) {
      sprintf(" (trim = %.2f)", x$trim)
    } else ""
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
#'   \item \strong{Group Distributions:} overlaid empirical kernel densities
#'         with the regularized group means marked by dashed vertical lines,
#'         annotated with the effect size estimate and confidence interval
#'         (if available).
#'   \item \strong{Regularized Quantile Maps:} both groups' fitted monotone
#'         polynomial quantile curves \eqn{X = f(Z)} on a shared latent normal
#'         axis, allowing direct visual comparison of location, scale, and
#'         shape asymmetry.
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
  col1 <- .hermite_pal$g1
  col2 <- .hermite_pal$g2

  # --- Panel 1: Empirical densities with regularized means -------------------
  d1 <- stats::density(x$x1)
  d2 <- stats::density(x$x2)

  plot(NA, xlim = range(d1$x, d2$x), ylim = range(0, d1$y, d2$y),
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
  graphics::legend("topright", legend = c(g1_name, g2_name),
                   col = c(col1, col2), lwd = 2.5, bty = "n", cex = 0.85)

  est_lab <- if (x$paired) "d_z" else "d_reg"
  est_txt <- sprintf("%s = %.3f", est_lab, x$estimate)
  if (!is.null(x$ci)) {
    est_txt <- sprintf("%s   [%.3f, %.3f]", est_txt, x$ci[1L], x$ci[2L])
  }
  graphics::mtext(est_txt, side = 3, line = 0.3, cex = 0.85,
                  col = .hermite_pal$annot)

  # --- Panel 2: Combined regularized quantile maps ----------------------------
  p1 <- .hermite_fit_grid(x$fit1)
  p2 <- .hermite_fit_grid(x$fit2)

  plot(NA, xlim = range(x$fit1$z, x$fit2$z), ylim = range(x$fit1$x, x$fit2$x),
       xlab = "Standard Normal Score (Z)", ylab = "Raw Value",
       main = "Regularized Quantile Maps", bty = "l",
       panel.first = graphics::grid(col = .hermite_pal$grid, lty = 1))
  graphics::points(x$fit1$z, x$fit1$x, pch = 19, cex = 0.7,
                   col = grDevices::adjustcolor(col1, alpha.f = 0.35))
  graphics::points(x$fit2$z, x$fit2$x, pch = 19, cex = 0.7,
                   col = grDevices::adjustcolor(col2, alpha.f = 0.35))
  graphics::lines(p1$z, p1$x, col = col1, lwd = 2.5)
  graphics::lines(p2$z, p2$x, col = col2, lwd = 2.5)
  graphics::legend("topleft", legend = c(g1_name, g2_name),
                   col = c(col1, col2), lwd = 2.5, bty = "n", cex = 0.85)

  invisible(x)
}

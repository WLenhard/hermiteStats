#' Diagnostic Visualizations for hermiteStats Objects
#'
#' @param x An object of class \code{"d_reg"}, \code{"cor_hermite"}, or \code{"hermite_fit"}.
#' @param ... Graphical parameters passed to base plotting functions.
#'
#' @export
plot.d_reg <- function(x, ...) {
  # Base R diagnostic plot showing empirical densities and modeled quantiles
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 2))

  # Plot Group 1
  plot(x$fit1, main = "Group 1 Quantile Map", ...)
  # Plot Group 2
  plot(x$fit2, main = "Group 2 Quantile Map", ...)
  invisible(x)
}

#' @export
plot.cor_hermite <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 2))

  # Latent Copula Scatter
  plot(x$fit_x$z, x$fit_y$z, pch = 19, col = grDevices::adjustcolor("steelblue", alpha.f = 0.6),
       xlab = "Normal Scores Z(X)", ylab = "Normal Scores Z(Y)",
       main = sprintf("Latent Copula (rho_z = %.3f)", x$rho_z), ...)
  graphics::abline(h = 0, v = 0, lty = 2, col = "gray60")
  graphics::abline(a = 0, b = x$rho_z, col = "darkblue", lwd = 2)

  # Manifest Scaled Scatter
  plot(x$x, x$y, pch = 19, col = grDevices::adjustcolor("coral", alpha.f = 0.6),
       xlab = "Raw X", ylab = "Raw Y",
       main = sprintf("Manifest r_HM = %.3f (A = %.2f)", x$r_hm, x$attenuation), ...)
  invisible(x)
}

#' @export
plot.hermite_fit <- function(x, main = "Fitted Quantile Map", ...) {
  z_grid <- seq(min(x$z), max(x$z), length.out = 200L)
  Zmat <- outer(z_grid, 0:x$degree, `^`)
  x_pred <- as.vector(Zmat %*% x$beta)

  plot(x$z, x$x, pch = 19, col = grDevices::adjustcolor("black", alpha.f = 0.4),
       xlab = "Standard Normal Score (Z)", ylab = "Raw Values (X)",
       main = paste0(main, " (Degree = ", x$degree, ")"), ...)
  graphics::lines(z_grid, x_pred, col = "firebrick", lwd = 2.5)
  invisible(x)
}

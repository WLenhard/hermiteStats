#' @title Distribution-Free Effect Size Estimation (d_reg)
#'
#' @description Computes the distribution-robust standardized mean difference (\eqn{d_{\mathrm{reg}}})
#' for independent or paired samples using monotone polynomial quantile modeling.
#'
#' @param x A numeric vector or a \code{formula} of the form \code{response ~ group}.
#' @param y A numeric vector of group 2 scores (if \code{x} is a vector).
#' @param data An optional data frame if a formula is supplied.
#' @param degree Polynomial degree for quantile fitting (default = 3).
#' @param monotonicity Monotonicity check: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param paired Logical; if \code{TRUE}, evaluates paired-samples designs (default = \code{FALSE}).
#' @param type Calculation type: \code{"regularized"} (default), \code{"hedges"}, \code{"glass"}, or \code{"combined"}.
#' @param conf_level Optional confidence level (e.g., 0.95), or \code{NULL}.
#' @param ci_method Method for CI: \code{"bootstrap"} (default) or \code{"nct"} (noncentral t).
#' @param B Number of bootstrap resamples (default = 1000).
#' @param ... Additional arguments.
#'
#' @return An S3 object of class \code{"d_reg"}.
#' @export
d_reg <- function(x, ...) {
  UseMethod("d_reg")
}

#' @rdname d_reg
#' @export
d_reg.formula <- function(formula, data = NULL, degree = 3L,
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  if (missing(formula) || (length(formula) != 3L)) {
    stop("Formula must be of the form 'response ~ group'.")
  }
  mf <- stats::model.frame(formula = formula, data = data)
  response <- mf[[1L]]
  group <- as.factor(mf[[2L]])

  levels_g <- levels(group)
  if (length(levels_g) != 2L) {
    stop("Grouping variable must have exactly two levels.")
  }

  x1 <- response[group == levels_g[1L]]
  x2 <- response[group == levels_g[2L]]

  res <- d_reg.default(x = x1, y = x2, degree = degree, monotonicity = monotonicity,
                       paired = paired, type = type, conf_level = conf_level,
                       ci_method = ci_method, B = B, ...)
  res$group_labels <- levels_g
  res
}

#' @rdname d_reg
#' @export
d_reg.default <- function(x, y = NULL, degree = 3L,
                          monotonicity = c("relaxed", "strict", "none"),
                          paired = FALSE,
                          type = c("regularized", "hedges", "glass", "combined"),
                          conf_level = NULL,
                          ci_method = c("bootstrap", "nct"),
                          B = 1000L, ...) {

  monotonicity <- match.arg(monotonicity)
  type         <- match.arg(type)
  ci_method    <- match.arg(ci_method)

  if (missing(y) || is.null(y)) stop("Vector 'y' must be supplied.")

  if (paired) {
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
    if (length(x1) != length(x2)) stop("Paired observations must have equal length.")
  } else {
    x1 <- x[is.finite(x)]
    x2 <- y[is.finite(y)]
  }

  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 3L || n2 < 3L) stop("Each group must contain at least 3 valid observations.")

  # Fit marginal polynomial quantile distributions
  fit1 <- hermite_fit(x1, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
  fit2 <- hermite_fit(x2, degree = degree, monotonicity = monotonicity, force_odd = TRUE)

  m1 <- hermite_moments(fit1)
  m2 <- hermite_moments(fit2)

  # Averaged standardizer for independent groups (Cohen 1988, Delacre 2021)
  sd_avg <- sqrt((m1$variance + m2$variance) / 2.0)

  # If paired, compute difference distribution and latent correlation
  if (paired) {
    diff_scores <- x2 - x1
    fit_diff <- hermite_fit(diff_scores, degree = degree, monotonicity = monotonicity, force_odd = TRUE)
    m_diff <- hermite_moments(fit_diff)
    r_hm_fit <- cor_hermite(x1, x2, poly_degree = degree, monotonicity = monotonicity)

    # d_z: standardized mean change; d_rm: repeated measures raw SMD
    d_z <- m_diff$mean / m_diff$sd
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  } else {
    fit_diff <- NULL; m_diff <- NULL; r_hm_fit <- NULL; d_z <- NULL
    d_reg_val <- (m2$mean - m1$mean) / sd_avg
  }

  # Classical benchmarks
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sd_pooled <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
  j_corr <- .J_correction(n1 + n2 - 2L)
  g_hedges <- ((mean(x2) - mean(x1)) / sd_pooled) * j_corr
  glass_delta <- (mean(x2) - mean(x1)) / s1

  # Type routing
  primary_d <- switch(
    type,
    "regularized" = if (paired) d_z else d_reg_val,
    "hedges"      = g_hedges,
    "glass"       = glass_delta,
    "combined"    = if (n1 > 50L && n2 > 50L && abs(d_reg_val) > 0.8) g_hedges else d_reg_val
  )

  res <- list(
    estimate = primary_d,
    d_reg = d_reg_val,
    d_z = d_z,
    hedges_g = g_hedges,
    glass_delta = glass_delta,
    type = type,
    paired = paired,
    n1 = n1, n2 = n2,
    group1 = list(mean = m1$mean, sd = m1$sd, variance = m1$variance, raw_mean = mean(x1), raw_sd = s1),
    group2 = list(mean = m2$mean, sd = m2$sd, variance = m2$variance, raw_mean = mean(x2), raw_sd = s2),
    diff_moments = m_diff,
    r_hm_paired = r_hm_fit,
    sd_standardizer = if (type == "glass") s1 else if (type == "hedges") sd_pooled else sd_avg,
    degrees = c(g1 = fit1$degree, g2 = fit2$degree),
    monotonicity = monotonicity,
    fit1 = fit1, fit2 = fit2,
    x1 = x1, x2 = x2
  )
  class(res) <- "d_reg"

  if (!is.null(conf_level)) {
    res$ci <- stats::confint(res, level = conf_level, method = ci_method, B = B)
    res$conf_level <- conf_level
    res$ci_method <- ci_method
  }
  res
}

#' Confidence Intervals for d_reg Objects
#' @export
confint.d_reg <- function(object, parm, level = 0.95,
                          method = c("bootstrap", "nct"), B = 1000L, ...) {
  method <- match.arg(method)
  n1 <- object$n1; n2 <- object$n2

  if (method == "bootstrap") {
    x1 <- object$x1; x2 <- object$x2
    boot_d <- numeric(B)
    for (b in seq_len(B)) {
      if (object$paired) {
        idx <- sample.int(n1, n1, replace = TRUE)
        fb <- d_reg(x1[idx], x2[idx], degree = max(object$degrees),
                    monotonicity = object$monotonicity, paired = TRUE, type = object$type)
      } else {
        s1 <- sample(x1, n1, replace = TRUE)
        s2 <- sample(x2, n2, replace = TRUE)
        fb <- d_reg(s1, s2, degree = max(object$degrees),
                    monotonicity = object$monotonicity, paired = FALSE, type = object$type)
      }
      boot_d[b] <- fb$estimate
    }
    alpha <- (1 - level) / 2
    ci <- stats::quantile(boot_d, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  } else {
    # Noncentral t inversion
    df <- if (object$paired) n1 - 1L else (if (object$type == "glass") n1 - 1L else n1 + n2 - 2L)
    n_tilde <- if (object$paired) n1 else (n1 * n2) / (n1 + n2)
    ci <- ci_nct(object$estimate, n1 = n1, n2 = n2, conf = level, df = df)
  }

  matrix(ci, nrow = 1L, dimnames = list("d_reg", c(paste0(100 * (1 - level)/2, " %"),
                                                   paste0(100 * (1 + level)/2, " %"))))
}

#' @export
#' @export
print.d_reg <- function(x, digits = 3L, ...) {
  cat("\n  Distribution-Free Effect Size Estimation (d_reg)\n")
  cat(strrep("-", 52), "\n", sep = "")

  if (x$paired) {
    cat(sprintf("  Effect Size (d_reg, raw scale)    :  %.*f\n", digits, x$d_reg))
    cat(sprintf("  Standardized Mean Change (d_z)    :  %.*f\n", digits, x$d_z))
    cat(sprintf("  Paired Hermite Correlation (r_HM) :  %.*f\n", digits, x$r_hm_paired$r_hm))
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n", digits, x$sd_standardizer))
    cat(sprintf("  Difference Model SD (sigma_diff)  :  %.*f\n", digits, x$diff_moments$sd))
  } else {
    cat(sprintf("  Effect Size (d_reg)               :  %.*f\n", digits, x$estimate))
    cat(sprintf("  Averaged Model SD (sigma_avg)     :  %.*f\n", digits, x$sd_standardizer))
  }

  cat(sprintf("  Hedges' g (Benchmark)             :  %.*f\n", digits, x$hedges_g))
  cat(sprintf("  Sample Sizes                      :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Polynomial Degrees                :  g1 = %d, g2 = %d\n", x$degrees["g1"], x$degrees["g2"]))
  cat(sprintf("  Monotonicity Check                :  %s\n", x$monotonicity))

  if (!is.null(x$ci)) {
    cat(sprintf("  %s CI (%s): [%.*f, %.*f]\n",
                paste0(round(x$conf_level * 100), "%"), x$ci_method,
                digits, x$ci[1L], digits, x$ci[2L]))
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.d_reg <- function(object, digits = 3L, ...) {
  print(object, digits = digits, ...)
  cat("  Group Distributional Moments (Polynomial Modeled):\n")
  cat(sprintf("    Group 1: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$group1$mean, digits, object$group1$sd, digits, object$group1$variance))
  cat(sprintf("    Group 2: Mean = %.*f, SD = %.*f, Var = %.*f\n",
              digits, object$group2$mean, digits, object$group2$sd, digits, object$group2$variance))
  cat("\n")
  invisible(object)
}

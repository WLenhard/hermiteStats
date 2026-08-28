#' @title Distribution-Robust Permutation Test for Mean Differences (Hermite t-Test)
#'
#' @description
#' Computes a permutation test for the difference in means between two
#' independent groups or paired (repeated) measurements, using the distribution-robust
#' polynomial-quantile-regularized effect size (\eqn{d_{\mathrm{reg}}} or \eqn{d_z})
#' from \code{\link{d_reg}} as the test statistic.
#'
#' @param x A numeric vector of observations for the first group (or baseline in a paired design),
#'   or a two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"t_hermite"}.
#'
#' @references
#' Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for multi-factorial analysis of variance. \emph{Journal of Statistical Computation and Simulation}, 73(2), 85–113. \doi{10.1080/00949650215733}
#'
#' Edgington, E. S., & Onghena, P. (2007). \emph{Randomization Tests} (4th ed.). Chapman & Hall/CRC.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. \emph{Statistical Applications in Genetics and Molecular Biology}, 9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' @seealso \code{\link{d_reg}}, \code{\link{shape_hermite}}, \code{\link{hermite_fit}}
#'
#' @examples
#' # 1. Independent groups (vectors)
#' set.seed(42)
#' ctrl <- rlnorm(25, meanlog = 2.0, sdlog = 0.5)
#' trt  <- rlnorm(25, meanlog = 2.3, sdlog = 0.5)
#' test1 <- t_hermite(ctrl, trt, nperm = 500)
#' print(test1)
#' plot(test1)
#'
#' # 2. Formula interface
#' df_test <- data.frame(
#'   score = c(ctrl, trt),
#'   group = factor(rep(c("Ctrl", "Trt"), each = 25))
#' )
#' test_form <- t_hermite(score ~ group, data = df_test, nperm = 500)
#' print(test_form)
#'
#' # 3. Paired / repeated-measures design
#' pre  <- rlnorm(20, meanlog = 3.0, sdlog = 0.4)
#' post <- pre + rnorm(20, mean = 3.0, sd = 1.5)
#' test_paired <- t_hermite(pre, post, paired = TRUE, nperm = 500)
#' print(test_paired)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
t_hermite <- function(x, ...) {
  UseMethod("t_hermite")
}

#' @rdname t_hermite
#' @param formula A two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param data An optional data frame, list, or environment containing variables in \code{formula}.
#' @export
t_hermite.formula <- function(formula, data = NULL, paired = FALSE, degree = 3L,
                              copula = c("none", "gaussian"),
                              monotonicity = c("relaxed", "strict", "none"),
                              type = c("regularized", "hedges", "glass", "combined"),
                              alternative = c("two.sided", "greater", "less"),
                              nperm = 2000L, min_success = 0.90, ...) {

  if (missing(formula) || length(formula) != 3L) {
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

  res <- t_hermite.default(x = x1, y = x2, paired = paired, degree = degree,
                           copula = copula, monotonicity = monotonicity,
                           type = type, alternative = alternative,
                           nperm = nperm, min_success = min_success, ...)
  res$group_labels <- levels_g
  res
}

#' @rdname t_hermite
#' @param y A numeric vector of observations for the second group (or post-measurement in paired designs).
#' @param paired Logical; \code{TRUE} for a paired-samples design. Default is \code{FALSE}.
#' @param degree Integer scalar; maximum polynomial degree used for quantile fits (default \code{3L}).
#' @param copula Character string; copula mode for paired designs (\code{"none"} or \code{"gaussian"}).
#' @param monotonicity Monotonicity constraint: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param type Effect-size variant: \code{"regularized"} (default), \code{"hedges"}, \code{"glass"}, or \code{"combined"}.
#' @param alternative Character string specifying the alternative hypothesis: \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param nperm Integer; number of permutation resamples (default \code{2000L}).
#' @param min_success Minimum proportion of successful permutation replicates required (default \code{0.90}).
#' @export
t_hermite.default <- function(x, y, paired = FALSE, degree = 3L,
                              copula = c("none", "gaussian"),
                              monotonicity = c("relaxed", "strict", "none"),
                              type = c("regularized", "hedges", "glass", "combined"),
                              alternative = c("two.sided", "greater", "less"),
                              nperm = 2000L, min_success = 0.90, ...) {

  copula       <- match.arg(copula)
  monotonicity <- match.arg(monotonicity)
  type         <- match.arg(type)
  alternative  <- match.arg(alternative)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (paired && length(x) != length(y)) {
    stop("'x' and 'y' must have identical length for a paired test.")
  }
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L) {
    stop("'nperm' must be a positive integer.")
  }

  obs_fit <- d_reg(x, y, degree = degree, copula = copula, monotonicity = monotonicity,
                   paired = paired, type = type, ...)

  n1 <- obs_fit$n1
  n2 <- obs_fit$n2
  n_tilde <- if (paired) n1 else (n1 * n2) / (n1 + n2)
  t_obs <- obs_fit$estimate * sqrt(n_tilde)

  x1 <- obs_fit$x1
  x2 <- obs_fit$x2

  perm_t <- rep(NA_real_, nperm)

  if (paired) {
    diffs <- x2 - x1
    for (b in seq_len(nperm)) {
      flips   <- sample(c(-1, 1), n1, replace = TRUE)
      x2_perm <- x1 + flips * diffs
      perm_t[b] <- tryCatch(
        d_reg(x1, x2_perm, degree = degree, copula = copula, monotonicity = monotonicity,
              paired = TRUE, type = type)$estimate * sqrt(n_tilde),
        error = function(e) NA_real_
      )
    }
  } else {
    combined <- c(x1, x2)
    N <- n1 + n2
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)
      perm_t[b] <- tryCatch(
        d_reg(combined[idx], combined[-idx], degree = degree, copula = copula,
              monotonicity = monotonicity, paired = FALSE, type = type)$estimate * sqrt(n_tilde),
        error = function(e) NA_real_
      )
    }
  }

  success      <- is.finite(perm_t)
  n_success    <- sum(success)
  success_rate <- n_success / nperm

  if (n_success < 1L) {
    warning("All permutation replicates failed; cannot compute a p-value.", call. = FALSE)
    p_val <- NA_real_
  } else {
    if (success_rate < min_success) {
      warning(sprintf("Only %d of %d permutation replicates succeeded (%.1f%%).",
                      n_success, nperm, 100 * success_rate), call. = FALSE)
    }
    perm_valid <- perm_t[success]
    tol <- 1e-10
    p_val <- switch(
      alternative,
      "two.sided" = (1 + sum(abs(perm_valid) >= (abs(t_obs) - tol))) / (n_success + 1),
      "greater"   = (1 + sum(perm_valid >= (t_obs - tol)))           / (n_success + 1),
      "less"      = (1 + sum(perm_valid <= (t_obs + tol)))           / (n_success + 1)
    )
  }

  res <- list(
    statistic         = t_obs,
    p_value           = p_val,
    estimate          = obs_fit$estimate,
    d_reg             = obs_fit$d_reg,
    d_z               = obs_fit$d_z,
    mean_diff         = obs_fit$group2$mean - obs_fit$group1$mean,
    n1                = n1,
    n2                = n2,
    paired            = paired,
    alternative       = alternative,
    nperm             = nperm,
    n_perm_success    = n_success,
    success_rate      = success_rate,
    perm_distribution = if (n_success > 0L) perm_t[success] else numeric(0),
    d_reg_fit         = obs_fit,
    method            = if (paired) "Paired Hermite Permutation Test" else "Two-Sample Hermite Permutation Test"
  )
  class(res) <- "t_hermite"
  res
}

#' @export
print.t_hermite <- function(x, digits = 3L, ...) {
  eff_lab <- if (x$paired) "d_z (standardized change)" else "d_reg"
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 58), "\n", sep = "")
  cat(sprintf("  t_Hermite Statistic        :  %.*f\n", digits, x$statistic))
  cat(sprintf("  Effect Size (%s) :  %.*f\n", eff_lab, digits, x$estimate))
  cat(sprintf("  Regularized Mean Difference:  %.*f\n", digits, x$mean_diff))
  cat(sprintf("  Sample Sizes               :  n1 = %d, n2 = %d\n", x$n1, x$n2))
  cat(sprintf("  Alternative Hypothesis     :  %s\n", x$alternative))
  cat(sprintf("  Permutations (successful)  :  %d / %d (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))

  p_txt <- if (is.na(x$p_value)) {
    "NA"
  } else if (x$p_value < 10^(-digits)) {
    sprintf("< %.*f", digits, 10^(-digits))
  } else {
    sprintf("%.*f", digits, x$p_value)
  }
  cat(sprintf("  p-value (Monte Carlo perm.):  %s\n", p_txt))

  if (!is.null(x$d_reg_fit$ci)) {
    ci_level <- round(x$d_reg_fit$conf_level * 100)
    cat(sprintf("  %d%% CI for Effect Size     : [%.*f, %.*f]\n",
                ci_level, digits, x$d_reg_fit$ci[1L], digits, x$d_reg_fit$ci[2L]))
  }
  cat("\n")
  invisible(x)
}


# =============================================================================
# shape_hermite
# =============================================================================

#' @title Distribution-Robust Permutation Test for Shape Differences (Variance, Skewness, Kurtosis)
#'
#' @description
#' Tests for differences in the regularized distributional shape
#' (variance, skewness, excess kurtosis) between two independent groups or paired measurements,
#' using closed-form Hermite moments and an exact permutation null distribution.
#'
#' @param x A numeric vector of observations, or a two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"shape_hermite"}.
#'
#' @references
#' O'Brien, P. C. (1984). Procedures for comparing samples with multiple endpoints. \emph{Biometrics}, 40(4), 1079–1087. \doi{10.2307/2531158}
#'
#' @seealso \code{\link{t_hermite}}, \code{\link{hermite_moments}}
#'
#' @examples
#' # 1. Independent groups
#' set.seed(1)
#' g1 <- rnorm(40)
#' g2 <- rgamma(40, shape = 2)
#' sres <- shape_hermite(g1, g2, nperm = 500)
#' print(sres)
#' plot(sres)
#'
#' # 2. Formula interface
#' df_shape <- data.frame(val = c(g1, g2), grp = factor(rep(c("A", "B"), each = 40)))
#' sres_form <- shape_hermite(val ~ grp, data = df_shape, nperm = 500)
#' print(sres_form)
#'
#' @author Wolfgang Lenhard, Alexandra Lenhard
#' @export
shape_hermite <- function(x, ...) {
  UseMethod("shape_hermite")
}

#' @rdname shape_hermite
#' @param formula A two-sided \code{\link[stats]{formula}} of the form \code{response ~ group}.
#' @param data An optional data frame containing variables in \code{formula}.
#' @export
shape_hermite.formula <- function(formula, data = NULL, paired = FALSE, degree = 3L,
                                  monotonicity = c("relaxed", "strict", "none"),
                                  ties_method = c("average", "random"),
                                  nperm = 2000L, min_success = 0.90, omnibus = TRUE, ...) {

  if (missing(formula) || length(formula) != 3L) {
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

  res <- shape_hermite.default(x = x1, y = x2, paired = paired, degree = degree,
                               monotonicity = monotonicity, ties_method = ties_method,
                               nperm = nperm, min_success = min_success, omnibus = omnibus, ...)
  res$group_labels <- levels_g
  res
}

#' @rdname shape_hermite
#' @param y A numeric vector of observations for the second group.
#' @param paired Logical; \code{TRUE} for a paired-samples design. Default is \code{FALSE}.
#' @param degree Integer scalar; maximum polynomial degree requested (default \code{3L}).
#' @param monotonicity Monotonicity constraint: \code{"relaxed"} (default), \code{"strict"}, or \code{"none"}.
#' @param ties_method Rank tie-handling method: \code{"average"} (default) or \code{"random"}.
#' @param nperm Integer; number of permutation resamples (default \code{2000L}).
#' @param min_success Minimum proportion of successful permutation replicates required (default \code{0.90}).
#' @param omnibus Logical; if \code{TRUE} (default), computes an omnibus sum-of-squares shape test.
#' @export
shape_hermite.default <- function(x, y, paired = FALSE, degree = 3L,
                                  monotonicity = c("relaxed", "strict", "none"),
                                  ties_method = c("average", "random"),
                                  nperm = 2000L, min_success = 0.90, omnibus = TRUE, ...) {

  monotonicity <- match.arg(monotonicity)
  ties_method  <- match.arg(ties_method)

  if (missing(y) || is.null(y)) stop("Argument 'y' must be supplied.")
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' and 'y' must be numeric vectors.")
  if (paired && length(x) != length(y)) {
    stop("'x' and 'y' must have identical length for a paired design.")
  }

  if (paired) {
    ok <- is.finite(x) & is.finite(y)
    x1 <- x[ok]; x2 <- y[ok]
  } else {
    x1 <- x[is.finite(x)]; x2 <- y[is.finite(y)]
  }
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 4L || n2 < 4L) stop("Each group must contain at least 4 valid observations.")

  moment_names <- c("variance", "skewness", "excess_kurtosis")

  fit_matched <- function(a, b) {
    cap <- degree
    repeat {
      fa <- hermite_fit(a, degree = cap, monotonicity = monotonicity,
                        ties_method = ties_method, force_odd = TRUE)
      fb <- hermite_fit(b, degree = cap, monotonicity = monotonicity,
                        ties_method = ties_method, force_odd = TRUE)
      common <- min(fa$degree, fb$degree)
      if (common >= cap) break
      cap <- common
    }
    list(fit1 = fa, fit2 = fb, degree = cap)
  }

  extract_diff <- function(fm) {
    m1 <- hermite_moments(fm$fit1)
    m2 <- hermite_moments(fm$fit2)
    c(variance        = m2$variance - m1$variance,
      skewness        = m2$skewness - m1$skewness,
      excess_kurtosis = m2$excess_kurtosis - m1$excess_kurtosis)
  }

  obs_fm   <- fit_matched(x1, x2)
  obs_diff <- extract_diff(obs_fm)

  perm_diff <- matrix(NA_real_, nrow = nperm, ncol = 3L,
                      dimnames = list(NULL, moment_names))

  if (paired) {
    for (b in seq_len(nperm)) {
      swap <- sample(c(FALSE, TRUE), n1, replace = TRUE)
      a  <- ifelse(swap, x2, x1)
      bb <- ifelse(swap, x1, x2)
      perm_diff[b, ] <- tryCatch(extract_diff(fit_matched(a, bb)),
                                 error = function(e) rep(NA_real_, 3L))
    }
  } else {
    combined <- c(x1, x2); N <- n1 + n2
    for (b in seq_len(nperm)) {
      idx <- sample.int(N, n1, replace = FALSE)
      perm_diff[b, ] <- tryCatch(extract_diff(fit_matched(combined[idx], combined[-idx])),
                                 error = function(e) rep(NA_real_, 3L))
    }
  }

  success      <- stats::complete.cases(perm_diff)
  n_success    <- sum(success)
  success_rate <- n_success / nperm

  if (n_success < 1L) stop("All permutation replicates failed; cannot compute p-values.")
  if (success_rate < min_success) {
    warning(sprintf("Only %d of %d permutation replicates succeeded (%.1f%%).",
                    n_success, nperm, 100 * success_rate), call. = FALSE)
  }
  perm_valid <- perm_diff[success, , drop = FALSE]

  tol <- 1e-10
  p_raw <- vapply(seq_along(moment_names), function(k) {
    (1 + sum(abs(perm_valid[, k]) >= (abs(obs_diff[k]) - tol))) / (n_success + 1)
  }, numeric(1L))
  names(p_raw) <- moment_names
  p_adj <- stats::p.adjust(p_raw, method = "holm")

  omnibus_res <- NULL
  if (omnibus) {
    ctr <- colMeans(perm_valid)
    scl <- apply(perm_valid, 2, stats::sd)
    ok_m <- is.finite(scl) & scl > 0
    if (!all(ok_m)) {
      warning("Degenerate permutation variance for: ",
              paste(moment_names[!ok_m], collapse = ", "),
              " -- excluded from the omnibus statistic.", call. = FALSE)
    }
    z_perm <- sweep(sweep(perm_valid[, ok_m, drop = FALSE], 2, ctr[ok_m], "-"), 2, scl[ok_m], "/")
    z_obs  <- (obs_diff[ok_m] - ctr[ok_m]) / scl[ok_m]
    stat_perm <- rowSums(z_perm^2)
    stat_obs  <- sum(z_obs^2)
    p_omni <- (1 + sum(stat_perm >= (stat_obs - tol))) / (length(stat_perm) + 1)
    omnibus_res <- list(statistic = stat_obs, p_value = p_omni)
  }

  res <- list(
    estimate          = obs_diff,
    p_value           = p_raw,
    p_adjusted        = p_adj,
    omnibus           = omnibus_res,
    n1                = n1,
    n2                = n2,
    degree_matched    = obs_fm$degree,
    degree_requested  = degree,
    paired            = paired,
    nperm             = nperm,
    n_perm_success    = n_success,
    success_rate      = success_rate,
    perm_distribution = perm_valid,
    fit1              = obs_fm$fit1,
    fit2              = obs_fm$fit2,
    method            = if (paired) "Paired Hermite Shape-Difference Permutation Test"
    else "Two-Sample Hermite Shape-Difference Permutation Test"
  )
  class(res) <- "shape_hermite"
  res
}

#' @export
print.shape_hermite <- function(x, digits = 3L, ...) {
  cat(sprintf("\n  %s\n", x$method))
  cat(strrep("-", 64), "\n", sep = "")
  cat(sprintf("  Sample Sizes              :  n1 = %d, n2 = %d\n", x$n1, x$n2))

  deg_txt <- if (x$degree_matched == x$degree_requested) {
    sprintf("%d", x$degree_matched)
  } else {
    sprintf("%d (reduced from requested: %d)", x$degree_matched, x$degree_requested)
  }
  cat(sprintf("  Matched Polynomial Degree :  %s\n", deg_txt))
  cat(sprintf("  Permutations (successful) :  %d / %d (%.1f%%)\n",
              x$n_perm_success, x$nperm, 100 * x$success_rate))

  m1 <- if (!is.null(x$fit1)) hermite_moments(x$fit1) else NULL
  m2 <- if (!is.null(x$fit2)) hermite_moments(x$fit2) else NULL

  cat("\n  Shape-Difference Profile (Group 2 - Group 1):\n")

  moment_labels <- c(
    variance        = "Variance",
    skewness        = "Skewness (g1)",
    excess_kurtosis = "Excess Kurtosis (g2)"
  )

  p_fmt <- function(p) {
    vapply(p, function(val) {
      if (is.na(val)) "NA"
      else if (val < 10^(-digits)) sprintf("< %.*f", digits, 10^(-digits))
      else sprintf("%.*f", digits, val)
    }, character(1L))
  }

  tab <- data.frame(
    Moment     = unname(moment_labels[names(x$estimate)]),
    Group_1    = if (!is.null(m1)) sprintf("%.*f", digits, c(m1$variance, m1$skewness, m1$excess_kurtosis)) else "-",
    Group_2    = if (!is.null(m2)) sprintf("%.*f", digits, c(m2$variance, m2$skewness, m2$excess_kurtosis)) else "-",
    Difference = sprintf("%+.*f", digits, unname(x$estimate)),
    p_value    = p_fmt(unname(x$p_value)),
    p_Holm     = p_fmt(unname(x$p_adjusted)),
    stringsAsFactors = FALSE
  )
  names(tab) <- c("Moment", "Group 1", "Group 2", "Difference", "p-value", "p (Holm)")
  print(tab, row.names = FALSE)

  if (!is.null(x$omnibus)) {
    omni_p <- if (x$omnibus$p_value < 10^(-digits)) {
      sprintf("< %.*f", digits, 10^(-digits))
    } else {
      sprintf("%.*f", digits, x$omnibus$p_value)
    }
    cat("\n  Omnibus Shape Test (Standardized Sum of Squares):\n")
    cat(sprintf("    Statistic = %.*f  |  p-value = %s\n",
                digits, x$omnibus$statistic, omni_p))
  }
  cat("\n")
  invisible(x)
}


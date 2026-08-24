#' Exact Small-Sample Correction Factor J(df)
#' @export
.J_correction <- function(df) {
  df <- max(df, 1.0001)
  exp(lgamma(df / 2.0) - 0.5 * log(df / 2.0) - lgamma((df - 1.0) / 2.0))
}

#' Standardized Mean Differences (Cohen's d and Hedges' g)
#'
#' @param x1 Numeric vector of group 1 scores.
#' @param x2 Numeric vector of group 2 scores.
#' @param type Standardizer: \code{"pooled"} (default), \code{"avg"}, or \code{"glass"}.
#' @param correct_bias Logical; apply Hedges' small-sample correction \eqn{J(df)}? Default is \code{TRUE}.
#'
#' @return Numeric scalar.
#' @export
d_cohen <- function(x1, x2, type = c("pooled", "avg", "glass"), correct_bias = TRUE) {
  type <- match.arg(type)
  n1 <- length(x1); n2 <- length(x2)
  if (n1 < 2L || n2 < 2L) return(NA_real_)
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  md <- mean(x2) - mean(x1)

  if (type == "pooled") {
    sdv <- sqrt(((n1 - 1L) * s1^2 + (n2 - 1L) * s2^2) / (n1 + n2 - 2L))
    df  <- n1 + n2 - 2L
  } else if (type == "avg") {
    sdv <- sqrt((s1^2 + s2^2) / 2.0)
    df  <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1 - 1L) + (s2^2/n2)^2/(n2 - 1L))
  } else {
    sdv <- s1
    df  <- n1 - 1L
  }
  if (sdv <= 0) return(0.0)
  d <- md / sdv
  if (correct_bias) d <- d * .J_correction(df)
  d
}

#' @rdname d_cohen
#' @export
hedges_g <- function(x1, x2) {
  d_cohen(x1, x2, type = "pooled", correct_bias = TRUE)
}

#' Noncentral-t (Lambda-Prime) Confidence Interval for SMD
#' @export
ci_nct <- function(d_point, n1, n2, conf = 0.95, df = n1 + n2 - 2L) {
  n_tilde <- (n1 * n2) / (n1 + n2)
  t_obs   <- d_point * sqrt(n_tilde)
  a       <- (1 - conf) / 2

  ncp_from_t <- function(t_val, df_val, p_val) {
    f <- function(ncp) suppressWarnings(stats::pt(t_val, df = df_val, ncp = ncp)) - p_val
    stats::uniroot(f, c(t_val - 6, t_val + 6), extendInt = "yes")$root
  }

  lo <- ncp_from_t(t_obs, df, 1 - a) / sqrt(n_tilde)
  hi <- ncp_from_t(t_obs, df, a) / sqrt(n_tilde)
  c(lower = lo, upper = hi)
}

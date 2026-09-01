# ==============================================================================
# sim_test_helper.R -- Tests on mean and median through Hermite polynomials
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(hermiteStats)
})
if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
library(pbapply)

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

# ------------------------------------------------------------------------------
# 1. Fast Yuen 20% Trimmed t-Test
# ------------------------------------------------------------------------------
yuen_test <- function(x, y, tr = 0.2) {
  trim_group <- function(v) {
    v <- sort(v[is.finite(v)])
    n <- length(v)
    k <- floor(tr * n)
    vals_trimmed <- v[(k + 1L):(n - k)]
    h <- length(vals_trimmed)
    m_trim <- mean(vals_trimmed)

    v_win <- v
    v_win[1:k] <- v[k + 1L]
    v_win[(n - k + 1L):n] <- v[n - k]
    var_win <- sum((v_win - mean(v_win))^2) / (n - 1L)

    list(m = m_trim, h = h, n = n, var_win = var_win, se2 = ((n - 1L) * var_win) / (h * (h - 1L)))
  }

  g1 <- trim_group(x)
  g2 <- trim_group(y)

  se_diff <- sqrt(g1$se2 + g2$se2)
  if (se_diff <= 1e-10) return(list(statistic = 0, p_value = 1.0))

  t_val <- (g2$m - g1$m) / se_diff
  df <- (g1$se2 + g2$se2)^2 / ((g1$se2^2 / (g1$h - 1L)) + (g2$se2^2 / (g2$h - 1L)))

  p_val <- 2 * stats::pt(-abs(t_val), df = df)
  list(statistic = t_val, p_value = p_val)
}

# ------------------------------------------------------------------------------
# 2. Synthetic Data Generator (Preserving Theoretical Population Nulls)
# ------------------------------------------------------------------------------
generate_location_data <- function(n1, n2, dist, param1, param2, param3 = NA,
                                   scale_ratio = 1.0, delta = 0.0) {

  # Generates zero-population-mean data with scaling
  gen_raw <- function(n, s_factor = 1.0) {
    if (dist == "norm") {
      rnorm(n, mean = 0, sd = param2 * s_factor)

    } else if (dist == "lnorm") {
      # Theoretical mean of Lognormal(0, s) is exp(s^2 / 2)
      pop_mean <- exp(param2^2 / 2)
      (rlnorm(n, meanlog = param1, sdlog = param2) - pop_mean) * s_factor

    } else if (dist == "t") {
      # Student-t with df > 1 has theoretical mean = 0
      rt(n, df = param1) * param2 * s_factor

    } else if (dist == "exp") {
      # Exponential with rate lambda has theoretical mean = 1/lambda
      pop_mean <- 1 / param1
      (rexp(n, rate = param1) - pop_mean) * s_factor

    } else if (dist == "outlier") {
      contam_prop <- if (is.finite(param3)) param3 else 0.02
      pop_mean <- contam_prop * param1
      contam <- runif(n) < contam_prop
      vals <- rnorm(n)
      vals[contam] <- vals[contam] + param1
      (vals - pop_mean) * s_factor

    } else if (dist == "1plirt") {
      n_items <- if (is.na(param3)) 20L else as.integer(param3)
      items <- seq(-2, 2, length.out = n_items)
      theta <- rnorm(n, sd = s_factor)
      pm <- stats::plogis(outer(theta, items, "-"))
      # Theoretical mean for symmetric items is n_items / 2
      rowSums(matrix(runif(n * n_items), n) < pm) - (n_items / 2.0)

    } else {
      stop("Unknown distribution family: ", dist)
    }
  }

  x1 <- gen_raw(n1, s_factor = 1.0)
  x2 <- gen_raw(n2, s_factor = scale_ratio)

  # Add standardized location shift to Group 2 if delta > 0
  if (delta != 0.0) {
    s1_pop <- if (dist == "lnorm") sqrt((exp(param2^2) - 1) * exp(param2^2))
    else if (dist == "t") sqrt(param1 / (param1 - 2)) * param2
    else if (dist == "exp") 1 / param1
    else if (dist == "norm") param2
    else if (dist == "outlier") sqrt(1 + 0.02 * param1^2 - (0.02 * param1)^2)
    else 2.5 # IRT approx SD
    s2_pop <- s1_pop * scale_ratio
    s_pooled <- sqrt((s1_pop^2 + s2_pop^2) / 2.0)
    x2 <- x2 + (delta * s_pooled)
  }

  list(x1 = x1, x2 = x2)
}

# ------------------------------------------------------------------------------
# 3. Ground Truth Pre-computation (Population Means & Medians)
# ------------------------------------------------------------------------------
.cond_key <- function(cond) {
  p3 <- if ("param3" %in% names(cond)) cond$param3 else NA
  sr <- if ("scale_ratio" %in% names(cond)) cond$scale_ratio else 1.0
  paste(cond$dist, cond$param1, cond$param2, p3, sr, cond$delta, sep = "_")
}

compute_location_truths <- function(conditions, pop_size = 5e5, seed = 1234) {
  keys <- vapply(seq_len(nrow(conditions)), function(i) .cond_key(conditions[i, ]), character(1L))
  uniq <- which(!duplicated(keys))

  cat(sprintf("\nPre-computing ground truths (%d unique distribution/effect cells, N = %g)...\n",
              length(uniq), pop_size))

  truth_cache <- list()
  pb <- txtProgressBar(min = 0, max = length(uniq), style = 3)

  for (j in seq_along(uniq)) {
    i <- uniq[j]
    cond <- conditions[i, ]
    p3 <- if ("param3" %in% names(cond)) cond$param3 else NA
    sr <- if ("scale_ratio" %in% names(cond)) cond$scale_ratio else 1.0

    set.seed(seed + j)
    pop <- generate_location_data(pop_size, pop_size, cond$dist, cond$param1, cond$param2, p3,
                                  scale_ratio = sr, delta = cond$delta)

    m1 <- mean(pop$x1); m2 <- mean(pop$x2)
    med1 <- stats::median(pop$x1); med2 <- stats::median(pop$x2)
    s1 <- stats::sd(pop$x1); s2 <- stats::sd(pop$x2)
    s_avg <- sqrt((s1^2 + s2^2) / 2.0)

    truth_cache[[keys[i]]] <- list(
      true_mean_diff   = m2 - m1,
      true_median_diff = med2 - med1,
      true_s_avg       = s_avg,
      true_d_mean      = (m2 - m1) / s_avg,
      true_d_median    = (med2 - med1) / s_avg
    )
    setTxtProgressBar(pb, j)
  }
  close(pb)
  truth_cache
}

# ------------------------------------------------------------------------------
# 4. Simulation Worker for a Single Condition Cell
# ------------------------------------------------------------------------------
simulate_location_cell <- function(cond_idx, conditions, truth_cache,
                                   n_reps, nperm, alpha = 0.05,
                                   seed_base = 20260101L) {
  cond  <- conditions[cond_idx, ]
  truth <- truth_cache[[.cond_key(cond)]]
  p3    <- if ("param3" %in% names(cond)) cond$param3 else NA
  sr    <- if ("scale_ratio" %in% names(cond)) cond$scale_ratio else 1.0

  n1    <- cond$n1
  n2    <- cond$n2
  delta <- cond$delta

  set.seed(as.integer(seed_base + cond_idx))

  p_herm_mean <- numeric(n_reps)
  p_herm_med  <- numeric(n_reps)
  p_welch     <- numeric(n_reps)
  p_student   <- numeric(n_reps)
  p_wilcox    <- numeric(n_reps)
  p_yuen      <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    dat <- generate_location_data(n1, n2, cond$dist, cond$param1, cond$param2, p3,
                                  scale_ratio = sr, delta = delta)
    x1 <- dat$x1
    x2 <- dat$x2

    # 1. Hermite Mean Test (t_hermite, permutation)
    res_tm <- tryCatch(
      t_hermite(x1, x2, method = "permutation", degree = 3L, nperm = nperm),
      error = function(e) list(p_value = NA_real_)
    )
    p_herm_mean[r] <- res_tm$p_value

    # 2. Hermite Median Test (median_hermite, permutation)
    res_med <- tryCatch(
      median_hermite(x1, x2, method = "permutation", degree = 3L, nperm = nperm),
      error = function(e) list(p_value = NA_real_)
    )
    p_herm_med[r] <- res_med$p_value

    # 3. Welch's t-test
    res_tw <- tryCatch(stats::t.test(x2, x1, var.equal = FALSE), error = function(e) list(p.value = NA_real_))
    p_welch[r] <- res_tw$p.value

    # 4. Student's t-test (Pooled variance)
    res_ts <- tryCatch(stats::t.test(x2, x1, var.equal = TRUE), error = function(e) list(p.value = NA_real_))
    p_student[r] <- res_ts$p.value

    # 5. Wilcoxon / Mann-Whitney U test
    res_wx <- tryCatch(stats::wilcox.test(x2, x1, exact = FALSE), error = function(e) list(p.value = NA_real_))
    p_wilcox[r] <- res_wx$p.value

    # 6. Yuen's 20% trimmed t-test
    res_yu <- tryCatch(yuen_test(x1, x2, tr = 0.2), error = function(e) list(p_value = NA_real_))
    p_yuen[r] <- res_yu$p_value
  }

  rej_rate <- function(p) safe_mean(p < alpha)

  data.frame(
    condition_id     = cond_idx,
    n1               = n1,
    n2               = n2,
    dist             = cond$dist,
    scale_ratio      = sr,
    delta            = delta,
    true_d_mean      = truth$true_d_mean,
    true_d_median    = truth$true_d_median,
    is_balanced      = (n1 == n2),
    is_homoscedastic = (sr == 1.0),

    # Rejection Rates (Alpha if delta = 0, Power if delta > 0)
    rate_t_hermite   = rej_rate(p_herm_mean),
    rate_med_hermite = rej_rate(p_herm_med),
    rate_welch       = rej_rate(p_welch),
    rate_student     = rej_rate(p_student),
    rate_wilcox      = rej_rate(p_wilcox),
    rate_yuen        = rej_rate(p_yuen)
  )
}

# ------------------------------------------------------------------------------
# 5. Simulation Runner
# ------------------------------------------------------------------------------
run_location_sim <- function(conditions, n_reps = 500, nperm = 1000L,
                             n_cores = 4, alpha = 0.05, seed = 20260101L,
                             truth_pop_size = 5e5) {
  start_sim <- Sys.time()

  truth_cache <- compute_location_truths(
    conditions = conditions,
    pop_size   = truth_pop_size,
    seed       = seed
  )

  n_conds <- nrow(conditions)
  cat(sprintf("Running %d conditions x %d replications (nperm = %d)...\n\n",
              n_conds, n_reps, nperm))

  pkg_root <- if (file.exists("DESCRIPTION")) {
    normalizePath(".", mustWork = TRUE)
  } else if (file.exists("../DESCRIPTION")) {
    normalizePath("..", mustWork = TRUE)
  } else {
    getwd()
  }

  cl <- NULL
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)

    parallel::clusterCall(cl, function(path) {
      devtools::load_all(path, quiet = TRUE)
      NULL
    }, path = pkg_root)

    parallel::clusterExport(
      cl,
      c("generate_location_data", ".cond_key", "yuen_test", "safe_mean",
        "simulate_location_cell"),
      envir = environment()
    )
  }

  res_list <- pblapply(
    seq_len(n_conds),
    simulate_location_cell,
    conditions  = conditions,
    truth_cache = truth_cache,
    n_reps      = n_reps,
    nperm       = nperm,
    alpha       = alpha,
    seed_base   = seed,
    cl          = cl
  )

  agg <- dplyr::bind_rows(res_list)
  cat("\n== Location simulation completed in:", format(Sys.time() - start_sim), "==\n\n")

  list(aggregated = agg, truth_cache = truth_cache)
}

# ------------------------------------------------------------------------------
# 6. Structured Summary & Win-Rate Diagnostics (Updated)
# ------------------------------------------------------------------------------
summary_location_sim <- function(agg, alpha = 0.05) {
  cat("\n========================================================================================\n")
  cat("              HERMITE LOCATION TESTING BENCHMARK: MEAN & MEDIAN DOMAINS\n")
  cat("========================================================================================\n\n")

  # --- 1. Homoscedastic Type I Error Control (delta = 0, scale_ratio = 1) -----
  cat("--- 1. TYPE I ERROR CONTROL (Standard Homoscedastic Conditions, delta = 0) ------------\n")
  t1_std <- agg %>%
    dplyr::filter(delta == 0.0, is_homoscedastic, is_balanced) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      t_Hermite   = safe_mean(rate_t_hermite),
      Med_Hermite = safe_mean(rate_med_hermite),
      Welch       = safe_mean(rate_welch),
      Student     = safe_mean(rate_student),
      Wilcoxon    = safe_mean(rate_wilcox),
      Yuen        = safe_mean(rate_yuen),
      .groups     = "drop"
    )
  print.data.frame(t1_std, digits = 4, row.names = FALSE)

  # --- 2. Behrens-Fisher Stress Test (Unequal n + Heteroscedasticity) ----------
  cat("\n--- 2. BEHRENS-FISHER STRESS TEST: Type I Error (delta = 0, Unequal n, s2 != s1) -----\n")
  cat("  [Direct: n1 < n2, s1 < s2  |  Inverse: n1 < n2, s1 > s2 (Student t breaks)]\n")
  t1_bf <- agg %>%
    dplyr::filter(delta == 0.0, (!is_homoscedastic | !is_balanced)) %>%
    dplyr::group_by(dist, n1, n2, scale_ratio) %>%
    dplyr::summarise(
      t_Hermite   = safe_mean(rate_t_hermite),
      Med_Hermite = safe_mean(rate_med_hermite),
      Welch       = safe_mean(rate_welch),
      Student     = safe_mean(rate_student),
      Wilcoxon    = safe_mean(rate_wilcox),
      Yuen        = safe_mean(rate_yuen),
      .groups     = "drop"
    )
  print.data.frame(t1_bf, digits = 4, row.names = FALSE)

  # --- 3. Domain A: Full Location Power Comparison (t_Hermite vs. Med_Hermite vs. Welch vs. Student) ---
  cat("\n--- 3. DOMAIN A: LOCATION POWER (t_Hermite vs. Med_Hermite vs. Welch vs. Student) ----\n")
  pwr_mean <- agg %>%
    dplyr::filter(delta > 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Power_t_Hermite   = safe_mean(rate_t_hermite),
      Power_Med_Hermite = safe_mean(rate_med_hermite),
      Power_Welch       = safe_mean(rate_welch),
      Power_Student     = safe_mean(rate_student),
      .groups           = "drop"
    )
  print.data.frame(pwr_mean, digits = 4, row.names = FALSE)

  # --- 4. Domain B: Median & Robust Location Power (median_hermite vs. Wilcoxon vs. Yuen) -
  cat("\n--- 4. DOMAIN B: MEDIAN & ROBUST LOCATION POWER (median_hermite vs. Wilcoxon vs. Yuen) -\n")
  pwr_med <- agg %>%
    dplyr::filter(delta > 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Power_Med_Hermite = safe_mean(rate_med_hermite),
      Power_Wilcoxon    = safe_mean(rate_wilcox),
      Power_Yuen        = safe_mean(rate_yuen),
      .groups           = "drop"
    )
  print.data.frame(pwr_med, digits = 4, row.names = FALSE)

  # --- 5. Head-to-Head Win-Rate Summary ---------------------------------------
  pwr_cells <- agg %>% dplyr::filter(delta > 0.0)
  total_cells <- nrow(pwr_cells)

  cat("\n--- 5. HEAD-TO-HEAD WIN-RATE SUMMARY (Cell-Level Win %) ------------------------------\n")
  cat("  [Mean & Parametric Benchmarks]\n")
  w_welch   <- sum(pwr_cells$rate_t_hermite >= pwr_cells$rate_welch)
  w_student <- sum(pwr_cells$rate_t_hermite >= pwr_cells$rate_student)
  w_med_t   <- sum(pwr_cells$rate_med_hermite >= pwr_cells$rate_t_hermite)
  cat(sprintf("    t_Hermite      >= Welch's t-test      : %2d / %2d cells (%5.1f%%)\n",
              w_welch, total_cells, 100 * w_welch / total_cells))
  cat(sprintf("    t_Hermite      >= Student's t-test    : %2d / %2d cells (%5.1f%%)\n",
              w_student, total_cells, 100 * w_student / total_cells))
  cat(sprintf("    median_hermite >= t_Hermite (Mean)    : %2d / %2d cells (%5.1f%%)\n",
              w_med_t, total_cells, 100 * w_med_t / total_cells))

  cat("\n  [Robust & Nonparametric Benchmarks]\n")
  w_wilcox <- sum(pwr_cells$rate_med_hermite >= pwr_cells$rate_wilcox)
  w_yuen   <- sum(pwr_cells$rate_med_hermite >= pwr_cells$rate_yuen)
  cat(sprintf("    median_hermite >= Wilcoxon Rank-Sum   : %2d / %2d cells (%5.1f%%)\n",
              w_wilcox, total_cells, 100 * w_wilcox / total_cells))
  cat(sprintf("    median_hermite >= Yuen Trimmed-Mean   : %2d / %2d cells (%5.1f%%)\n",
              w_yuen, total_cells, 100 * w_yuen / total_cells))

  cat("========================================================================================\n\n")
}

# ------------------------------------------------------------------------------
# 7. Diagnostic Visualizations
# ------------------------------------------------------------------------------

#' Visualize Behrens-Fisher Robustness & Standard Type I Error
#' @export
plot_type1_error_comparison <- function(agg, nominal = 0.05) {
  t1 <- agg %>% dplyr::filter(delta == 0.0)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, 2), mar = c(5, 4.5, 3.5, 1.2), font.main = 1)

  # Panel A: Homoscedastic balanced Type I error
  std_t1 <- t1 %>%
    dplyr::filter(is_homoscedastic, is_balanced) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      t_Hermite   = safe_mean(rate_t_hermite),
      Med_Hermite = safe_mean(rate_med_hermite),
      Welch       = safe_mean(rate_welch),
      Student     = safe_mean(rate_student),
      Wilcoxon    = safe_mean(rate_wilcox),
      .groups     = "drop"
    )

  mat1 <- as.matrix(std_t1[, -1])
  rownames(mat1) <- std_t1$dist

  barplot(t(mat1), beside = TRUE,
          col = c("blue3", "cyan3", "orchid", "gray60", "darkorange"),
          ylim = c(0, 0.10), ylab = "Empirical Type I Error",
          main = "A. Standard Homoscedastic Nulls (delta = 0)",
          las = 2)
  abline(h = nominal, lty = 2, col = "red", lwd = 2)
  abline(h = c(0.025, 0.075), lty = 3, col = "gray40")

  # Panel B: Behrens-Fisher Stress Conditions
  bf_t1 <- t1 %>%
    dplyr::filter(!is_homoscedastic | !is_balanced) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      t_Hermite   = safe_mean(rate_t_hermite),
      Med_Hermite = safe_mean(rate_med_hermite),
      Welch       = safe_mean(rate_welch),
      Student     = safe_mean(rate_student),
      Wilcoxon    = safe_mean(rate_wilcox),
      .groups     = "drop"
    )

  mat2 <- as.matrix(bf_t1[, -1])
  rownames(mat2) <- bf_t1$dist

  barplot(t(mat2), beside = TRUE,
          col = c("blue3", "cyan3", "orchid", "gray60", "darkorange"),
          ylim = c(0, 0.25), ylab = "Empirical Type I Error",
          main = "B. Behrens-Fisher Stress Nulls (s1 != s2, n1 != n2)",
          las = 2)
  abline(h = nominal, lty = 2, col = "red", lwd = 2)
  abline(h = c(0.025, 0.075), lty = 3, col = "gray40")

  legend("topright",
         legend = c("t_Hermite (Mean)", "median_Hermite", "Welch t", "Student t", "Wilcoxon"),
         fill = c("blue3", "cyan3", "orchid", "gray60", "darkorange"),
         bty = "n", cex = 0.75)
}

#' Multi-Panel Power Comparison Across Domains
#' @export
plot_location_power_curves <- function(agg, dist_focus = "lnorm") {
  sub <- agg %>% dplyr::filter(dist == dist_focus, delta > 0.0, is_homoscedastic, is_balanced)
  deltas <- sort(unique(sub$delta))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, length(deltas)), mar = c(4.5, 4.2, 3, 1), font.main = 1)

  for (d in deltas) {
    dat_d <- sub %>% dplyr::filter(delta == d) %>% dplyr::arrange(n1)

    plot(dat_d$n1, dat_d$rate_med_hermite, type = "b", pch = 19, col = "cyan3", lwd = 2.5,
         ylim = c(0, 1), log = "x", xlab = "Sample Size (n per group)", ylab = "Statistical Power",
         main = sprintf("%s (delta = %.1f)", tools::toTitleCase(dist_focus), d), bty = "l",
         panel.first = grid(col = "gray90", lty = 1))

    lines(dat_d$n1, dat_d$rate_wilcox,    type = "b", pch = 15, col = "darkorange", lwd = 2, lty = 2)
    lines(dat_d$n1, dat_d$rate_yuen,      type = "b", pch = 4,  col = "darkcyan",   lwd = 2, lty = 3)
    lines(dat_d$n1, dat_d$rate_t_hermite, type = "b", pch = 17, col = "blue3",      lwd = 2, lty = 4)
    lines(dat_d$n1, dat_d$rate_welch,     type = "b", pch = 18, col = "orchid",     lwd = 1.8, lty = 5)

    if (d == deltas[1L]) {
      legend("bottomright",
             legend = c("median_Hermite", "Wilcoxon", "Yuen 20%", "t_Hermite (Mean)", "Welch t"),
             col = c("cyan3", "darkorange", "darkcyan", "blue3", "orchid"),
             pch = c(19, 15, 4, 17, 18), lty = c(1, 2, 3, 4, 5), lwd = 2,
             bty = "n", cex = 0.75)
    }
  }
}

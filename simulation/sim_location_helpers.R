# ==============================================================================
# sim_location_helpers.R -- Simulation Infrastructure for Mean & Median Tests
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
})
if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
library(pbapply)

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

# ------------------------------------------------------------------------------
# 1. Raw Sample Median Permutation Test
# ------------------------------------------------------------------------------
raw_median_perm_test <- function(x1, x2, nperm = 499L) {
  obs_diff <- stats::median(x2) - stats::median(x1)

  # Location-aligned pool under H0
  r1 <- x1 - stats::median(x1)
  r2 <- x2 - stats::median(x2)
  pool <- c(r1, r2)
  N <- length(pool); n1 <- length(x1)

  perm_diffs <- numeric(nperm)
  for (b in seq_len(nperm)) {
    idx <- sample.int(N, n1, replace = FALSE)
    perm_diffs[b] <- stats::median(pool[-idx]) - stats::median(pool[idx])
  }

  tol <- 1e-10
  (1.0 + sum(abs(perm_diffs) >= (abs(obs_diff) - tol))) / (nperm + 1.0)
}

# ------------------------------------------------------------------------------
# 2. Data Generator for Location Shifts
# ------------------------------------------------------------------------------
generate_location_data <- function(n1, n2, dist, delta = 0.0) {
  gen_raw <- function(n) {
    if (dist == "norm") {
      rnorm(n, mean = 0, sd = 1)
    } else if (dist == "lnorm") {
      rlnorm(n, meanlog = 0, sdlog = 0.8) # Right-skewed
    } else if (dist == "t3") {
      rt(n, df = 3)                       # Heavy tails
    } else if (dist == "exp") {
      rexp(n, rate = 1)                   # Exponential decay
    } else if (dist == "outlier") {
      vals <- rnorm(n)
      contam <- runif(n) < 0.02
      vals[contam] <- vals[contam] + 4.0  # 2% severe outliers
      vals
    } else {
      stop("Unknown dist: ", dist)
    }
  }

  x1 <- gen_raw(n1)
  x2 <- gen_raw(n2)

  # Apply standardized location shift delta
  if (delta != 0.0) {
    scale_factor <- if (dist == "lnorm") sqrt((exp(0.8^2) - 1) * exp(0.8^2))
    else if (dist == "t3") sqrt(3 / (3 - 2))
    else if (dist == "exp") 1.0
    else 1.0
    x2 <- x2 + delta * scale_factor
  }

  list(x1 = x1, x2 = x2)
}

# ------------------------------------------------------------------------------
# 3. Single Condition Worker
# ------------------------------------------------------------------------------
simulate_single_location_cell <- function(cond_idx, conditions, n_reps, nperm, alpha = 0.05) {
  cond <- conditions[cond_idx, ]
  n <- cond$n_sample
  dist <- cond$dist
  delta <- cond$delta

  p_mean_anal <- numeric(n_reps)
  p_mean_perm <- numeric(n_reps)
  p_welch     <- numeric(n_reps)
  p_student   <- numeric(n_reps)

  p_med_anal  <- numeric(n_reps)
  p_med_perm  <- numeric(n_reps)
  p_med_raw   <- numeric(n_reps)
  p_wilcox    <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    dat <- generate_location_data(n, n, dist = dist, delta = delta)
    x1 <- dat$x1; x2 <- dat$x2

    # --- A. Mean Tests ---
    # 1. Hermite Mean (Analytical closed-form)
    m_an <- tryCatch(t_hermite(x1, x2, method = "analytical", degree = 3L),
                     error = function(e) list(p_value = NA_real_))
    p_mean_anal[r] <- m_an$p_value

    # 2. Hermite Mean (Permutation)
    m_pm <- tryCatch(t_hermite(x1, x2, method = "permutation", degree = 3L, nperm = nperm),
                     error = function(e) list(p_value = NA_real_))
    p_mean_perm[r] <- m_pm$p_value

    # 3. Welch t-test
    tw <- tryCatch(stats::t.test(x2, x1, var.equal = FALSE), error = function(e) list(p.value = NA_real_))
    p_welch[r] <- tw$p.value

    # 4. Student t-test
    ts <- tryCatch(stats::t.test(x2, x1, var.equal = TRUE), error = function(e) list(p.value = NA_real_))
    p_student[r] <- ts$p.value

    # --- B. Median Tests ---
    # 5. Hermite Median (Analytical closed-form)
    med_an <- tryCatch(median_hermite(x1, x2, method = "analytical", degree = 3L),
                       error = function(e) list(p_value = NA_real_))
    p_med_anal[r] <- med_an$p_value

    # 6. Hermite Median (Permutation)
    med_pm <- tryCatch(median_hermite(x1, x2, method = "permutation", degree = 3L, nperm = nperm),
                       error = function(e) list(p_value = NA_real_))
    p_med_perm[r] <- med_pm$p_value

    # 7. Raw Sample Median Permutation
    p_med_raw[r] <- raw_median_perm_test(x1, x2, nperm = nperm)

    # 8. Wilcoxon Rank-Sum Test
    wx <- tryCatch(stats::wilcox.test(x2, x1, exact = FALSE), error = function(e) list(p.value = NA_real_))
    p_wilcox[r] <- wx$p.value
  }

  rej_rate <- function(p) safe_mean(p < alpha)

  data.frame(
    condition_id     = cond_idx,
    n_sample         = n,
    dist             = dist,
    delta            = delta,

    # Mean rejection rates
    rate_mean_anal   = rej_rate(p_mean_anal),
    rate_mean_perm   = rej_rate(p_mean_perm),
    rate_welch       = rej_rate(p_welch),
    rate_student     = rej_rate(p_student),

    # Median rejection rates
    rate_med_anal    = rej_rate(p_med_anal),
    rate_med_perm    = rej_rate(p_med_perm),
    rate_med_raw     = rej_rate(p_med_raw),
    rate_wilcox      = rej_rate(p_wilcox)
  )
}

# ------------------------------------------------------------------------------
# 4. Simulation Runner
# ------------------------------------------------------------------------------
run_location_sim <- function(conditions, n_reps = 500L, nperm = 499L,
                             n_cores = 4L, alpha = 0.05, seed = 20260301L) {
  start_sim <- Sys.time()
  n_conds <- nrow(conditions)

  cat(sprintf("\nRunning %d location test conditions x %d replications (nperm = %d)...\n\n",
              n_conds, n_reps, nperm))

  pkg_root <- if (file.exists("DESCRIPTION")) normalizePath(".") else normalizePath("..")

  cl <- NULL
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)

    parallel::clusterCall(cl, function(path) {
      devtools::load_all(path, quiet = TRUE)
      NULL
    }, path = pkg_root)

    parallel::clusterExport(
      cl,
      c("generate_location_data", "raw_median_perm_test", "safe_mean",
        "simulate_single_location_cell"),
      envir = .GlobalEnv
    )
    clusterSetRNGStream(cl, seed)
  } else {
    set.seed(seed)
  }

  res_list <- pblapply(
    seq_len(n_conds),
    simulate_single_location_cell,
    conditions = conditions,
    n_reps     = n_reps,
    nperm      = nperm,
    alpha      = alpha,
    cl         = cl
  )

  agg <- dplyr::bind_rows(res_list)
  cat("\n== Location simulation completed in:", format(Sys.time() - start_sim), "==\n\n")
  agg
}

# ------------------------------------------------------------------------------
# 5. Summary & Win-Rate Tables
# ------------------------------------------------------------------------------
summary.sim_location <- function(agg) {
  cat("\n========================================================================================\n")
  cat("             HERMITE MEAN & MEDIAN TESTS SIMULATION BENCHMARK SUMMARY\n")
  cat("========================================================================================\n\n")

  # 1. Type I Error Control (delta = 0)
  t1 <- agg %>%
    dplyr::filter(delta == 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Mean_Herm_Anal = safe_mean(rate_mean_anal),
      Mean_Herm_Perm = safe_mean(rate_mean_perm),
      Mean_Welch     = safe_mean(rate_welch),
      Med_Herm_Anal  = safe_mean(rate_med_anal),
      Med_Herm_Perm  = safe_mean(rate_med_perm),
      Med_Raw_Perm   = safe_mean(rate_med_raw),
      Wilcoxon       = safe_mean(rate_wilcox),
      .groups = "drop"
    )

  cat("--- 1. TYPE I ERROR CONTROL (Nominal alpha = 0.05, delta = 0.0) -----------------------\n")
  print.data.frame(t1, digits = 4, row.names = FALSE)

  # 2. Statistical Power (delta > 0)
  pwr <- agg %>%
    dplyr::filter(delta > 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Mean_Herm_Anal = safe_mean(rate_mean_anal),
      Mean_Herm_Perm = safe_mean(rate_mean_perm),
      Mean_Welch     = safe_mean(rate_welch),
      Med_Herm_Anal  = safe_mean(rate_med_anal),
      Med_Herm_Perm  = safe_mean(rate_med_perm),
      Med_Raw_Perm   = safe_mean(rate_med_raw),
      Wilcoxon       = safe_mean(rate_wilcox),
      .groups = "drop"
    )

  cat("\n--- 2. AVERAGE STATISTICAL POWER (delta > 0) -------------------------------------------\n")
  print.data.frame(pwr, digits = 4, row.names = FALSE)

  # 3. Head-to-Head Win Rates
  pwr_conds <- agg %>% dplyr::filter(delta > 0.0)
  total_cells <- nrow(pwr_conds)

  cat("\n--- 3. HEAD-TO-HEAD POWER COMPARISONS --------------------------------------------------\n")
  # Mean tests
  win_m_anal_w <- sum(pwr_conds$rate_mean_anal >= pwr_conds$rate_welch)
  win_m_perm_w <- sum(pwr_conds$rate_mean_perm >= pwr_conds$rate_welch)
  cat(sprintf("  Mean: Hermite (Analytical)  >= Welch's t-test    : %2d / %2d cells (%5.1f%%)\n",
              win_m_anal_w, total_cells, 100 * win_m_anal_w / total_cells))
  cat(sprintf("  Mean: Hermite (Permutation) >= Welch's t-test    : %2d / %2d cells (%5.1f%%)\n",
              win_m_perm_w, total_cells, 100 * win_m_perm_w / total_cells))

  # Median tests
  win_med_anal_raw <- sum(pwr_conds$rate_med_anal >= pwr_conds$rate_med_raw)
  win_med_perm_raw <- sum(pwr_conds$rate_med_perm >= pwr_conds$rate_med_raw)
  win_med_anal_wlx <- sum(pwr_conds$rate_med_anal >= pwr_conds$rate_wilcox)
  cat(sprintf("  Median: Hermite (Analytical) >= Raw Median Perm   : %2d / %2d cells (%5.1f%%)\n",
              win_med_anal_raw, total_cells, 100 * win_med_anal_raw / total_cells))
  cat(sprintf("  Median: Hermite (Permutation)>= Raw Median Perm   : %2d / %2d cells (%5.1f%%)\n",
              win_med_perm_raw, total_cells, 100 * win_med_perm_raw / total_cells))
  cat(sprintf("  Median: Hermite (Analytical) >= Wilcoxon (Rank)   : %2d / %2d cells (%5.1f%%)\n",
              win_med_anal_wlx, total_cells, 100 * win_med_anal_wlx / total_cells))
  cat("========================================================================================\n\n")
}

# ------------------------------------------------------------------------------
# 6. Visualizations
# ------------------------------------------------------------------------------
plot_location_comparison <- function(agg, dist_focus = "lnorm") {
  sub <- agg %>% dplyr::filter(dist == dist_focus, delta > 0.0)
  deltas <- sort(unique(sub$delta))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, length(deltas)), mar = c(4.5, 4.2, 3.2, 1.2), font.main = 1)

  for (d in deltas) {
    dat <- sub %>% dplyr::filter(delta == d) %>% dplyr::arrange(n_sample)

    plot(dat$n_sample, dat$rate_med_anal, type = "b", pch = 19, col = "blue3", lwd = 2.5,
         ylim = c(0, 1), log = "x", xlab = "Sample Size (n per group)", ylab = "Statistical Power",
         main = sprintf("%s (delta = %.2f)", tools::toTitleCase(dist_focus), d), bty = "l",
         panel.first = grid(col = "gray90", lty = 1))

    lines(dat$n_sample, dat$rate_mean_anal, type = "b", pch = 17, col = "skyblue3",   lwd = 2, lty = 2)
    lines(dat$n_sample, dat$rate_welch,     type = "b", pch = 18, col = "orchid",     lwd = 2, lty = 3)
    lines(dat$n_sample, dat$rate_med_raw,   type = "b", pch = 4,  col = "darkgreen",  lwd = 2, lty = 4)
    lines(dat$n_sample, dat$rate_wilcox,    type = "b", pch = 15, col = "darkorange", lwd = 2, lty = 5)

    if (d == deltas[1L]) {
      legend("bottomright",
             legend = c("Hermite Median (Analyt)", "Hermite Mean (Analyt)", "Welch Mean t", "Raw Median Perm", "Wilcoxon Rank"),
             col = c("blue3", "skyblue3", "orchid", "darkgreen", "darkorange"),
             pch = c(19, 17, 18, 4, 15), lty = c(1, 2, 3, 4, 5), lwd = 2, bty = "n", cex = 0.75)
    }
  }
}

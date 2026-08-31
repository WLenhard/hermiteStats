# ==============================================================================
# sim_test_helper.R -- Comprehensive Simulation Infrastructure for t_hermite
# ==============================================================================
#
# Benchmark Comparators:
#   1. Hermite (Permutation): Regularized Permutation Test (t_Hermite)
#   2. Welch t-test:          Heteroscedastic Parametric t-Test
#   3. Student t-test:        Classical Pooled Equal-Variance t-Test
#   4. Wilcoxon:              Mann-Whitney U Rank-Sum Test
#   5. Yuen t-test:           Robust 20% Trimmed-Mean t-Test
#
# ==============================================================================

library(parallel)
library(dplyr)
if (!requireNamespace("pbapply", quietly = TRUE)) {
  install.packages("pbapply")
}
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
# 2. Data Generation
# ------------------------------------------------------------------------------
generate_group_data <- function(n1, n2, dist, param1, param2, param3 = NA, delta = 0.0) {
  gen_raw <- function(n) {
    if (dist == "norm") {
      rnorm(n, mean = param1, sd = param2)
    } else if (dist == "lnorm") {
      rlnorm(n, meanlog = param1, sdlog = param2)
    } else if (dist == "t") {
      rt(n, df = param1) * param2
    } else if (dist == "exp") {
      rexp(n, rate = param1)
    } else if (dist == "outlier") {
      contam_prop <- if (is.finite(param3)) param3 else 0.02
      contam <- runif(n) < contam_prop
      vals <- rnorm(n)
      vals[contam] <- vals[contam] + param1
      vals
    } else if (dist == "1plirt") {
      items <- seq(-2, 2, length.out = if (is.na(param3)) 20L else as.integer(param3))
      theta <- rnorm(n)
      pm <- plogis(outer(theta, items, "-"))
      as.integer(rowSums(matrix(runif(n * length(items)), n) < pm))
    } else {
      stop("Unknown distribution: ", dist)
    }
  }

  x1 <- gen_raw(n1)
  x2 <- gen_raw(n2)

  if (delta != 0.0) {
    shift_val <- delta * (if (dist == "lnorm") sqrt((exp(param2^2) - 1) * exp(2 * param1 + param2^2))
                          else if (dist == "t") sqrt(param1 / (param1 - 2)) * param2
                          else if (dist == "exp") 1 / param1
                          else if (dist == "norm") param2
                          else stats::sd(x1))
    x2 <- x2 + shift_val
  }
  list(x1 = x1, x2 = x2)
}

# ------------------------------------------------------------------------------
# 3. Ground Truth Precomputation
# ------------------------------------------------------------------------------
.cond_key_test <- function(cond) {
  p3 <- if ("param3" %in% names(cond)) cond$param3 else NA
  paste(cond$dist, cond$param1, cond$param2, p3, cond$delta, sep = "_")
}

compute_test_population_truths <- function(conditions, pop_size = 5e5, seed = 1234) {
  keys <- vapply(seq_len(nrow(conditions)), function(i) .cond_key_test(conditions[i, ]), character(1L))
  uniq <- which(!duplicated(keys))

  cat(sprintf("\nPre-computing reference ground truths (%d unique distributions x effect sizes, N = %g)...\n",
              length(uniq), pop_size))

  truth_cache <- list()
  pb <- txtProgressBar(min = 0, max = length(uniq), style = 3)

  for (j in seq_along(uniq)) {
    i <- uniq[j]
    cond <- conditions[i, ]
    p3 <- if ("param3" %in% names(cond)) cond$param3 else NA

    set.seed(seed + j)
    pop <- generate_group_data(pop_size, pop_size, cond$dist, cond$param1, cond$param2, p3, delta = cond$delta)

    m1 <- mean(pop$x1); m2 <- mean(pop$x2)
    s1 <- sd(pop$x1);   s2 <- sd(pop$x2)
    s_avg <- sqrt((s1^2 + s2^2) / 2.0)
    true_d <- (m2 - m1) / s_avg

    truth_cache[[keys[i]]] <- list(
      mean_diff  = m2 - m1,
      true_s_avg = s_avg,
      true_d     = true_d
    )
    setTxtProgressBar(pb, j)
  }
  close(pb)
  truth_cache
}

# ------------------------------------------------------------------------------
# 4. Per-Condition Simulation Worker
# ------------------------------------------------------------------------------
simulate_single_test_condition <- function(cond_idx, conditions, truth_cache,
                                           n_reps, nperm, alpha = 0.05,
                                           seed_base = 20260101L) {
  cond <- conditions[cond_idx, ]
  truth <- truth_cache[[.cond_key_test(cond)]]
  p3 <- if ("param3" %in% names(cond)) cond$param3 else NA

  n <- cond$n_sample
  delta <- cond$delta
  true_d <- truth$true_d

  set.seed(as.integer(seed_base + cond_idx))

  p_hermite_perm <- numeric(n_reps)
  p_welch        <- numeric(n_reps)
  p_student      <- numeric(n_reps)
  p_wilcox       <- numeric(n_reps)
  p_yuen         <- numeric(n_reps)

  est_d_reg    <- numeric(n_reps)
  est_hedges_g <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    dat <- generate_group_data(n, n, cond$dist, cond$param1, cond$param2, p3, delta = delta)
    x1 <- dat$x1; x2 <- dat$x2

    # 1. Hermite Permutation Test
    res_perm <- tryCatch(
      t_hermite(x1, x2, method = "permutation", degree = 3L, nperm = nperm),
      error = function(e) list(p_value = NA_real_, estimate = NA_real_, d_reg_fit = list(hedges_g = NA_real_))
    )
    p_hermite_perm[r] <- res_perm$p_value
    est_d_reg[r]      <- res_perm$estimate
    est_hedges_g[r]   <- if (!is.null(res_perm$d_reg_fit$hedges_g)) res_perm$d_reg_fit$hedges_g else NA_real_

    # 2. Welch t-test
    res_tw <- tryCatch(stats::t.test(x2, x1, var.equal = FALSE), error = function(e) list(p.value = NA_real_))
    p_welch[r] <- res_tw$p.value

    # 3. Student t-test
    res_ts <- tryCatch(stats::t.test(x2, x1, var.equal = TRUE), error = function(e) list(p.value = NA_real_))
    p_student[r] <- res_ts$p.value

    # 4. Wilcoxon / Mann-Whitney
    res_wx <- tryCatch(stats::wilcox.test(x2, x1, exact = FALSE), error = function(e) list(p.value = NA_real_))
    p_wilcox[r] <- res_wx$p.value

    # 5. Yuen 20% Trimmed t-test
    res_yu <- tryCatch(yuen_test(x1, x2, tr = 0.2), error = function(e) list(p_value = NA_real_))
    p_yuen[r] <- res_yu$p_value
  }

  rej_rate <- function(p) safe_mean(p < alpha)

  data.frame(
    condition_id      = cond_idx,
    n_sample          = n,
    dist              = cond$dist,
    delta             = delta,
    true_d            = true_d,

    # Rejection rates (Alpha if delta=0, Power if delta>0)
    rate_hermite_perm = rej_rate(p_hermite_perm),
    rate_welch        = rej_rate(p_welch),
    rate_student      = rej_rate(p_student),
    rate_wilcox       = rej_rate(p_wilcox),
    rate_yuen         = rej_rate(p_yuen),

    # Effect size recovery (d_reg vs Hedges' g)
    bias_d_reg        = safe_mean(est_d_reg - true_d),
    var_d_reg         = var(est_d_reg, na.rm = TRUE),
    mse_d_reg         = safe_mean((est_d_reg - true_d)^2),

    bias_hedges_g     = safe_mean(est_hedges_g - true_d),
    var_hedges_g      = var(est_hedges_g, na.rm = TRUE),
    mse_hedges_g      = safe_mean((est_hedges_g - true_d)^2)
  )
}

# ------------------------------------------------------------------------------
# 5. Simulation Runner
# ------------------------------------------------------------------------------
run_test_sim <- function(conditions, n_reps = 500, nperm = 1000L,
                         n_cores = 4, alpha = 0.05, seed = 20260101L,
                         truth_pop_size = 5e5) {
  start_sim <- Sys.time()

  truth_cache <- compute_test_population_truths(
    conditions = conditions,
    pop_size   = truth_pop_size,
    seed       = seed
  )

  n_conds <- nrow(conditions)
  cat(sprintf("Running %d test conditions x %d Monte Carlo replications (nperm = %d)...\n\n",
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
      c("generate_group_data", ".cond_key_test", "yuen_test", "safe_mean",
        "simulate_single_test_condition"),
      envir = .GlobalEnv
    )
  }

  res_list <- pblapply(
    seq_len(n_conds),
    simulate_single_test_condition,
    conditions  = conditions,
    truth_cache = truth_cache,
    n_reps      = n_reps,
    nperm       = nperm,
    alpha       = alpha,
    seed_base   = seed,
    cl          = cl
  )

  agg <- dplyr::bind_rows(res_list)
  cat("\n== Test simulation finished in:", format(Sys.time() - start_sim), "==\n\n")

  list(aggregated = agg, truth_cache = truth_cache)
}

# ------------------------------------------------------------------------------
# 6. Summary & Head-to-Head Win Rates
# ------------------------------------------------------------------------------
summary.sim_test <- function(agg) {
  cat("\n========================================================================================\n")
  cat("                    HERMITE t-TEST SIMULATION BENCHMARK SUMMARY\n")
  cat("========================================================================================\n\n")

  # 1. Type I Error Control (delta = 0)
  t1 <- agg %>%
    dplyr::filter(delta == 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Alpha_t_Hermite = safe_mean(rate_hermite_perm),
      Alpha_Welch     = safe_mean(rate_welch),
      Alpha_Student   = safe_mean(rate_student),
      Alpha_Wilcoxon  = safe_mean(rate_wilcox),
      Alpha_Yuen      = safe_mean(rate_yuen),
      .groups = "drop"
    )

  cat("--- 1. TYPE I ERROR CONTROL (Nominal alpha = 0.05) ------------------------------------\n")
  print.data.frame(t1, digits = 4, row.names = FALSE)

  # 2. Statistical Power (delta > 0)
  pwr <- agg %>%
    dplyr::filter(delta > 0.0) %>%
    dplyr::group_by(dist) %>%
    dplyr::summarise(
      Power_t_Hermite = safe_mean(rate_hermite_perm),
      Power_Welch     = safe_mean(rate_welch),
      Power_Student   = safe_mean(rate_student),
      Power_Wilcoxon  = safe_mean(rate_wilcox),
      Power_Yuen      = safe_mean(rate_yuen),
      .groups = "drop"
    )

  cat("\n--- 2. AVERAGE STATISTICAL POWER (delta > 0) -------------------------------------------\n")
  print.data.frame(pwr, digits = 4, row.names = FALSE)

  # 3. Comprehensive Head-to-Head Power Win Rates
  pwr_conds <- agg %>% dplyr::filter(delta > 0.0)
  total_cells <- nrow(pwr_conds)

  comparators <- list(
    list(col = "rate_welch",   name = "Welch's t-test"),
    list(col = "rate_student", name = "Student's t-test"),
    list(col = "rate_wilcox",  name = "Wilcoxon (Mann-Whitney U)"),
    list(col = "rate_yuen",    name = "Yuen's trimmed t-test")
  )

  cat("\n--- 3. HEAD-TO-HEAD POWER COMPARISONS (t_Hermite vs. Benchmarks) ----------------------\n")
  for (comp in comparators) {
    wins <- sum(pwr_conds$rate_hermite_perm >= pwr_conds[[comp$col]])
    cat(sprintf("  t_Hermite >= %-26s : %2d / %2d condition cells (%5.1f%%)\n",
                comp$name, wins, total_cells, 100 * wins / total_cells))
  }
  cat("========================================================================================\n\n")
}

# ------------------------------------------------------------------------------
# 7. Diagnostic Visualizations
# ------------------------------------------------------------------------------

#' Multi-Panel Type I Error Rate Visualization
#' @export
plot_type1_error <- function(agg, by = c("sample_size", "distribution"), nominal = 0.05) {
  by <- match.arg(by)
  t1 <- agg %>% dplyr::filter(delta == 0.0)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (by == "sample_size") {
    # 6-Panel Grid: Type I error across sample sizes for each distribution
    dists <- sort(unique(t1$dist))
    par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3.2, 1.2), font.main = 1)

    for (d in dists) {
      sub_d <- t1 %>% dplyr::filter(dist == d) %>% dplyr::arrange(n_sample)

      plot(sub_d$n_sample, sub_d$rate_hermite_perm, type = "b", pch = 19, col = "blue3", lwd = 2.5,
           ylim = c(0.01, 0.09), log = "x",
           xlab = "Sample Size (n per group)", ylab = "Empirical Alpha",
           main = tools::toTitleCase(d), bty = "l",
           panel.first = grid(col = "gray90", lty = 1))

      lines(sub_d$n_sample, sub_d$rate_welch,   type = "b", pch = 17, col = "orchid",     lwd = 1.8, lty = 2)
      lines(sub_d$n_sample, sub_d$rate_student, type = "b", pch = 18, col = "gray50",    lwd = 1.8, lty = 3)
      lines(sub_d$n_sample, sub_d$rate_wilcox,  type = "b", pch = 15, col = "darkorange", lwd = 1.8, lty = 4)
      lines(sub_d$n_sample, sub_d$rate_yuen,    type = "b", pch = 4,  col = "darkcyan",   lwd = 1.8, lty = 5)

      abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
      abline(h = c(0.025, 0.075), col = "gray40", lty = 3) # Bradley's bounds

      if (d == dists[1L]) {
        legend("topright",
               legend = c("t_Hermite", "Welch t", "Student t", "Wilcoxon", "Yuen"),
               col = c("blue3", "orchid", "gray50", "darkorange", "darkcyan"),
               pch = c(19, 17, 18, 15, 4), lty = c(1, 2, 3, 4, 5), lwd = 1.8,
               bty = "n", cex = 0.75)
      }
    }
  } else {
    # Grouped Summary Barplot across distributions
    par(mar = c(6, 4.5, 3.5, 1.5))
    by_dist <- t1 %>%
      dplyr::group_by(dist) %>%
      dplyr::summarise(
        Hermite  = safe_mean(rate_hermite_perm),
        Welch    = safe_mean(rate_welch),
        Student  = safe_mean(rate_student),
        Wilcoxon = safe_mean(rate_wilcox),
        Yuen     = safe_mean(rate_yuen),
        .groups = "drop"
      )

    mat <- as.matrix(by_dist[, -1])
    rownames(mat) <- by_dist$dist

    barplot(t(mat), beside = TRUE,
            col = c("blue3", "orchid", "gray60", "darkorange", "darkcyan"),
            ylim = c(0, 0.10), ylab = "Empirical Type I Error",
            main = "Type I Error Rate Control Across Distributions (alpha = 0.05)",
            las = 2)
    abline(h = nominal, lty = 2, col = "red", lwd = 2)
    abline(h = c(0.025, 0.075), lty = 3, col = "gray40")

    legend("topright", legend = c("t_Hermite", "Welch t", "Student t", "Wilcoxon", "Yuen"),
           fill = c("blue3", "orchid", "gray60", "darkorange", "darkcyan"), bty = "n", cex = 0.85)
  }
}

#' Multi-Panel Statistical Power Visualization
#' @export
plot_power_curves <- function(agg, dist_focus = "lnorm") {
  sub <- agg %>% dplyr::filter(dist == dist_focus, delta > 0.0)
  deltas <- sort(unique(sub$delta))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, length(deltas)), mar = c(4.5, 4.2, 3, 1), font.main = 1)

  for (d in deltas) {
    dat_d <- sub %>% dplyr::filter(delta == d) %>% dplyr::arrange(n_sample)

    plot(dat_d$n_sample, dat_d$rate_hermite_perm, type = "b", pch = 19, col = "blue3", lwd = 2.5,
         ylim = c(0, 1), log = "x", xlab = "Sample Size (n per group)", ylab = "Statistical Power",
         main = sprintf("%s (delta = %.1f)", tools::toTitleCase(dist_focus), d), bty = "l",
         panel.first = grid(col = "gray90", lty = 1))

    lines(dat_d$n_sample, dat_d$rate_welch,   type = "b", pch = 17, col = "orchid",     lwd = 2, lty = 2)
    lines(dat_d$n_sample, dat_d$rate_student, type = "b", pch = 18, col = "gray50",    lwd = 2, lty = 3)
    lines(dat_d$n_sample, dat_d$rate_wilcox,  type = "b", pch = 15, col = "darkorange", lwd = 2, lty = 4)
    lines(dat_d$n_sample, dat_d$rate_yuen,    type = "b", pch = 4,  col = "darkcyan",   lwd = 2, lty = 5)

    if (d == deltas[1L]) {
      legend("bottomright",
             legend = c("t_Hermite", "Welch t", "Student t", "Wilcoxon", "Yuen"),
             col = c("blue3", "orchid", "gray50", "darkorange", "darkcyan"),
             pch = c(19, 17, 18, 15, 4), lty = c(1, 2, 3, 4, 5), lwd = 2,
             bty = "n", cex = 0.8)
    }
  }
}

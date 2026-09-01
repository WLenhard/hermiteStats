# ==============================================================================
# sim_shape_helper.R -- Simulation infrastructure for shape_hermite
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
  if (!length(x)) NA_real_ else mean(x)
}

# ------------------------------------------------------------------------------
# 1. Classical and Non-Parametric Benchmark Statistics
# ------------------------------------------------------------------------------

g1_raw <- function(v) {
  vc <- v - mean(v)
  m2 <- mean(vc^2)
  if (m2 <= 1e-12) 0 else mean(vc^3) / m2^1.5
}

g2_raw <- function(v) {
  vc <- v - mean(v)
  m2 <- mean(vc^2)
  if (m2 <= 1e-12) 0 else mean(vc^4) / m2^2 - 3
}

logvar <- function(v) log(max(stats::var(v), 1e-12))

ks_stat <- function(a, b) {
  na <- length(a)
  nb <- length(b)
  o  <- order(c(a, b))
  lab <- c(rep(1, na), rep(0, nb))[o]
  max(abs(cumsum(lab) / na - cumsum(1 - lab) / nb))
}

# Pettitt form of the two-sample Anderson-Darling statistic
ad_stat <- function(a, b) {
  na <- length(a)
  nb <- length(b)
  N  <- na + nb
  o  <- order(c(a, b))
  lab <- c(rep(1L, na), rep(0L, nb))[o]
  M  <- cumsum(lab)[seq_len(N - 1L)]
  i  <- seq_len(N - 1L)
  sum((M * N - na * i)^2 / (i * (N - i))) / (na * nb)
}

bf_p <- function(a, b) {
  z1 <- abs(a - stats::median(a))
  z2 <- abs(b - stats::median(b))
  tryCatch(stats::t.test(z1, z2, var.equal = TRUE)$p.value, error = function(e) NA_real_)
}

fk_p <- function(a, b) {
  tryCatch(stats::fligner.test(list(a, b))$p.value, error = function(e) NA_real_)
}

f_p <- function(a, b) {
  tryCatch(stats::var.test(a, b)$p.value, error = function(e) NA_real_)
}

# ------------------------------------------------------------------------------
# 2. Shared Permutation Engine for Resampling Benchmarks
# ------------------------------------------------------------------------------
perm_bank <- function(x1, x2, perm_idx, stats_spec) {
  n1 <- length(x1)
  n2 <- length(x2)
  pools <- list(
    raw = c(x1, x2),
    loc = c(x1 - stats::median(x1), x2 - stats::median(x2)),
    std = c((x1 - stats::median(x1)) / stats::mad(x1),
            (x2 - stats::median(x2)) / stats::mad(x2))
  )
  out <- numeric(length(stats_spec))
  names(out) <- names(stats_spec)
  B <- nrow(perm_idx)

  for (nm in names(stats_spec)) {
    sp   <- stats_spec[[nm]]
    pool <- pools[[sp$align]]
    if (any(!is.finite(pool))) { out[nm] <- NA_real_; next }
    obs <- if (sp$align != "raw") sp$fun(x1, x2) else sp$fun(pool[seq_len(n1)], pool[n1 + seq_len(n2)])
    pv <- numeric(B)
    for (b in seq_len(B)) {
      idx <- perm_idx[b, ]
      pv[b] <- sp$fun(pool[idx], pool[-idx])
    }
    ctr <- mean(pv, na.rm = TRUE)
    out[nm] <- (1 + sum(abs(pv - ctr) >= abs(obs - ctr) - 1e-10, na.rm = TRUE)) / (B + 1)
  }
  out
}

# ------------------------------------------------------------------------------
# 3. Data Generating Scenarios
# ------------------------------------------------------------------------------
std <- function(v, m, s) (v - m) / s

SCENARIOS <- list(
  # ---- Complete nulls: F1 = F2 -----------------------------------------------
  null_norm  = list(kind = "null", target = "-", ref = NA,
                    gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rnorm(n2))),
  null_lnorm = list(kind = "null", target = "-", ref = NA,
                    gen = function(n1, n2) list(x1 = rlnorm(n1, 0, .6), x2 = rlnorm(n2, 0, .6))),
  null_t3    = list(kind = "null", target = "-", ref = NA,
                    gen = function(n1, n2) list(x1 = rt(n1, 3), x2 = rt(n2, 3))),
  null_unif  = list(kind = "null", target = "-", ref = NA,
                    gen = function(n1, n2) list(x1 = runif(n1), x2 = runif(n2))),

  # ---- Partial nulls: tested aspect is null, others differ -------------------
  pnull_loc  = list(kind = "pnull", target = "all", ref = NA,
                    gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rnorm(n2, mean = 1))),
  pnull_loc_skew = list(kind = "pnull", target = "all", ref = NA,
                        gen = function(n1, n2) list(x1 = rlnorm(n1, 0, .6), x2 = rlnorm(n2, 0, .6) + 2)),
  pnull_scale = list(kind = "pnull", target = "shape", ref = NA,
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rnorm(n2, sd = 1.5))),
  pnull_shape = list(kind = "pnull", target = "scale", ref = NA,
                     gen = function(n1, n2) list(x1 = rnorm(n1),
                                                 x2 = std(rgamma(n2, shape = 2), 2, sqrt(2)))),

  # ---- Scale alternatives ---------------------------------------------------
  scale_norm  = list(kind = "power", target = "scale", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rnorm(n2, sd = 1.5))),
  scale_t3    = list(kind = "power", target = "scale", ref = "null_t3",
                     gen = function(n1, n2) list(x1 = rt(n1, 3), x2 = 1.5 * rt(n2, 3))),
  scale_lnorm = list(kind = "power", target = "scale", ref = "null_lnorm",
                     gen = function(n1, n2) list(x1 = rlnorm(n1, 0, .6), x2 = 1.5 * rlnorm(n2, 0, .6))),
  scale_contam= list(kind = "power", target = "scale", ref = "null_norm",
                     gen = function(n1, n2) {
                       f <- function(n) { v <- rnorm(n); k <- runif(n) < .05; v[k] <- v[k] * 4; v }
                       list(x1 = f(n1), x2 = 1.5 * f(n2))
                     }),
  scale_shift = list(kind = "power", target = "scale", ref = "pnull_loc",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rnorm(n2, mean = 0.8, sd = 1.5))),

  # ---- Asymmetry alternatives -----------------------------------------------
  asym_gamma  = list(kind = "power", target = "asymmetry", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = std(rgamma(n2, shape = 2), 2, sqrt(2)))),
  asym_gamma8 = list(kind = "power", target = "asymmetry", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = std(rgamma(n2, shape = 8), 8, sqrt(8)))),
  asym_lnorm  = list(kind = "power", target = "asymmetry", ref = "null_norm",
                     gen = function(n1, n2) {
                       s <- .4; m <- exp(s^2 / 2); v <- sqrt((exp(s^2) - 1) * exp(s^2))
                       list(x1 = rnorm(n1), x2 = std(rlnorm(n2, 0, s), m, v))
                     }),

  # ---- Tail-weight alternatives ---------------------------------------------
  tail_t5     = list(kind = "power", target = "tailweight", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rt(n2, 5) * sqrt(3 / 5))),
  tail_t3     = list(kind = "power", target = "tailweight", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = rt(n2, 3) / sqrt(3))),
  tail_unif   = list(kind = "power", target = "tailweight", ref = "null_norm",
                     gen = function(n1, n2) list(x1 = rnorm(n1), x2 = std(runif(n2), .5, sqrt(1 / 12)))),
  tail_contam = list(kind = "power", target = "tailweight", ref = "null_norm",
                     gen = function(n1, n2) {
                       v <- rnorm(n2); k <- runif(n2) < .05; v[k] <- v[k] * 4
                       list(x1 = rnorm(n1), x2 = std(v, mean(v), stats::sd(v)))
                     })
)

generate_shape_data <- function(n1, n2, scenario) SCENARIOS[[scenario]]$gen(n1, n2)

# ------------------------------------------------------------------------------
# 4. Population Projection of Hermite Indices (Estimand Transparency)
# ------------------------------------------------------------------------------
shape_population_truth <- function(scenarios = names(SCENARIOS), N = 2e5, seed = 99) {
  set.seed(seed)
  do.call(rbind, lapply(scenarios, function(s) {
    d  <- SCENARIOS[[s]]$gen(N, N)
    p1 <- hermiteStats:::.hermite_profile(d$x1, degree = 3L)
    p2 <- hermiteStats:::.hermite_profile(d$x2, degree = 3L)

    if (is.null(p1) || is.null(p2)) {
      stop(sprintf("Quantile projection failed for scenario '%s'.", s))
    }

    data.frame(
      scenario  = s,
      kind      = SCENARIOS[[s]]$kind,
      target    = SCENARIOS[[s]]$target,
      d_scale   = p2$log_sd - p1$log_sd,
      d_asym    = p2$c2     - p1$c2,
      d_tail    = p2$c3     - p1$c3,
      true_g1_1 = g1_raw(d$x1),
      true_g1_2 = g1_raw(d$x2),
      true_g2_1 = g2_raw(d$x1),
      true_g2_2 = g2_raw(d$x2),
      stringsAsFactors = FALSE
    )
  }))
}

# ------------------------------------------------------------------------------
# 5. Worker: Cell Simulator
# ------------------------------------------------------------------------------
simulate_shape_cell <- function(cond_idx, conditions, n_reps, nperm, ...) {
  cond     <- conditions[cond_idx, ]
  n1       <- cond$n1
  n2       <- cond$n2
  scenario <- as.character(cond$scenario)
  N        <- n1 + n2

  cols <- c("herm_scale", "herm_asym", "herm_tail",
            "herm_scale_wy", "herm_asym_wy", "herm_tail_wy",
            "herm_omni_minP", "herm_omni_maxT",
            "herm_scale_ex", "herm_asym_ex", "herm_tail_ex",
            "f_test", "bf", "fligner",
            "raw_logvar_perm", "raw_g1_perm", "raw_g2_perm",
            "ks_perm", "ad_perm")

  P <- matrix(NA_real_, nrow = n_reps, ncol = length(cols), dimnames = list(NULL, cols))
  deg_ok <- logical(n_reps)

  spec <- list(
    raw_logvar_perm = list(fun = function(a, b) logvar(b) - logvar(a), align = "loc"),
    raw_g1_perm     = list(fun = function(a, b) g1_raw(b) - g1_raw(a), align = "std"),
    raw_g2_perm     = list(fun = function(a, b) g2_raw(b) - g2_raw(a), align = "std"),
    ks_perm         = list(fun = function(a, b) ks_stat(a, b),          align = "raw"),
    ad_perm         = list(fun = function(a, b) ad_stat(a, b),          align = "raw")
  )

  for (r in seq_len(n_reps)) {
    d  <- generate_shape_data(n1, n2, scenario)
    x1 <- d$x1
    x2 <- d$x2

    # Shared permutation indices across non-parametric procedures
    pidx <- t(vapply(seq_len(nperm), function(i) sample.int(N, n1, replace = FALSE), integer(n1)))

    # Location-aligned Hermite shape test (canonical)
    sh <- tryCatch(
      shape_hermite(x1, x2, degree = 3L, align = "location",
                    contrasts = c("asymmetry", "tailweight", "scale"),
                    omnibus = "minP", nperm = nperm),
      error = function(e) NULL
    )

    if (!is.null(sh)) {
      P[r, "herm_scale"]    <- sh$p_value[["scale"]]
      P[r, "herm_asym"]     <- sh$p_value[["asymmetry"]]
      P[r, "herm_tail"]     <- sh$p_value[["tailweight"]]
      P[r, "herm_scale_wy"] <- sh$p_adjusted[["scale"]]
      P[r, "herm_asym_wy"]  <- sh$p_adjusted[["asymmetry"]]
      P[r, "herm_tail_wy"]  <- sh$p_adjusted[["tailweight"]]
      if (!is.null(sh$omnibus)) P[r, "herm_omni_minP"] <- sh$omnibus$p_value
      deg_ok[r] <- TRUE
    }

    # MaxT omnibus test comparison
    sh_maxT <- tryCatch(
      shape_hermite(x1, x2, degree = 3L, align = "location",
                    contrasts = c("asymmetry", "tailweight", "scale"),
                    omnibus = "maxT", nperm = nperm),
      error = function(e) NULL
    )
    if (!is.null(sh_maxT) && !is.null(sh_maxT$omnibus)) {
      P[r, "herm_omni_maxT"] <- sh_maxT$omnibus$p_value
    }

    # Exchangeable Hermite shape test (align = "none") for partial-null benchmarking
    sh_ex <- tryCatch(
      shape_hermite(x1, x2, degree = 3L, align = "none",
                    contrasts = c("asymmetry", "tailweight", "scale"),
                    omnibus = "none", nperm = nperm),
      error = function(e) NULL
    )
    if (!is.null(sh_ex)) {
      P[r, "herm_scale_ex"] <- sh_ex$p_value[["scale"]]
      P[r, "herm_asym_ex"]  <- sh_ex$p_value[["asymmetry"]]
      P[r, "herm_tail_ex"]  <- sh_ex$p_value[["tailweight"]]
    }

    # Classical tests
    P[r, "f_test"]  <- f_p(x1, x2)
    P[r, "bf"]      <- bf_p(x1, x2)
    P[r, "fligner"] <- fk_p(x1, x2)

    # Resampling benchmark statistics
    P[r, names(spec)] <- perm_bank(x1, x2, pidx, spec)
  }

  cbind(
    data.frame(condition_id = cond_idx, n1 = n1, n2 = n2, scenario = scenario,
               kind = SCENARIOS[[scenario]]$kind,
               target = SCENARIOS[[scenario]]$target,
               ref = SCENARIOS[[scenario]]$ref,
               rep = seq_len(n_reps), fit_ok = deg_ok,
               stringsAsFactors = FALSE),
    as.data.frame(P)
  )
}

# ------------------------------------------------------------------------------
# 6. Cluster Runner (L'Ecuyer-CMRG)
# ------------------------------------------------------------------------------
run_shape_sim <- function(conditions, n_reps = 1000L, nperm = 499L,
                          n_cores = 4L, seed = 20260201L) {
  t0 <- Sys.time()
  cat(sprintf("\nRunning %d cells x %d replications (nperm = %d)\n\n",
              nrow(conditions), n_reps, nperm))

  pkg_root <- if (file.exists("DESCRIPTION")) normalizePath(".") else normalizePath("..")

  cl <- NULL
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterCall(cl, function(p) { devtools::load_all(p, quiet = TRUE); NULL }, p = pkg_root)
    clusterExport(cl, c("SCENARIOS", "generate_shape_data", "perm_bank", "g1_raw", "g2_raw",
                        "logvar", "ks_stat", "ad_stat", "bf_p", "fk_p", "f_p", "std",
                        "safe_mean", "simulate_shape_cell"), envir = environment())
    clusterSetRNGStream(cl, seed)
  } else {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  res <- pblapply(seq_len(nrow(conditions)), simulate_shape_cell,
                  conditions = conditions, n_reps = n_reps, nperm = nperm, cl = cl)
  out <- dplyr::bind_rows(res)
  cat("\nCompleted in", format(Sys.time() - t0), "\n\n")
  out
}

# ------------------------------------------------------------------------------
# 7. Analysis: Rates, Size-Adjustment, Paired Contrasts & Win Rates
# ------------------------------------------------------------------------------
TESTS <- c("herm_scale", "herm_asym", "herm_tail", "herm_scale_wy", "herm_asym_wy",
           "herm_tail_wy", "herm_omni_minP", "herm_omni_maxT", "herm_scale_ex",
           "herm_asym_ex", "herm_tail_ex", "f_test", "bf", "fligner",
           "raw_logvar_perm", "raw_g1_perm", "raw_g2_perm", "ks_perm", "ad_perm")

rates_table <- function(pv, alpha = 0.05) {
  pv %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::group_by(scenario, kind, target, n1, n2, test) %>%
    dplyr::summarise(
      R    = sum(is.finite(p)),
      rate = safe_mean(p < alpha),
      mcse = sqrt(rate * (1 - rate) / R),
      .groups = "drop"
    )
}

size_adjusted <- function(pv, alpha = 0.05) {
  crit <- pv %>%
    dplyr::filter(kind == "null") %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::group_by(ref_scen = scenario, n1, n2, test) %>%
    dplyr::summarise(crit = stats::quantile(p, alpha, na.rm = TRUE), .groups = "drop")

  pv %>%
    dplyr::filter(kind == "power", !is.na(ref)) %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::left_join(crit, by = c("ref" = "ref_scen", "n1", "n2", "test")) %>%
    dplyr::group_by(scenario, target, n1, n2, test) %>%
    dplyr::summarise(power_adj = safe_mean(p <= crit), .groups = "drop")
}

paired_compare <- function(pv, a, b, alpha = 0.05) {
  pv %>%
    dplyr::filter(kind == "power") %>%
    dplyr::mutate(da = as.numeric(.data[[a]] < alpha) - as.numeric(.data[[b]] < alpha)) %>%
    dplyr::group_by(scenario, target, n1, n2) %>%
    dplyr::summarise(
      diff = mean(da, na.rm = TRUE),
      se   = stats::sd(da, na.rm = TRUE) / sqrt(sum(is.finite(da))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      A   = a,
      B   = b,
      sig = ifelse(abs(diff) > 1.96 * se, ifelse(diff > 0, "+", "-"), ".")
    )
}

#' Compute Global and Domain-Specific Win Rates for shape_hermite
win_rate_summary <- function(pv, alpha = 0.05) {
  pairs <- list(
    # Scale contrasts
    list(target = "scale", herm = "herm_scale", comp = "f_test", label = "Scale: Hermite vs F-test"),
    list(target = "scale", herm = "herm_scale", comp = "bf",     label = "Scale: Hermite vs Brown-Forsythe"),
    list(target = "scale", herm = "herm_scale", comp = "fligner",label = "Scale: Hermite vs Fligner-Killeen"),
    list(target = "scale", herm = "herm_scale", comp = "raw_logvar_perm", label = "Scale: Hermite vs Perm-LogVar"),
    # Asymmetry contrasts
    list(target = "asymmetry", herm = "herm_asym", comp = "raw_g1_perm", label = "Asym: Hermite vs Perm-g1"),
    list(target = "asymmetry", herm = "herm_asym", comp = "ks_perm",     label = "Asym: Hermite vs Perm-KS"),
    list(target = "asymmetry", herm = "herm_asym", comp = "ad_perm",     label = "Asym: Hermite vs Perm-AD"),
    # Tail-weight contrasts
    list(target = "tailweight", herm = "herm_tail", comp = "raw_g2_perm", label = "Tail: Hermite vs Perm-g2"),
    list(target = "tailweight", herm = "herm_tail", comp = "ks_perm",     label = "Tail: Hermite vs Perm-KS"),
    list(target = "tailweight", herm = "herm_tail", comp = "ad_perm",     label = "Tail: Hermite vs Perm-AD"),
    # Omnibus comparisons
    list(target = "all", herm = "herm_omni_minP", comp = "ks_perm", label = "Omnibus: minP vs Perm-KS"),
    list(target = "all", herm = "herm_omni_minP", comp = "ad_perm", label = "Omnibus: minP vs Perm-AD")
  )

  results <- lapply(pairs, function(p) {
    sub_pv <- if (p$target == "all") {
      pv %>% dplyr::filter(kind == "power")
    } else {
      pv %>% dplyr::filter(kind == "power", target == p$target)
    }

    rej_h <- sub_pv[[p$herm]] < alpha
    rej_c <- sub_pv[[p$comp]] < alpha

    wins_h <- sum(rej_h & !rej_c, na.rm = TRUE)
    wins_c <- sum(!rej_h & rej_c, na.rm = TRUE)
    ties   <- sum(rej_h == rej_c, na.rm = TRUE)
    disagreements <- wins_h + wins_c

    win_rate_cond <- if (disagreements > 0) wins_h / disagreements else 0.5
    win_rate_all  <- (wins_h + 0.5 * ties) / nrow(sub_pv)

    # Cell-level statistical significance
    pc <- paired_compare(sub_pv, p$herm, p$comp, alpha = alpha)
    cell_wins   <- sum(pc$sig == "+")
    cell_losses <- sum(pc$sig == "-")
    cell_ties   <- sum(pc$sig == ".")

    data.frame(
      Comparison       = p$label,
      Domain           = p$target,
      Hermite_Test     = p$herm,
      Comparator       = p$comp,
      Mean_Power_Diff  = mean(rej_h, na.rm = TRUE) - mean(rej_c, na.rm = TRUE),
      Rep_Win_Rate_Pct = win_rate_cond * 100,
      Cell_Wins        = cell_wins,
      Cell_Losses      = cell_losses,
      Cell_Ties        = cell_ties,
      Total_Cells      = nrow(pc),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(results)
}

summary_shape <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha)

  cat("\n=== 1. TYPE I ERROR, COMPLETE NULLS (F1 = F2) =============================\n")
  print(
    rt %>%
      dplyr::filter(kind == "null",
                    test %in% c("herm_scale", "herm_asym", "herm_tail", "herm_omni_minP",
                                "f_test", "bf", "fligner", "raw_g1_perm", "raw_g2_perm",
                                "ks_perm", "ad_perm")) %>%
      dplyr::select(scenario, n1, test, rate) %>%
      tidyr::pivot_wider(names_from = test, values_from = rate, id_cols = c(scenario, n1)) %>%
      as.data.frame(),
    digits = 3
  )

  cat("\n=== 2. ROBUSTNESS UNDER PARTIAL NULLS =====================================\n")
  cat("   (Valid tests must hold nominal alpha for their target contrast)\n")
  print(
    rt %>%
      dplyr::filter(kind == "pnull") %>%
      dplyr::select(scenario, target, n1, test, rate) %>%
      tidyr::pivot_wider(names_from = test, values_from = rate) %>%
      dplyr::select(scenario, target, n1, herm_scale, herm_asym, herm_tail,
                    herm_scale_ex, herm_asym_ex, herm_tail_ex,
                    f_test, bf, raw_g1_perm, raw_g2_perm) %>%
      as.data.frame(),
    digits = 3
  )

  cat("\n=== 3. NOMINAL POWER BY SCENARIO ==========================================\n")
  print(
    rt %>%
      dplyr::filter(kind == "power") %>%
      dplyr::select(scenario, target, n1, test, rate) %>%
      tidyr::pivot_wider(names_from = test, values_from = rate) %>%
      as.data.frame(),
    digits = 3
  )

  cat("\n=== 4. SIZE-ADJUSTED POWER ================================================\n")
  print(
    size_adjusted(pv, alpha) %>%
      tidyr::pivot_wider(names_from = test, values_from = power_adj) %>%
      as.data.frame(),
    digits = 3
  )

  cat("\n=== 5. WIN-RATE BENCHMARK SUMMARY (HERMITE vs COMPETITORS) =================\n")
  win_summary <- win_rate_summary(pv, alpha = alpha)
  print(win_summary %>%
          dplyr::select(Comparison, Mean_Power_Diff, Rep_Win_Rate_Pct, Cell_Wins, Cell_Losses, Total_Cells) %>%
          as.data.frame(),
        digits = 3)

  invisible(rt)
}

# ------------------------------------------------------------------------------
# 8. Visualizations
# ------------------------------------------------------------------------------
plot_shape_alpha <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha) %>% dplyr::filter(kind %in% c("null", "pnull"))
  sc <- unique(rt$scenario)
  op <- par(mfrow = c(2, ceiling(length(sc) / 2)), mar = c(4, 4, 3, 1), font.main = 1)
  on.exit(par(op))

  keys <- c("herm_scale", "herm_asym", "herm_tail", "herm_omni_minP", "f_test", "bf")
  cols <- c("blue3", "#d95f02", "#7570b3", "black", "grey50", "darkgreen")

  for (s in sc) {
    d <- rt %>% dplyr::filter(scenario == s)
    plot(NA, xlim = range(d$n1), ylim = c(0, 0.20), log = "x",
         xlab = "n per group", ylab = "Empirical Type I Error", main = s, bty = "l")
    abline(h = alpha, col = "red", lty = 2)
    abline(h = c(.025, .075), col = "grey60", lty = 3)
    for (i in seq_along(keys)) {
      dd <- d %>% dplyr::filter(test == keys[i]) %>% dplyr::arrange(n1)
      if (nrow(dd) > 0) lines(dd$n1, dd$rate, type = "b", pch = 19, col = cols[i], lwd = 2)
    }
    if (s == sc[1]) legend("topright", keys, col = cols, lwd = 2, bty = "n", cex = .7)
  }
}

plot_shape_power <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha) %>% dplyr::filter(kind == "power")
  sc <- unique(rt$scenario)
  op <- par(mfrow = c(3, ceiling(length(sc) / 3)), mar = c(4, 4, 3, 1), font.main = 1)
  on.exit(par(op))

  for (s in sc) {
    tgt  <- SCENARIOS[[s]]$target
    keys <- switch(tgt,
                   scale      = c("herm_scale", "f_test", "bf", "fligner", "ad_perm"),
                   asymmetry  = c("herm_asym", "raw_g1_perm", "herm_omni_minP", "ks_perm", "ad_perm"),
                   tailweight = c("herm_tail", "raw_g2_perm", "herm_omni_minP", "ks_perm", "ad_perm"))
    cols <- c("blue3", "grey40", "darkgreen", "darkcyan", "orchid")
    d    <- rt %>% dplyr::filter(scenario == s)
    plot(NA, xlim = range(d$n1), ylim = c(0, 1), log = "x",
         xlab = "n per group", ylab = "Power", main = s, bty = "l")
    for (i in seq_along(keys)) {
      dd <- d %>% dplyr::filter(test == keys[i]) %>% dplyr::arrange(n1)
      if (nrow(dd) > 0) lines(dd$n1, dd$rate, type = "b", pch = 19, col = cols[i], lwd = 2)
    }
    legend("bottomright", keys, col = cols, lwd = 2, bty = "n", cex = .65)
  }
}

# ==============================================================================
# sim_shape_helper.R -- Simulation infrastructure for shape_hermite (v2)
# ==============================================================================
#
# What is different from v1:
#   * per-replication p-values are stored -> paired comparisons + size-adjusted power
#   * one shared permutation stream per replication -> all permutation tests use
#     identical resamples (fair, paired, and ~4x faster)
#   * comparators: F-test, Brown-Forsythe, Fligner-Killeen (scale);
#     aligned raw-g1 / raw-g2 permutation (shape);
#     permutation-KS, permutation Anderson-Darling (omnibus)
#   * scenario set includes non-normal scale contrasts, partial nulls, and a
#     location-shift null (the realistic case v1 omitted)
#   * n = 10, 15 included; unequal-n cells included
#   * L'Ecuyer-CMRG streams
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel); library(dplyr); library(tidyr)
})
if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
library(pbapply)

safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }

# ------------------------------------------------------------------------------
# 1. Comparator statistics
# ------------------------------------------------------------------------------

g1_raw <- function(v) { vc <- v - mean(v); m2 <- mean(vc^2)
if (m2 <= 1e-12) 0 else mean(vc^3) / m2^1.5 }
g2_raw <- function(v) { vc <- v - mean(v); m2 <- mean(vc^2)
if (m2 <= 1e-12) 0 else mean(vc^4) / m2^2 - 3 }
logvar <- function(v) log(max(stats::var(v), 1e-12))

ks_stat <- function(a, b) {
  na <- length(a); nb <- length(b)
  o  <- order(c(a, b)); lab <- c(rep(1, na), rep(0, nb))[o]
  max(abs(cumsum(lab) / na - cumsum(1 - lab) / nb))
}

# Pettitt form of the two-sample Anderson-Darling statistic
ad_stat <- function(a, b) {
  na <- length(a); nb <- length(b); N <- na + nb
  o  <- order(c(a, b)); lab <- c(rep(1L, na), rep(0L, nb))[o]
  M  <- cumsum(lab)[seq_len(N - 1L)]; i <- seq_len(N - 1L)
  sum((M * N - na * i)^2 / (i * (N - i))) / (na * nb)
}

bf_p <- function(a, b) {
  z1 <- abs(a - stats::median(a)); z2 <- abs(b - stats::median(b))
  tryCatch(stats::t.test(z1, z2, var.equal = TRUE)$p.value, error = function(e) NA_real_)
}
fk_p <- function(a, b) {
  tryCatch(stats::fligner.test(list(a, b))$p.value, error = function(e) NA_real_)
}
f_p <- function(a, b) tryCatch(stats::var.test(a, b)$p.value, error = function(e) NA_real_)

# ------------------------------------------------------------------------------
# 2. Shared permutation engine for all resampling comparators
# ------------------------------------------------------------------------------
# Returns two-sided permutation p-values for a list of statistic functions,
# each with its own alignment mode:
#   "raw"  - permute pooled raw data (tests F1 = F2)
#   "loc"  - permute location-aligned residuals (tests scale)
#   "std"  - permute standardized residuals    (tests shape)
perm_bank <- function(x1, x2, perm_idx, stats_spec) {
  n1 <- length(x1); n2 <- length(x2)
  pools <- list(
    raw = c(x1, x2),
    loc = c(x1 - stats::median(x1), x2 - stats::median(x2)),
    std = c((x1 - stats::median(x1)) / stats::mad(x1),
            (x2 - stats::median(x2)) / stats::mad(x2))
  )
  out <- numeric(length(stats_spec)); names(out) <- names(stats_spec)
  B <- nrow(perm_idx)

  for (nm in names(stats_spec)) {
    sp   <- stats_spec[[nm]]
    pool <- pools[[sp$align]]
    if (any(!is.finite(pool))) { out[nm] <- NA_real_; next }
    obs  <- sp$fun(pool[seq_len(n1)], pool[n1 + seq_len(n2)])
    if (sp$align != "raw") obs <- sp$fun(x1, x2)   # invariant statistics: use raw
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
# 3. Scenarios
# ------------------------------------------------------------------------------
# Each entry: gen(n1, n2) -> list(x1, x2); ref = null scenario used for
# size-adjusted critical values; kind = "null" | "pnull" | "power"
# target = which contrast the alternative lives in ("-" for nulls)

std <- function(v, m, s) (v - m) / s

SCENARIOS <- list(
  # ---- complete nulls: F1 = F2 ------------------------------------------------
  null_norm  = list(kind="null", target="-", ref=NA,
                    gen=function(n1,n2) list(x1=rnorm(n1), x2=rnorm(n2))),
  null_lnorm = list(kind="null", target="-", ref=NA,
                    gen=function(n1,n2) list(x1=rlnorm(n1,0,.6), x2=rlnorm(n2,0,.6))),
  null_t3    = list(kind="null", target="-", ref=NA,
                    gen=function(n1,n2) list(x1=rt(n1,3), x2=rt(n2,3))),
  null_unif  = list(kind="null", target="-", ref=NA,
                    gen=function(n1,n2) list(x1=runif(n1), x2=runif(n2))),

  # ---- partial nulls: the tested contrast is null, others are not ------------
  pnull_loc  = list(kind="pnull", target="all", ref=NA,   # pure location shift
                    gen=function(n1,n2) list(x1=rnorm(n1), x2=rnorm(n2, mean=1))),
  pnull_loc_skew = list(kind="pnull", target="all", ref=NA,
                        gen=function(n1,n2) list(x1=rlnorm(n1,0,.6), x2=rlnorm(n2,0,.6)+2)),
  pnull_scale = list(kind="pnull", target="shape", ref=NA, # scale differs, shape null
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=rnorm(n2, sd=1.5))),
  pnull_shape = list(kind="pnull", target="scale", ref=NA, # shape differs, scale null
                     gen=function(n1,n2) list(x1=rnorm(n1),
                                              x2=std(rgamma(n2, shape=2), 2, sqrt(2)))),

  # ---- scale alternatives ----------------------------------------------------
  scale_norm  = list(kind="power", target="scale", ref="null_norm",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=rnorm(n2, sd=1.5))),
  scale_t3    = list(kind="power", target="scale", ref="null_t3",
                     gen=function(n1,n2) list(x1=rt(n1,3), x2=1.5*rt(n2,3))),
  scale_lnorm = list(kind="power", target="scale", ref="null_lnorm",
                     gen=function(n1,n2) list(x1=rlnorm(n1,0,.6), x2=1.5*rlnorm(n2,0,.6))),
  scale_contam= list(kind="power", target="scale", ref="null_norm",
                     gen=function(n1,n2){ f<-function(n){v<-rnorm(n); k<-runif(n)<.05; v[k]<-v[k]*4; v}
                     list(x1=f(n1), x2=1.5*f(n2))}),
  scale_shift = list(kind="power", target="scale", ref="pnull_loc",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=rnorm(n2, mean=0.8, sd=1.5))),

  # ---- asymmetry alternatives ------------------------------------------------
  asym_gamma  = list(kind="power", target="asymmetry", ref="null_norm",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=std(rgamma(n2,shape=2), 2, sqrt(2)))),
  asym_gamma8 = list(kind="power", target="asymmetry", ref="null_norm",   # milder
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=std(rgamma(n2,shape=8), 8, sqrt(8)))),
  asym_lnorm  = list(kind="power", target="asymmetry", ref="null_norm",
                     gen=function(n1,n2){ s<-.4; m<-exp(s^2/2); v<-sqrt((exp(s^2)-1)*exp(s^2))
                     list(x1=rnorm(n1), x2=std(rlnorm(n2,0,s), m, v))}),

  # ---- tail-weight alternatives ----------------------------------------------
  tail_t5     = list(kind="power", target="tailweight", ref="null_norm",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=rt(n2,5)*sqrt(3/5))),
  tail_t3     = list(kind="power", target="tailweight", ref="null_norm",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=rt(n2,3)/sqrt(3))),
  tail_unif   = list(kind="power", target="tailweight", ref="null_norm",
                     gen=function(n1,n2) list(x1=rnorm(n1), x2=std(runif(n2), .5, sqrt(1/12)))),
  tail_contam = list(kind="power", target="tailweight", ref="null_norm",
                     gen=function(n1,n2){ v<-rnorm(n2); k<-runif(n2)<.05; v[k]<-v[k]*4
                     list(x1=rnorm(n1), x2=std(v, mean(v), stats::sd(v)))})
)

generate_shape_data <- function(n1, n2, scenario) SCENARIOS[[scenario]]$gen(n1, n2)

# ------------------------------------------------------------------------------
# 4. Population projection of the Hermite indices (estimand transparency)
# ------------------------------------------------------------------------------
shape_population_truth <- function(scenarios = names(SCENARIOS), N = 2e5, seed = 99) {
  set.seed(seed)
  do.call(rbind, lapply(scenarios, function(s) {
    d  <- SCENARIOS[[s]]$gen(N, N)
    p1 <- hermiteStats:::.shape_vector(d$x1, 3L)
    p2 <- hermiteStats:::.shape_vector(d$x2, 3L)
    data.frame(scenario = s, kind = SCENARIOS[[s]]$kind, target = SCENARIOS[[s]]$target,
               d_scale = p2["log_scale"] - p1["log_scale"],
               d_asym  = p2["asymmetry"] - p1["asymmetry"],
               d_tail  = p2["tailweight"] - p1["tailweight"],
               true_g1_1 = g1_raw(d$x1), true_g1_2 = g1_raw(d$x2),
               true_g2_1 = g2_raw(d$x1), true_g2_2 = g2_raw(d$x2),
               row.names = NULL)
  }))
}

# ------------------------------------------------------------------------------
# 5. Worker: one condition cell -> per-replication p-value matrix
# ------------------------------------------------------------------------------
simulate_shape_cell <- function(cond_idx, conditions, n_reps, nperm, ...) {
  cond <- conditions[cond_idx, ]
  n1 <- cond$n1; n2 <- cond$n2; scenario <- as.character(cond$scenario)
  N  <- n1 + n2

  cols <- c("herm_scale", "herm_asym", "herm_tail",
            "herm_scale_wy", "herm_asym_wy", "herm_tail_wy",
            "herm_omni_minP", "herm_omni_maxT",
            "herm_scale_ex", "herm_asym_ex", "herm_tail_ex",   # align = "none"
            "f_test", "bf", "fligner",
            "raw_logvar_perm", "raw_g1_perm", "raw_g2_perm",
            "ks_perm", "ad_perm")
  P <- matrix(NA_real_, nrow = n_reps, ncol = length(cols), dimnames = list(NULL, cols))
  deg_ok <- logical(n_reps)

  spec <- list(
    raw_logvar_perm = list(fun = function(a,b) logvar(b) - logvar(a), align = "loc"),
    raw_g1_perm     = list(fun = function(a,b) g1_raw(b) - g1_raw(a), align = "std"),
    raw_g2_perm     = list(fun = function(a,b) g2_raw(b) - g2_raw(a), align = "std"),
    ks_perm         = list(fun = function(a,b) ks_stat(a,b),          align = "raw"),
    ad_perm         = list(fun = function(a,b) ad_stat(a,b),          align = "raw")
  )

  for (r in seq_len(n_reps)) {
    d <- generate_shape_data(n1, n2, scenario)
    x1 <- d$x1; x2 <- d$x2

    # one shared permutation stream for every resampling test
    pidx <- t(vapply(seq_len(nperm), function(i) sample.int(N, n1), integer(n1)))

    sh <- tryCatch(shape_hermite(x1, x2, degree = 3L, align = "moment",
                                 omnibus = "minP", nperm = nperm, perm_idx = pidx),
                   error = function(e) NULL)
    if (!is.null(sh)) {
      P[r, c("herm_scale","herm_asym","herm_tail")]          <- sh$p_value
      P[r, c("herm_scale_wy","herm_asym_wy","herm_tail_wy")] <- sh$p_adjusted
      P[r, "herm_omni_minP"] <- sh$omnibus$p_value
      deg_ok[r] <- TRUE
    }
    sh2 <- tryCatch(shape_hermite(x1, x2, degree = 3L, align = "moment",
                                  omnibus = "maxT", nperm = nperm, perm_idx = pidx),
                    error = function(e) NULL)
    if (!is.null(sh2)) P[r, "herm_omni_maxT"] <- sh2$omnibus$p_value

    shE <- tryCatch(shape_hermite(x1, x2, degree = 3L, align = "none",
                                  omnibus = "none", nperm = nperm, perm_idx = pidx),
                    error = function(e) NULL)
    if (!is.null(shE)) P[r, c("herm_scale_ex","herm_asym_ex","herm_tail_ex")] <- shE$p_value

    P[r, "f_test"]  <- f_p(x1, x2)
    P[r, "bf"]      <- bf_p(x1, x2)
    P[r, "fligner"] <- fk_p(x1, x2)
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
# 6. Runner (L'Ecuyer streams)
# ------------------------------------------------------------------------------
run_shape_sim <- function(conditions, n_reps = 2000L, nperm = 999L,
                          n_cores = 4L, seed = 20260201L) {
  t0 <- Sys.time()
  cat(sprintf("\n%d cells x %d replications (nperm = %d)\n\n",
              nrow(conditions), n_reps, nperm))

  pkg_root <- if (file.exists("DESCRIPTION")) normalizePath(".") else normalizePath("..")

  cl <- NULL
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores); on.exit(stopCluster(cl), add = TRUE)
    clusterCall(cl, function(p) { devtools::load_all(p, quiet = TRUE); NULL }, p = pkg_root)
    clusterExport(cl, c("SCENARIOS","generate_shape_data","perm_bank","g1_raw","g2_raw",
                        "logvar","ks_stat","ad_stat","bf_p","fk_p","f_p","std",
                        "safe_mean","simulate_shape_cell"), envir = .GlobalEnv)
    clusterSetRNGStream(cl, seed)
  } else {
    RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  }

  res <- pblapply(seq_len(nrow(conditions)), simulate_shape_cell,
                  conditions = conditions, n_reps = n_reps, nperm = nperm, cl = cl)
  out <- dplyr::bind_rows(res)
  cat("\n== finished in", format(Sys.time() - t0), "==\n\n")
  out
}

# ------------------------------------------------------------------------------
# 7. Analysis: nominal power, size-adjusted power, paired comparisons
# ------------------------------------------------------------------------------
TESTS <- c("herm_scale","herm_asym","herm_tail","herm_scale_wy","herm_asym_wy",
           "herm_tail_wy","herm_omni_minP","herm_omni_maxT","herm_scale_ex",
           "herm_asym_ex","herm_tail_ex","f_test","bf","fligner",
           "raw_logvar_perm","raw_g1_perm","raw_g2_perm","ks_perm","ad_perm")

rates_table <- function(pv, alpha = 0.05) {
  pv %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::group_by(scenario, kind, target, n1, n2, test) %>%
    dplyr::summarise(R = sum(is.finite(p)),
                     rate = safe_mean(p < alpha),
                     mcse = sqrt(rate * (1 - rate) / R), .groups = "drop")
}

#' Size-adjusted power using the matched null scenario's empirical 5% quantile
size_adjusted <- function(pv, alpha = 0.05) {
  crit <- pv %>% dplyr::filter(kind == "null") %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::group_by(ref_scen = scenario, n1, n2, test) %>%
    dplyr::summarise(crit = stats::quantile(p, alpha, na.rm = TRUE), .groups = "drop")

  pv %>% dplyr::filter(kind == "power", !is.na(ref)) %>%
    tidyr::pivot_longer(dplyr::all_of(TESTS), names_to = "test", values_to = "p") %>%
    dplyr::left_join(crit, by = c("ref" = "ref_scen", "n1", "n2", "test")) %>%
    dplyr::group_by(scenario, target, n1, n2, test) %>%
    dplyr::summarise(power_adj = safe_mean(p <= crit), .groups = "drop")
}

#' Paired head-to-head: mean difference in rejection with paired MCSE
paired_compare <- function(pv, a, b, alpha = 0.05) {
  pv %>% dplyr::filter(kind == "power") %>%
    dplyr::mutate(da = as.numeric(.data[[a]] < alpha) - as.numeric(.data[[b]] < alpha)) %>%
    dplyr::group_by(scenario, target, n1, n2) %>%
    dplyr::summarise(diff = mean(da, na.rm = TRUE),
                     se   = stats::sd(da, na.rm = TRUE) / sqrt(sum(is.finite(da))),
                     .groups = "drop") %>%
    dplyr::mutate(A = a, B = b,
                  sig = ifelse(abs(diff) > 1.96 * se, ifelse(diff > 0, "+", "-"), "."))
}

summary_shape <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha)
  ci <- function(r, R) sprintf("%.3f", r)

  cat("\n=== 1. TYPE I ERROR, complete nulls (F1 = F2) ==============================\n")
  print(rt %>% dplyr::filter(kind == "null",
                             test %in% c("herm_scale","herm_asym","herm_tail","herm_omni_minP",
                                         "f_test","bf","fligner","raw_g1_perm","raw_g2_perm",
                                         "ks_perm","ad_perm")) %>%
          dplyr::select(scenario, n1, test, rate, mcse) %>%
          tidyr::pivot_wider(names_from = test, values_from = rate,
                             id_cols = c(scenario, n1)) %>% as.data.frame(), digits = 3)

  cat("\n=== 2. VALIDITY UNDER PARTIAL NULLS =======================================\n")
  cat("   (target = which contrast should be at nominal alpha)\n")
  print(rt %>% dplyr::filter(kind == "pnull") %>%
          dplyr::select(scenario, target, n1, test, rate) %>%
          tidyr::pivot_wider(names_from = test, values_from = rate) %>%
          dplyr::select(scenario, target, n1, herm_scale, herm_asym, herm_tail,
                        herm_scale_ex, herm_asym_ex, herm_tail_ex,
                        f_test, bf, raw_g1_perm, raw_g2_perm) %>%
          as.data.frame(), digits = 3)

  cat("\n=== 3. NOMINAL POWER ======================================================\n")
  print(rt %>% dplyr::filter(kind == "power") %>%
          dplyr::select(scenario, target, n1, test, rate) %>%
          tidyr::pivot_wider(names_from = test, values_from = rate) %>%
          as.data.frame(), digits = 3)

  cat("\n=== 4. SIZE-ADJUSTED POWER ================================================\n")
  print(size_adjusted(pv, alpha) %>%
          tidyr::pivot_wider(names_from = test, values_from = power_adj) %>%
          as.data.frame(), digits = 3)

  cat("\n=== 5. PAIRED HEAD-TO-HEAD (diff in rejection, +/- = sig. at 5%) ==========\n")
  pairs <- list(c("herm_scale","f_test"), c("herm_scale","bf"), c("herm_scale","fligner"),
                c("herm_asym","raw_g1_perm"), c("herm_tail","raw_g2_perm"),
                c("herm_omni_minP","ks_perm"), c("herm_omni_minP","ad_perm"),
                c("herm_omni_minP","herm_omni_maxT"))
  hh <- dplyr::bind_rows(lapply(pairs, function(p) paired_compare(pv, p[1], p[2], alpha)))
  print(hh %>% dplyr::group_by(A, B, target) %>%
          dplyr::summarise(mean_diff = mean(diff),
                           se = sqrt(mean(se^2) / dplyr::n()),
                           n_cells = dplyr::n(),
                           wins = sum(sig == "+"), losses = sum(sig == "-"),
                           .groups = "drop") %>% as.data.frame(), digits = 3)
  invisible(rt)
}

# ------------------------------------------------------------------------------
# 8. Plots
# ------------------------------------------------------------------------------
plot_shape_alpha <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha) %>% dplyr::filter(kind %in% c("null","pnull"))
  sc <- unique(rt$scenario)
  op <- par(mfrow = c(2, ceiling(length(sc)/2)), mar = c(4,4,3,1), font.main = 1)
  on.exit(par(op))
  keys <- c("herm_scale","herm_asym","herm_tail","herm_omni_minP","f_test","bf")
  cols <- c("blue3","#d95f02","#7570b3","black","grey50","darkgreen")
  for (s in sc) {
    d <- rt %>% dplyr::filter(scenario == s)
    plot(NA, xlim = range(d$n1), ylim = c(0, 0.20), log = "x",
         xlab = "n per group", ylab = "empirical alpha", main = s, bty = "l")
    abline(h = alpha, col = "red", lty = 2); abline(h = c(.025,.075), col="grey60", lty=3)
    for (i in seq_along(keys)) {
      dd <- d %>% dplyr::filter(test == keys[i]) %>% dplyr::arrange(n1)
      lines(dd$n1, dd$rate, type = "b", pch = 19, col = cols[i], lwd = 2)
    }
    if (s == sc[1]) legend("topright", keys, col = cols, lwd = 2, bty = "n", cex = .7)
  }
}

plot_shape_power <- function(pv, alpha = 0.05) {
  rt <- rates_table(pv, alpha) %>% dplyr::filter(kind == "power")
  sc <- unique(rt$scenario)
  op <- par(mfrow = c(3, ceiling(length(sc)/3)), mar = c(4,4,3,1), font.main = 1)
  on.exit(par(op))
  for (s in sc) {
    tgt  <- SCENARIOS[[s]]$target
    keys <- switch(tgt,
                   scale      = c("herm_scale","f_test","bf","fligner","ad_perm"),
                   asymmetry  = c("herm_asym","raw_g1_perm","herm_omni_minP","ks_perm","ad_perm"),
                   tailweight = c("herm_tail","raw_g2_perm","herm_omni_minP","ks_perm","ad_perm"))
    cols <- c("blue3","grey40","darkgreen","darkcyan","orchid")
    d <- rt %>% dplyr::filter(scenario == s)
    plot(NA, xlim = range(d$n1), ylim = c(0,1), log = "x",
         xlab = "n per group", ylab = "power", main = s, bty = "l")
    for (i in seq_along(keys)) {
      dd <- d %>% dplyr::filter(test == keys[i]) %>% dplyr::arrange(n1)
      lines(dd$n1, dd$rate, type = "b", pch = 19, col = cols[i], lwd = 2)
    }
    legend("bottomright", keys, col = cols, lwd = 2, bty = "n", cex = .65)
  }
}

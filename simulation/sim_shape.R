# ==============================================================================
# sim_shape.R -- Benchmark Driver for shape_hermite
# ==============================================================================

source("simulation/sim_shape_helper.R")


# ---- 1. Simulation Parameters ------------------------------------------------
REPS_NULL  <- 2500L          # MCSE = 0.0044 at alpha = 0.05 (Tight 95% CI: [0.041, 0.059])
REPS_POWER <- 1000L          # Paired stream yields high power precision
NPERM      <- 999L           # Clean 0.001 resolution (Publication standard)
SEED       <- 20260201L
CORES      <- max(1L, parallel::detectCores() - 1L)

# ---- 2. Design Grid ----------------------------------------------------------
balanced   <- data.frame(n1 = c(10, 15, 20, 30, 50, 100))
balanced$n2 <- balanced$n1
unbalanced <- data.frame(n1 = c(10, 15, 25), n2 = c(30, 45, 75))
sizes      <- rbind(balanced, unbalanced)

grid_for <- function(kinds) {
  scen <- names(SCENARIOS)[vapply(SCENARIOS, function(s) s$kind %in% kinds, logical(1))]
  merge(sizes, data.frame(scenario = scen, stringsAsFactors = FALSE))
}

cond_null  <- grid_for(c("null", "pnull"))
cond_power <- grid_for("power")

cat(sprintf("Configured Cells -> Null / Partial-Null: %d | Power: %d | Total: %d\n",
            nrow(cond_null), nrow(cond_power), nrow(cond_null) + nrow(cond_power)))

# ---- 3. Population Estimand Projection ---------------------------------------
truth <- shape_population_truth(N = 2e5)
cat("\n--- Population Hermite Contrasts (Degree-3 Projection) ---\n")
print(truth, digits = 3, row.names = FALSE)
write.csv(truth, "sim_shape_population_truth.csv", row.names = FALSE)

# ---- 4. Execution ------------------------------------------------------------
pv_null  <- run_shape_sim(cond_null,  n_reps = REPS_NULL,  nperm = NPERM,
                          n_cores = CORES, seed = SEED)
pv_power <- run_shape_sim(cond_power, n_reps = REPS_POWER, nperm = NPERM,
                          n_cores = CORES, seed = SEED + 7L)

pv <- dplyr::bind_rows(pv_null, pv_power)
saveRDS(pv, "sim_shape_pvalues.rds")

# ---- 5. Analysis, Summaries & Win-Rate Output --------------------------------
rt <- summary_shape(pv, alpha = 0.05)
win_summary <- win_rate_summary(pv, alpha = 0.05)

write.csv(rt, "sim_shape_rates.csv", row.names = FALSE)
write.csv(size_adjusted(pv), "sim_shape_power_sizeadjusted.csv", row.names = FALSE)
write.csv(win_summary, "sim_shape_win_rates.csv", row.names = FALSE)

# ---- 6. Diagnostic Visualizations --------------------------------------------
pdf("sim_shape_alpha.pdf", width = 12, height = 7)
plot_shape_alpha(pv)
dev.off()

pdf("sim_shape_power.pdf", width = 13, height = 9)
plot_shape_power(pv)
dev.off()

# ---- 7. Fit Failure Checks ---------------------------------------------------
cat("\n--- Numerical Fit Failure Rates by Cell (Target: 0) ---\n")
fails <- pv %>%
  dplyr::group_by(scenario, n1, n2) %>%
  dplyr::summarise(fail_rate = 1 - mean(fit_ok), .groups = "drop") %>%
  dplyr::filter(fail_rate > 0)

if (nrow(fails) > 0) {
  print(as.data.frame(fails))
} else {
  cat("All quantile fits completed without degeneracy across all replications.\n")
}

cat("\nBenchmark run completed successfully.\n")

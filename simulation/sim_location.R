# ==============================================================================
# sim_location.R -- Driver Script for Benchmarking Mean & Median Tests
# ==============================================================================

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".")
  source("simulations/sim_location_helpers.R")
} else {
  devtools::load_all("..")
  source("sim_location_helpers.R")
}

# ------------------------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------------------------
repetitions <- 500L          # Replications per condition cell
n_perms     <- 499L          # Permutations for resampling tests (Analytical tests run instantly)
SEED        <- 20260301L

# ------------------------------------------------------------------------------
# 2. Design Grid
# ------------------------------------------------------------------------------
sample_sizes <- c(15, 25, 50, 100)
effect_sizes <- c(0.0, 0.3, 0.6)     # 0.0 = Type I error, 0.3 = small-med, 0.6 = medium-large
distributions <- c("norm", "lnorm", "t3", "exp", "outlier")

sim_conditions <- expand.grid(
  n_sample = sample_sizes,
  delta    = effect_sizes,
  dist     = distributions,
  stringsAsFactors = FALSE
)

cat(sprintf("Configured %d condition cells for Mean & Median testing.\n", nrow(sim_conditions)))

# ------------------------------------------------------------------------------
# 3. Execute Simulation
# ------------------------------------------------------------------------------
sim_results <- run_location_sim(
  conditions = sim_conditions,
  n_reps     = repetitions,
  nperm      = n_perms,
  n_cores    = max(1L, parallel::detectCores() - 1L),
  alpha      = 0.05,
  seed       = SEED
)

save(sim_results, file = "sim_location_results.RData")

# ------------------------------------------------------------------------------
# 4. Summary Tables & Visualizations
# ------------------------------------------------------------------------------
summary.sim_location(sim_results)

# Plot Power curves across effect sizes for skewed and heavy-tailed data
plot_location_comparison(sim_results, dist_focus = "lnorm")
plot_location_comparison(sim_results, dist_focus = "exp")
plot_location_comparison(sim_results, dist_focus = "t3")

write.csv(sim_results, file = "sim_location_aggregated_results.csv", row.names = FALSE)
cat("Simulation complete. Results exported to 'sim_location_aggregated_results.csv'.\n")

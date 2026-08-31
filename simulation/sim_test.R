# ==============================================================================
# sim_test.R -- Driver Script for Benchmarking t_hermite
# ==============================================================================

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".")
  source("simulations/sim_test_helper.R")
} else {
  devtools::load_all("..")
  source("sim_test_helper.R")
}

# ------------------------------------------------------------------------------
# 1. Simulation Parameters
# ------------------------------------------------------------------------------
repetitions <- 500L          # Replications per condition cell
n_perms     <- 1000L         # Permutations per test (Publication Grade)
SEED        <- 20260101L

# ------------------------------------------------------------------------------
# 2. Design Grid
# ------------------------------------------------------------------------------
sample_sizes <- c(10, 15, 25, 50, 100)
effect_sizes <- c(0.0, 0.2, 0.5, 0.8)

base_conditions <- expand.grid(
  n_sample = sample_sizes,
  delta    = effect_sizes
)

dist_params <- rbind(
  data.frame(dist = "norm",    param1 = 0,  param2 = 1,   param3 = NA),   # Normal
  data.frame(dist = "lnorm",   param1 = 0,  param2 = 0.8, param3 = NA),   # Right skew
  data.frame(dist = "t",       param1 = 3,  param2 = 1,   param3 = NA),   # Heavy tails / t3
  data.frame(dist = "exp",     param1 = 1,  param2 = NA,  param3 = NA),   # Exponential decay
  data.frame(dist = "outlier", param1 = 4,  param2 = NA,  param3 = 0.02), # 2% contamination
  data.frame(dist = "1plirt",  param1 = 0,  param2 = 1,   param3 = 20)    # Discrete IRT
)

sim_conditions <- merge(base_conditions, dist_params, all = TRUE)

# ------------------------------------------------------------------------------
# 3. Run Simulation
# ------------------------------------------------------------------------------
sim_results <- run_test_sim(
  conditions     = sim_conditions,
  n_reps         = repetitions,
  nperm          = n_perms,
  n_cores        = max(1, parallel::detectCores() - 1),
  alpha          = 0.05,
  seed           = SEED,
  truth_pop_size = 1e6
)

agg <- sim_results$aggregated
save(sim_results, agg, file = "sim_test_results.RData")

# ------------------------------------------------------------------------------
# 4. Summary & Visualizations
# ------------------------------------------------------------------------------
summary.sim_test(agg)

# Multi-panel Type I error across sample sizes
plot_type1_error(agg, by = "sample_size")

# Power curves across sample sizes
plot_power_curves(agg, dist_focus = "lnorm")
plot_power_curves(agg, dist_focus = "t")
plot_power_curves(agg, dist_focus = "exp")
plot_power_curves(agg, dist_focus = "norm")
plot_power_curves(agg, dist_focus = "1plirt")
plot_power_curves(agg, dist_focus = "outlier")

write.csv(agg, file = "sim_test_aggregated_results.csv", row.names = FALSE)

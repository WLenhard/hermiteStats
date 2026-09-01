# TODO
# 1. Rename Wilcoxon to Mann-Whitney-U test (independent sample)
# 2. Add Mood's Median test as a benchmark (-> median_hermite)
# 3. Draft dependent sample simulation, which draws an cor_hermite, Gaussian copula?
# 4. Probably integrate with the shape test simulation






# ==============================================================================
# sim_test.R -- Driver for Location Benchmarking (Mean & Median Domains)
# ==============================================================================

source("simulation/sim_test_helper.R")


# ------------------------------------------------------------------------------
# 1. Simulation Parameters
# ------------------------------------------------------------------------------
repetitions <- 1e3                # Monte Carlo replications per cell
n_perms     <- 1e3 - 1L           # Permutation replicates per test
population  <- 1e6                # Population size
SEED        <- 20260101L
CORES       <- max(1L, parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 2. Design Grid (Balanced + Behrens-Fisher Unbalanced Cells)
# ------------------------------------------------------------------------------

# A. Standard balanced sample sizes
balanced_sizes <- data.frame(
  n1 = c(15, 25, 50, 100),
  n2 = c(15, 25, 50, 100)
)

# B. Behrens-Fisher unbalanced sample sizes (1:3 ratio)
unbalanced_sizes <- data.frame(
  n1 = c(15, 20),
  n2 = c(45, 60)
)

all_sizes <- rbind(balanced_sizes, unbalanced_sizes)
effect_sizes <- c(0.0, 0.2, 0.5, 0.8)

# Core distribution specifications
dist_specs <- rbind(
  data.frame(dist = "norm",    param1 = 0,  param2 = 1,   param3 = NA, scale_ratio = 1.0), # Gaussian
  data.frame(dist = "lnorm",   param1 = 0,  param2 = 0.8, param3 = NA, scale_ratio = 1.0), # Skewed
  data.frame(dist = "t",       param1 = 3,  param2 = 1,   param3 = NA, scale_ratio = 1.0), # Heavy-tailed (t3)
  data.frame(dist = "exp",     param1 = 1,  param2 = NA,  param3 = NA, scale_ratio = 1.0), # Exponential
  data.frame(dist = "outlier", param1 = 4,  param2 = NA,  param3 = 0.02, scale_ratio = 1.0), # 2% contamination
  data.frame(dist = "1plirt",  param1 = 0,  param2 = 1,   param3 = 20, scale_ratio = 1.0), # Psychometric IRT

  # Behrens-Fisher stress conditions: variance ratio 1:4 (sigma2 / sigma1 = 2.0 or 0.5)
  data.frame(dist = "norm",    param1 = 0,  param2 = 1,   param3 = NA, scale_ratio = 2.0), # Direct heteroscedastic
  data.frame(dist = "norm",    param1 = 0,  param2 = 1,   param3 = NA, scale_ratio = 0.5), # Inverse heteroscedastic
  data.frame(dist = "lnorm",   param1 = 0,  param2 = 0.8, param3 = NA, scale_ratio = 1.5)  # Skewed + heteroscedastic
)

# Construct full factorial grid
sim_conditions <- merge(
  expand.grid(
    size_idx = seq_len(nrow(all_sizes)),
    delta    = effect_sizes,
    stringsAsFactors = FALSE
  ),
  cbind(size_idx = seq_len(nrow(all_sizes)), all_sizes),
  by = "size_idx"
)
sim_conditions$size_idx <- NULL

sim_conditions <- merge(sim_conditions, dist_specs, all = TRUE)

# Filter out redundant combinations (only run heteroscedastic scales on relevant size cells)
sim_conditions <- sim_conditions %>%
  dplyr::filter(
    scale_ratio == 1.0 | (n1 != n2) | (n1 == 25 & n2 == 25)
  )

cat(sprintf("Configured Location Simulation Grid: %d condition cells\n", nrow(sim_conditions)))

# ------------------------------------------------------------------------------
# 3. Execution
# ------------------------------------------------------------------------------
sim_output <- run_location_sim(
  conditions     = sim_conditions,
  n_reps         = repetitions,
  nperm          = n_perms,
  n_cores        = CORES,
  alpha          = 0.05,
  seed           = SEED,
  truth_pop_size = POPULATION
)

agg <- sim_output$aggregated

# ------------------------------------------------------------------------------
# 4. Summary Tables & Publication Figures
# ------------------------------------------------------------------------------
summary_location_sim(agg, alpha = 0.05)

# Visualizations
plot_type1_error_comparison(agg, nominal = 0.05)
plot_location_power_curves(agg, dist_focus = "norm")
plot_location_power_curves(agg, dist_focus = "lnorm")
plot_location_power_curves(agg, dist_focus = "t")
plot_location_power_curves(agg, dist_focus = "exp")
plot_location_power_curves(agg, dist_focus = "1plirt")
plot_location_power_curves(agg, dist_focus = "outlier")

# write.csv(agg, "sim_location_aggregated.csv", row.names = FALSE)
# saveRDS(sim_output, "sim_location_results.rds")

cat("\nLocation simulation pipeline completed successfully.\n")

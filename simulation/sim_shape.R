# ==============================================================================
# sim_shape.R -- Driver for the shape_hermite benchmark (v2)
# ==============================================================================

if (file.exists("DESCRIPTION")) {
  devtools::load_all("."); source("simulations/sim_shape_helper.R")
} else {
  devtools::load_all(".."); source("sim_shape_helper.R")
}

# ---- 1. parameters -----------------------------------------------------------
REPS_NULL  <- 1000L    # nulls need precision: MCSE ~ 0.003
REPS_POWER <- 1000L
NPERM      <- 499L
SEED       <- 20260201L
CORES      <- max(1L, parallel::detectCores() - 1L)

# ---- 2. design grid ----------------------------------------------------------
balanced <- data.frame(n1 = c(10, 15, 20, 30, 50, 100))
balanced$n2 <- balanced$n1
unbalanced <- data.frame(n1 = c(10, 15, 25), n2 = c(30, 45, 75))
sizes <- rbind(balanced, unbalanced)

grid_for <- function(kinds) {
  scen <- names(SCENARIOS)[vapply(SCENARIOS, function(s) s$kind %in% kinds, logical(1))]
  merge(sizes, data.frame(scenario = scen, stringsAsFactors = FALSE))
}

cond_null  <- grid_for(c("null", "pnull"))
cond_power <- grid_for("power")

cat(sprintf("Null/partial-null cells: %d   Power cells: %d\n",
            nrow(cond_null), nrow(cond_power)))

# ---- 3. estimand transparency -----------------------------------------------
truth <- shape_population_truth(N = 2e5)
cat("\n--- Population Hermite contrasts (degree-3 projection) ---\n")
print(truth, digits = 3, row.names = FALSE)
write.csv(truth, "sim_shape_population_truth.csv", row.names = FALSE)

# ---- 4. run ------------------------------------------------------------------
pv_null  <- run_shape_sim(cond_null,  n_reps = REPS_NULL,  nperm = NPERM,
                          n_cores = CORES, seed = SEED)
pv_power <- run_shape_sim(cond_power, n_reps = REPS_POWER, nperm = NPERM,
                          n_cores = CORES, seed = SEED + 7L)

pv <- dplyr::bind_rows(pv_null, pv_power)
saveRDS(pv, "sim_shape_pvalues.rds")

# ---- 5. analysis -------------------------------------------------------------
rt <- summary_shape(pv, alpha = 0.05)

write.csv(rt, "sim_shape_rates.csv", row.names = FALSE)
write.csv(size_adjusted(pv), "sim_shape_power_sizeadjusted.csv", row.names = FALSE)

# ---- 6. plots ----------------------------------------------------------------
pdf("sim_shape_alpha.pdf", width = 12, height = 7); plot_shape_alpha(pv);  dev.off()
pdf("sim_shape_power.pdf", width = 13, height = 9); plot_shape_power(pv); dev.off()

# ---- 7. diagnostics ----------------------------------------------------------
cat("\n--- Fit failure rate by cell (should be ~0) ---\n")
print(pv %>% dplyr::group_by(scenario, n1) %>%
        dplyr::summarise(fail = 1 - mean(fit_ok), .groups = "drop") %>%
        dplyr::filter(fail > 0) %>% as.data.frame())

cat("\nDone.\n")

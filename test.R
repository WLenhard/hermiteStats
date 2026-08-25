# Simulate non-normal data (Lognormal)
set.seed(42)
x <- rlnorm(30, meanlog = 0, sdlog = 0.8)
y <- rlnorm(30, meanlog = 0.5, sdlog = 0.8)

# 1. Hermite-Mehler Correlation
r_fit <- cor_hermite(x, y, conf_level = 0.95)
print(r_fit)
plot(r_fit)

# 2. Distribution-Free Effect Size (Independent)
d_fit <- d_reg(x, y, conf_level = 0.95, ci_method="nct")
print(d_fit)
plot(d_fit)

# 3. Distribution-Free Effect Size (Paired)
d_paired <- d_reg(x, y, paired = TRUE, conf_level = 0.95)
print(d_paired)
plot(d_paired)

# 4. Diagnostic Plots
plot(r_fit)
plot(d_fit)

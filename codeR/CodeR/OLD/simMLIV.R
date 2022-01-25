
# Setup ----------------------------------------------------------------------------------
#  - y = β0 + β1 x1 + β2 x2 + ε
#  - x1 and x2 are correlated
#  - x2 is unobserved (as is ε)
#  - z is a valid instrument for x1
#      - z affects x1
#      - z is uncorrelated with x2 and ε
#  - z and x2^2 are correlated

# Setup ----------------------------------------------------------------------------------
  # Options
  options(stringsAsFactors = F)
  # Packages
  library(pacman)
  p_load(
    ggplot2, ggthemes, latex2exp, Cairo,
    party, estimatr, lfe,
    tidyverse, haven, data.table, lubridate, magrittr, parallel
  )

# Function: One iteration ----------------------------------------------------------------
  one_iter <- function(iter, n) {
    # Set parameters
    β0 <- 1
    β1 <- 1
    β2 <- 1
    # Generate the instrument z
    i_dt <- data.table(
      z = runif(n = n, min = -1, max = 1),
      e_common = rnorm(n)
    )
    # z affects x1 (plus correlated noise)
    i_dt[, x1 := 1 + z + e_common]
    # z^2 and affects x2 (plus correlated noise)
    i_dt[, x2 := 1 + 10 * z^2 + e_common]
    # Calculate y
    i_dt[, y := β0 + β1 * x1 + β2 * x2 + runif(n, min = -1, max = 1)]
    # OLS first stage
    i_dt[, x1_hat_ols := predict(lm(x1 ~ z, data = i_dt))]
    # Machine learn the first stage
    i_dt[, x1_hat_ml := predict(ctree(x1 ~ z, data = i_dt))]
    bind_rows(
      # OLS: DGP
      lm(y ~ x1 + x2, data = i_dt) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "OLS: DGP"),
      # OLS: OVB
      lm(y ~ x1, data = i_dt) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "OLS: OVB"),
      # Second stage: OLS
      lm(y ~ x1_hat_ols, data = i_dt) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "Second stage: OLS"),
      # Second stage: ML
      lm(y ~ x1_hat_ml, data = i_dt) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "Second stage: ML")
    ) %>% mutate(iter = iter) %>% data.table()
  }

# Function: Run iterations ---------------------------------------------------------------
  run_sim <- function(n, n_sims, n_cores = 4, seed = 12345) {
    # Set the seed
    set.seed(seed)
    # Run a parallelized simulation
    mclapply(X = 1:n_sims, FUN = one_iter, n = n, mc.cores = n_cores) %>% rbindlist()
  }

# Run simulation -------------------------------------------------------------------------
  sim_dt <- run_sim(n = 1e4, n_sims = 1e4)

# Plot simulation ------------------------------------------------------------------------
  ggplot(
    data = sim_dt %>% filter(grepl("Second", model)),
    aes(x = estimate, fill = model)
  ) +
  geom_density(color = NA, alpha = 0.65) +
  geom_vline(xintercept = 1) +
  theme_pander() +
  scale_fill_viridis_d("Model", begin = 0.1, end = 0.85, option = "B")

sim_dt[, mean(estimate, na.rm = T), by = model]

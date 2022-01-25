
# Setup ----------------------------------------------------------------------------------
#  - y = β0 + β1 x1 + β2 x2 + ε 
#  - x1 and x2 are correlated
#  - x2 is unobserved (as is ε)
#  - z is a valid instrument for x1
#      - z affects x1
#      - z is uncorrelated with x2 and ε
#  - z and x2 are correlated

# Setup ----------------------------------------------------------------------------------
  # Options
  options(stringsAsFactors = F)
  # Packages
  library(pacman)
  p_load(
    ggplot2, ggthemes, latex2exp, Cairo,
    party, estimatr, lfe, randomForest,
    tidyverse, haven, data.table, lubridate, magrittr, parallel
  )

# Function: One iteration ----------------------------------------------------------------
  one_iter <- function(iter, n, ntree) {
    # Set parameters
     b0=1;  b1=1;  b2=1;  b3=1;  b4=1; b5=1; 
    b11=0; b22=0; b33=0; b44=0; b55=0; 
    b11=0; b13=0; b14=0; b15=0; 
    b23=0; b24=0; b25=0; 
    b34=0; b35=0; 
    b45=0; 
    # Generate the instrument z
    i_dt <- data.table(
      z = runif(n = n, min = -1, max = 1),
      e_common = rnorm(n, 0, 1)
    )
    # z affects x1 (plus correlated noise)
    i_dt[, x1 := 1 + z + z^2 + e_common]
    # z^2 and affects x2 (plus correlated noise)
    i_dt[, x2 := 1 + 1 * z^2 + e_common]
    # other covariates
    i_dt[, x3 := rnorm(n, 0, 1)]
    i_dt[, x4 := rnorm(n, 0, 2)]
    i_dt[, x5 := rnorm(n, 0, 3)]
    # Calculate y
    i_dt[, y := b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 +
           b11*x1*x1 + b12*x1*x2 + b13*x1*x3 + b14*x1*x4 + b15*x1*x5 +
           b22*x2*x2 + b23*x2*x3 + b24*x2*x4 + b25*x2*x5 +
           b33*x3*x3 + b34*x3*x4 + b35*x3*x5 + 
           b44*x4*x4 + b45*x4*x5 + 
           rnorm(n, 0, 1)] # + runif(n, min = -1, max = 1)]
    # OLS first stage
    i_dt[, x1_hat_ols := predict(lm(x1 ~ z, data = i_dt))]
    # Machine learn the first stage
    i_dt[, x1_hat_ml := predict(ctree(x1 ~ z, data = i_dt))]
    # Machine learn the first stage
    i_dt[, x1_hat_rf := predict(randomForest(x1 ~ z, data = i_dt, ntree=ntree))]
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
        mutate(model = "Second stage: ML"),
      # Second stage: RF
      lm(y ~ x1_hat_rf, data = i_dt) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "Second stage: RF")
    ) %>% mutate(iter = iter) %>% data.table()
  }

  
  
# Function: Run iterations ---------------------------------------------------------------
  run_sim <- function(n, n_sims, ntree, n_cores = 4, seed = 12345) {
    # Set the seed
    set.seed(seed)
    # Run a parallelized simulation
    mclapply(X = 1:n_sims, FUN = one_iter, n = n, ntree = ntree, mc.cores = n_cores) %>% rbindlist()
  }

  
  
# Run simulation -------------------------------------------------------------------------
  sim_dt <- run_sim(n = 1e4, n_sims = 1e2, ntree = 10)


  
  
# Plot simulation ------------------------------------------------------------------------

df <- merge(
  sim_dt %>% filter(grepl("Second stage: ML", model)),
  sim_dt %>% filter(grepl("Second stage: OLS", model)),
  by="iter")
df <- merge(
  df,
  sim_dt %>% filter(grepl("Second stage: RF", model)),
  by="iter")
df$MLgeOLS <- ifelse( df$estimate.x > df$estimate.y, 1, 0)
df$MLlessOLS <- df$estimate.x - df$estimate.y
setDT(df)
df[, mean(MLgeOLS, na.rm = T)]


  ggplot(
    data = sim_dt %>% filter(grepl("Second", model)),
    aes(x = estimate, fill = model)
  ) +
  geom_density(color = NA, alpha = 0.65) +
  geom_vline(xintercept = 1) +
  theme_pander() +
  theme(plot.caption = element_text(face = "italic", size=6)) +
  scale_fill_viridis_d("Model", begin = 0.1, end = 0.85, option = "B") +
  labs(caption=paste("Fraction in which ML produces the larger estimate: ", 
                  df[, mean(MLgeOLS, na.rm = T)],
                  sep=""), size=5)

  sim_dt[, mean(estimate, na.rm = T), by = model] %>% filter(!grepl("SUMMARY", model))



  
n=1000
# Set parameters
    b0 <- 1
    b1 <- 1
    b2 <- 1
    # Generate the instrument z
    i_dt <- data.table(
      z = runif(n = n, min = -1, max = 1),
      e_common = rnorm(n, 0, 1)
    )
    # z affects x1 (plus correlated noise)
    i_dt[, x1 := 1 + z + e_common]
    # z^2 and affects x2 (plus correlated noise)
    i_dt[, x2 := 1 + 10 * z^2 + e_common]
    # Calculate y
    i_dt[, y := b0 + b1 * x1 + b2 * x2 + runif(n, min = -1, max = 1)]
#    i_dt[, y := β0 + β1 * x1 + β2 * x2 + rnorm(n, 0, 1)]
    # OLS first stage
    i_dt[, x1_hat_ols := predict(lm(x1 ~ z, data = i_dt))]
    # Machine learn the first stage
    i_dt[, x1_hat_ml := predict(ctree(x1 ~ z, data = i_dt))]
    # Machine learn the first stage
    i_dt[, x1_hat_rf := predict(randomforest(x1 ~ z, data = i_dt, ntree=500))]
    
  

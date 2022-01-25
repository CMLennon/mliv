

# Notes ----------------------------------------------------------------------------------


# Setup ----------------------------------------------------------------------------------
  # Load packages
  library(pacman)
  p_load(
    tidyverse, scico, patchwork, data.table, fst,
    sf, geodist, 
    tidymodels, fixest, 
    parallel, magrittr, here
  )


# Generate data --------------------------------------------------------------------------
  # Bounds of the city
  xm = 10
  ym = 10
  # Generate obsevations
  gen_dt = CJ(
    x = seq(-xm, xm, by = 0.1),
    y = seq(-ym, ym, by = 0.1)
  )
  # Add an ID
  gen_dt[, id := 1:.N]
  # Set column order
  setcolorder(gen_dt, c("id"))
  # Calculate distance to college (located at 0,0)
  gen_dt[, dist := sqrt(x^2 + y^2)]
  # Generate education
  # - River runs east-west at y = 0
  river = 0
  # - North of river: 1/(2 * dist^0.3)
  # - South of river: 0.2 - 0.1 * dist
  gen_dt[y > river, educ := 1.2/(1.5 * dist^0.3)]
  gen_dt[y <= river, educ := 1.2/(1.5 * dist^0.3) - 0.3]
  # Cap Education at 1
  gen_dt[educ > 1, educ := 1]
  # Education at univeristy is '1'
  gen_dt[x == 0 & y == 0, educ := 1]
  # Generate wealth (confounder with neighborhoods)
  # - North of river, default: low wealth (0)
  # - North of river, near univ ([-3,3] x [-3x3]): high wealth (1)
  # - South of river: wealth increases linearly from 0
  gen_dt[y > river, wealth := 0]
  gen_dt[y > river & between(x, -3, 3) & between(y, -3, 3), wealth := 1]
  gen_dt[y <= river, wealth := 0 +  0.059 * dist]
  # Generate income = 30 + 20 * education + 10 * wealth + N(0,5)
  gen_dt[, income := 30 + 20 * educ + 10 * wealth + rnorm(.N, sd = 5)]


# Empirics: OLS predictions --------------------------------------------------------------
  # First-stage predictions
  gen_dt[, educ_pred_ols := feols(educ ~ dist, gen_dt)$fitted.values]


# Empirics: Random forest ----------------------------------------------------------------
# Goal: Train random forest to predict education as a function of distance to uni
  # Define the recipe
  the_recipe = recipe(
    educ ~ dist,
    data = gen_dt
  )
  # Define the hyperparameter grid
# NOTE: I've tried min_n down to 1 and trees from 100 to 5,000.
  rf_grid = CJ(
    mtry = 1,
    min_n = seq(1e3, 1e4, by = 1e3),
    trees = c(1000)
  )
  # Write a random forest model function (so we can try with OOB and in parallel)
  rf_i = function(i) {
    # Define the random forest model for model i
    model_i = rand_forest(
      mode = "regression",
      mtry = rf_grid$mtry[i],
      trees = rf_grid$trees[i],
      min_n = rf_grid$min_n[i]
    ) %>% set_engine(
      engine = "ranger"
    )
    # Define workflow
    wf_i = workflow() %>% add_model(model_i) %>% add_recipe(the_recipe)
    # Fit the workflow
    fit_i = wf_i %>% fit(gen_dt)
    # Return DF w/ OOB error and the hyperparameters
    data.table(
      mtry = rf_grid$mtry[i],
      min_n = rf_grid$min_n[i],
      trees = rf_grid$trees[i],
      # Note: OOB error is a bit buried
      error_oob = fit_i$fit$fit$fit$prediction.error,
      r2_oob = fit_i$fit$fit$fit$r.squared
    )  
  }
  # Train the forests in parallel
  rf_dt = mclapply(
    X = 1:nrow(rf_grid),
    FUN = rf_i,
    mc.cores = min(c(10, nrow(rf_grid)))
  ) %>% rbindlist(use.names = T, fill = T)
  # Arrange by OOB error
  setorder(rf_dt, error_oob)
 # Split the full dataset into two halves
  set.seed(12345)
  n = gen_dt[,.N]
  split_assign = sample(x = 1:n, size = n, replace = F)
  split1 = split_assign[1:floor(n/2)] %>% gen_dt[.,]
  split2 = split_assign[(floor(n/2)+1):n] %>% gen_dt[.,]
  # Define chosen model
  chosen_model = rand_forest(
    "regression",
    mtry = rf_dt$mtry[1],
    trees = rf_dt$trees[1],
    min_n = rf_dt$min_n[1]
  ) %>% set_engine(engine = "ranger")
  # Fit 'chosen' models on split 1
  fit1 = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      educ ~ dist,
      data = split1
    )
  ) %>% fit(split1)
  # Fit 'chosen' models on split 2
  fit2 = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      educ ~ dist,
      data = split2
    )
  ) %>% fit(split2)
  # Use fitted models to predict onto the other splits
  split1[, educ_pred_rfss := predict(fit2, new_data = split1)$.pred]
  split2[, educ_pred_rfss := predict(fit1, new_data = split2)$.pred]
  # Combine splits
  gen_dt = list(split1, split2) %>% rbindlist(use.names = T, fill = T)
  setorder(gen_dt, id)
  # Regress share_degree on educ_pred_rf AND tract fixed effects
  gen_dt[, `:=`(
    educ_pred_rfss_linearized = feols(
      educ ~ educ_pred_rfss,
      data = gen_dt
    )$fitted.values
  )]
  # Fit chosen model on full dataset
  fit_full = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      educ ~ dist,
      data = gen_dt
    )
  ) %>% fit(gen_dt)
  # Full (non-split-sample) predictions
  gen_dt[, educ_pred_rf := predict(fit_full, new_data = gen_dt)]


# Save results ---------------------------------------------------------------------------
  # Save 'gen_dt'
  write_fst(
    x = gen_dt,
    path = here("data-spatial-sim", "results.fst"),
    compress = 100
  )


# Load results ---------------------------------------------------------------------------
  # Load the dataset
  gen_dt = here("data-spatial-sim", "results.fst") %>% read_fst(as.data.table = T)


# Results: OLS and 2SLS ------------------------------------------------------------------
  # Linear exclusion restriction: Fine
  reg1 = feols(wealth ~ dist, data = gen_dt) %>% summary(se = "hetero")
  # Endogeneity concern: Wealth and education are correlated (negatively)
  reg2 = feols(educ ~ wealth, data = gen_dt) %>% summary(se = "hetero")
  # Naïve OLS: Biased (downward)
  reg3 = feols(income ~ educ, data = gen_dt) %>% summary(se = "hetero")
  # First stage: Strong
  reg4 = feols(income ~ 1 | 1 | educ ~ dist, data = gen_dt) %>% summary(se = "hetero", stage = 1)
  # 2SLS: Near to truth
  reg5 = feols(income ~ 1 | 1 | educ ~ dist, data = gen_dt) %>% summary(se = "hetero", stage = 2)


# Results: RF results --------------------------------------------------------------------
  # 2SLS with RF in the first stage (not MLSS)
  reg6 = feols(income ~ educ_pred_rf, data = gen_dt) %>% 
    summary(se = "hetero")
  # MLSS: Second stage
  reg7 = feols(income ~ 1 | 1 | educ ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero", stage = 2)
  # MLSS: First stage
  reg8 = feols(income ~ 1 | 1 | educ ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero", stage = 1)
  # MLSS: Reduced form
  reg9 = feols(income ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero")


# Make table(s) --------------------------------------------------------------------------
  # Dictionary
  setFixest_dict(c(
    income = "Income",
    educ = "Education",
    wealth = "Wealth",
    dist = "Distance to univ.",
    educ_pred_rf = "RF-pred. education",
    educ_pred_rfss = "SS RF-pred. education",
    educ_pred_rfss_linearized = "Linearized SS RF-pred. education"
  ))
  # Least-squares regressions
  etable(
    "Linear excl. rest." = reg1,
    "Endog. concern" = reg2,
    "Naïve OLS" = reg3,
    "2SLS, stage 1" = reg4,
    "2SLS, stage 2" = reg5,
    fitstat = "n"
  )
  # 2SLS and MLIV results
  etable(
    "2SLS" = reg5,
    "2SLS w/ RF" = reg6,
    "MLSS w/ RF" = reg7,
    "MLSS w/ RF" = reg8,
    subtitles = c("Stage 2", "Stage 2", "Stage 2", "Stage 1"),
    fitstat = "n"
  )

# Figures: Spatial figures: Main variables -----------------------------------------------
  # Spatial plot of education
  spatial_educ = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ)) +
  scale_fill_viridis_c("Education", option = "magma") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of wealth
  spatial_wealth = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = wealth)) +
  scale_fill_viridis_c("Wealth", option = "magma") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of income
  spatial_income = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = income)) +
  scale_fill_viridis_c("Income", option = "magma") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Combine spatial plots
  spatial_educ + spatial_wealth + spatial_income


# Figures: Spatial figures: Predicted education ------------------------------------------
  # Spatial plot of wealth
  spatial_wealth = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = wealth)) +
  scale_fill_scico("Wealth (endog.)", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of education
  spatial_educ = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ)) +
  scale_fill_scico("Education", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of OLS predictions
  spatial_ols = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ_pred_ols)) +
  scale_fill_scico("Pred. educ., OLS", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RF 
  spatial_rf = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ_pred_rf)) +
  scale_fill_scico("Pred. educ., RF", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RFSS
  spatial_rfss = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ_pred_rfss)) +
  scale_fill_scico("Pred. educ., RFSS", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RFSS, linearized
  spatial_rfssl = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_raster(aes(fill = educ_pred_rfss_linearized)) +
  scale_fill_scico("Pred. educ., RFSS line.", palette = "roma", limits = c(0,1), direction = -1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()
  # Combine spatial plots
  (spatial_wealth + spatial_educ + spatial_ols) / (spatial_rf + spatial_rfss + spatial_rfssl)


# Figures: Variables as functions of income ----------------------------------------------
  # Education as a function of distance to uni
  ggplot(data = gen_dt, aes(x = dist, y = educ)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Education level") +
  theme_minimal()
  # Wealth as a function of distance to uni
  ggplot(data = gen_dt, aes(x = dist, y = wealth)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Wealth") +
  theme_minimal()
  # Income as a function of distance to uni
  ggplot(data = gen_dt, aes(x = dist, y = income)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Income") +
  theme_minimal()


# Figure: Education and predictions as a function of distance ----------------------------
  # (Predicted) Education as a function of distance
  ggplot(data = gen_dt, aes(x = dist)) +
  geom_point(aes(color = letters[1], y = educ), size = 1/4) +
  geom_point(aes(color = letters[2], y = educ_pred_ols), size = 1/4) +
  geom_point(aes(color = letters[3], y = educ_pred_rf), size = 1/4) +
  geom_point(aes(color = letters[4], y = educ_pred_rfss), size = 1/4) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Eduation level") +
  scale_color_manual(
    "Education:",
    labels = c("Actual", "OLS pred.", "RF pred.", "SS RF pred."),
    values = c("black", viridis::magma(3, begin = 0.5, end = 0.92))
  ) +
  theme_minimal()

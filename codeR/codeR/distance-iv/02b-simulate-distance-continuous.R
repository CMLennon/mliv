

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
  # Set a seed
  set.seed(1)
  # Bounds of the city
  dist = seq(from = 0, to = 10, by = 0.1)
  deg = seq(from = 0, to = 359, by = 1)
  # Generate obsevations
  gen_dt = CJ(dist, deg)
  # Add an ID
  gen_dt[, id := 1:.N]
  # Set column order
  setcolorder(gen_dt, c("id"))
  # Generate education
  gen_dt[dist < quantile(dist, probs = 0.5), educ := 1.2/(1.5 * dist^0.3) - 0.3]
  gen_dt[dist >= quantile(dist, probs = 0.5), educ := 0.05]
  # Education at univeristy is '1'
  gen_dt[dist == 0, educ := 1]
  # Generate wealth (nonlinear confounder)
  wb = gen_dt[, quantile(dist, probs = c(1/3, 2/3))]
  gen_dt[dist < wb[1], wealth := 0.75]
  gen_dt[between(dist, wb[1], wb[2]), wealth := 0.25]
  gen_dt[dist > wb[2], wealth := 0.75]
  # Add noise to education and distance
  gen_dt[, `:=`(
    educ = educ + rnorm(.N, mean = 0, sd = 0.1),
    wealth = wealth + rnorm(.N, mean = 0, sd = 0.1)
  )]
  # Bound education and wealth between 0 and 1
  gen_dt[educ > 1, educ := 1]
  gen_dt[wealth > 1, wealth := 1]
  gen_dt[educ < 0, educ := 0]
  gen_dt[wealth < 0, wealth := 0]
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
  rf_grid = CJ(
    mtry = 1,
    min_n = c(50, 100, 1000),
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
  # Regress education on SS RF predictions
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


# Empirics: XGBoost ----------------------------------------------------------------------
# Goal: Train XGBoost to predict education as a function of distance to uni
  # Define the recipe
  the_recipe = recipe(
    educ ~ dist,
    data = gen_dt
  )
  # Define the hyperparameter grid
  xg_grid = CJ(
    tree_depth = 1:5,
    learn_rate = 10^(-4:-1),
    trees = c(1000),
    sample_size = seq(0.6, 1, 0.1),
    stop_iter = 5
  )
  # Define CV splits
  xg_cv = gen_dt %>% vfold_cv(v = 5)
  # Write a random forest model function (so we can try with OOB and in parallel)
  xg_i = function(i) {
    # Define the random forest model for model i
    model_i = boost_tree(
      mode = "regression",
      tree_depth = xg_grid$tree_depth[i],
      learn_rate = xg_grid$learn_rate[i],
      trees = xg_grid$trees[i],
      sample_size = xg_grid$sample_size[i],
      stop_iter = xg_grid$stop_iter[i]
    ) %>% set_engine(
      engine = "xgboost"
    )
    # Define workflow
    wf_i = workflow() %>% add_model(model_i) %>% add_recipe(the_recipe)
    # Fit the workflow
    fit_i = wf_i %>% fit_resamples(xg_cv)
    # Return DF w/ OOB error and the hyperparameters
    data.table(
      tree_depth = xg_grid$tree_depth[i],
      learn_rate = xg_grid$learn_rate[i],
      trees = xg_grid$trees[i],
      sample_size = xg_grid$sample_size[i],
      stop_iter = xg_grid$stop_iter[i],
      rmse = collect_metrics(fit_i)$mean[1],
      rsq = collect_metrics(fit_i)$mean[2]
    )  
  }
  # Train the xgboost model in parallel
  xg_dt = mclapply(
    X = 1:nrow(xg_grid),
    FUN = xg_i,
    mc.cores = min(c(10, nrow(xg_grid)))
  ) %>% rbindlist(use.names = T, fill = T)
  # Arrange by OOB error
  setorder(xg_dt, rmse)
  # Split the full dataset into two halves
  set.seed(12345)
  n = gen_dt[,.N]
  split_assign = sample(x = 1:n, size = n, replace = F)
  split1 = split_assign[1:floor(n/2)] %>% gen_dt[.,]
  split2 = split_assign[(floor(n/2)+1):n] %>% gen_dt[.,]
  # Define chosen model
  chosen_model = boost_tree(
    "regression",
    tree_depth = xg_dt$tree_depth[1],
    learn_rate = xg_dt$learn_rate[1],
    trees = xg_dt$trees[1],
    sample_size = xg_dt$sample_size[1],
    stop_iter = xg_dt$stop_iter[1]
  ) %>% set_engine(engine = "xgboost")
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
  split1[, educ_pred_xgss := predict(fit2, new_data = split1)$.pred]
  split2[, educ_pred_xgss := predict(fit1, new_data = split2)$.pred]
  # Combine splits
  gen_dt = list(split1, split2) %>% rbindlist(use.names = T, fill = T)
  setorder(gen_dt, id)
  # Regress education on SS RF predictions
  gen_dt[, `:=`(
    educ_pred_xgss_linearized = feols(
      educ ~ educ_pred_xgss,
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
  gen_dt[, educ_pred_xg := predict(fit_full, new_data = gen_dt)]


# Save results ---------------------------------------------------------------------------
  # Save 'gen_dt'
  write_fst(
    x = gen_dt,
    path = here("data-spatial-sim", "simple-results.fst"),
    compress = 100
  )


# Load results ---------------------------------------------------------------------------
  # Load the dataset
  gen_dt = here("data-spatial-sim", "simple-results.fst") %>% read_fst(as.data.table = T)
  # Convert polar coords to cartesian (for plotting)
  gen_dt[, `:=`(
    x = dist * cospi(deg/180),
    y = dist * sinpi(deg/180)
  )]


# Results: OLS and 2SLS ------------------------------------------------------------------
  # Linear exclusion restriction: Fine
  reg1 = feols(wealth ~ dist, data = gen_dt) %>% summary(se = "hetero")
  # Endogeneity concern: Wealth and education are correlated (negatively)
  reg2 = feols(educ ~ wealth, data = gen_dt) %>% summary(se = "hetero")
  # Naïve OLS: Biased
  reg3 = feols(income ~ educ, data = gen_dt) %>% summary(se = "hetero")
  # First stage: Strong
  reg4 = feols(income ~ 1 | 1 | educ ~ dist, data = gen_dt) %>% summary(se = "hetero", stage = 1)
  # 2SLS: Near to truth
  reg5 = feols(income ~ 1 | 1 | educ ~ dist, data = gen_dt) %>% summary(se = "hetero", stage = 2)


# Results: RF results --------------------------------------------------------------------
  # 2SLS with RF in the first stage (not MLSS)
  reg6 = feols(income ~ educ_pred_rf, data = gen_dt) %>% 
    summary(se = "hetero")
  # MLSS, RF: Second stage
  reg7 = feols(income ~ 1 | 1 | educ ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero", stage = 2)
  # MLSS, RF: First stage
  reg8 = feols(income ~ 1 | 1 | educ ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero", stage = 1)
  # MLSS, RF: Reduced form
  reg9 = feols(income ~ educ_pred_rfss, data = gen_dt) %>%
    summary(se = "hetero")
  # 2SLS with XGBoost in the first stage (not MLSS)
  reg10 = feols(income ~ educ_pred_xg, data = gen_dt) %>% 
    summary(se = "hetero")
  # MLSS, XGB: Second stage
  reg11 = feols(income ~ 1 | 1 | educ ~ educ_pred_xgss, data = gen_dt) %>%
    summary(se = "hetero", stage = 2)
  # MLSS, XGB: First stage
  reg12 = feols(income ~ 1 | 1 | educ ~ educ_pred_xgss, data = gen_dt) %>%
    summary(se = "hetero", stage = 1)
  # MLSS, XGB: Reduced form
  reg13 = feols(income ~ educ_pred_xgss, data = gen_dt) %>%
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
    educ_pred_rfss_linearized = "Linearized SS RF-pred. education",
    educ_pred_xg = "XGB-pred. education",
    educ_pred_xgss = "SS XGB-pred. education",
    educ_pred_xgss_linearized = "Linearized SS XGB-pred. education"
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
  # 2SLS and MLIV results: Random forest
  etable(
    "2SLS w/ RF" = reg6,
    "MLSS w/ RF" = reg7,
    "MLSS w/ RF" = reg8,
    "MLSS w/ RF" = reg9,
    subtitles = c("Stage 2", "Stage 2", "Stage 1", "Reduced form"),
    fitstat = "n"
  )
  # 2SLS and MLIV results: XGBoost
  etable(
    "2SLS w/ RF" = reg10,
    "MLSS w/ RF" = reg11,
    "MLSS w/ RF" = reg12,
    "MLSS w/ RF" = reg13,
    subtitles = c("Stage 2", "Stage 2", "Stage 1", "Reduced form"),
    fitstat = "n"
  )


# Figures: Spatial figures: Main variables -----------------------------------------------
  # Spatial plot of education
  spatial_educ = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ), size = 0.5) +
  scale_color_viridis_c("Education", option = "magma") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "left") +
  coord_equal()
  # Spatial plot of wealth
  spatial_wealth = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = wealth), size = 0.5) +
  scale_color_viridis_c("Wealth", option = "magma") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "left") +
  coord_equal()
  # Spatial plot of income
  spatial_income = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = income), size = 0.5) +
  scale_color_viridis_c("Income", option = "magma") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "left") +
  coord_equal()
  # Combine spatial plots
  spatial_educ / spatial_wealth / spatial_income


# Figures: Variables as functions of distances -------------------------------------------
  # Education as a function of distance to uni
  dist_educ = ggplot(data = gen_dt, aes(x = dist, y = educ)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Education level") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14)
  # Wealth as a function of distance to uni
  dist_wealth = ggplot(data = gen_dt, aes(x = dist, y = wealth)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Wealth") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14)
  # Income as a function of distance to uni
  dist_income = ggplot(data = gen_dt, aes(x = dist, y = income)) +
  geom_point(size = 1/4) +
  geom_smooth(method = lm, se = F) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Income") +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14)
  dist_educ / dist_wealth / dist_income


# Figurs: Variables in space and as a function of distance -------------------------------
  # Combine already-made figures
  (spatial_educ + dist_educ) / 
  (spatial_wealth + dist_wealth) / 
  (spatial_income + dist_income)


# Figures: Spatial figures: Predicted education ------------------------------------------
  # Find limits of education and its predictions
  min_e = gen_dt %>% select(starts_with("educ")) %>% min() %>% round(1)
  max_e = gen_dt %>% select(starts_with("educ")) %>% max() %>% round(1)
  lim_e = c(min_e,max_e)
  # Spatial plot of wealth
  spatial_wealth = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = wealth), size = 0.5) +
  scale_color_viridis_c("Wealth (endog.)", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of education
  spatial_educ = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ), size = 0.5) +
  scale_color_viridis_c("Education", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of OLS predictions
  spatial_ols = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ_pred_ols), size = 0.5) +
  scale_color_viridis_c("Pred. educ., OLS", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RF 
  spatial_rf = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ_pred_rf), size = 0.5) +
  scale_color_viridis_c("Pred. educ., RF", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RFSS
  spatial_rfss = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ_pred_rfss), size = 0.5) +
  scale_color_viridis_c("Pred. educ., RFSS", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Spatial plot of RFSS, linearized
  spatial_rfssl = ggplot(data = gen_dt, aes(x = x, y = y)) +
  geom_point(aes(color = educ_pred_rfss_linearized), size = 0.5) +
  scale_color_viridis_c("Pred. educ., RFSS line.", option = "magma", limits = lim_e, direction = -1) +
  theme_minimal(base_family = "Fira Sans Book", base_size = 14) +
  theme(legend.position = "bottom") +
  coord_equal()
  # Combine spatial plots
  (spatial_wealth + spatial_educ) / (spatial_ols + spatial_rf) / (spatial_rfss + spatial_rfssl)


# Figure: Education and predictions as a function of distance ----------------------------
  # (Predicted) Education as a function of distance
  ggplot(data = gen_dt, aes(x = dist)) +
  geom_hline(yintercept = 0, size = 1/4) +
  geom_point(aes(color = letters[1], y = educ), shape = 4, size = 1/5, alpha = 0.05) +
  geom_point(aes(color = letters[2], y = educ_pred_ols), size = 1.5) +
  # geom_point(aes(color = letters[3], y = educ_pred_rf), size = 1.5) +
  geom_point(aes(color = letters[3], y = educ_pred_xgss), size = 1.5) +
  geom_point(aes(color = letters[4], y = educ_pred_rfss), size = 1.5) +
  scale_x_continuous("Distance to university") +
  scale_y_continuous("Education level") +
  scale_color_manual(
    "Education:",
    labels = c("Actual", "OLS pred.", "RF pred.", "SS RF pred."),
    values = c("black", viridis::magma(3, begin = 0.5, end = 0.9))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
  
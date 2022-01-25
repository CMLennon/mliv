
# Notes ----------------------------------------------------------------------------------


# Setup ----------------------------------------------------------------------------------
  # Load packages
  library(pacman)
  p_load(
    tidyverse, qs, data.table,
    sf, geodist, 
    tidymodels, fixest, 
    parallel, magrittr, here
  )

# Load data: Locations of schools --------------------------------------------------------
  # Load dataset
  school_sf = here(
    "data-raw", "colleges-universities", "colleges-universities.shp"
  ) %>% st_read()
  # Names to lowercase
  school_sf %<>% janitor::clean_names()
  # Reproject to WGS84
  school_sf %<>% st_transform(crs = 4326)
  # Drop schools outside of contiguous US
  school_sf %<>% filter(
    !(state %in% c("AK", "HI", "PR", "GU", "PW", "MP", "MH", "VI", "AS", "FM"))
  )
  # To data table
  school_dt = copy(school_sf)
  school_dt %<>% cbind(st_coordinates(school_dt))
  setDT(school_dt)
  setnames(school_dt, old = c("X", "Y"), new = c("lon", "lat"))
  school_dt %<>% .[,.(
    school_id = objectid, 
    school = name,
    state, 
    fips = countyfips,
    lat, lon,
    type, 
    naics_code, naics_desc
  )]

# Load data: ACS data --------------------------------------------------------------------
  # Load the dataset's descriptions
  acs_info = here(
    "data-raw", "nhgis-acs", "nhgis-acs.csv"
  ) %>% fread(nrows = 1)  
  # Load the dataset
  acs_dt = here(
    "data-raw", "nhgis-acs", "nhgis-acs.csv"
  ) %>% fread(header = F, skip = 2)
  # Add names
  setnames(acs_dt, names(acs_info))
  # Transform 'info' dataset
  acs_info = data.table(
    name = names(acs_info),
    desc = acs_info[1,] %>% unlist()
  )

# Data work: Clean ACS dataset -----------------------------------------------------------
  # Grab and rename desired variables
  acs_clean = acs_dt[, .(
    gisjoin = GISJOIN,
    cbg = paste0(
      str_pad(STATEA, 2, "left", 0),
      str_pad(COUNTYA, 3, "left", 0),
      str_pad(TRACTA, 6, "left", 0),
      str_pad(BLKGRPA, 1, "left", 0)
    ),
    year = YEAR,
    state = STATE,
    county = paste0(
      str_pad(STATEA, 2, "left", 0),
      str_pad(COUNTYA, 3, "left", 0)
    ),
    tract = paste0(
      str_pad(STATEA, 2, "left", 0),
      str_pad(COUNTYA, 3, "left", 0),
      str_pad(TRACTA, 6, "left", 0)
    ),
    pop = JMAE001,
    pop_white = JMBE002,
    pop_black = JMBE003,
    pop_ai = JMBE004,
    pop_asian = JMBE005,
    pop_hispanic = JMJE012,
    pop_male = JNXE002,
    pop_femlae = JNXE011,
    pop_male_married = JNXE004,
    pop_female_married = JNXE013,
    pop_school = JNYE003 + JNYE027,
    educ_nohs = 
      JN9E003 + JN9E004 + JN9E005 + JN9E006 + JN9E007 + JN9E008 + JN9E009 + JN9E010 +
      JN9E020 + JN9E021 + JN9E022 + JN9E023 + JN9E024 + JN9E025 + JN9E026 + JN9E027,
    educ_hs = JN9E011 + JN9E028,
  # NOTE: Some college, no degree 
    educ_college = JN9E012 + JN9E013 + JN9E029 + JN9E030,
    educ_associate = JN9E014 + JN9E031,
    educ_bachelor = JN9E015 + JN9E032,
    educ_masters = JN9E016 + JN9E033,
    educ_pro = JN9E017 + JN9E034,
    educ_doc = JN9E018 + JN9E035,
  # NOTE: Population living within 2 times poverty rate
    pop_poverty = JOCE002 + JOCE003 + JOCE004 + JOCE005 + JOCE006 + JOCE007 + JOCE008,
    income_median_hh = JOIE001,
    pop_foreign = JT2E005
  )]
  # Add a few additional
  acs_clean[, `:=`(
    educ_college_any = educ_college + educ_associate + educ_bachelor + educ_masters + educ_pro + educ_doc,
    educ_college_degree = educ_associate + educ_bachelor + educ_masters + educ_pro + educ_doc
  )]
  
# Load data: CBG data --------------------------------------------------------------------
  # Load CBG centroids
  cbg_sf = here("data-raw", "cbg-centroids.qs") %>% qread()
  # Find CBG's state
  cbg_sf %<>% mutate(state = str_sub(GEOID, 1, 2))
  # Drop CBGs outside of contiguous US
  cbg_sf %<>% filter(!(state %in% c("78", "72", "69", "66", "60", "15", "02")))
  # Drop 'state'
  cbg_sf %<>% select(-state)
  # Define CRS
  st_crs(cbg_sf) = 4326

# Data work: Calculate CBG-to-school distances -------------------------------------------
  # Calculate distances: Nearest public university
  dist_uni_pub = geodist(
    st_coordinates(cbg_sf),
    school_dt[type == 1 & naics_code == 611310,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Calculate distances: Nearest private university
  dist_uni_priv = geodist(
    st_coordinates(cbg_sf),
    school_dt[type == 2 & naics_code == 611310,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Calculate distances: Nearest public junior college
  dist_juco_pub = geodist(
    st_coordinates(cbg_sf),
    school_dt[type == 1 & naics_code == 611210,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Calculate distances: Nearest private junior college
  dist_juco_priv = geodist(
    st_coordinates(cbg_sf),
    school_dt[type == 2 & naics_code == 611210,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Calculate distances: Nearest trade school
  dist_trade = geodist(
    st_coordinates(cbg_sf),
    school_dt[naics_code == 611519,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Calculate distances: Nearest cosmo/barber school
  dist_cosmo = geodist(
    st_coordinates(cbg_sf),
    school_dt[naics_code == 611511,.(lon,lat)],
    measure = "cheap"
  ) %>% apply(MARGIN = 1, FUN = min) %>% divide_by(1e3)
  # Create a distance dataset
  cbg_dist = data.table(
    cbg = cbg_sf$GEOID,
    dist_uni_pub,
    dist_uni_priv,
    dist_juco_pub,
    dist_juco_priv,
    dist_trade,
    dist_cosmo
  )
  # Add a few more
  cbg_dist[, `:=`(
    dist_uni = pmin(dist_uni_pub, dist_uni_priv),
    dist_juco = pmin(dist_juco_pub, dist_juco_priv)
  )]
  cbg_dist[, `:=`(
    dist_uni_juco = pmin(dist_uni, dist_juco)
  )]

# Data work: Join CBG datasets -----------------------------------------------------------
  # Join distances to ACS CBG data
  cbg_dt = merge(
    x = acs_clean,
    y = cbg_dist,
    by = "cbg",
    all.x = F,
    all.y = F
  )
  # Drop zero-population CBGs
  cbg_dt %<>% .[pop > 0]

# Regressions: Education on distance -----------------------------------------------------
  # Degrees per 1,000 on distance to nearest uni/college/juco
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni_juco |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Degrees per 1,000 on distance to nearest uni/college
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Degrees per 1,000 on distance to nearest uni/college and nearest juco
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni + dist_juco |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Degrees per 1,000 on distance to nearest uni/college and nearest juco w/ public/private
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # County FEs
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    county,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Tract FEs
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    tract,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Tract FEs, adding in trade schools
  feols(
    I(educ_college_degree/pop * 1e3) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv + dist_trade + dist_cosmo |
    tract,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])

# Regressions: Potential confounders on distance -----------------------------------------
  # Variables to check: 
  #   - pop
  #   - pop_white/pop
  #   - pop_black/pop
  #   - pop_foreign/pop
  #   - pop_poverty/pop
  #   - income_median_hh
# NOTE: Most confounders make distance look bad at all levels... except maybe tract.
#       No correlation with poverty rate (or median income) when we control for tract.
  # Distance to nearest uni/college/juco
  feols(
    I(pop_black/pop) ~ 
    dist_uni_juco |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Distance to nearest uni/college
  feols(
    I(pop_black/pop) ~ 
    dist_uni |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Distance to nearest uni/college and nearest juco
  feols(
    I(pop_black/pop) ~ 
    dist_uni + dist_juco |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Distance to nearest uni/college and nearest juco w/ public/private
  feols(
    I(pop_black/pop) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    state,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # County FEs
  feols(
    I(pop_black/pop) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    county,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])
  # Tract FEs
  feols(
    I(pop_black/pop) ~ 
    dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
    tract,
    data = cbg_dt
  ) %>% summary(cluster = cbg_dt[,state])

# Predict share of population with a degree using education, linear reg. -----------------
  # Add share of population with a degree
  cbg_dt[, share_degree := educ_college_degree/pop]
  # Run regression(s)
  cbg_dt[, `:=`(
    share_degree_hat = feols(
      share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv |
      tract,
      data = cbg_dt
    )$fitted.values
  )]

# Predict share of population with a degree using education, random forest ---------------
  # Define the recipe
  the_recipe = recipe(
    share_degree ~ dist_uni_pub,
    # share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
    data = cbg_dt
  )
  # Define the hyperparameter grid
  rf_grid = CJ(
    mtry = 1,
    # mtry = 1:4,
    # min_n = c(1, 5, 10, 25, 50),
    min_n = c(10, 25, 50), # NOTE: Did not choose 1 or 5 
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
    fit_i = wf_i %>% fit(cbg_dt)
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
    mc.cores = min(c(12, nrow(rf_grid)))
  ) %>% rbindlist(use.names = T, fill = T)
  # Arrange by OOB error
  setorder(rf_dt, error_oob)
  # Split the full dataset into two halves
  set.seed(12345)
  n = cbg_dt[,.N]
  split_assign = sample(x = 1:n, size = n, replace = F)
  split1 = split_assign[1:floor(n/2)] %>% cbg_dt[.,]
  split2 = split_assign[(floor(n/2)+1):n] %>% cbg_dt[.,]
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
      share_degree ~ dist_uni_pub,
      # share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
      data = split1
    )
  ) %>% fit(split1)
  # Fit 'chosen' models on split 2
  fit2 = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      share_degree ~ dist_uni_pub,
      # share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
      data = split2
    )
  ) %>% fit(split2)
  # Use fitted models to predict onto the other splits
  split1[, share_degree_rf := predict(fit2, new_data = split1)$.pred]
  split2[, share_degree_rf := predict(fit1, new_data = split2)$.pred]
  # Combine splits
  cbg_dt = list(split1, split2) %>% rbindlist(use.names = T, fill = T)
  setorder(cbg_dt, cbg)
  # Regress share_degree on share_degree_rf AND tract fixed effects
  cbg_dt[, `:=`(
    share_degree_rf_hat = feols(
      share_degree ~ share_degree_rf |
      tract,
      data = cbg_dt
    )$fitted.values
  )]

# In-sample RF predictions for first stage -----------------------------------------------
  # Fit 'chosen' model on full dataset
  fit_full = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
      data = cbg_dt
    )
  ) %>% fit(cbg_dt)
  # Use fitted model (trained on all observations) to predict onto the full dataset
  cbg_dt[, share_degree_rf_insample := predict(fit_full, new_data = cbg_dt)$.pred]

# Predict share of population with a degree using education, random forest with FEs ------
  # Create new dataset with distance and degree shares residualized on Census tracts
  cbg_resid = cbg_dt[, .(
    cbg,
    share_degree = feols(share_degree ~ 1 | tract, data = cbg_dt)$residuals,
    dist_uni_pub = feols(dist_uni_pub ~ 1 | tract, data = cbg_dt)$residuals,
    dist_uni_priv = feols(dist_uni_priv ~ 1 | tract, data = cbg_dt)$residuals,
    dist_juco_pub = feols(dist_juco_pub ~ 1 | tract, data = cbg_dt)$residuals,
    dist_juco_priv = feols(dist_juco_priv ~ 1 | tract, data = cbg_dt)$residuals
  )]
  # Define the recipe
  the_recipe = recipe(
    # share_degree ~ dist_uni_pub,
    share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
    data = cbg_resid
  )
  # Define the hyperparameter grid
  rf_grid = CJ(
    # mtry = 1,
    mtry = 1:4,
    min_n = c(1, 5, 10, 25, 50, 100),
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
    fit_i = wf_i %>% fit(cbg_resid)
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
  rf_resid_dt = mclapply(
    X = 1:nrow(rf_grid),
    FUN = rf_i,
    mc.cores = min(c(12, nrow(rf_grid)))
  ) %>% rbindlist(use.names = T, fill = T)
  # Arrange by OOB error
  setorder(rf_resid_dt, error_oob)
  # Split the full dataset into two halves
  set.seed(12345)
  n = cbg_resid[,.N]
  split_assign = sample(x = 1:n, size = n, replace = F)
  split1 = split_assign[1:floor(n/2)] %>% cbg_resid[.,]
  split2 = split_assign[(floor(n/2)+1):n] %>% cbg_resid[.,]
  # Define chosen model
  chosen_model = rand_forest(
    "regression",
    mtry = rf_resid_dt$mtry[1],
    trees = rf_resid_dt$trees[1],
    min_n = rf_resid_dt$min_n[1]
  ) %>% set_engine(engine = "ranger")
  # Fit 'chosen' models on split 1
  fit1 = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
      # share_degree ~ dist_uni_pub,
      data = split1
    )
  ) %>% fit(split1)
  # Fit 'chosen' models on split 2
  fit2 = workflow(
  ) %>% add_model(
    chosen_model
  ) %>% add_recipe(
    recipe(
      share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
      # share_degree ~ dist_uni_pub,
      data = split2
    )
  ) %>% fit(split2)
  # Use fitted models to predict onto the other splits
  split1[, share_degree_rf := predict(fit2, new_data = split1)$.pred]
  split2[, share_degree_rf := predict(fit1, new_data = split2)$.pred]
  # Combine splits
  cbg_resid = list(split1, split2) %>% rbindlist(use.names = T, fill = T)
  setorder(cbg_resid, cbg)
  setnames(cbg_resid, old = "share_degree_rf", new = "share_degree_rf_resid")
  # Join predictions to main dataset
  cbg_dt %<>% merge(
    y = cbg_resid[, .(cbg, share_degree_rf_resid)],
    by = "cbg",
    all = T
  )
  # Regress share_degree on share_degree_rf AND tract fixed effects
  cbg_dt[, `:=`(
    share_degree_rf_resid_hat = feols(
      share_degree ~ share_degree_rf_resid |
      tract,
      data = cbg_dt
    )$fitted.values
  )]

# Check correlation of predictions with potential confounders ----------------------------
  # Summary function: Cluster at state level
  s = . %>% summary(cluster = cbg_dt[,state])
  # OLS predictions: Correlations with potential confounders
  feols(pop ~ share_degree_hat | tract, data = cbg_dt) %>% s()
  feols(pop_white/pop ~ share_degree_hat | tract, data = cbg_dt) %>% s()
  feols(pop_black/pop ~ share_degree_hat | tract, data = cbg_dt) %>% s()
  feols(pop_foreign/pop ~ share_degree_hat | tract, data = cbg_dt) %>% s()
  feols(pop_poverty/pop ~ share_degree_hat | tract, data = cbg_dt) %>% s()
  # RF predictions: Correlations with potential confounders
  feols(pop ~ share_degree_rf_insample | tract, data = cbg_dt) %>% s()
  feols(pop_white/pop ~ share_degree_rf_insample | tract, data = cbg_dt) %>% s()
  feols(pop_black/pop ~ share_degree_rf_insample | tract, data = cbg_dt) %>% s()
  feols(pop_foreign/pop ~ share_degree_rf_insample | tract, data = cbg_dt) %>% s()
  feols(pop_poverty/pop ~ share_degree_rf_insample | tract, data = cbg_dt) %>% s()
  # RF-based instrument predictions: Correlations with potential confounders
  feols(pop ~ share_degree_rf_hat | tract, data = cbg_dt) %>% s()
  feols(pop_white/pop ~ share_degree_rf_hat | tract, data = cbg_dt) %>% s()
  feols(pop_black/pop ~ share_degree_rf_hat | tract, data = cbg_dt) %>% s()
  feols(pop_foreign/pop ~ share_degree_rf_hat | tract, data = cbg_dt) %>% s()
  feols(pop_poverty/pop ~ share_degree_rf_hat | tract, data = cbg_dt) %>% s()
  # RF-based resid. instrument predictions: Correlations with potential confounders
  feols(pop ~ share_degree_rf_resid_hat | tract, data = cbg_dt) %>% s()
  feols(pop_white/pop ~ share_degree_rf_resid_hat | tract, data = cbg_dt) %>% s()
  feols(pop_black/pop ~ share_degree_rf_resid_hat | tract, data = cbg_dt) %>% s()
  feols(pop_foreign/pop ~ share_degree_rf_resid_hat | tract, data = cbg_dt) %>% s()
  feols(pop_poverty/pop ~ share_degree_rf_resid_hat | tract, data = cbg_dt) %>% s()

# Regressions: 2SLS results --------------------------------------------------------------
  # Option 0: Naïve OLS (still with tract FEs)
  result0 = feols(
    income_median_hh ~
    share_degree |
    tract,
    data = cbg_dt[dist_uni_pub < Inf]
  ) %>% summary(cluster = cbg_dt[dist_uni_pub < Inf,state])
  # Option 1: Classic OLS-based 2SLS
  result1 = feols(
    income_median_hh ~ 
    1 |
    tract |
    # share_degree ~ dist_uni_pub + dist_uni_priv + dist_juco_pub + dist_juco_priv,
    share_degree ~ dist_uni_pub,
    data = cbg_dt[dist_uni_pub < Inf]
  ) %>% summary(cluster = cbg_dt[dist_uni_pub < Inf,state])
  # Option 2: RF in the first stage
  result2 = feols(
    income_median_hh ~
    share_degree_rf_insample |
    tract,
    data = cbg_dt[dist_uni_pub < Inf]
  ) %>% summary(cluster = cbg_dt[dist_uni_pub < Inf,state])
  # Option 3: MLSS w/ RF-synthesized instrument
  result3 = feols(
    income_median_hh ~
    1 |
    tract |
    share_degree ~ share_degree_rf,
    data = cbg_dt[dist_uni_pub < Inf]
  ) %>% summary(cluster = cbg_dt[dist_uni_pub < Inf,state])
  # Option 4: MLSS w/ residualized RF-synthesized instrument
  result4 = feols(
    income_median_hh ~
    1 |
    tract |
    share_degree ~ share_degree_rf_resid,
    data = cbg_dt[dist_uni_pub < Inf]
  ) %>% summary(cluster = cbg_dt[dist_uni_pub < Inf,state])
  # Results tables
  setFixest_dict(c(
    income_median_hh = "Median HH income",
    share_degree = "Pop. share with college degree",
    share_degree_rf = "RF-pred. share with college degree",
    share_degree_rf_resid = "Residualized RF-pred. share with college degree",
    tract = "Census tract",
    dist_uni_pub = "Dist. to public univ. (km)",
    dist_uni_priv = "Dist. to private univ. (km)",
    dist_juco_pub = "Dist. to public junior coll. (km)",
    dist_juco_priv = "Dist. to private junior coll. (km)"
  ))
  etable(
     "Naive OLS" = result0,
     "2SLS" = result1,
     "RF in 2SLS" = result2,
     "MLSS w/ RF" = result3,
     "Resid. MLSS w/ RF" = result4,
     cluster = cbg_dt[dist_uni_pub < Inf, state],
     stage = 2
  )
  etable(
     "Naive OLS" = result0,
     "2SLS" = result1,
     "RF in 2SLS" = result2,
     "MLSS w/ RF" = result3,
     "Resid. MLSS w/ RF" = result4,
     cluster = cbg_dt[dist_uni_pub < Inf, state],
     stage = 1
  )

# Figures --------------------------------------------------------------------------------
  # 
  ggplot(
    data = cbg_dt,
    aes(x = dist_uni_pub, y = educ_college_any/pop * 1e3)
  ) + 
  geom_point(size = 0.1) +
  geom_smooth(method = lm, color = "orange") +
  geom_smooth() +
  scale_x_continuous("Distance to nearest public univeristy (km)") +
  scale_y_continuous("College degrees per 1,000 residents (CBG)", labels = scales::comma) +
  theme_minimal()
  # 
  ggplot(
    data = cbg_dt,
    aes(x = dist_uni_pub, y = share_degree_rf_resid * 1e3)
  ) + 
  geom_point(size = 0.1) +
  geom_smooth(method = lm, color = "orange") +
  geom_smooth() +
  scale_x_continuous("Distance to nearest public univeristy (km)") +
  scale_y_continuous("College degrees per 1,000 residents (CBG)", labels = scales::comma) +
  theme_minimal()

# Todo -----------------------------------------------------------------------------------
#   - Add ML 'best practice' to first stage (i.e., CV): Does it reduce bias?
#   - done. It does - kind of. It makes the bias less predictible. It's possible    a RF with CV may avoid this issue with a bit of luck.
#   - added this feature - had conversation with Ed 10/20. Will add the real intent today, 10/21
#
#     NOTE: Add (not replace)
#   - Compare to PCA (maybe think of correlating the z's?)
#   - added as of 1/19/20. Reduces variance compared to 2SLS and is unbiased. Likely a function of correlation between instruments.

#     Added both lasso and post-lasso procedures to the model set.
#     TODO: Add a regression using the error term for y and e for each model.

# Setup ----------------------------------------------------------------------------------
# Packages
library(pacman)
p_load(
  tidyverse, broom, data.table,
  lfe, fixest, estimatr, party,
  parallel, magrittr, MASS, partykit, rsample, furrr, future.apply, adabag,
  mboost, promises, elasticnet, coefplot, glmnet, tidymodels, here
)

# Function: Generate Xs and Ys ---------------------------------------------------------------
generate_system = function(i_dt, n){
  #Utility function to generate our system of equations so 
  #that it can be done for multiple data sets
  β0 = 1
  β1 = 1
  β2 = 1
  num_zs = 7
  corr = .2
  # z affects x1 (plus correlated noise)
  i_dt[, x1 := 1 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + e_common]
  # z^2 and affects x2 (plus correlated noise)
  # i_dt[, x2 := 1 + 10 * z^2 + y_err_full]
  i_dt[, x2 := 1 + e_common]
  # Calculate y
  i_dt[, y_err := runif(n, min = -1, max = 1)]
  i_dt[, y := β0 + β1 * x1 + β2 * x2 + y_err]
  i_dt[, y_err_full := y_err + e_common]
  return(i_dt)
}

original_data = function(n){
  #  Original DGP - replaced with a multivariate process to allow
  # correlation  between instruments.
  
  #Takes in n, which is the number of observations in output
  
  # Generate the instruments
  i_dt = data.table(
    z1 = runif(n = n, min = -1, max = 1),
    z2 = runif(n = n, min = -1, max = 1),
    z3 = runif(n = n, min = -1, max = 1),
    z4 = runif(n = n, min = -1, max = 1),
    z5 = runif(n = n, min = -1, max = 1),
    z6 = runif(n = n, min = -1, max = 1),
    z7 = runif(n = n, min = -1, max = 1),
    e_common = rnorm(n)
  )
  
  #use generate system ancillary function to generate system of data
  i_dt = generate_system(i_dt, n)
  return(i_dt)
}

# Functions for K-fold and Cross Validation--------------------------------------------

#Caret implementations of lasso/enet and adaboost
adamod = function(dt, k_folds) {
  adaGrid <-  expand.grid(mstop = c(30:90), maxdepth = c(1:10))
  adafit <- train(as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), 
                  method = 'blackboost',
                  trControl = trainControl(sample, method = 'cv', number = k_folds),
                  data = dt,
                  tuneGrid = adaGrid
  )
  
  return(adafit)
}

lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), 
                       nfolds = k_folds,
                       alpha = 1,
                       lambda = seq(from = 0, to = 1000, length.out = 10000)
                       )
  return(lassofit)
}

post_lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), 
                       nfolds = k_folds,
                       alpha = 1,
                       lambda = seq(from = 0, to = 1000, length.out = 10000)
  )
  pl = rownames(coef(lassofit, s = 'lambda.min'))[coef(lassofit, 
                                                  s = 'lambda.min')[,1]!= 0]
  reg = lm(as.formula(paste0('x1 ~ ', paste(tail(pl,-1), collapse = " + "))), data = dt)
  return(reg)
}

forest_fcn = function(dta, num_zs = 7, znames){
  
  #Function that builds and finalizes a tidymodels randomforest, using
  #ranger as its engine. It takes data, a formula to estimate, and returns
  #a finalized workflow object
  
  ### Inputs:
  
  # dta: a data-frame like object. Must be compatible with tidymodels
  
  # num_zs: the number of instruments used
  
  #znames: the vector of instrument variable names, as character values
  
  #lossfcn: a string containing 'mse' or 'mae' to determine how the model
  #will select hyperparameters. Default is mse.
  #dta = i_dt_corr_test
  #set recipipe = recipe(as.formula(paste0('x1 ~ ', paste0(znames, collapse = '+'))), 
  forest_recipe = recipe(as.formula(paste0('x1 ~ ', paste0(znames, collapse = '+'))), 
                         data = dta)
  
  #build splits (5-fold cv)
  resamps = vfold_cv(dta, v =5)
  
  #set model, make the models tunable on mtry and smallest leaf
  forest_model = rand_forest() %>% set_engine('ranger') %>% 
    set_args(trees = tune(), mtry = tune(), min_n = tune()) %>% set_mode('regression')
  
  forest_tuner = expand_grid(trees = seq(from = 50, to = 300, by = 50), mtry = 2:(num_zs-1), min_n = 1:8)
  
  #build workflow object
  forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)
  
  #find best CVd parameters
  best = tune_grid(forest_wf, grid = forest_tuner, 
                   resamples = resamps, metrics = metric_set(rmse, rsq)) %>% 
    select_best(metric = 'rmse', maximize = F)
  
  forest_wf = forest_wf %>% finalize_workflow(best)
  
  return(forest_wf)
}

# Function: One iteration ----------------------------------------------------------------
one_iter = function(iter, n) {
  # Set parameters
  β0 = 1
  β1 = 1
  β2 = 1
  num_zs = 7
  corr = .7
  #set up a covar mat with off-diagonals = corr^(|j-h|) where j is own-z (ie, for z2 = 2) and h is the correlated variable's z. Eg. z1 and z3 have corr = corr^2 and z1 has corr = corr^0=1.
  #cvmat = matrix(c(rep(c(1,rep(corr, num_zs)), num_zs-1)), nrow = num_zs, byrow = TRUE) {archived to have more variation in correlation levels across zs}
  cvmat = matrix(0,num_zs,num_zs)
  for (j in (1:num_zs)) {
    for (h in (1:num_zs)) {
      cvmat[c(j),c(h)] = corr^abs(j-h)
    }
  }
  
  i_dt_corr_test = data.table(
    mvrnorm(n = n, mu = rep(0, num_zs), Sigma = cvmat),
    rnorm(n))
  
  znames = sapply(1:num_zs, FUN = function(x) paste0('z_',as.character(x)))
 
  #replace names in multivariate matrix with z-names.
  i_dt_corr_test = setnames(i_dt_corr, c(znames, "e_common"))
  
  #replicate Ed's code using the new utility function (which is really just Ed's code)
  
  #generate system for correlated system, then split for cv process.
  i_dt_corr_test = generate_system(i_dt_corr, n)
  
  # OLS first stage
  i_dt_corr_test[, x1_hat_ols := predict(lm(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test
  ))]
  
  #error from OLS first stage
  i_dt_corr_test[, err_ols := x1 - x1_hat_ols]
  
  # Random forest, no CV
  i_dt_corr_test[, x1_hat_rf := predict(party::cforest(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test
  ))]
  
  #Error from error, random forest
  i_dt_corr_test[, err_rf := x1 - x1_hat_rf]
  
  #PCA for meta-instruments (dimensionality reduction)
  i_dt.pca = prcomp(i_dt_corr_test[,names(i_dt_corr) %in% znames, with = FALSE], center = TRUE,scale. = TRUE)
  
  #PCA, OLS process for our regression
  i_dt_corr_test[, x1_hat_pca := predict(lm(
    x1 ~ i_dt.pca$x[,1] + i_dt.pca$x[,2] + i_dt.pca$x[,3],
    data = i_dt_corr_test
  ))]
  i_dt_corr_test[, err_pca := x1 - x1_hat_pca]
  
  #Generate tidymodels workflow 
  rfobj = forest_fcn(dta = i_dt_corr_test, znames = znames, num_zs = 7)
  
  #fit workflow
  rffit = parsnip::fit(rfobj, data = i_dt_corr_test)
  
  #find in-sample fit values
  i_dt_corr_test[, x1_hat_rf_cv := rffit %>% predict(new_data = i_dt_corr_test),
                 ]
  
  #error for the cross-validated rf
  i_dt_corr_test[, err_rf_cv := x1-x1_hat_rf_cv]
  
  #Lasso model to select variables for first stage (post lasso)
  plasso_mod = post_lasso(i_dt_corr_test, k = 5, znames)
  
  i_dt_corr_test[, x1_hat_plasso_cv := predict(plasso_mod,
                                              newdata = i_dt_corr_test)]
  
  
  i_dt_corr_test[, err_plasso := x1 - x1_hat_plasso_cv]

  #Lasso model
  lasso = lasso(i_dt_corr_test, k = 5, znames)
  
  i_dt_corr_test[, x1_hat_lasso_cv := predict(lasso,
                                               newx = i_dt_corr_test[,znames, with = FALSE] %>% 
                                                as.matrix())]
  
  i_dt_corr_test[, err_lasso := x1 - x1_hat_lasso_cv]
  
  #Boosted Trees model for first stage
  adamodel = adamod(i_dt_corr_test, k = 5)
  
  i_dt_corr_test[, x1_hat_boost_cv := predict(adamodel,
                                              newdata = i_dt_corr_test)]
  
  i_dt_corr_test[, err_boost_cv := x1 - x1_hat_boost_cv]
  
  #Run second regressions
  bind_rows(
    # OLS: DGP
    lm(y ~ x1 + x2, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "OLS: DGP"),
    
    # OLS: OVB
    lm(y ~ x1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "OLS: OVB"),
    # OLS: OVB error coef
    lm(y ~ x1 + y_err_full, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: OLS: OVB"),
    
    # Second stage: OLS
    lm(y ~ x1_hat_ols, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: OLS"),
    #2S-LS Y error
    lm(y ~ x1 + err_ols, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: OLS "),
    
    # Second stage: Random Forest, no CV
    lm(y ~ x1_hat_rf, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Random Forest"),
    # Errors: Random Forest, no CV
    lm(y ~ x1 + err_rf, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: Random Forest"),
    
    # Second stage: Random Forest with CV
    lm(y ~ x1_hat_rf_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Random Forest, CV"),
    # Errors: Random Forest with CV
    lm(y ~ x1 + err_rf_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: Random Forest, CV"),
    
    # Second Stage: PCA instruments
    lm(y ~ x1_hat_pca, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: PCA"),
    #Errors: PCA instruments
    lm(y ~ x1 + err_pca, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: PCA"),
    
    #Second Stage: Lasso Selection (post-lasso)
    lm(y ~ x1_hat_plasso_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: post-LASSO selection"),
    #Errors: Post-lasso
    lm(y ~ x1 + err_plasso, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: post-LASSO selection"),
    
    #Second Stage: Pure Lasso
    lm(y ~ x1_hat_lasso_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: LASSO selection"),
    #Errors: Pure Lasso
    lm(y ~ x1 + err_lasso, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: LASSO selection"),
    
    #Second Stage: Boosted Trees
    lm(y ~ x1_hat_boost_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: Boosted Trees"),
    #Errors: Boosted Trees
    lm(y ~ x1 + err_boost_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: Boosted Trees")
    
    #important to add a CV'd PCA?
  ) %>% mutate(iter = iter) %>% data.table()
}

# Function: Run iterations ---------------------------------------------------------------
run_sim <-function(n, n_sims, seed = 12345,future_cores = 12) {
  future::plan("multicore", workers = future_cores)
  # Set the seed
  set.seed(seed)
  # Run a parallelized simulation
  future_lapply(X = 1:n_sims, FUN = one_iter, n = n) %>% rbindlist()
}
#stopCluster

# Run simulation -------------------------------------------------------------------------
sim_dt <- run_sim(n = 1e3, n_sims = 5)

#write out data
#write_csv(sim_dt, here('Plots:Resources/SimulationDT.csv'))


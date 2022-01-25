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


################## OUTPUT PATH
write_path = '/Users/connor/Dropbox/MLIV/Plots:Resources/SimulationDT.csv'
##################

library(pacman)
p_load(
  tidyverse, broom, data.table,
  lfe, fixest, estimatr, party,
  parallel, magrittr, MASS, partykit, rsample, furrr, future.apply, adabag,
  mboost, promises, elasticnet, coefplot, glmnet, tidymodels, rlang
)

# Function: Generate Xs and Ys ---------------------------------------------------------------
generate_system = function(i_dt, n){
  
  #Utility function to generate our system of equations so that it can be done for multiple data sets
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

# Functions for K-fold and Cross Validation---------------------------------------------------------------
#function to produce squared errors on holdout fold
num_zs = 7

#function to create mse metric in tidymodels, from documentation
mse_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  
  mse_impl <- function(truth, estimate) {
    mean((truth - estimate) ^ 2)
  }
  
  metric_vec_template(
    metric_impl = mse_impl,
    truth = truth, 
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric",
    ...
  )
  
}

mse <- function(data, ...) {
  UseMethod("mse")
}

mse.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  
  metric_summarizer(
    metric_nm = "mse",
    metric_fn = mse_vec,
    data = data,
    truth = !! enquo(truth),
    estimate = !! enquo(estimate), 
    na_rm = na_rm,
    ...
  )
  
}

class(mse) <- c("numeric_metric", class(mse))
class(mse_vec) <- c('numeric_metric', class(mse_vec))

forest_fcn = function(dta, num_zs = 7, znames){
  
  #Function that builds and finalizes a tidymodels randomforest, using
  #ranger as its engine. It takes data, a formula to estimate, and returns
  #a list of two items, a finalized workflow object and the chosen hyper-
  #parameters
  
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
  
  forest_tuner = expand_grid(mtry = 2:num_zs, min_n = 3:20, trees = seq(from = 50, to = 350, by = 50))
  
  #build workflow object
  forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)

  #find best CVd parameters
  best = tune_grid(forest_wf, grid = forest_tuner, 
                   resamples = resamps, metrics = metric_set(rmse, rsq)) %>% 
    select_best(metric = 'rmse', maximize = F)
  
  forest_wf = forest_wf %>% finalize_workflow(best)
  
  deliverable = list(forest_wf, best)
  
  return(deliverable)
  }

# Function: One iteration ----------------------------------------------------------------
one_iter = function(iter, n) {
  # Set parameters
  β0 = 1
  β1 = 1
  β2 = 1
  num_zs = 7
  corr = .2
  # Generate the instrument
  
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
  
  #set up a covar mat with off-diagonals = corr^(|j-h|) where j is own-z (ie, for z2 = 2) and h is the correlated variable's z. Eg. z1 and z3 have corr = corr^2 and z1 has corr = corr^0=1.
  #cvmat = matrix(c(rep(c(1,rep(corr, num_zs)), num_zs-1)), nrow = num_zs, byrow = TRUE) {archived to have more variation in correlation levels across zs}
  cvmat = matrix(0,num_zs,num_zs)
  corr = .5
  for (j in (1:num_zs)) {
    for (h in (1:num_zs)) {
      cvmat[c(j),c(h)] = corr^abs(j-h)
    }
  }
  
  i_dt_corr = data.table(
    mvrnorm(n = n, mu = rep(0, num_zs), Sigma = cvmat),
    rnorm(n))
  
  znames = sapply(1:num_zs, FUN = function(x) paste0('z',as.character(x)))
  
  #replace names with z-names.
  i_dt_corr = setnames(i_dt_corr, c(znames, "e_common"))
  
  #replicate Ed's code using the new utility function (which is really just Ed's code)
  i_dt = generate_system(i_dt, n)
  
  #generate system for correlated system, then split for cv process.
  i_dt_corr_test = generate_system(i_dt_corr, n)
  
  #Machine learn with cv on the first stage
  #find optimal oos control parameters
 
  #build control parameter object with CV'd optimal parameters
  rfobj = forest_fcn(dta = i_dt_corr, znames = znames, num_zs = 7)
  
  best = rfobj[[2]]
  rfwf= rfobj[[1]]
  rfobj = parsnip::fit(rfwf, data = i_dt_corr_test)
  
  #Build model with test data
  i_dt_corr_test[, x1_hat_rf_cv := rfobj %>% predict(new_data = i_dt_corr_test),
    ]
  
  
  #Run second regression
  bind_rows(
    #check tuning parameter variation
    lm(y ~ x1_hat_rf_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Random Forest, CV", mtry = best[[1]], 
             min_n = best[[2]], trees = best[[3]]),
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
sim_dt <- run_sim(n = 1e3, n_sims = 90)
write_csv(sim_dt, write_path)
# Plot simulation ------------------------------------------------------------------------
sim_dt[, median(mtry, na.rm = TRUE), by = model]
sim_dt[, median(min_n, na.rm = TRUE), by = model]
sim_dt[, median(estimate, na.rm = TRUE), by = model]


# The beginning of the R script (run-sim.R) ----------------------------------------------

#read in command-line arguments
command_args = commandArgs(trailingOnly = T)
# Note:
#   - i[1] determines data generation: 'F' for simple and 'T' for Belloni
#   - i[2] gives type of Belloni β pattern. These follow the paper's
# "exponential" design, which is β1 = .7, β2 = .7^2... We modify this to test
# PCA's functionality more thoroughly
#       - "1" is a shuffled β pattern
#       - "2" is standard Belloni paper approach
#       - "3" is set so the strongest β is β50

#   - i[3] gives number of observations per iteration
#   - i[4] gives number of simulations to run in total


#   - tail(i, -4) gives methods to NOT run
#       - "dgp" : run a model that uses x without shared error term
#       - 'ovb' : run a naive regression y ~ x1
#       - "2sls": run two stage least squares
#       - "lasso": run lasso on the first stage, with predicted results in 2nd stage
#       - "plasso": as above, but run regression using all non-zero coefficients from a lasso model (post lasso method)
#       - "rf": run a non-crossvalidated model in ranger (using ranger defaults) for first stage
#       - "rfcv" : cross-validate first stage random forest.
#       - "pca": build a set of synthetic instruments using priciple component analysis, cross-validated against x-hat/x mse.
#       - "gboost": xgboost algorithm to estimate the first stage, cross-validated.
#       - "split": implemented, split-sample IV
#       - "LIML": Limited information maximum likelihood. No variance or covariance can be returned given how this estimate is generated
#       - "JIVE": Jack-knife Instrumental Variables Estimator. No variance or covariance can be returned given how this estimate is generated



#load necessary packages
library(pacman)
p_load(tidymodels, tidyverse, data.table, here, 
       MASS, glmnet, xgboost, ranger, future, future.apply, plotly, SteinIV, ivmodel)

# All methods
methods_all = c("dgp", 'ovb', "2sls", "lasso", "plasso", "rf", "rfcv", "pca", "gboost", "split", "LIML", 'JIVE')

# Methods to run
methods_run = setdiff(methods_all, tolower(tail(command_args, -4)))

#set parameters from 1st command argument
if (command_args[1] == "F") {
  # Source the function that generates the simple data
  source(here("codeR", "simple_data_function.R"))
  k = 7
  complex = F
  # Define the type of beta distribution (irrelevant for 'simple')
  shuffle = command_args[2]
} else {
  # Source the function that generates the Belloni data
  source(here("codeR", "Belloni_aux_function_shuffle.R"))
  k = 100
  complex = T
  # Define the type of beta distribution
  shuffle = command_args[2]
}

#pull argument 3 (number of observations) into the simulation
numobs = as.numeric(command_args[3])
#pull argument 4 (number of simulations) into full workspace
numsims = as.numeric(command_args[4])

one_iter = function(iter, n, k, shuffle, methods_run, complex) {
  #Purpose: generate data, generate 5-fold set to do consistent
  #cross-validation, run methods according to command-argument 3.
  
  #-----------STRATEGY-------------
  
  #toggle models from methods_run vector, reading in auxiliary functions
  #that use tidymodels crossvalidation/evaluation to produce meaningful
  #hyperparameters. These functions return differing versions of x-hat
  #for a first stage.
  
  #They take (plus or minus certain elements specific to models): a dataset
  #a pre-specified rsample object, and a randomized seed when necessary
  
  #-----------INPUTS---------------
  
  #iter: simulation number for within data performance comparisons
  
  #n: number of datapoints in our dataset
  
  #k: number of instruments
  
  #shuffle: **read from command line** shuffle strategy for beta
  #pattern. "1": fully shuffled, "2": original. "3": 50th instrument is
  #strongest.
  
  #methods_run: **partially read from command line** a vector of strings
  #for what 1st stage method used to estimate x1.
  
  #-----------OUTPUTS----------------
  
  #iter_dt : a data.table object with 1 row per item in methods_run

  #-----------------------------------------
  #ITER_DT Contents:

  #For each sub-model in the simulation:

  #model: model used in first stage (character)
  #iter: iteration number - identifies simulation number
  #term: beta_IV estimate
  #true_var/true_cov = variance of x1, covariance between u_i and x_i, if u_i = all of y != x
  #var/cov = variance of x1_hat, covariance of x1_hat and u_i


  #-----------------------------------------
  
  #generate data based on "complex" argument from command line (T/F)
  i_dt_corr = data_gen(n,k,shuffle)
  
  #folds object for cross validation
  folds = vfold_cv(i_dt_corr, v = 5)
  
  #initialize a data table to create new rows inside of. This holds second
  #stage beta results and an iter number to track intra-iteration variation
  
  iter_dt = data.table()
    # OLS: Oracle Model
  if('dgp' %in% methods_run){
    if(complex){
      vncoef = summary(lm(y ~ x1 + vn1, data = i_dt_corr))$coefficients[3,1]
      iter_dt = lm(y ~ x1 + vn1, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model")  %>%mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1)) %>% mutate(cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1 - vncoef*i_dt_corr$vn1)) %>% bind_rows(.,iter_dt)
    }
    else{
      iter_dt = lm(y ~ x1 + x2, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model") %>%  mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% mutate(true_var = var(i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1)) %>% mutate(cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1 - i_dt_corr$x2)) %>% bind_rows(.,iter_dt)
    }
  }

    
    # OLS: Naive Regression
  if('ovb' %in% methods_run){
      iter_dt = lm(y ~ x1, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First Stage: Naive OLS") %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% mutate(true_var = var(i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1)) %>% mutate(cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% data.table()
      }
    
    # Second stage: 2SLS
    if('2sls' %in% methods_run){
    zs = grep("z", names(i_dt_corr), value = "T")
    
    formula = as.formula(paste0("x1 ~ ", paste0(zs, collapse = '+')))
    
    i_dt_corr[,x1_hat_ols := lm(formula, data = i_dt_corr) %>% predict(new_data = i_dt_corr)]  
    
    iter_dt = lm(y ~ x1_hat_ols , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: 2SLS") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_ols)) %>% mutate(cov = cov(i_dt_corr$x1_hat_ols, i_dt_corr$y - i_dt_corr$x1)) %>%  bind_rows(.,iter_dt) %>% 
        data.table()
    }
    
    # Second stage: Random Forest, no CV instrument
    if('rf' %in% methods_run){
    source(here("codeR", "RandomForestCVParsnip.R"))
    i_dt_corr[,x1_hat_rf := forest_fcn(datatable = i_dt_corr,
                                      folds = folds,
                                      cv = F)]
    
    iter_dt = lm(y ~ x1_hat_rf, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: Random Forest")   %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_rf)) %>% mutate(cov = cov(i_dt_corr$x1_hat_rf, i_dt_corr$y - i_dt_corr$x1))%>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
    
    # Second stage: Random Forest with CV instrument
    if('rfcv' %in% methods_run){
      source(here("codeR", "RandomForestCVParsnip.R"))
      i_dt_corr[,x1_hat_rf_cv := forest_fcn(datatable = i_dt_corr,
                                        folds = folds,
                                        cv = T)]
      
    iter_dt = lm(y ~ x1_hat_rf_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: Random Forest, CV") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_rf_cv)) %>% mutate(cov = cov(i_dt_corr$x1_hat_rf_cv, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
    
    # Second stage: PCA synthetic instruments
    if('pca' %in% methods_run){
      source(here("codeR", "PCAParsnipCV.R"))
      i_dt_corr[,x1_hat_pca := pca_fit_cv(datatable = i_dt_corr,
                                        folds = folds,
                                        cv = T)]
      
    iter_dt = lm(y ~ x1_hat_pca , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: PCA") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_pca)) %>% mutate(cov = cov(i_dt_corr$x1_hat_pca, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
    
    #Second Stage: Lasso Instruments
    if('lasso' %in% methods_run){
    source(here("codeR", "LassoCVParsnip.R"))
    i_dt_corr[,x1_hat_lasso_cv := lasso_fit_cv(datatable = i_dt_corr,
                                         folds = folds)]
    
    iter_dt = lm(y ~ x1_hat_lasso_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: LASSO selection") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_lasso_cv)) %>% mutate(cov = cov(i_dt_corr$x1_hat_lasso_cv, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
    
    #Second Stage: Post-Lasso Instruments
    if('plasso' %in% methods_run){
    source(here("codeR", "PostLassoCVParsnip.R"))
    i_dt_corr[,x1_hat_plasso_cv := plasso_fit_cv(datatable = i_dt_corr,
                                         folds = folds)]
    
    
    iter_dt = lm(y ~ x1_hat_plasso_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: post-LASSO selection") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_plasso_cv)) %>% mutate(cov = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
  
    #Second Stage: Boosted Trees Instrument
    if('gboost' %in% methods_run){
      source(here("codeR", "BoostedTreesParsnipCV.R"))
    i_dt_corr[,x1_hat_boost_cv := xgbst_fit(datatable = i_dt_corr,
                                               folds = folds)]
    
    iter_dt = lm(y ~ x1_hat_boost_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: Boosted Trees")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_boost_cv)) %>% mutate(cov = cov(i_dt_corr$x1_hat_boost_cv, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    }
  
  if('split' %in% methods_run){
    source(here('codeR', "SplitSampleIV.R"))
    i_dt_corr[,x1_hat_ussiv := ssiv_fcn(datatable = i_dt_corr,
                                        folds = folds,
                                        unbiased = T)]
    iter_dt = lm(y ~ x1_hat_ussiv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: USSIV")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_ussiv)) %>% mutate(cov = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
  }
  if('LIML' %in% methods_run){
    #LIML estimation here is using the IVModels package. There is a handwritten version as well that we can use should we want it,
    # especially if we need some form of X-hat.
    source(here('codeR', "LIMLest.R"))
    i_dt_corr[,x1_hat_liml := NA]
    LIMLobj = LIML_fit(datatable = i_dt_corr, folds = folds, iter = iter)
    iter_dt = data.frame(term = 'x1_LIML') %>% mutate(estimate = LIMLobj$point.est[1], std.error = LIMLobj$std.err[1], 
                                                      statistic =LIMLobj$test.stat[1], p.value = LIMLobj$p.value[1]) %>%
      mutate(model = "First stage: LIML (Fuller)") %>% mutate(var = NA, cov = NA, true_var = var(i_dt_corr$x1), true_cov =  cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      bind_rows(.,iter_dt) %>% data.table()
  }
  if('JIVE' %in% methods_run){
    source(here('codeR', 'JIVEest.R'))
    i_dt_corr[,x1_hat_jive := NA]
    JIVEobj = JIVE_fit(datatable = i_dt_corr, folds = folds, iter = iter)
    iter_dt = data.frame(term = 'x1_JIVE') %>% mutate(estimate = JIVEobj$est[1], std.error = JIVEobj$se[1], 
                                                      statistic =JIVEobj$test_statisitc, p.value = JIVEobj$p.value) %>%
      mutate(model = "First stage: JIVE") %>% mutate(var = NA, cov = NA, true_var = var(i_dt_corr$x1), true_cov =  cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      bind_rows(.,iter_dt) %>% data.table()
  }
    
  iter_dt %>% mutate(iter = iter)
  return(iter_dt)
}

run_sim <-function(n, n_sims, seed = 12345,
                   future_cores = availableCores(),
                   methods_run,
                   shuffle, k, complex) {
  #---------------------------------------------------------------
  # Using future, we run a parallelized simulation of all models in 
  # methods_run vector, across k instruments, utilizing β-pattern strat
  # specified in shuffle. Defaults to utilizing all available cores,
  # and running all models. Returns:
  
  # sim_dt, a datatable containing an iteration-labeled set of second
  # stage estimates of β-hat.
  #--------------------------------------------------------------
  
  future::plan("multicore", workers = future_cores)
  # Set the seed
  set.seed(seed)
  # Run a parallelized simulation
  future_lapply(X = 1:n_sims, FUN = one_iter, n = n, 
                methods_run = methods_run,
                shuffle=shuffle, k = k, complex) %>% rbindlist()
}


#run simulation
sim_dt <- run_sim(n = parent.frame()$numobs, n_sims = parent.frame()$numsims, methods_run = parent.frame()$methods_run,
                  shuffle = parent.frame()$shuffle, k = parent.frame()$k, complex = complex)

#print median beta values to console
sim_dt[, mean(estimate, na.rm = TRUE), by = model]

#write to a file with the following file-naming conventions:

# sim-results-(complex dgp T/F)-(shuffle_strategy == (1,2,3))-(models-run).csv

# Save results ---------------------------------------------------------------------------
  # Source post-processing function from codeR/postprocess-results.R
  source(here("codeR", "postprocess-results.R"))
  # Save method-specific results files (run the post-processing function)
  post_process(
    sim_dt = sim_dt,
    is_belloni = as.logical(command_args[1]),
    n_dgp = command_args[2]
  )
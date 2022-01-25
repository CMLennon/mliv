library(pacman)
p_load(tidymodels, tidyverse, data.table, here, 
       MASS, glmnet, xgboost, ranger, future, future.apply, caret)

command_args = commandArgs(trailingOnly = T)
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

# All methods
methods_all = c("dgp", 'ovb', "2sls", "lasso", "plasso", "rf", "rfcv", "pca", "gboost", "split", "classo")

# Methods to run
methods_run = setdiff(methods_all, tolower(tail(command_args, -4)))

#pull argument 3 (number of observations) into the simulation
numobs = as.numeric(command_args[3])
numohbs = 1000
#pull argument 4 (number of simulations) into full workspace
numsims = as.numeric(command_args[4])
numsims = 1000
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
      iter_dt = lm(y ~ x1 + vn1, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model")  %>%mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1 + i_dt_corr$vn1)) %>% mutate(cov = cov(i_dt_corr$x1 + i_dt_corr$vn1, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt)
    }
    else{
      iter_dt = lm(y ~ x1 + x2, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model") %>%  mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% mutate(true_var = var(i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1 + i_dt_corr$x2)) %>% mutate(cov = i_dt_corr$x1 + i_dt_corr$x2, i_dt_corr$y - i_dt_corr$x1) %>% bind_rows(.,iter_dt)
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

    if('classo' %in% methods_run){
    source(here("codeR", "LassoCVCaret.R"))
    i_dt_corr[,x1_hat_lasso_cv_caret := lasso_fit_cv_c(dt = i_dt_corr,
                                         folds = folds)]
    
    iter_dt = lm(y ~ x1_hat_lasso_cv_caret , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: LASSO selection caret") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_hat_lasso_cv_caret)) %>% mutate(cov = cov(i_dt_corr$x1_hat_lasso_cv_caret, i_dt_corr$y - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
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
#print(sim_dt)
#print median beta values to console
sim_dt[, mean(estimate, na.rm = TRUE), by = model]

write_csv(sim_dt, '/Users/connor/Downloads/test.csv')
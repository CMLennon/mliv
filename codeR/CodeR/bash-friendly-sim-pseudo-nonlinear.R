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

#New Argument:
#-----------------------------------------
#   - i[5] is T (True)/F (False) and instructs script to build/rebuild data from scratch or not
#-----------------------------------------

#   - tail(i, -5) gives methods to NOT run
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
#       - "nnetf": Neural net, full
#       - "nnetw": Neural net, wide
#       - "nnetd": Neural net, deep

#load necessary packages
library(pacman)
p_load(tidymodels, tidyverse, data.table, here, 
       MASS, glmnet, xgboost, ranger, future, future.apply, plotly, SteinIV, ivmodel, pbapply, R.utils, magrittr, tensorflow)

#install_tensorflow(version = 'nightly-gpu')

# Not in function
`%nin%` = negate(`%in%`)

# All methods
methods_run = c('nnet1', 'nnet2', 'nnet3')

# Methods to run
#methods_run = setdiff(tolower(methods_all), tolower(tail(command_args, -5)))

#set parameters from 1st command argument
if (command_args[1] == "F") {
  # Source the function that generates the simple data
  #source(here("codeR", "simple_data_function.R"))
  k = 7
  complex = F
  # Define the type of beta distribution (irrelevant for 'simple')
  shuffle = command_args[2]
} else {
  # Source the function that generates the Belloni data
  #source(here("codeR", "Belloni_aux_function_shuffle.R"))
  k = 100
  complex = T
  # Define the type of beta distribution
  shuffle = command_args[2]
}

#build data toggle
build_data = FALSE
if(command_args[5] == "T"|command_args[4]=='True'|command_args[5] == "t"){
    print('Building datasets to run simulations')
    build_data = TRUE
} else if(command_args[5] == 'F'|command_args[5] == 'False'){
  print('Using prebuilt datasets to run simulations')
    build_data = FALSE
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
  #data_gen(n,k,shuffle=shuffle) #old command

  if(complex){
    complex = 'T'
    filename = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '.csv')
    filepath = here('dgp_folder',filename)
    filepath = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath, fixed = TRUE), fixed = TRUE),fixed = TRUE)
    i_dt_corr = fread(filepath)
  } else{
    complex = 'F'
    filename = paste0('dgp_data_', complex, '-', iter, '.csv')
    filepath = here('dgp_folder',filename)
    i_dt_corr = fread(filepath)
  }
  
  
  #folds object for cross validation
  folds = vfold_cv(i_dt_corr, v = 5)
  
  #initialize a data table to create new rows inside of. This holds second
  #stage beta results and an iter number to track intra-iteration variation
  
  iter_dt = data.table()
 if('nnet1' %in% methods_run){
   source(here('codeR', "nnet_final.R"))
   
   filename1 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_1', '.csv')
   filename2 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_1oos', '.csv')
   filepath1 = here('dgp_folder_weird',filename1)
   filepath1 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath1, fixed = TRUE), fixed = TRUE),fixed = TRUE)
   filepath2 = here('dgp_folder_weird',filename2)
   filepath2 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath2, fixed = TRUE), fixed = TRUE),fixed = TRUE)
   

  #  nnet_top5 = best_x_hyp(shuffle = shuffle, complex = complex, deep = FALSE, top_x = 2, numfolds = 5, iter = iter)

  #  nn_best = finetune_cv(top_x_obj = nnet_top5, datatable = i_dt_corr, shuffle = shuffle, iter = iter, folds = folds, complex = complex, numfolds = 5)


   i_dt_corr[,x1_hat_nnetw := final_nnet(datatable = i_dt_corr, folds = folds, 
      width = 32, dropout = .1, depth = 6, 
      epoch_num = 40, batch_size = 10, 
      activationfcn = 'relu', oos = FALSE, foldnum = 5)]
   
   i_dt_corr$fst_stage_preds_is_1 = predict(lm(x1 ~ x1_hat_nnetw, data = i_dt_corr))
   
   fwrite(i_dt_corr, filepath1)


   data_oosw = final_nnet(datatable = i_dt_corr, folds = folds, 
                          width = 32, dropout = .1, depth = 6, 
                         epoch_num = 40, batch_size = 10, 
                         activationfcn = 'relu', oos = TRUE, foldnum = 5) %>% mutate(x1_hat_nnetw_oos = nnet_x1_hat)
   
   data_oosw$fst_stage_preds_oos_1 = predict(lm(x1 ~ x1_hat_nnetw_oos, data = data_oosw))
   
   fwrite(data_oosw, filepath2)
   

  iter_dt = lm(y ~ fst_stage_preds_is_1 , data = i_dt_corr) %>% tidy(quick = T) %>%
    filter(grepl('fst', term)) %>%
    mutate(model = "First stage: Deep Neural Net norm")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
    mutate(var = var(i_dt_corr$fst_stage_preds_is_1)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds_is_1, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds_is_1, i_dt_corr$x1 - i_dt_corr$fst_stage_preds_is_1)) %>% bind_rows(.,iter_dt) %>% 
    data.table()

   iter_dt = lm(y ~ fst_stage_preds_oos_1, data = data_oosw) %>% tidy(quick = T) %>%
     filter(grepl('fst', term)) %>%
     mutate(model = "First stage: Deep Neural Net norm OOS")%>% mutate(true_var = var(data_oosw$x1)) %>% mutate(true_cov = cov(data_oosw$x1, data_oosw$y - data_oosw$x1)) %>%
     mutate(var = var(data_oosw$fst_stage_preds_oos_1)) %>% mutate(cov_u = cov(data_oosw$fst_stage_preds_oos_1, data_oosw$y - data_oosw$x1), cov_e = cov(data_oosw$fst_stage_preds_oos_1, data_oosw$x1 - data_oosw$fst_stage_preds_oos_1)) %>% bind_rows(.,iter_dt) %>% 
     data.table()
 }
  
   if('nnet2' %in% methods_run){
   source(here('codeR', "nnet_final.R"))
     
     filename1 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_2', '.csv')
     filename2 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_2oos', '.csv')
     filepath1 = here('dgp_folder_weird',filename1)
     filepath1 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath1, fixed = TRUE), fixed = TRUE),fixed = TRUE)
     filepath2 = here('dgp_folder_weird',filename2)
     filepath2 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath2, fixed = TRUE), fixed = TRUE),fixed = TRUE)

     i_dt_corr[,x1_hat_nnetd := final_nnet(datatable = i_dt_corr, folds = folds, 
                                           width = 32, dropout = .1, depth = 6, 
                                           epoch_num = 40, batch_size = 10, 
                                           activationfcn = 'relu', oos = FALSE, foldnum = 5, z_prob = TRUE)]
     i_dt_corr = i_dt_corr %>% mutate(z_prob= ifelse(abs(sin(z_1)) > .995, round(sin(z_1)/10,2), 0))
     
     i_dt_corr %<>% mutate(x1 = x1 + z_prob, y = y + 2*z_prob)
     i_dt_corr$fst_stage_preds_is_2 = predict(lm(x1 ~ x1_hat_nnetw, data = i_dt_corr))
     
     data_oosd = final_nnet(datatable = i_dt_corr, folds = folds, 
                            width = 32, dropout = .1, depth = 6, 
                            epoch_num = 40, batch_size = 10, 
                            activationfcn = 'relu', oos = TRUE, foldnum = 5, z_prob = TRUE) %>% mutate(x1_hat_nnetd_oos = nnet_x1_hat)
     
     data_oosd$fst_stage_preds_oos_2 = predict(lm(x1 ~ x1_hat_nnetd_oos, data = data_oosd))
     
     fwrite(i_dt_corr, filepath1)
     fwrite(data_oosd, filepath2)
     
     iter_dt = lm(y ~ fst_stage_preds_is_2, data = i_dt_corr) %>% tidy(quick = T) %>%
       filter(grepl('fst', term)) %>%
       mutate(model = "First stage: Deep Neural Net norm")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
       mutate(var = var(i_dt_corr$fst_stage_preds_is_2)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds_is_2, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds_is_2, i_dt_corr$x1 - i_dt_corr$fst_stage_preds_is_2)) %>% bind_rows(.,iter_dt) %>% 
       data.table()
     
     iter_dt = lm(y ~ fst_stage_preds_oos_2, data = data_oosd) %>% tidy(quick = T) %>%
       filter(grepl('fst', term)) %>%
       mutate(model = "First stage: Deep Neural Net norm OOS")%>% mutate(true_var = var(data_oosd$x1)) %>% mutate(true_cov = cov(data_oosd$x1, data_oosw$y - data_oosd$x1)) %>%
       mutate(var = var(data_oosd$fst_stage_preds_oos_2)) %>% mutate(cov_u = cov(data_oosd$fst_stage_preds_oos_2, data_oosw$y - data_oosw$x1), cov_e = cov(data_oosd$fst_stage_preds_oos_2, data_oosd$x1 - data_oosd$fst_stage_preds_oos_2)) %>% bind_rows(.,iter_dt) %>% 
       data.table()
     
     i_dt_corr[,x1_oracle := true_x - vn1]
     vncoef = summary(lm(y ~ true_x + vn1 + z_prob, data = data_oosd))$coefficients[3,1]
     
     iter_dt = lm(y ~ true_x + vn1 + z_prob, data = data_oosd) %>% tidy(quick = T) %>%
       filter(grepl("x1", term)) %>%
       mutate(model = "First Stage: Oracle Model")  %>%mutate(true_var = var(data_oosd$true_x)) %>% mutate(true_cov = cov(data_oosd$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
       mutate(var = var(i_dt_corr$x1_oracle)) %>% mutate(cov_u = cov(i_dt_corr$x1_oracle, i_dt_corr$y - i_dt_corr$x1 - vncoef*i_dt_corr$vn1), cov_e = cov(i_dt_corr$x1_oracle, i_dt_corr$x1 - i_dt_corr$x1_oracle)) %>% bind_rows(.,iter_dt)
     
     source(here('codeR', "SplitSampleIV.R"))
     i_dt_corr[,x1_hat_ussiv := ssiv_fcn(datatable = i_dt_corr,
                                         folds = folds,
                                         unbiased = T)]
     iter_dt = lm(y ~ x1_hat_ussiv , data = i_dt_corr) %>% tidy(quick = T) %>%
       filter(grepl('x1', term)) %>%
       mutate(model = "First stage: USSIV-z2")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
       mutate(var = var(i_dt_corr$x1_hat_ussiv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$x1 - i_dt_corr$x1_hat_ussiv)) %>% bind_rows(.,iter_dt) %>% 
       data.table()
     
     source(here("codeR", "PostLassoCVParsnip.R"))
     i_dt_corr[,x1_hat_plasso_cv := plasso_fit_cv(datatable = i_dt_corr,
                                                  folds = folds)]
     
     
     iter_dt = lm(y ~ x1_hat_plasso_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
       filter(grepl('x1', term)) %>%
       mutate(model = "First stage: post-LASSO selection-z2") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
       mutate(var = var(i_dt_corr$x1_hat_plasso_cv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$x1 - i_dt_corr$x1_hat_plasso_cv)) %>% bind_rows(.,iter_dt) %>% 
       data.table()
   }
  
  if('nnet3' %in% methods_run){
    source(here('codeR', "nnet_final.R"))
    
    source(here('codeR', 'Belloni_aux_function_shuffle.R'))
    
    filename1 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_3', '.csv')
    filename2 = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '_3oos', '.csv')
    filepath1 = here('dgp_folder_weird',filename1)
    filepath1 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath1, fixed = TRUE), fixed = TRUE),fixed = TRUE)
    filepath2 = here('dgp_folder_weird',filename2)
    filepath2 = gsub("\\s+", "\\\\ ", gsub("\\(","\\\\(",gsub("\\)", "\\\\)", filepath2, fixed = TRUE), fixed = TRUE),fixed = TRUE)
    
    i_dt_corr = data_gen(n = 1000, k = 100, shuffle = '1', z_prob = TRUE)
    
    folds_z = vfold_cv(i_dt_corr, v = 5)
    
    i_dt_corr[,x1_hat_nnetd := final_nnet(datatable = i_dt_corr, folds = folds_z, 
                                          width = 32, dropout = .1, depth = 6, 
                                          epoch_num = 40, batch_size = 10, 
                                          activationfcn = 'relu', oos = FALSE, foldnum = 5, z_prob = FALSE)]
    
    i_dt_corr$fst_stage_preds_is_2 = predict(lm(x1 ~ x1_hat_nnetd, data = i_dt_corr))
    
    data_oosd = final_nnet(datatable = i_dt_corr, folds = folds_z, 
                           width = 32, dropout = .1, depth = 6, 
                           epoch_num = 40, batch_size = 10, 
                           activationfcn = 'relu', oos = TRUE, foldnum = 5, z_prob = FALSE) %>% mutate(x1_hat_nnetd_oos = nnet_x1_hat)
    
    data_oosd$fst_stage_preds_oos_2 = predict(lm(x1 ~ x1_hat_nnetd_oos, data = data_oosd))
    
    fwrite(i_dt_corr, filepath1)
    fwrite(data_oosd, filepath2)
    
    iter_dt = lm(y ~ fst_stage_preds_is_2, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('fst', term)) %>%
      mutate(model = "First stage: Deep Neural Net norm")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$fst_stage_preds_is_2)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds_is_2, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds_is_2, i_dt_corr$x1 - i_dt_corr$fst_stage_preds_is_2)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    
    iter_dt = lm(y ~ fst_stage_preds_oos_2, data = data_oosd) %>% tidy(quick = T) %>%
      filter(grepl('fst', term)) %>%
      mutate(model = "First stage: Deep Neural Net norm OOS")%>% mutate(true_var = var(data_oosd$x1)) %>% mutate(true_cov = cov(data_oosd$x1, data_oosw$y - data_oosd$x1)) %>%
      mutate(var = var(data_oosd$fst_stage_preds_oos_2)) %>% mutate(cov_u = cov(data_oosd$fst_stage_preds_oos_2, data_oosw$y - data_oosw$x1), cov_e = cov(data_oosd$fst_stage_preds_oos_2, data_oosd$x1 - data_oosd$fst_stage_preds_oos_2)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    
    i_dt_corr[,x1_oracle := x1 - vn1]
    vncoef = summary(lm(y ~ x1 + vn1, data = i_dt_corr))$coefficients[3,1]
    
    iter_dt = lm(y ~ x1 + vn1, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First Stage: Oracle Model-z3")  %>%mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_oracle)) %>% mutate(cov_u = cov(i_dt_corr$x1_oracle, i_dt_corr$y - i_dt_corr$x1 - vncoef*i_dt_corr$vn1), cov_e = cov(i_dt_corr$x1_oracle, i_dt_corr$x1 - i_dt_corr$x1_oracle)) %>% bind_rows(.,iter_dt)
    
    source(here('codeR', "SplitSampleIV.R"))
    i_dt_corr[,x1_hat_ussiv := ssiv_fcn(datatable = i_dt_corr,
                                        folds = folds_z,
                                        unbiased = T)]
    iter_dt = lm(y ~ x1_hat_ussiv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: USSIV-z3")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_ussiv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$x1 - i_dt_corr$x1_hat_ussiv)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    
    source(here("codeR", "PostLassoCVParsnip.R"))
    i_dt_corr[,x1_hat_plasso_cv := plasso_fit_cv(datatable = i_dt_corr,
                                                 folds = folds_z)]
    
    
    iter_dt = lm(y ~ x1_hat_plasso_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: post-LASSO selection-z3") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_plasso_cv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$x1 - i_dt_corr$x1_hat_plasso_cv)) %>% bind_rows(.,iter_dt) %>% 
      data.table()

  }
    
  iter_dt %>% mutate(iter = iter)
  return(iter_dt)
  }

run_sim <-function(n, n_sims, seed = 12345,
                   future_cores = availableCores(),
                   methods_run,
                   shuffle, k, complex, set_seed = TRUE, build_data = FALSE) {
  #---------------------------------------------------------------
  # Using future, we run a parallelized simulation of all models in 
  # methods_run vector, across k instruments, utilizing β-pattern strat
  # specified in shuffle. Defaults to utilizing all available cores,
  # and running all models. Returns:
  
  # sim_dt, a datatable containing an iteration-labeled set of second
  # stage estimates of β-hat.
  #--------------------------------------------------------------

  #-----------------------------------------
  ## IMPORTANT ##
  #-----------------------------------------
  #set to sequential to avoid vmemory issues with nnet methods
  plan_char = 'sequential'
  #-----------------------------------------
  # Set the seed if set_seed command is true
  if(set_seed){
    set.seed(seed)
  }
  # Run a parallelized simulation

  #set up future_seed to follow directions from set_seed
  if(set_seed){
    fseed =seed
  } else{fseed = NULL}
  #-----------------------------------------
  # Build data in advance
  #-----------------------------------------

  if(build_data){
  # Source the function that generates simple and all three complex datasets
  
  dgp_writer = function(iter, n = n){
  #-----------------------------------------
    # Utility function to write a set of dgps to disk for reading by run-sim.

    #Keeps data identical across methods for a given seed.
  #-----------------------------------------
    source(here("codeR", "simple_data_function.R"))
    k = 7
    complex = 'F'
    i_dt_corr = data_gen(n,k,shuffle='1')
    filename = paste0('dgp_data_', complex, '-', iter, '.csv')
    fwrite(i_dt_corr,here('dgp_folder',filename))
    
    source(here("codeR", "Belloni_aux_function_shuffle.R"))

    k = 100
    complex = 'T'
    shuffle = '1'

    i_dt_corr = data_gen(n,k,shuffle=shuffle)
    filename = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '.csv')
    fwrite(i_dt_corr,here('dgp_folder',filename))

    shuffle = '2'

    i_dt_corr = data_gen(n,k,shuffle=shuffle)
    filename = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '.csv')
    fwrite(i_dt_corr,here('dgp_folder',filename))

    shuffle = '3'

    i_dt_corr = data_gen(n,k,shuffle=shuffle)
    filename = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '.csv')
    fwrite(i_dt_corr,here('dgp_folder',filename))
 
}
  if(n_sims <= 1000){
    print('Progress on generating new data:')
    pblapply(1:1000, dgp_writer, n = n)
  } else{
    print('Progress on generating new data:')
    pblapply(1:n_sims, dgp_writer, n = n)
  }
}
  #-----------------------------------------


  future::plan(plan_char)


  future_lapply(X = 1:n_sims, future.seed = fseed, FUN = one_iter, n = n, 
                methods_run = methods_run,
                shuffle=shuffle, k = k, complex) %>% rbindlist()
}


#run simulation
sim_dt <- run_sim(n = 1000, n_sims =1000, methods_run = methods_run,
                  shuffle = shuffle, k = k, complex = complex, build_data = build_data)

#print median beta values to console
sim_dt[, mean(estimate, na.rm = TRUE), by = model]

#write to a file with the following file-naming conventions:

# sim-results-(complex dgp T/F)-(shuffle_strategy == (1,2,3))-(models-run).csv

# Save results ---------------------------------------------------------------------------
  # Source post-processing function from codeR/postprocess-results.R
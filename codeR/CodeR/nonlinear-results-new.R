# Load all results -----------------------------------------------------------------------
# Read in the results, grouping DGPs into separate objects
library(pacman)
p_load(tidymodels, tidyverse, data.table, here, 
       MASS, glmnet, xgboost, ranger, future, future.apply, plotly, SteinIV, ivmodel, pbapply, R.utils, magrittr, matrixStats)

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

#New Arguments:
#-----------------------------------------
#   - i[5] is T (True)/F (False) and instructs script to build/rebuild data from scratch or not
#   - i[6] is 1-4. Determines type of data-disruption to do.
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


# Not in function
`%nin%` = negate(`%in%`)

# All methods
methods_all = c("dgp", 'ovb', "2sls", "lasso", "plasso","rfcv", "pca", "gboost", "split", "liml", 'jive', 'nnetw')

# Methods to run
methods_run = setdiff(tolower(methods_all), tolower(tail(command_args, -6)))

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
#complex = T
#build data toggle
build_data = FALSE
if(command_args[5] == "T"|command_args[5]=='True'|command_args[5] == "t"){
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

#pull argument 5 (data_transform) into full workspace
data_transform = as.numeric(command_args[6])
#data_transform = 1

one_iter = function(iter,shuffle,n,k, data_transform = data_transform, complex, threshold = .95, pow_degree = NULL, int_degree = NULL, methods_run) {
 #This function reads maimaiin a dataset, runs regressions using IV, SSIV, oracle model, neural network predictions
 #and OOS NN predictions, alongside the SSML approach. These get put into a dataframe
  #----------------------------------------
  ## INPUTS ##
  #Complex - if true for this function, maps to Belloni et al.
             #if FALSE - 7 instrument case is used.We use false for data_transform 2 and 3.
  
  #shuffle - 1 or 2. Takes a string that changes how the strength of instruments corresponds to their power
    # 1. - shuffled strength of coefficients (shuffled beta pattern)
    # 2. - strongest instrument at 1 and declining to 100.
  
  #iter - gives iteration of simulation we are on
  
  #data_transform - this can take on one of three values
    #1. - Slight learnable disturbance. This adds the noise term if(abs(sin(z_1)) > .995), sign(sin(z_1))*.1, else, 0 twice over to y and once to x
    #2. - Quadratic form of bad-variation. Requires a 'pow_degree' value as well
    #3. - Interaction form of bad-variation. Requires a 'int_degree' value as well
  
    #5. - Interaction and subinteractions from bad interaction
  
    #numeric - this maps to the dgp where x = f(z) + e => z_1 + z_2 + z_3 +... z_i^(data_transform) + e and y = x + z_i^data_transform + v
   
  if(complex=='t'|complex =='T'|complex == TRUE){
    complex = 'T'
    filename = paste0('dgp_data_', complex, '-', shuffle, '-', iter, '.csv')
    filepath = here('dgp_folder',filename)
    i_dt_corr = fread(filepath)
  } else {
    complex = 'F'
    filename = paste0('dgp_data_', complex, '-', iter, '.csv')
    filepath = here('dgp_folder',filename)
    i_dt_corr = fread(filepath)
  }
  
  zs = grep("z", names(i_dt_corr), value = "T")
  if(data_transform == 1){
    i_dt_corr = i_dt_corr %>% mutate(e_prob = ifelse(abs(sin(6*z_1)) > threshold, .1*round(abs(sin(2*z_1))), 0), e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                     x1 = x1 + e_prob, y = y + 4*abs(e_prob))
    
  } 
  if(data_transform == 2){
    if(pow_degree == 0){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1^pow_degree,
                                       x1 = x1 + e_prob, y = y + 2*e_prob)
    }
    else{
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1^pow_degree, e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 + e_prob, y = y + 2*e_prob)
    }
  }
  
  if(data_transform == 3){
    coef = 3
    if(int_degree == 1){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1, e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 + coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 2){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2, e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 + coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 3){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2*z_3, e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 +coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 4){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2*z_3*z_4,e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 + coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 5){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2*z_3*z_4*z_5,e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 +coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 6){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2*z_3*z_4*z_5*z_6,e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 +coef*e_prob, y = y + 2*coef*e_prob)
    }
    if(int_degree == 7){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = z_1*z_2*z_3*z_4*z_5*z_6*z_7,e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 +coef*e_prob, y = y + 2*coef*e_prob)
    }
  }
    
    if(data_transform == 4){
      i_dt_corr = i_dt_corr %>% mutate(e_prob = ifelse(abs(sin(z_1 + z_2 + z_3 + z_4)) > threshold, .1*round(sin(z_1 + z_2 + z_3 + z_4)), 0),
                                       e_prob = (e_prob - mean(e_prob))/sd(e_prob),
                                       x1 = x1 + e_prob, y = y + 2*e_prob)
    }
  
  #interactions up to order provided in iterations
  if(data_transform == 5){
    
    #separate out the instruments from the main observations, preserving ordering.
    i_dt_corr %<>% data.frame(i_dt_corr)
    z = i_dt_corr[,colnames(i_dt_corr)[grepl(pattern = 'z', colnames(i_dt_corr))]] %>% as.matrix()
    
    product_matrix = lapply(X =2:7, FUN = function(n_inter){
      #all n-iter combos
      combinations = combn(1:7, m =n_inter) %>% t()
      #transform into vectors
      
      combinations = lapply(1:nrow(combinations), . %>% combinations[.,])
      
      #find numeric value of the row product
      
      interaction_terms = lapply(X = seq_along(combinations),
                                 FUN= function(i) z[,combinations[[i]]] %>% matrixStats::rowProds()
      ) %>% unlist() %>%matrix(ncol =length(combinations)) %>% rowSums()
      # Return
      return(interaction_terms)
    }) %>% unlist() %>% matrix(ncol = 6)
    
    #normalize the matrix
   norm_product_mat = lapply(1:6, function(k) {
     
     #need to except the 1 case because 1 col of a matrix is a vector in R
     if(k ==1){
       k_mat = product_matrix[,1:k]
     } else{
       k_mat = product_matrix[,1:k] %>% rowSums()
     }
     #actually normalize
     k_mat = (k_mat - mean(k_mat))/sd(k_mat)
     # return obj
   }) %>% unlist() %>% matrix(ncol = 6)
   
   #add relevant column to i_dt_corr dgp matrix
   i_dt_corr$e_prob = norm_product_mat[,(int_degree -1)]
   
   i_dt_corr %<>% mutate(x1 = x1 + e_prob, y = y + 2*e_prob) %>% data.table()
    
  }
  
  folds = vfold_cv(i_dt_corr, 5)
  
  #data_table_os = dir(here("dgp_folder_weird"),paste0("T-",shuffle, '-', iter, '_', data_transform,'oos.csv'), full.names = T) %>% 
  #  map_dfr(fread)
  
  iter_dt = data.table()
  #print(paste0('i_dt_corr is found? ', !is.null(i_dt_corr)))
  if('dgp' %in% methods_run){
    if(complex==TRUE|complex == 'T'|complex == 't'){
      vncoef = summary(lm(y ~ x1 + vn1 + e_prob, data = i_dt_corr))$coefficients[3,1]
      i_dt_corr[,x1_oracle := true_x]
      iter_dt = lm(y ~ true_x, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model")  %>%mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$x1_oracle)) %>% mutate(cov_u = cov(i_dt_corr$x1_oracle, i_dt_corr$y - i_dt_corr$x1 - vncoef*i_dt_corr$vn1), cov_e = cov(i_dt_corr$x1_oracle, i_dt_corr$x1 - i_dt_corr$x1_oracle)) %>% bind_rows(.,iter_dt) %>% data.table()
    }
    else{
      iter_dt = lm(y ~ x1 + x2 +e_prob, data = i_dt_corr) %>% tidy(quick = T) %>%
        filter(grepl("x1", term)) %>%
        mutate(model = "First Stage: Oracle Model") %>%  mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% mutate(true_var = var(i_dt_corr$x1)) %>%
        mutate(var = var(i_dt_corr$true_x)) %>% mutate(cov_u = cov(i_dt_corr$true_x, i_dt_corr$y - i_dt_corr$x1 - i_dt_corr$x2 -i_dt_corr$e_prob), cov_e = cov(i_dt_corr$x1, i_dt_corr$x1 - i_dt_corr$true_x)) %>% bind_rows(.,iter_dt)%>% data.table()
    }
  }
  
  if('nnetw' %in% methods_run){
    #print('starting nnet')
    source(here('codeR', "nnet_final.R"))
    
    i_dt_corr$x1_hat_nnetf_oos = final_nnet(datatable = i_dt_corr, folds = folds, 
                                            width = 64, dropout = .1, depth = 7, 
                                            epoch_num = 100, batch_size = 1, 
                                            activationfcn = 'relu', oos = TRUE, foldnum = 5)
    
    i_dt_corr$fst_stage_preds = lm(x1 ~ x1_hat_nnetf_oos, data= i_dt_corr) %>% predict(newdata = i_dt_corr)
    
    iter_dt = lm(y ~ x1_hat_nnetf_oos, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: Unrestricted Neural Net OOS")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_nnetf_oos)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_nnetf_oos, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_nnetf, i_dt_corr$x1 - i_dt_corr$x1_hat_nnetf)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
    
    iter_dt = lm(y ~ fst_stage_preds, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('fst', term)) %>%
      mutate(model = "First stage: Unrestricted Neural Net norm OOS")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$fst_stage_preds)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds, i_dt_corr$x1 - i_dt_corr$fst_stage_preds)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  if('ovb' %in% methods_run){
    #print('starting ovb')
    iter_dt = lm(y ~ x1, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First Stage: Naive OLS") %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>% mutate(true_var = var(i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1)) %>% mutate(cov_u = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1, i_dt_corr$x1 - i_dt_corr$x1)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  if('2sls' %in% methods_run){
    #print('starting 2sls')
    zs = grep("z", names(i_dt_corr), value = "T")
    
    formula = as.formula(paste0("x1 ~ ", paste0(zs, collapse = '+')))
    
    i_dt_corr[,x1_hat_ols := lm(formula, data = i_dt_corr) %>% predict(new_data = i_dt_corr)]  
    
    iter_dt = lm(y ~ x1_hat_ols , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: 2SLS") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_ols)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_ols, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_ols, i_dt_corr$x1 - i_dt_corr$x1_hat_ols)) %>%  bind_rows(.,iter_dt) %>% data.table()
  }
  
  if('rfcv' %in% methods_run){
    #print('starting randforest')
    source(here("codeR", "RandomForestCVParsnip.R"))
    i_dt_corr$x1_hat_rfcv = forest_fcn(datatable = i_dt_corr,
                                          folds = folds,
                                          cv = F,
                                          oos = T)
    
    iter_dt = lm(y ~ x1_hat_rfcv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: Random Forest, CV") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov( i_dt_corr$x1,  i_dt_corr$y -  i_dt_corr$x1)) %>%
      mutate(var = var( i_dt_corr$x1_hat_rfcv)) %>% mutate(cov_u = cov( i_dt_corr$x1_hat_rfcv,  i_dt_corr$y -  i_dt_corr$x1), cov_e = cov( i_dt_corr$x1_hat_rfcv,  i_dt_corr$x1 -  i_dt_corr$x1_hat_rfcv)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
    
    i_dt_corr$fst_stage_preds_rfcv = lm(data = i_dt_corr, x1 ~ x1_hat_rfcv) %>% predict(new_data = i_dt_corr)
    
    iter_dt = lm(y ~ fst_stage_preds_rfcv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("fst", term)) %>%
      mutate(model = "First stage: Random Forest norm, CV") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov( i_dt_corr$x1,  i_dt_corr$y -  i_dt_corr$x1)) %>%
      mutate(var = var( i_dt_corr$fst_stage_preds_rfcv)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds_rfcv,  i_dt_corr$y -  i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds_rfcv,  i_dt_corr$x1 -  i_dt_corr$fst_stage_preds_rfcv)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  if('gboost' %in% methods_run){
    source(here('codeR', 'BoostedTreesParsnip.R'))
    i_dt_corr[,x1_hat_boost_cv := xgbst_fit(datatable = i_dt_corr,
                                            folds = folds)]
    
    iter_dt = lm(y ~ x1_hat_boost_cv, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: Boosted Trees")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_boost_cv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_boost_cv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_boost_cv, i_dt_corr$x1 - i_dt_corr$x1_hat_boost_cv)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
    
    i_dt_corr$fst_stage_preds_boost = lm(data = i_dt_corr, x1 ~ x1_hat_boost_cv) %>% predict(new_data = i_dt_corr)
    
    iter_dt = lm(y ~ fst_stage_preds_boost, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("fst", term)) %>%
      mutate(model = "First stage: Random Forest norm, CV") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov( i_dt_corr$x1,  i_dt_corr$y -  i_dt_corr$x1)) %>%
      mutate(var = var( i_dt_corr$fst_stage_preds_boost)) %>% mutate(cov_u = cov( i_dt_corr$fst_stage_preds_boost,  i_dt_corr$y -  i_dt_corr$x1), cov_e = cov( i_dt_corr$fst_stage_preds_boost,  i_dt_corr$x1 -  i_dt_corr$fst_stage_preds_boost)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  if('pca' %in% methods_run){
    #print('starting pca')
    source(here("codeR", "PCAParsnipCV.R"))
    i_dt_corr[,x1_hat_pca := pca_fit_cv(datatable = i_dt_corr,
                                        folds = folds,
                                        cv = T)]
    
    iter_dt = lm(y ~ x1_hat_pca , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: PCA") %>%  mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_pca)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_pca, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_pca, i_dt_corr$x1 - i_dt_corr$x1_hat_pca)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  if('lasso' %in% methods_run){
    #print('starting lasso')
    source(here("codeR", "LassoCVParsnip.R"))
    i_dt_corr$x1_hat_lasso = lasso_fit_cv(datatable = i_dt_corr,
                                               folds = folds, oos = T)
    
    iter_dt = lm(y ~ x1_hat_lasso , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: LASSO selection") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_lasso)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_lasso, i_dt_corr$y -i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_lasso, i_dt_corr$x1 - i_dt_corr$x1_hat_lasso)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
    
    i_dt_corr$fst_stage_preds_lasso = lm(data = i_dt_corr, x1 ~ x1_hat_lasso) %>% predict(new_data = i_dt_corr)
    
    iter_dt = lm(y ~ fst_stage_preds_lasso, data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('fst', term)) %>%
      mutate(model = "First stage: norm LASSO selection") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$fst_stage_preds_lasso)) %>% mutate(cov_u = cov(i_dt_corr$fst_stage_preds_lasso, i_dt_corr$y -i_dt_corr$x1), cov_e = cov(i_dt_corr$fst_stage_preds_lasso, i_dt_corr$x1 - i_dt_corr$fst_stage_preds_lasso)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  
  #Second Stage: Post-Lasso Instruments
  if('plasso' %in% methods_run){
    #print('starting postlasso')
    source(here("codeR", "PostLassoCVParsnip.R"))
    i_dt_corr[,x1_hat_plasso_cv := plasso_fit_cv(datatable = i_dt_corr,
                                                 folds = folds, oos = T, plugin = F)]
    
    
    iter_dt = lm(y ~ x1_hat_plasso_cv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: post-LASSO selection") %>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_plasso_cv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_plasso_cv, i_dt_corr$x1 - i_dt_corr$x1_hat_plasso_cv)) %>% bind_rows(.,iter_dt) %>% 
      data.table()
  }
  
  if('split' %in% methods_run){
    #print('starting ssiv')
    source(here('codeR', "SplitSampleIV.R"))
    i_dt_corr[,x1_hat_ussiv := ssiv_fcn(datatable = i_dt_corr,
                                        folds = folds,
                                        unbiased = T)]
    iter_dt = lm(y ~ x1_hat_ussiv , data = i_dt_corr) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: USSIV")%>% mutate(true_var = var(i_dt_corr$x1)) %>% mutate(true_cov = cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      mutate(var = var(i_dt_corr$x1_hat_ussiv)) %>% mutate(cov_u = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$y - i_dt_corr$x1), cov_e = cov(i_dt_corr$x1_hat_ussiv, i_dt_corr$x1 - i_dt_corr$x1_hat_ussiv)) %>% bind_rows(.,iter_dt) %>% 
      data.table() 
  }
  if('liml' %in% methods_run){
    #print('starting liml')
    #LIML estimation here is using the IVModels package. There is a handwritten version as well that we can use should we want it,
    # especially if we need some form of X-hat.
    source(here('codeR', "LIMLest.R"))
    i_dt_corr[,x1_hat_liml := NA]
    LIMLobj = LIML_fit(datatable = i_dt_corr, folds = folds, iter = iter)
    iter_dt = data.frame(term = 'x1_LIML') %>% mutate(estimate = LIMLobj$point.est[1], std.error = LIMLobj$std.err[1], 
                                                      statistic =LIMLobj$test.stat[1], p.value = LIMLobj$p.value[1]) %>%
      mutate(model = "First stage: LIML (Fuller)") %>% mutate(var = NA, cov_u = NA,cov_e = NA, true_var = var(i_dt_corr$x1), true_cov =  cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1)) %>%
      bind_rows(.,iter_dt) %>% data.table() 
  }
  if('jive' %in% methods_run){
    #print('starting jive')
    source(here('codeR', 'new_jive.R'))
    #grep instruments
    zs = grep("z", names(i_dt_corr), value = "T")
    JIVEobj = new_jive(y = i_dt_corr$y %>% as.matrix(), X = i_dt_corr$x1 %>% as.matrix(), Z = i_dt_corr[,..zs] %>% as.matrix())
    
    i_dt_corr[,x1_hat_jive := JIVEobj$x_hat]
    iter_dt = data.frame(term = 'x1_JIVE') %>% mutate(estimate = JIVEobj$beta[1]) %>%
      mutate(model = "First stage: JIVE") %>% mutate(var = var(i_dt_corr$x1_hat_jive), 
                                                     cov_u = cov(i_dt_corr$x1_hat_jive, i_dt_corr$y - i_dt_corr$x1), true_var = var(i_dt_corr$x1), 
                                                     true_cov =  cov(i_dt_corr$x1, i_dt_corr$y - i_dt_corr$x1),
                                                     cov_e = cov(i_dt_corr$x1_hat_jive, i_dt_corr$x1 - i_dt_corr$x1_hat_jive)) %>%
      bind_rows(.,iter_dt) %>% data.table() 
  }
  #print('starting return')
  return(iter_dt %>% mutate(iter = iter, pow_degree = pow_degree, int_degree = int_degree, data_transform = data_transform) %>% data.table())
}


run_sim <-function(n, n_sims, seed = 12345,
                   future_cores = availableCores(),
                   methods_run,
                   shuffle, k, threshold, data_transform, complex, set_seed = TRUE, build_data = FALSE,int_degree=1,pow_degree = 0
                   ) {
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
  if(any(grepl('nnet', methods_run))){
    plan_char = 'sequential'
    future::plan(plan_char)
  } else{
    plan_char = 'multicore'
    future::plan(plan_char, workers = future_cores)
  }
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
  
  if(data_transform ==1){
  
  return(future_lapply(X = 1:n_sims,  FUN = one_iter,future.seed = TRUE, n = n, 
                methods_run = methods_run,
                shuffle=shuffle, k = k, complex=complex,data_transform = data_transform, 
                pow_degree = pow_degree, int_degree = int_degree) %>% rbindlist() %>% data.table())
    
  }
  if(data_transform == 2){
    pow_grid = expand_grid(iter = 1:n_sims,pow_degree = seq(-4,4,by =1))
    return(future_mapply(one_iter, iter = pow_grid$iter, pow_degree = pow_grid$pow_degree,MoreArgs = list(n = n, 
                  methods_run = methods_run,
                  shuffle=shuffle, k = k, complex=complex,data_transform = data_transform, int_degree = int_degree), SIMPLIFY = FALSE,future.seed = TRUE) %>% rbindlist()%>% data.table())
  }
  
  if(data_transform == 3){
    int_grid = expand_grid(iter = 1:n_sims,int_degree = seq(1,7,by = 1))
    return(future_mapply(FUN = one_iter,iter = int_grid$iter, int_degree = int_grid$int_degree, MoreArgs=list(n = n, 
                  methods_run = methods_run,
                  shuffle=shuffle, k = k, complex=complex,data_transform = data_transform, 
                  pow_degree = pow_degree), SIMPLIFY = FALSE, future.seed = TRUE) %>% rbindlist()%>% data.table())
  }
  
  if(data_transform ==4){
    return(future_lapply(X = 1:n_sims, future.seed = TRUE, FUN = one_iter, n = n, 
                  methods_run = methods_run,
                  shuffle=shuffle, k = k, complex=complex,data_transform = data_transform, 
                  pow_degree = pow_degree, int_degree = int_degree) %>% rbindlist()%>% data.table())
    
  }
  
  if(data_transform == 5){
    int_grid = expand_grid(iter = 1:n_sims,int_degree = seq(2,7,by = 1))
    return(future_mapply(FUN = one_iter,iter = int_grid$iter, int_degree = int_grid$int_degree, MoreArgs=list(n = n, 
                           methods_run = methods_run,
                           shuffle=shuffle, k = k, complex=complex,data_transform = data_transform, 
                           pow_degree = pow_degree), SIMPLIFY = FALSE, future.seed = TRUE) %>% rbindlist()%>% data.table())
    
  }
}


#run simulation
sim_dt <- run_sim(n = parent.frame()$numobs, n_sims = parent.frame()$numsims, methods_run = parent.frame()$methods_run,
                  shuffle = parent.frame()$shuffle, k = parent.frame()$k, complex = complex, build_data = parent.frame()$build_data, 
                  data_transform =parent.frame()$data_transform) %>% data.table()

#print median beta values to console
if(any(grepl('nnet', methods_run))){
  nnet_run = 'nnet'
  fpath = '/home/connor/Desktop/data-out'
} else{
  nnet_run = 'non-nnet'
  fpath = here('Resources-Nonlinear')
}

fname = paste0('sim-results-nonlin-', complex, '-',shuffle,'-datatrans-',data_transform,'-','cf3', '-', paste0(methods_run, collapse = '-'),'.csv')
fullpath = paste0(fpath, '/',fname)

fwrite(x = sim_dt, paste0(fpath, '/', fname))

sim_dt[, median(estimate, na.rm = TRUE), by = model]


# Change names
sim_dt[, model := str_replace_all(model, "First Stage", "First stage")]

sim_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]

sim_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]

sim_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]

sim_dt[, model := str_replace_all(model, "Random Forest", "Random forest")]

sim_dt[, model := str_replace_all(model, "LASSO", "Lasso")]

sim_dt[, model := str_replace_all(model, "post-", "Post-")]

sim_dt[, model := str_replace_all(model, "Trees", "trees")]

sim_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]

sim_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]

sim_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]

# Rename neural nets, OOS
sim_dt[str_detect(model, "Unrestricted") & !str_detect(model, "norm"), model := "First stage: Neural net OOS, unrestricted"]
sim_dt[str_detect(model, "Unrestricted") & str_detect(model, "norm"), model := "First stage: Neural net OOS-linearized, unrestricted"]

# Abbreviation for unrestricted
sim_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]

# Order the models
all_models = c(
  "First stage: Random forest, CV",
  "First stage: Boosted trees",
  "First stage: Neural net OOS, unrest.",
  "First stage: Neural net OOS, narrow",
  "First stage: Neural net OOS, shallow",
  "First stage: Neural net, unrest.",
  "First stage: Neural net, narrow",
  "First stage: Neural net, shallow",
  "First stage: Lasso selection",
  "First stage: Post-Lasso selection",
  "First stage: PCA",
  "Jackknife IV (JIVE)",
  "Split-sample IV",
  "LIML (Fuller)",
  "First stage: OLS",
  "Naive OLS",
  "\"Oracle\" model"
)

sim_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]


# Data work: Drop unwanted columns -------------------------------------------------------
sim_dt[, c("std.error", "statistic", "p.value") := NULL]


# Data work: Covariances -----------------------------------------------------------------
# Rename covariance
setnames(sim_dt, old = "cov_u", new = "cov_xhat_u")

# Transfer over 'cov_u' to 'cov_xhat_u'
sim_dt[is.na(cov_xhat_u), cov_xhat_u := cov_u]

# Drop 'cov_u'
sim_dt[, cov_u := NULL]

# Calculate the covariance between e and 
beta = 1
sim_dt[, cov_xhat_e := (var * (estimate - beta) - cov_xhat_u) / beta]

# Drop 'cov_e' (missing for most models; calculations match output)
sim_dt[, cov_e := NULL]

# Calculate the bias (error for beta)
sim_dt[, beta_bias := estimate - beta]







# # Figures: Coefficient densities ---------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = estimate, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab("Estimate") +
#       ylab("Density") +
#       scale_fill_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$estimate %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )

#----------------------------------------------
#produce basic figures.
# Figures: Coefficient densities with 'ggridges' -----------------------------------------
# Create the individual figures

if(data_transform == 1 | data_transform == 4){
  
X= ggplot(
    data = sim_dt,
    aes(x = estimate, y = model)
  ) +
  geom_density_ridges(
    aes(fill = model, color = model, vline_color = model),
    size = 0.2,
    alpha = 0.85, 
    rel_min_height = 0.0005,
    scale = 5,
    vline_size = 0.8,
    quantile_lines = T,
    quantiles = 2
  ) +
  geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
  xlab("Estimate") +
  ylab("") +
  scale_fill_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
  scale_color_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
  scale_discrete_manual(
    "vline_color",
    values = magma(sim_dt$model %>% uniqueN(), end = 0.95, direction = 1),
    name = "Method"
  ) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(legend.position = "none")

ggsave(
  plot = X,
  path = here("Figures"),
  filename = paste0("sim-density-ridge-dt-",data_transform, '-', complex, '-', shuffle, '-',nnet-run, ".pdf"),
  device = cairo_pdf,
  width = 4.5,
  height = 2.5
)
}

#----------------------------------------------
#produce basic figures.
# Figures: Bias plot-line with exponential term in the x-axis-----------------------------------------
# Create the individual figures

if(data_transform == 2){
  
  X= ggplot(
    data = results,
    aes(x = pow_degree, group = model, y = estimate, color = model)
  ) +
    geom_smooth(
      aes(color = model)
    ) +
    geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
    geom_vline(xintercept = 0, size = 0.2)+
    xlab("Degree of Polynomial Noise") +
    ylab("Bias") +
    scale_color_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
    scale_discrete_manual(
      "vline_color",
      values = magma(n = sim_dt$model %>% unique() %>% length(), end = 0.95, direction = 1),
      name = "Method"
    ) +
    theme_minimal(
      base_size = 7.5,
      base_family = "Charter"
    ) +
    theme(legend.position = "none")
  
  ggsave(
    plot = X,
    path = here("Figures"),
    filename = paste0("sim-line-dt-",data_transform, '-', complex, '-', shuffle, '-', 'nnet_run', ".pdf"),
    device = cairo_pdf,
    width = 4.5,
    height = 2.5
  )
}

if(data_transform == 3){
  
  X= ggplot(
    data = sim_dt %>% group_by(model, int_degree) %>% summarize(estimate = mean(estimate-1)),
    aes(x = int_degree, color = model, y = estimate)
  ) +
    geom_line(
      aes(color = model)
    ) +
    xlab("Size of Interaction Term in Noise") +
    ylab("Bias") +
    scale_color_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
    scale_discrete_manual(
      "vline_color",
      values = magma(sim_dt$model %>% uniqueN(), end = 0.95, direction = 1),
      name = "Method"
    ) +
    theme_minimal(
      base_size = 7.5,
      base_family = "Charter"
    )
  
  ggsave(
    plot = X,
    path = here("Figures"),
    filename = paste0("sim-line-dt-",data_transform, '-', complex, '-', shuffle, '-',nnet-run, ".pdf"),
    device = cairo_pdf,
    width = 4.5,
    height = 2.5
  )
}

deprecated == FALSE 
if(deprecated == TRUE){
  lapply(
    X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
    FUN = function(d) { 
      ggplot(
        data = d,
        aes(x = estimate, y = model)
      ) +
        geom_density_ridges(
          aes(fill = model, color = model, vline_color = model),
          size = 0.2,
          alpha = 0.85, 
          rel_min_height = 0.0005,
          scale = 5,
          vline_size = 0.8,
          quantile_lines = T,
          quantiles = 2
        ) +
        geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
        xlab("Estimate") +
        ylab("") +
        scale_fill_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
        scale_color_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
        scale_discrete_manual(
          "vline_color",
          values = magma(d[,model] %>% uniqueN(), end = 0.95, direction = 1),
          name = "Method"
        ) +
        theme_minimal(
          base_size = 7.5,
          base_family = "Charter"
        ) +
        theme(legend.position = "none")
    }
  )
  
  for (i in seq_along(figure_list)) {
    ggsave(
      plot = figure_list[[i]],
      path = here("Figures"),
      filename = paste0("sim-density-ridge-", names(figure_list)[i], ".pdf"),
      device = cairo_pdf,
      width = 4.5,
      height = 2.5
    )
  }
}

# Save individual plots



# # Two-plot figure (stacked): Comparing non-Belloni and Belloni
# two_gg = (
#   (figure_list$f1 + ggtitle("'Low-complexity' case: Seven strong instruments")) /
#   (figure_list$t2 + ggtitle("'High-complexity' case: 100 mixed-strength instruments"))
# ) +
# # plot_annotation(tag_levels = "A", tag_suffix = ".") +
# plot_layout(guides = "collect") & 
# theme(legend.position = "none") &
# guides(fill = guide_legend(ncol = 3, byrow = T))
# # Three-plot figure (stacked): Comparing Belloni
# three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
# plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
# plot_layout(guides = "collect") & 
# theme(legend.position = "none") &
# xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$estimate %>% range()) &
# guides(fill = guide_legend(ncol = 3, byrow = T))
# # Save the two-plot figure
# ggsave(
#   plot = two_gg,
#   path = here("Figures"),
#   filename = "sim-density-ridge-two.pdf",
#   device = cairo_pdf,
#   height = 7,
#   width = 6.5
# )
# # Save the three-plot figure
# ggsave(
#   plot = three_gg,
#   path = here("Figures"),
#   filename = "sim-density-ridge-three.pdf",
#   device = cairo_pdf,
#   height = 7,
#   width = 6.5
# )


# # Figures: Covariance densities with u ---------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = cov_xhat_u, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 0, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Cov($\\widehat{x}$,u)")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         # values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#         values = magma(d[,model] %>% uniqueN(), end = 0.95, direction = 1)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$cov %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-u-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-u-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Covariance densities with e ---------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = cov_xhat_e, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 0, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Cov($\\widehat{x}$,e)")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme() +
#       coord_cartesian(ylim = c(0,40))
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$cov %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Variance densities ------------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d[str_detect(model, "Oracle", negate = T)],
#         aes(x = var, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Var(\\widehat{x})")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$var %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Bias-share densities ----------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d[(
#           str_detect(model, "Lasso|forest|trees") & str_detect(model, "Post", negate = T)
#         )],
#         aes(x = cov_xhat_e / (cov_xhat_e + beta * cov_xhat_u), fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Var(\\widehat{x})")) +
#       ylab("Density") +
#       # scale_fill_manual(
#       #   "Method",
#       #   values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       # ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$var %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )

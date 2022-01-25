# Todo -----------------------------------------------------------------------------------
#   - Add ML 'best practice' to first stage (i.e., CV): Does it reduce bias?
#   - done. It does - kind of. It makes the bias less predictible. It's possible    a RF with CV may avoid this issue with a bit of luck.
#   - added this feature - had conversation with Ed 10/20. Will add the real intent today, 10/21
#
#     NOTE: Add (not replace)
#   - Compare to PCA (maybe think of correlating the z's?)
#   - added as of 1/19/20. Reduces variance compared to 2SLS and is unbiased. Likely a function of correlation between instruments.

# Setup ----------------------------------------------------------------------------------
# Packages
library(pacman)
p_load(
  tidyverse, broom, data.table,
  lfe, fixest, estimatr, party,
  parallel, magrittr, MASS, partykit, rsample, furrr, future.apply, adabag,
  mboost, promises, elasticnet, tidymodels, glmnet
)
hypchoice <- data.table(
  maxdepth = c(),
  minprob = c(),
  iter = c()
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
  # i_dt[, x2 := 1 + 10 * z^2 + e_common]
  i_dt[, x2 := 1 + e_common]
  # Calculate y
  i_dt[, y := β0 + β1 * x1 + β2 * x2 + runif(n, min = -1, max = 1)]
  return(i_dt)
}

# Functions for K-fold and Cross Validation---------------------------------------------------------------
#function to produce squared errors on holdout fold
num_zs = 7

holdout_results = function(splits, ...) {
  # Fit the model to the majority
  treemod <- party::cforest(..., data = analysis(splits))
  # Save the extra fold
  holdout <- assessment(splits)
  # bind results
  preds = predict(treemod, newdata = holdout)
  res = cbind(holdout, preds = preds)
  # save x1 vals, run preds, calculate sq error
  res$sqerr <- (res$preds - holdout$x1)^2
  # Return the assessment data set with the additional columns
  return(mean(res$sqerr%>% unlist()))
}


#function to run above for all folds, inserting control parameters to ctree as they are searched in grid search

#temporarily until I find a workaround - found. I wasn't using the syntax correctly.
belloni = function(n, k){
set_C <- function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
    #solving for a C value that changes behavior of lasso selectors.
    C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
    return(C)
  }
  #Generate two errors according to Bellloni (pg 26)
  
  #first need covariance and variance set, to produce covar matrix
  musq = 5000
  var_er = 1
  cor_er = .5
  var_inst = 1
  
  #produce Sigma Z. variance is 1, correlation (thus cov) is equal to .5^(|j-h|) for z_j and z_h
  Sig_z = matrix(0,k,k)
  for (j in (1:k)) {
    for (h in (1:k)) {
      Sig_z[c(j),c(h)] = .5^abs(j-h)
    }
  }
  z_n1 = mvrnorm(n = n, mu = rep(0,k), Sigma = Sig_z, empirical = TRUE)
  
  #using exponential
  beta_pattern <-  unlist(lapply(c(0:(k-1)), function(x) (.9)^x), use.names = FALSE)
  #remove one here
  
  
  Sig_z <- as.matrix(Sig_z)
  beta_pattern <- as.vector(beta_pattern)
  
  print("first lin alg step")
  C = set_C(musq=musq, n = n, beta_pattern = beta_pattern,Sig_z = Sig_z)
  C = c(C)
  
  #now that C is known, set variance for error term on X
  
  print("lin alg step 2")
  sigma_v = abs(1 - C*crossprod(beta_pattern,Sig_z) %*% beta_pattern)
  
  #We can also set our full pi matrix to find our true instrument coefficients
  betas <- as.vector(C*beta_pattern)
  cov_er<-cor_er*sqrt(sigma_v)
  
  #use multivariate distribution given matrix above to get a matrix of values for sample
  #size of 'n1'
  covarmat <- matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)
  
  #we now need to get and separate our error terms
  evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
  errors = as.data.frame(evmat_n1)
  names(errors) <- c("e", "v")
  
  #separate errors
  en1 = evmat_n1[,1]
  vn1 = evmat_n1[,2]
  
  #let's construct x's and y's and add them to our z's to create a full data matrix
  

  print("x being created")
  X = z_n1 %*% betas + vn1
  print('x created')
  x_true = X - vn1
  print('x_true')
  Y = X + en1
  print('y created')
  znames = paste(rep("z", k), c(1:k), sep = "_")
  z = as.data.frame(z_n1)
  colnames(z) <- znames
  
  #X (or d, as in Belloni paper) name
  X = as.data.frame(X)
  print('x to df')
  colnames(X) <- "X"
  print('x renamed')
  
  #Target value Y name
  Y = as.data.frame(Y, names = "Y")
  colnames(Y) <- "Y"
  
  #final matrix
  #true_x = X - vn1
  product<-cbind(z,X[,1],Y,errors, x_true)
  product = data.table(product)
  product = setnames(product, c(znames, "x1", "y", "en1", 'vn1', 'true_x'))
  print('product created')
  return(product)
}
bigfolds = function(mtry,ntrees, data, fnum, mod_form = form) {
  #Generate folds:
  
  folds <- data %>% vfold_cv(v=fnum, repeats = 1)
  
  #build a ctree_control object with updated parameters
  control <- cforest_control(mtry = mtry, ntree = ntrees)
  
  #map the list of splits to the holdout reults function
  folds$mse <- folds$splits %>% future_map(~ holdout_results(.x, mod_form, control = control))
  
  #return the set of mses for each fold.
  
  #No longer needed, to avoid cpu overhead moved this to holdout_results. As per future's documentation for efficiency.
  
  #folds$mse <- future_map_dbl(folds$results, function(x) mean(x$sqerr))
  
  #return average mse across all folds
  mse = mean(folds$mse%>%unlist())
  return(mse)
}

#gridsearch function to check across numerous parameter values for minprob (minimum proportion of data in final leaf) and maximum tree depth.
gridsearch = function(dt, fnum, mod_form = form, max_trees = 600, mtry = seq(from = 1, to = 100, by = 20)) {
  #build grid of plausible values for maxdepth and minprob
  paramgrid <- list(mtry = mtry,
                    ntrees = seq(from = 100, to = max_trees, by = 50)) %>%
    cross_df()
  #attach fold count to grid
  paramgrid$fnum = fnum
  #map values of grid to bigfolds function above, return as new variable mse
  paramgrid = paramgrid %>% mutate(mse = future_pmap(paramgrid, bigfolds, data = dt, mod_form = mod_form) %>% unlist())
  #unlist mse
  #paramgrid$mse %<>% unlist()
  #sort by mse
  paramgrid %<>% arrange(mse) %>% head(1)
  return(paramgrid[,1:2])
}

#Caret implementations of lasso/enet and adaboost
adamod = function(dt, k_folds, form) {
  adaGrid <-  expand.grid(mstop = seq(from = 30, to = 1000, by = 200), maxdepth = seq(from = 1, to =100, by= 20))
  adafit <- train(form, 
                  method = 'blackboost',
                  trControl = trainControl(sample, method = 'cv', number = k_folds),
                  data = dt,
                  tuneGrid = adaGrid
  )
  
  return(adafit)
}

enet = function(dt, k_folds, form) {
  enetfit <- train(form, 
                   method = 'glmnet',
                   trControl = trainControl(sample,method = 'cv', number = k_folds),
                   data = dt,
                   tuneGrid = expand.grid(alpha = 1, lambda = seq(from = 0, to = 1, by = .01))
  )
  return(enetfit)
}

ldamod = function(dt, k_folds, form) {
  ldafit <- train(form, 
                  method = 'lda',
                  trControl = trainControl(sample,method = 'cv', number = k_folds),
                  data = dt
  )
  
  return(ldafit)
}

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
  
  forest_tuner = expand_grid(trees = seq(from = 50, to = 300, by = 50), mtry = 2:6, min_n = 2:8)
  
  #build workflow object
  forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)
  
  #find best CVd parameters
  best = tune_grid(forest_wf, grid = forest_tuner, 
                   resamples = resamps, metrics = metric_set(rmse, rsq)) %>% 
    select_best(metric = 'rmse', maximize = F)
  
  forest_wf = forest_wf %>% finalize_workflow(best)
  
  return(forest_wf)
}

post_lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames %>% unlist(), with = F] %>% as.matrix(), 
                       nfolds = k_folds
  )
  pl = rownames(coef(lassofit, s = 'lambda.min'))[coef(lassofit, 
                                                       s = 'lambda.min')[,1]!= 0]
  reg = lm(as.formula(paste0('x1 ~ ', paste(tail(pl,-1), collapse = " + "))), data = dt)
  return(reg)
}

# Function: One iteration ----------------------------------------------------------------
one_iter = function(iter, n, original = FALSE) {
  # Set parameters
  β0 = 1
  β1 = 1
  β2 = 1
  num_zs = 7
  corr = .2
  # Generate the instrument
if (original == TRUE){
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
  for (j in (1:num_zs)) {
    for (h in (1:num_zs)) {
      cvmat[c(j),c(h)] = corr^abs(j-h)
    }
  }
  
  i_dt_corr = data.table(
    mvrnorm(n = 2*n, mu = rep(0, num_zs), Sigma = cvmat),
    rnorm(2*n))
  
  znames = sapply(1:num_zs, FUN = function(x) paste0('z',as.character(x)))
  
  #replace names with z-names.
  i_dt_corr = setnames(i_dt_corr, c(znames, "e_common"))
  
  #replicate Ed's code using the new utility function (which is really just Ed's code)
  i_dt = generate_system(i_dt, n)
  
  #generate system for correlated system, then split for cv process.
  i_dt_corr = generate_system(i_dt_corr, 2*n)
  i_dt_corr_train <- i_dt_corr %>% sample_frac(.5)
  i_dt_corr_test = setdiff(i_dt_corr, i_dt_corr_train)
  i_dt_corr_test[, true_x := x1 + x2]
  form = as.formula(paste0('x1 ~ ', paste(znames, collapse = '+')))
} else{
  i_dt_corr_test = belloni(n, 100)
  i_dt_corr_test = data.table(i_dt_corr_test)
  znames = lapply(1:100, FUN = function(x) paste0('z_',as.character(x)))
  zform = paste(znames, collapse = '+')
  form = as.formula(paste0('x1 ~ ', zform))
}
  
  # OLS first stage
  i_dt_corr_test[, x1_hat_ols := predict(lm(
    formula = form,
    data = i_dt_corr_test
  ))]
  
  # Machine learn the first stage
  i_dt_corr_test[, x1_hat_rf := predict(party::cforest(
  formula = form,
    data = i_dt_corr_test
  ))]
  
  #Use PCA to create meta-instrument then run linear first-stage regression
  i_dt.pca = prcomp(i_dt_corr_test[,names(i_dt_corr_test) %in% znames, with = FALSE], center = TRUE,scale. = TRUE)
  
  i_dt_corr_test[, x1_hat_pca := predict(lm(
    x1 ~ i_dt.pca$x[,1] + i_dt.pca$x[,2] + i_dt.pca$x[,3],
    data = i_dt_corr_test
  ))]
  
  #Machine learn with cv on the first stage
  
  #Build model with test data
  rfobj = forest_fcn(dta = i_dt_corr_test, znames = znames, num_zs = 100)
  rffit = parsnip::fit(rfobj, data = i_dt_corr_test)
  i_dt_corr_test[, x1_hat_rf_cv := rffit %>% predict(new_data = i_dt_corr_test),
                 ]
  
  #Lasso model to select variables for first stag
  enetmod = enet(i_dt_corr_test, k = 5, form = form)
  i_dt_corr_test[, x1_hat_lasso_cv := predict(enetmod,
                                              newdata = i_dt_corr_test)]

  #Boosted Trees model for first stage
  adamodel = adamod(i_dt_corr_test, k = 5, form = form)
  i_dt_corr_test[, x1_hat_boost_cv := predict(adamodel,
                                              newdata = i_dt_corr_test)]
  
  plasso_mod = post_lasso(dt = i_dt_corr_test, 5, znames = znames)
  i_dt_corr_test[, x1_hat_plasso_cv := predict(plasso_mod,
                                              newdata = i_dt_corr_test)]
  
  #LDA, as per Ed's idea. We'd need a binary predictee I think. We could come
  #up with an LDA that returns nonclass results.
  #ldamodel = ldamod(i_dt_corr_test, k = 5)
  #i_dt_corr_test[, x1_hat_lda_cv := predict(ldamodel,
  #                                            newdata = i_dt_corr_test)]
  i_dt_corr_test = na.omit(i_dt_corr_test)
  #Run second regression
  bind_rows(
    # OLS: DGP
    lm(y ~ true_x + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("true_x", term)) %>%
      mutate(model = "First Stage: Oracle Model"),
    # OLS: OVB
    lm(y ~ x1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First Stage: Naive OLS"),
    # Second stage: OLS
    lm(y ~ x1_hat_ols + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: 2SLS"),
    # Second stage: ML
    lm(y ~ x1_hat_rf + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: Random Forest"),
    # Second stage: ML with CV
    lm(y ~ x1_hat_rf_cv + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: Random Forest, CV"),
    lm(y ~ x1_hat_pca + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "First stage: PCA"),
    lm(y ~ x1_hat_lasso_cv + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: LASSO selection"),
    lm(y ~ x1_hat_plasso_cv + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: LASSO selection"),
    lm(y ~ x1_hat_boost_cv + en1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "First stage: Boosted Trees")
    #important to add a CV'd PCA?
  ) %>% mutate(iter = iter) %>% data.table()
}

#big_data = fread('/Users/connor/Desktop/BelloniData_n_1e+06_typeof_exponential_mu_180_Nne_final.csv')
#one_iter(n= 1000, iter = 1, big_data = big_data, original = FALSE)

# Function: Run iterations ---------------------------------------------------------------
run_sim <-function(n, n_sims, seed = 12345,future_cores = 12) {
  future::plan("multicore", workers = future_cores)
  # Set the seed
  set.seed(seed)
  # Run a parallelized simulation
  future_lapply(X = 1:n_sims, FUN = one_iter, n = n, future.seed = seed) %>% rbindlist()
}
#stopCluster

# Run simulation -------------------------------------------------------------------------
#one_iter(n = 1000, iter =1)
sim_dt <- run_sim(n = 1e3, n_sims = 500)

# Plot simulation ------------------------------------------------------------------------
sim_dt[, median(estimate, na.rm = TRUE), by = model]
plotter2<-ggplot(
  #data = sim_dt %>% filter(grepl("Second", model)),
  data = sim_dt,
  aes(x = estimate, fill = model)
) +
  geom_density(color = NA, alpha = 0.5) +
  geom_vline(xintercept = 1) +
  theme_minimal() +
  scale_fill_viridis_d("Model", begin = 0.1, end = 0.85, option = "A")

#write out data and plot
fwrite(sim_dt, file = '/Users/connor/Dropbox/MLIV/Plots:Resources/SimulationDT_Belloni.csv')
ggsave('/Users/connor/Dropbox/MLIV/Plots:Resources/PCA.CVML.Full.Belloni.png',plotter2, device = png(),scale = 1.8, dpi = 320)

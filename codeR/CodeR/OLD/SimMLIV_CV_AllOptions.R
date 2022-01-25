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
  mboost, promises, elasticnet, coefplot, glmnet
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

bigfolds = function(mtry,ntrees, data, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7)) {
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
gridsearch = function(dt, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), max_trees = 300, mtry = c(1:7)) {
  #build grid of plausible values for maxdepth and minprob
  paramgrid <- list(mtry = mtry,
                    ntrees = seq(from = 50, to = max_trees, by = 80)) %>%
    cross_df()
  #attach fold count to grid
  paramgrid$fnum = fnum
  #map values of grid to bigfolds function above, return as new variable mse
  paramgrid = paramgrid %>% mutate(mse = future_pmap(paramgrid, bigfolds, data = dt) %>% unlist())
  #unlist mse
  #paramgrid$mse %<>% unlist()
  #sort by mse
  paramgrid %<>% arrange(mse) %>% head(1)
  return(paramgrid[,1:2])
}

#Caret implementations of lasso/enet and adaboost
adamod = function(dt, k_folds) {
  adaGrid <-  expand.grid(mstop = c(30:90), maxdepth = c(1:10))
  adafit %<-% train(as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), 
                  method = 'blackboost',
                  trControl = trainControl(sample, method = 'cv', number = k_folds),
                  data = dt,
                  tuneGrid = adaGrid
  )
  
  return(adafit)
}

lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), 
                       nfolds = k_folds
                       )
  return(lassofit)
}

post_lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), 
                       nfolds = k_folds
  )
  pl = rownames(coef(lassofit, s = 'lambda.min'))[coef(lassofit, 
                                                  s = 'lambda.min')[,1]!= 0]
  reg = lm(as.formula(paste0('x1 ~ ', paste(tail(pl,-1), collapse = " + "))), data = dt)
  return(reg)
}

ldamod = function(dt, k_folds) {
  ldafit %<-% train(as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), 
                  method = 'lda',
                  trControl = trainControl(sample,method = 'cv', number = k_folds),
                  data = dt
  )
  
  return(ldafit)
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
  
  # OLS first stage
  i_dt_corr_test[, x1_hat_ols := predict(lm(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test
  ))]
  
  i_dt_corr_test[, err_ols := x1 - x1_hat_ols]
  
  # Machine learn the first stage
  i_dt_corr_test[, x1_hat_rf := predict(party::cforest(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test
  ))]
  
  i_dt_corr_test[, err_rf := x1 - x1_hat_rf]
  
  #Use PCA to create meta-instrument then run linear first-stage regression
  i_dt.pca = prcomp(i_dt_corr_test[,names(i_dt_corr) %in% znames, with = FALSE], center = TRUE,scale. = TRUE)
  
  i_dt_corr_test[, x1_hat_pca := predict(lm(
    x1 ~ i_dt.pca$x[,1] + i_dt.pca$x[,2] + i_dt.pca$x[,3],
    data = i_dt_corr_test
  ))]
  i_dt_corr_test[, err_pca := x1 - x1_hat_pca]
  #Machine learn with cv on the first stage
  #find optimal oos control parameters
  cvparams = gridsearch(dt = i_dt_corr_test, fnum = 5)
  #build control parameter object with CV'd optimal parameters
  ctrl = cforest_control(mtry= cvparams[1,1] %>% unlist(), ntree = cvparams[1,2] %>% unlist())
  #Build model with test data
  i_dt_corr_test[, x1_hat_rf_cv := predict(party::cforest(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test,
    control = ctrl
  ), newdata = i_dt_corr_test)]
  
  i_dt_corr_test[, err_rf_cv := x1-x1_hat_rf_cv]
  #Lasso model to select variables for first stage (post lasso)
  enetmod = post_lasso(i_dt_corr_test, k = 5, znames)
  i_dt_corr_test[, x1_hat_plasso_cv := predict(enetmod,
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
  #LDA, as per Ed's idea. We'd need a binary predictee I think. We could come
  #up with an LDA that returns nonclass results.
  #ldamodel = ldamod(i_dt_corr_test, k = 5)
  #i_dt_corr_test[, x1_hat_lda_cv := predict(ldamodel,
  #                                            newdata = i_dt_corr_test)]
  
  #Run second regression
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
    
    # Second stage: Random Forest
    lm(y ~ x1_hat_rf, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Random Forest"),
    # Errors: Random Forest, no CV
    lm(y ~ x1 + err_rf, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: Random Forest"),
    
    # Second stage: ML with CV
    lm(y ~ x1_hat_rf_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Random Forest, CV"),
    # Errors: ML with CV
    lm(y ~ x1 + err_rf_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: Random Forest, CV"),
    
    # Second Stage: PCA, 3 components
    lm(y ~ x1_hat_pca, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: PCA"),
    #Errors: PCA
    lm(y ~ x1 + err_pca, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: PCA"),
    
    #Second Stage: Lasso Selection (post-lasso)
    lm(y ~ x1_hat_plasso_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: post-LASSO selection"),
    #Errors: post-lasso
    lm(y ~ x1 + err_plasso, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: post-LASSO selection"),
    
    #Second Stage: Lasso
    lm(y ~ x1_hat_lasso_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: LASSO selection"),
    #Errors: lasso
    lm(y ~ x1 + err_lasso, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('err', term)) %>%
      mutate(model = "err: Second stage: LASSO selection"),
    
    #Second Stage: boosted trees
    lm(y ~ x1_hat_boost_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl('x1', term)) %>%
      mutate(model = "Second stage: Boosted Trees"),
    #Errors: Boosted trees
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

one_iter(n = 1000, iter = 1)

#write out data and plot
#fwrite(sim_dt, file = '/Users/connor/Dropbox/MLIV/Plots:Resources/SimulationDT.csv')
#ggsave('/Users/connor/Dropbox/MLIV/Plots:Resources/PCA.CVML.Full.png',plotter2, device = png(),scale = 1.8, dpi = 320)

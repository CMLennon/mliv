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
  parallel, magrittr, MASS, partykit, rsample, furrr, future.apply
)

future::plan('multiprocess')
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

holdout_results <- function(splits, ...) {
  # Fit the model to the majority
  treemod <- ctree(..., data = analysis(splits))
  # Save the extra fold
  holdout <- assessment(splits)
  # bind results
  preds = predict(treemod, newdata = holdout)
  res = cbind(holdout, preds = preds)
  # save x1 vals, run preds, calculate sq error
  res$sqerr <- (res$preds - holdout$x1)^2
  # Return the assessment data set with the additional columns
  return(mean(res$sqerr))
}


#function to run above for all folds, inserting control parameters to ctree as they are searched in grid search

#temporarily until I find a workaround - found. I wasn't using the syntax correctly.

bigfolds = function(maxdepth,minprob, data, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7)) {
                          #Generate folds:
  
  folds = data %>% vfold_cv(v=fnum, repeats = 1)
  
  #build a ctree_control object with updated parameters
  control = ctree_control(maxdepth = maxdepth, minprob = minprob)

  #map the list of splits to the holdout reults function
  folds$mse <- folds$splits %>%
                      future_map(~ holdout_results(.x, mod_form, control = control))
    
    #return the set of mses for each fold
  #folds$mse <- future_map_dbl(folds$results, function(x) mean(x$sqerr))
  #return average mse across all folds
  mse = mean(folds$mse)%>%unlist()
  return(mse)
}
#gridsearch function to check across numerous parameter values for minprob (minimum proportion of data in final leaf) and maximum tree depth.
gridsearch = function(dt, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), max_proportion = 8, maxdepth = round(2*num_zs)) {
    #build grid of plausible values for maxdepth and minprob
  paramgrid <- list(maxdepth = c(2:(maxdepth)),
  minprob = c(1:max_proportion)/100) %>%
  cross_df()
  #attach fold count to grid
  paramgrid$fnum = fnum
  #map values of grid to bigfolds function above, return as new variable mse
  paramgrid = paramgrid %>% mutate(mse = future_pmap(paramgrid, bigfolds, data = dt))
  #unlist mse
  paramgrid$mse %<>% unlist()
  #sort by mse
  paramgrid %<>% arrange(mse) %>% head(1)
  return(paramgrid[,1:2] %>% head(1))
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
  
  # Machine learn the first stage
  i_dt_corr_test[, x1_hat_ml := predict(ctree(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7, control = ctree_control(),
    data = i_dt_corr_test
  ))]
  
  #Use PCA to create meta-instrument then run linear first-stage regression
  i_dt.pca = prcomp(i_dt_corr_test[,names(i_dt_corr) %in% znames, with = FALSE], center = TRUE,scale. = TRUE)
  
  i_dt_corr_test[, x1_hat_pca := predict(lm(
    x1 ~ i_dt.pca$x[,1] + i_dt.pca$x[,2] + i_dt.pca$x[,3],
    data = i_dt_corr_test
  ))]
  
  #Machine learn with cv on the first stage

  #find optimal oos control parameters
  cvparams = gridsearch(dt = i_dt_corr_test, fnum = 5)
  #print(cvparams)
  #build control parameter object with CV'd optimal parameters
  ctrl = ctree_control(maxdepth= cvparams[,1] %>% unlist(), minprob = cvparams[,2] %>% unlist())
  
  #build model with training data, predicted on test data (still 50% split at this point.)
  i_dt_corr_test[, x1_hat_ml_cv := predict(ctree(
    x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7,
    data = i_dt_corr_test,
    control = ctrl
  ), newdata = i_dt_corr_test)]
  #itercounter<<-itercounter + 1
  #print(itercounter)
  bind_rows(
    # OLS: DGP
    lm(y ~ x1 + x2, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "OLS: DGP"),
    # OLS: OVB
    lm(y ~ x1, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "OLS: OVB"),
    # Second stage: OLS
    lm(y ~ x1_hat_ols, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: OLS"),
    # Second stage: ML
    lm(y ~ x1_hat_ml, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: ML"),
    # Second stage: ML with CV
    lm(y ~ x1_hat_ml_cv, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: Ctree, CV"),
    lm(y ~ x1_hat_pca, data = i_dt_corr_test) %>% tidy(quick = T) %>%
      filter(grepl("x1", term)) %>%
      mutate(model = "Second stage: PCA")
    
    #important to add a CV'd PCA?
  ) %>% mutate(iter = iter) %>% data.table()
}

# Function: Run iterations ---------------------------------------------------------------
run_sim = function(n, n_sims, seed = 12345) {
  # Set the seed
  set.seed(seed)
  # Run a parallelized simulation
  future_lapply(X = 1:n_sims, FUN = one_iter, n = n) %>% rbindlist()
}

# Run simulation -------------------------------------------------------------------------
sim_dt = run_sim(n = 1e3, n_sims = 500)


 # Plot simulation ------------------------------------------------------------------------
sim_dt[, sd(estimate, na.rm = TRUE), by = model]

plotter2<-ggplot(
  #data = sim_dt %>% filter(grepl("Second", model)),
  data = sim_dt,
  aes(x = estimate, fill = model)
) +
  geom_density(color = NA, alpha = 0.65) +
  geom_vline(xintercept = 1) +
  theme_minimal() +
  scale_fill_viridis_d("Model", begin = 0.1, end = 0.85, option = "B")
plotter

sim_dt[, mean(estimate, na.rm = T), by = model]

#one_iter(1,1000)

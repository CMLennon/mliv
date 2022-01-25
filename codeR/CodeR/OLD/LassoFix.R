
library(pacman)
p_load(
  tidyverse, broom, data.table,
  lfe, fixest, estimatr, party,
  parallel, magrittr, MASS, partykit, rsample, furrr, future.apply, adabag,
  mboost, promises, elasticnet, coefplot, glmnet, tidymodels
)

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

lasso = function(dt, k_folds, znames) {
  lassofit = cv.glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), 
                       nfolds = k_folds,
                       alpha = 1,
                       lambda = seq(from = 0, to = 1000, length.out = 10000)
  )
  return(glmnet(y = dt$x1, x = dt[, znames, with = FALSE] %>% as.matrix(), lambda = lassofit$lambda.min, alpha = 1))
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





one_iter = function(iter, n) {
  
  num_zs = 7
  corr = .6
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

enetmod = post_lasso(i_dt_corr_test, k = 5, znames)
i_dt_corr_test[, x1_hat_plasso_cv := predict(enetmod,
                                             newdata = i_dt_corr_test)]
i_dt_corr_test[, err_plasso := x1 - x1_hat_plasso_cv]

lasso = lasso(i_dt_corr_test, k = 5, znames)
i_dt_corr_test[, x1_hat_lasso_cv := predict(lasso,
                                            newx = i_dt_corr_test[,znames, with = FALSE] %>% 
                                              as.matrix())]
i_dt_corr_test[, err_lasso := x1 - x1_hat_lasso_cv]

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
  #2S-LS Y error
  lm(y ~ x1 + err_ols, data = i_dt_corr_test) %>% tidy(quick = T) %>%
    filter(grepl('err', term)) %>%
    mutate(model = "err: Second stage: OLS "),
  
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

  
  #important to add a CV'd PCA?
) %>% mutate(iter = iter) %>% data.table()
}

run_sim <-function(n, n_sims, seed = 12345,future_cores = 12) {
  future::plan("multicore", workers = future_cores)
  # Set the seed
  set.seed(seed)
  # Run a parallelized simulation
  future_lapply(X = 1:n_sims, FUN = one_iter, n = n) %>% rbindlist()
}
#stopCluster

# Run simulation -------------------------------------------------------------------------
sim_dt <- run_sim(n = 1e3, n_sims = 1000)

write_csv(sim_dt, path = '/Users/connor/Dropbox/MLIV/Plots:Resources/LassoCheck.csv')
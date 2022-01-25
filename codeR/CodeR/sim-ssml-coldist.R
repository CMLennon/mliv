library(pacman)

p_load(AER, tidyverse, tidymodels, data.table, future.apply, future, magrittr)

data("CollegeDistance")
fwrite(CollegeDistance, '/Users/connor/Desktop/GithubProjects/IV-test/cdist.csv')

cdist = fread('/Users/connor/Desktop/GithubProjects/IV-test/cdist.csv')
cdist[, i_black := 1*(ethnicity == 'afam')]

cdist %<>% filter(region != 'west')




one_iter= function(iter, data){
  
  
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
  mse <- new_numeric_metric(mse_vec, direction = "minimize")
  
  mse <- function(data, ...) {
    UseMethod("mse")
  }
  
  mse <- new_numeric_metric(mse, direction = "minimize")
  
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
  
  set = metric_set(mse)
  data_used = sample_frac(data, size = .90, replace = T)
  init_split = vfold_cv(data_used, 5)

ss_xgboost = function(split, fol_num = 1,
                            metric = set){
    
  ssformula = 'education ~ tuition'
    train_dat = analysis(split)
    test_dat = assessment(split)
  
  bt_recipe = recipe(as.formula(ssformula), data = train_dat)
  
  bt_model = boost_tree(tree_depth = tune(),
                        learn_rate = tune(),
                        sample_size = tune(),
                        trees = 1000) %>% set_engine("xgboost")
  
  wf = workflow() %>% add_recipe(bt_recipe) %>% add_model(bt_model)
  
  resamp = vfold_cv(train_dat, 3)
  
  cvres = tune_grid(wf, resamples = resamp, metrics = metric, grid = 50) %>% select_best()

  reg_res = wf  %>% finalize_workflow(cvres) %>% fit(train_dat)
  test_dat$rf_instrument = predict(reg_res, new_data = test_dat)$.pred
  test_dat$ss_iv_instrument = predict(lm(data = train_dat, as.formula(ssformula)), newdata = test_dat)
  return(test_dat)
  
}

dt_tuit = lapply(init_split$splits, FUN = ss_xgboost, metric = set) %>% rbindlist()

dt[, first_stage_rf := (lm(education ~ rf_instrument, dt)$fitted.values)]
dt[, first_stage_ols := (lm(education ~ ss_iv_instrument, dt)$fitted.values)]

iter_dt = data.table()

iter_dt = lm(i_black ~ first_stage_rf, data = dt) %>% tidy(quick = T) %>%
  filter(grepl("rf", term)) %>%
  mutate(model = "First stage: Random Forest SSML")   %>%  mutate(true_var = var(dt$education)) %>%
  mutate(var = var(dt$first_stage_rf)) %>% bind_rows(.,iter_dt) %>% 
  data.table()

iter_dt = lm(i_black ~ first_stage_ols, data = dt) %>% tidy(quick = T) %>%
  filter(grepl("ols", term)) %>%
  mutate(model = "First stage: OLS")   %>%  mutate(true_var = var(dt$education)) %>% mutate(true_cov = cov(dt$education, dt$wage - dt$education)) %>%
  mutate(var = var(dt$first_stage_ols)) %>% mutate(cov_u = cov(dt$first_stage_rf, dt$wage - dt$education), cov_e = cov(dt$first_stage_rf, dt$education - dt$first_stage_ols))%>% bind_rows(.,iter_dt) %>% 
  data.table()

iter_dt %>% mutate(iter = iter)
return(iter_dt)
}

one_iter(data = cdist)
run_sim = function(iters = 300, func = one_iter, data = cdist){
  
  
  future::plan('multicore', workers = availableCores())
  future_seed = 42
  sim_res = future.apply::future_lapply(1:iters ,FUN = func, data = data, future.seed = future_seed) %>% rbindlist()
  fwrite(sim_res, 'random_forest_ssml_res.csv')
}

run_sim_nopar = function(iters = 300, func = one_iter, data = cdist){
  sim_res = lapply(1:iters ,FUN = func, data = data) %>% rbindlist()
  #fwrite(sim_res, 'random_forest_ssml_res_samp.csv')
  return(sim_res)
  
}

results_dist = run_sim_nopar()

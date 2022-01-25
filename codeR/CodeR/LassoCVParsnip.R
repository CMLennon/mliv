#Utility function to run a reproducible lasso function

library(pacman)

p_load(tidymodels)

lasso_fit_cv = function(seed, datatable, folds, oos = TRUE){
  #take seed, data and formula and estimate/crossvalidate a lasso model
  
  #seed: take seedno.
  #datatable: pass in data
  # formula: No longer needed - formula built based on data fed in, assumes all instruments start with z_
  #Folds: rsample Rsplits object for crossvalidation. This is passesd separately
  #to avoid any irregularities due to being trained on different datasets
  
  #returns x-hat on datatable
  
  zs = grep("z", names(datatable), value = "T")
  num_zs = length(zs)
  formula = paste0("x1 ~ ", paste0(zs, collapse = '+')) %>% as.formula()
  
  lswf = workflow()
  lsrecipe = recipe(formula, data = datatable)
  search_grid = expand.grid(penalty = seq(.01, 1, length.out = 100))
  lasso_model = linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet")
  
  lswf = lswf %>% add_recipe(lsrecipe) %>% add_model(lasso_model)
  
  best = tune_grid(lswf, grid = search_grid, resamples = folds, metrics = metric_set(rmse)) %>%
    select_best(metric='rmse', maximize = F)
  
  lswf = lswf %>% finalize_workflow(best)
  
  if(oos){
    produce_oos_preds = function(foldnum,folds, lswf){
      train_dt = folds$splits[[foldnum]] %>% analysis()
      test_dt = folds$splits[[foldnum]] %>% assessment()
      lasso_fit = parsnip::fit(lswf, data = train_dt)
      test_dt$out = predict(lasso_fit, new_data = test_dt)
      return(test_dt)
    }
    
    out = lapply(1:5, produce_oos_preds, folds = folds, lswf = lswf) %>% rbindlist() %>% arrange(id)
    return(out$out)
  }
  else{
  x1_xhat_lasso_cv = parsnip::fit(lswf, data = datatable)
  out = predict(x1_xhat_lasso_cv, new_data = datatable)
  return(out)
  }
}



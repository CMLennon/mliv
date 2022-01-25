pca_fit_cv = function(seed, datatable, folds, cv = F){
  #takes in a seed (fix randomness in order)
  #datatable is data to train/predict on
  #formula is created based on datatable, reads any variable with a z_* prefix as an instrument
  #folds is an rsample rsplit object of some kind
  
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)
  
  form = paste0("x1 ~ ", paste0(zs, collapse = '+')) %>% as.formula()
  
  pcawf = workflow()
  
  pcarecipe = recipe(form, data = datatable) %>% step_pca(all_predictors(), num_comp = tune())
  
  if(cv == T){
    search_grid = expand.grid(num_comp = c(1,2,3, seq(4,round((2*num_zs)/3), by = round(num_zs/5))))
  }
  else{
    search_grid = expand.grid(num_comp = 3)
  }
  
  pca_mod = linear_reg() %>% set_engine('lm')
  
  pcawf = pcawf %>% add_recipe(pcarecipe) %>% add_model(pca_mod)
  
  best = tune_grid(pcawf, resamples = folds, grid = search_grid,
                   metrics = metric_set(rmse)) %>% select_best(metric = 'rmse', maximize = F)
  
  pcawf = pcawf %>% finalize_workflow(best)
  
  modfit = parsnip::fit(pcawf, data = datatable)
  
  out = predict(modfit, new_data = datatable)
  
  return(out)
}

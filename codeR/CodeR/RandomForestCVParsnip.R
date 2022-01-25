#Random Forest Parsnip utility function

forest_fcn = function(seed, datatable, folds, cv = T, oos = T){
  
  #Function that builds and finalizes a tidymodels randomforest, using
  #ranger as its engine. It takes data, a formula to estimate, and returns
  #a finalized workflow object
  
  ### Inputs:
  
  # datatable: a data-frame like object. Must be compatible with tidymodels
  
  # num_zs: the number of instruments used
  
  #znames: the vector of instrument variable names, as character values
  
  #lossfcn: a string containing 'mse' or 'mae' to determine how the model
  #will select hyperparameters. Default is mse.
  #datatable = i_dt_corr_test
  #set recipipe = recipe(as.formula(paste0('x1 ~ ', paste0(znames, collapse = '+'))), 
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)
  
  form = paste0("x1 ~ ", paste0(zs, collapse = '+')) %>% as.formula()
  
  forest_recipe = recipe(form, 
                         data = datatable)
  
  #set model, make the models tunable onnumber of trees, mtry, smallest leaf
  forest_model = rand_forest() %>% set_engine('ranger') %>% 
    set_args(trees = tune(), mtry = tune(), min_n = tune()) %>% set_mode('regression')
  
  forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)
  
  if(cv == T){
    #choose optimal hps by grid search
    forest_tuner = expand_grid(trees = seq(from = 50, to = 800, by = 100), 
                               mtry = seq(2,(num_zs-1), by = (num_zs/3) %>% floor()), 
                               min_n = seq(1,8, by = 2))
    best = tune_grid(forest_wf, grid = forest_tuner, 
                     resamples = folds, metrics = metric_set(rmse, rsq)) %>% 
      select_best(metric = 'rmse', maximize = F)
  }
  else{
    #replace with ranger default values
    best = expand_grid(trees = 500, 
                               mtry = sqrt(num_zs) %>% floor(), 
                               min_n = 5)
  }
  #build workflow object
  
  
  #find best CVd parameters
  
  
  forest_wf = forest_wf %>% finalize_workflow(best)
  
  if(oos){
    produce_oos_preds = function(foldnum,folds, forest_wf){
      train_dt = folds$splits[[foldnum]] %>% analysis()
      test_dt = folds$splits[[foldnum]] %>% assessment()
      forest_fit = parsnip::fit(forest_wf, data = train_dt)
      test_dt$out = predict(forest_fit, new_data = test_dt)
      return(test_dt)
    }
    
    out = lapply(1:5, produce_oos_preds, folds = folds, forest_wf = forest_wf) %>% rbindlist() %>% arrange(id)
    return(out$out)
  }
  else{
  forest_fit = parsnip::fit(forest_wf, data = datatable)
  out = predict(forest_fit, new_data = datatable)
  
  return(out)
  }
}

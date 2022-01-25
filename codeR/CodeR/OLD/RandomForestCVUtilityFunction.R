#Random Forest Parsnip utility function

forest_fcn = function(seed, datable, form, folds){
  
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
  forest_recipe = recipe(form, 
                         data = datatable)
  
  #set model, make the models tunable on mtry and smallest leaf
  forest_model = rand_forest() %>% set_engine('ranger') %>% 
    set_args(trees = tune(), mtry = tune(), min_n = tune()) %>% set_mode('regression')
  
  forest_tuner = expand_grid(trees = seq(from = 50, to = 300, by = 50), mtry = 2:(num_zs-1), min_n = 1:8)
  
  #build workflow object
  forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)
  
  #find best CVd parameters
  best = tune_grid(forest_wf, grid = forest_tuner, 
                   resamples = folds, metrics = metric_set(rmse, rsq)) %>% 
    select_best(metric = 'rmse', maximize = F)
  
  forest_wf = forest_wf %>% finalize_workflow(best)
  
  out = predict(forest_wf, new_data = datatable)
  
  return(forest_wf)
}
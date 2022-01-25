#Boosted trees parsnip specification

xgbst_fit = function(seed, datatable, folds, oos = T){
  #takes in a seed (fix randomness in order)
  #datatable is data to train/predict on
  #formula is created based on datatable, reads any variable with a z_* prefix as an instrument
  #folds is an rsample rsplit object of some kind
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)
  
  form = paste0("x1 ~ ", paste0(zs, collapse = '+'))
  
  bstwf = workflow()
  
  bstrecipe = recipe(as.formula(form), data = datatable)
  
  search_grid = expand.grid(trees = 400)
  
  bst_mod = boost_tree(mode = 'regression', trees = tune()) %>% set_engine('xgboost')
  
  bstwf = bstwf %>% add_recipe(bstrecipe) %>% add_model(bst_mod)
  
  best = tune_grid(bstwf, resamples = folds, grid = search_grid,
                   metrics = metric_set(rmse)) %>% select_best(metric = 'rmse', maximize = F)
  
  bstwf = bstwf %>% finalize_workflow(best)
  
  if(oos){
    produce_oos_preds = function(foldnum,folds, bstwf){
      train_dt = folds$splits[[foldnum]] %>% analysis()
      test_dt = folds$splits[[foldnum]] %>% assessment()
      bstfit = parsnip::fit(bstwf, data = train_dt)
      test_dt$out = predict(bstfit, new_data = test_dt)
      return(test_dt)
    }
    
    out = lapply(1:5, produce_oos_preds, folds = folds, bstwf = bstwf) %>% rbindlist() %>% arrange(id)
    return(out$out)
  }
  out = parsnip::fit(bstwf, data = datatable)
  return(out)
}

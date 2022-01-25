#Utility function to run a reproducible JIVE estimation

#Uses ivmodels for this implementation - I really don't like ivmodels, I would prefer a self-implementation,
#but time.

library(pacman)

p_load(tidymodels, glmnet, tidyverse, ivmodel, data.table, SteinIV)

JIVE_fit = function(seed, datatable, folds, iter){
  #take seed, data and formula and estimate/crossvalidate a post-lasso model
  
  #seed: take seedno.

  #datatable: pass in data

  #dynamically reads instruments from the dataset. Assumes instruments begin with 'z_'

  #Folds: rsample Rsplits object for crossvalidation. This is passesd separately
  #to avoid any irregularities due to being trained on different datasets
  
  #returns x-hat on datatable
  
  n = nrow(datatable)
  
  foldmapper = function(i, folds){
    fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    data.table(fold)
    fold[,fold.id := i]
    return(fold)
  }
  
  datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()

  datatable = datatable %>% arrange(id)
  
  Z_vars = grep("z", names(datatable), value = "T")
  
  num_zs = length(Z_vars)
  
  X_vars = c('x1')

  y_vars = c('y')
  
  mod = jive.est(y = datatable[, ..y_vars] %>% as.matrix(), X = datatable[, ..X_vars] %>% as.matrix(), Z = datatable[, ..Z_vars] %>% as.matrix(), SE = TRUE, n.bt = n/10)

  mod$test_statistic = mod$est[1]/mod$se[1]

  mod$p.value = pnorm(mod$test_statistic, mean = 0, sd = 1)
  
  return(mod)
}


#Utility function to run a reproducible LIML estimation

#Uses ivmodels for this implementation - I really don't like ivmodels, I would prefer a self-implementation,
#but time.

#Maybe someone else knows - why would LIML do so badly with our data? It often produces estimates WAY lower
#than it should. Like .86 when using 100 instruments.

library(pacman)

p_load(tidymodels, glmnet, tidyverse, ivmodel, data.table)

LIML_fit = function(seed, datatable, folds, iter, plugin = F){
  #take seed, data and formula and estimate/crossvalidate a post-lasso model
  
  #seed: take seedno.
  #datatable: pass in data
  #formula: now read dynamically from the dataset. Assumes instruments begin with 'z_'
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
  
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)
  
  form = paste0("y ~ x1 |", paste0(zs, collapse = '+'))
  
  form = as.formula(form)
  
  mod = ivmodelFormula(formula = form, data = datatable)
  
  Fuller_fit = Fuller(mod, beta0 = 0, alpha = .05)
  
  return(Fuller_fit)
}


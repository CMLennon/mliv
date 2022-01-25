#Utility function to run a reproducible SVA, as described in phillipines health-stock paper

library(pacman)

p_load(tidymodels, glmnet, tidyverse, ivmodel, data.table)

SVA_fit = function(seed, datatable, folds, iter, plugin = F){
  #take seed, data and formula and estimate/crossvalidate a post-lasso model
  
  #seed: take seedno.
  #datatable: pass in data
  #formula: now read dynamically from the dataset. Assumes instruments begin with 'z_'
  #Folds: rsample Rsplits object for crossvalidation. This is passesd separately
  #to avoid any irregularities due to being trained on different datasets
  
  #Because this method uses 'OLS' to fit on a modified y, this returns an estimate of beta rather
  #than X-hat.
  
  n = nrow(datatable)
  
  foldmapper = function(i, folds){
    fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    data.table(fold)
    fold[,fold.id := i]
    return(fold)
  }
  foldlist = lapply(c(1:5), foldmapper, folds)
  datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
  
  datatable = datatable %>% arrange(id)
  
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)
  
  form = paste0("y ~ x1 |", paste0(zs, collapse = '+'))
  
  form = as.formula(form)
  
  SVA = function(fold, foldlist, unbiased, zs){
    trainlist = setdiff(foldlist, fold)
    SVAOnce = function(train, fold, unbiased, zs){
      SVD = svd(fold[,..zs])
      d = SVD[[1]]
      U = SVD[[2]]
      Vt = SVD[[3]]
      S = t(U) %*% as.matrix(fold[,..zs]) %*% t(Vt)
      g = t(U) %*% fold$x1
      
      newdat = data.table(
        g = g,
        S
      )
      SVAObj = lm(data = newdat, g ~ .)
      fold$x_iv = predict(init_iv, newdata = fold)
      return(fold$x_iv)
    }
  
  mod = ivmodelFormula(formula = form, data = datatable)
  mod
  Fuller_fit = Fuller(mod, beta0 = 0, alpha = .05)
  LIML_fit$point.est
  return(predict(plsfit, new_data = datatable)$.pred)
}


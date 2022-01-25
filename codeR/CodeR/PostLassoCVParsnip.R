#Utility function to run a reproducible post-lasso function

library(pacman)

p_load(tidymodels, glmnet, tidyverse)

plasso_fit_cv = function(seed, datatable, folds, iter, plugin = T, oos = F){
  #take seed, data and formula and estimate/crossvalidate a post-lasso model
  
  #seed: take seedno.
  #datatable: pass in data
  #formula: now read dynamically from the dataset. Assumes instruments begin with 'z_'
  #Folds: rsample Rsplits object for crossvalidation. This is passesd separately
  #to avoid any irregularities due to being trained on different datasets
  
  #returns x-hat on datatable
  zs = grep("z", names(datatable), value = T)
  
  n = nrow(datatable)
  
  foldmapper = function(i, folds){
    fold = assessment(folds$splits[[i]]) %>% data.table()
    data.table(fold)
    fold[,fold.id := i]
    return(fold)
  }
  
  datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
  
  #weirdly, using the double pipe will break this command...
  datatable = datatable %>% arrange(id)
 
  
  
  num_zs = length(zs)
  
  form = paste0("x1 ~ ", paste0(zs, collapse = '+'))
  
  form = as.formula(form)
  if(plugin == F){
    #crossvalidate lambda
    lassofit = cv.glmnet(y = datatable$x1, x = datatable[, ..zs] %>% as.matrix(), 
                       foldid = datatable$fold.id)
  }
  else{
    #Calculate and use plug-in penalty from Belloni 2013
    #gamma constant from Belloni2012 - .1/log(p v n) - translates to maximum of n or num of instr.
    #replace num_zs with new variable if non-structured dataset breaks this
    gam = .1/log(max(num_zs, n)) #as per Belloni 2012.
    c = 1.1 #reccomendation from Belloni 2013
    k_e = 1 #number of endogenous variables. Set this differently if k_e changes.
    p = num_zs
    lmd = c*2*sqrt(n)*qnorm(1 - (gam)/(2*k_e*p))
    lassofit = glmnet::glmnet(y = datatable$x1, x = datatable[, ..zs] %>% as.matrix(),
                              lambda = lmd)
  }
  
  #find variables such that coefficients are equal to 0
  regvars = rownames(coef(lassofit, s = 'lambda.min'))[coef(lassofit, s = 'lambda.min')[,1]!= 0]
  
  post_lasso = linear_reg() %>% set_engine('lm')
  
  if(regvars != "(Intercept)"){
    #print('made it')
    formula2_string = paste0('x1 ~ ', paste0(grep('z_',regvars,value = 'T'), collapse = '+'))
    formula2 = as.formula(formula2_string)
  
  }
  
  else{
    
    formula2 = as.formula('x1 ~ z_1')
  }
  
  plsrecipe = recipe(formula2, data = datatable)
  
  plswf = workflow() %>% add_recipe(plsrecipe) %>% add_model(post_lasso)
  
  if(oos){
    produce_oos_preds = function(foldnum,folds, plswf){
      train_dt = folds$splits[[foldnum]] %>% analysis()
      test_dt = folds$splits[[foldnum]] %>% assessment()
      plsfit = parsnip::fit(plswf, data = train_dt)
      test_dt$out = predict(plsfit, new_data = test_dt)$.pred
      return(test_dt)
    }
    
    out = lapply(1:5, produce_oos_preds, folds = folds, plswf=plswf) %>% rbindlist() %>% arrange(id)
    return(out$out)
  }
  
  plsfit = parsnip::fit(plswf, data = datatable)
  
  return(predict(plsfit, new_data = datatable)$.pred)
}


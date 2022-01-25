#Split Sample IV Utility Function

ssiv_fcn = function(seed, datatable, folds, unbiased = T, cv = T){
  
  #Function estimates an unbiased split-sample IV procedure, using
  #lm as its engine. It takes a data table object, and returns a set of Xs.
  
  ### Inputs:
  
  # datatable: a data-frame like object. Must be compatible with tidymodels
  
  # seed: random seed to maintain constant randomness

  # folds: a rsample/rsets object
  
  # cv: Traditional SSIV uses two samples. However, this isn't a necessary requirement for the unbiased version. For traditional, set cv to false.
  
  # Unbiased: logical toggling the beta correction on or off. Functionally, takes beta_IV and multiplies it by Theta inverse, where Theta is X_{1,2}(X_{1,2} inv(X_{1,2}))X_1

	# Note: Classic SSIV, as described in Angrist 1995 would be parametrized by ssiv_fcn(seed, datatable, folds, unbiased = F, cv = F)

 foldmapper = function(i, folds){
    fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    data.table(fold)
    return(fold)
  }

  foldlist = lapply(c(1:5), foldmapper, folds)
  datatable = foldlist %>% rbindlist() %>% data.table()

  datatable %>% arrange(id)

  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)

  if(cv){
  ssiv = function(fold, foldlist, unbiased){
  	trainlist = setdiff(foldlist, fold)
  	ssivonce = function(train, fold, unbiased){
  		
  		init_iv = lm(data = train, as.formula(paste0('x1 ~ ', paste0(zs, collapse = '+'))))
  		fold$x_iv = predict(init_iv, newdata = fold)
  	return(fold$x_iv)
  }
  x_iv = lapply(trainlist, ssivonce, fold, unbiased) %>% bind_cols() %>% rowMeans()
  fold$x_iv = x_iv
  
  if(unbiased){
  	fold$x_iv = predict(lm(data = fold, x1 ~ x_iv))
}
  return(fold)
  	}
  	dta = lapply(foldlist, ssiv, foldlist, unbiased) %>% rbindlist() %>% arrange(id)
  		return(dta$x_iv)
  }

  else{
  	ssiv = function(fold, datatable, unbiased){
  	train = setdiff(datatable, fold)

  	init_iv = lm(data = train, as.formula(paste0('x1 ~ ', paste0(zs, collapse = '+'))))
  	fold$x_iv = predict(init_iv, newdata = fold)
  		
  	if(unbiased){
  		fold$x_iv = predict(lm(data = fold, x1 ~ x_iv))
  	}
  	return(fold)
  }
  dta = lapply(foldlist, ssiv, datatable, unbiased) %>% rbindlist() %>% arrange(id)
  dta$x_iv_2 = lm(data = dta, x1 ~ x_iv) %>% predict(new_data = dta)
  		return(dta$x_iv_2)
  }
}
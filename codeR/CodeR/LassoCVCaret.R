## Testing Caret 

## Function that returns predictions using caret

lasso_fit_cv_c = function(dt, k_folds, folds) {
	

  foldmapper = function(i, folds){
    fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    data.table(fold)
    fold[,fold.id := i]
    return(fold)
  }

  zs = grep("z", names(dt), value = "T")
  form = paste0("x1 ~ ", paste0(zs, collapse = '+'))
  k_folds = 5
  form = as.formula(form)
  dt = lapply(c(1:k_folds), foldmapper, folds) %>% rbindlist() %>% data.table()
  dt = dt %>% arrange(id)

  train.control <- trainControl( 
    method = "cv",
    number = 5,
    savePredictions = T 
    )

  enetfit <- train(form, 
                   method = 'glmnet',
                   trControl = train.control,
                   data = dt,
                   tuneGrid = expand.grid(alpha = 1, lambda = seq(from = 0, to = 1, by = .01))
  )
  x1_hat_lasso_cv = predict(enetfit, newdata = dt)
  return(x1_hat_lasso_cv)
}
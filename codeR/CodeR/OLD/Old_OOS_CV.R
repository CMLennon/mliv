
holdout_results = function(splits, ...) {
  # Fit the model to the majority
  treemod <- party::cforest(..., data = analysis(splits))
  # Save the extra fold
  holdout <- assessment(splits)
  # bind results
  preds = predict(treemod, newdata = holdout)
  res = cbind(holdout, preds = preds)
  # save x1 vals, run preds, calculate sq error
  res$sqerr <- (res$preds - holdout$x1)^2
  # Return the assessment data set with the additional columns
  return(mean(res$sqerr%>% unlist()))
}


#function to run above for all folds, inserting control parameters to ctree as they are searched in grid search

#temporarily until I find a workaround - found. I wasn't using the syntax correctly.

bigfolds = function(mtry,ntrees, data, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7)) {
  #Generate folds:
  
  folds <- data %>% vfold_cv(v=fnum, repeats = 1)
  
  #build a ctree_control object with updated parameters
  control <- cforest_control(mtry = mtry, ntree = ntrees)
  
  #map the list of splits to the holdout reults function
  folds$mse <- folds$splits %>% future_map(~ holdout_results(.x, mod_form, control = control))
  
  #return the set of mses for each fold.
  
  #No longer needed, to avoid cpu overhead moved this to holdout_results. As per future's documentation for efficiency.
  
  #folds$mse <- future_map_dbl(folds$results, function(x) mean(x$sqerr))
  
  #return average mse across all folds
  mse = mean(folds$mse%>%unlist())
  return(mse)
}

#gridsearch function to check across numerous parameter values for minprob (minimum proportion of data in final leaf) and maximum tree depth.
gridsearch = function(dt, fnum, mod_form = as.formula(x1 ~ z1 + z2 + z3 + z4 + z5 + z6 + z7), max_trees = 300, mtry = c(1:7)) {
  #build grid of plausible values for maxdepth and minprob
  paramgrid <- list(mtry = mtry,
                    ntrees = seq(from = 50, to = max_trees, by = 80)) %>%
    cross_df()
  #attach fold count to grid
  paramgrid$fnum = fnum
  #map values of grid to bigfolds function above, return as new variable mse
  paramgrid = paramgrid %>% mutate(mse = future_pmap(paramgrid, bigfolds, data = dt) %>% unlist())
  #unlist mse
  #paramgrid$mse %<>% unlist()
  #sort by mse
  paramgrid %<>% arrange(mse) %>% head(1)
  return(paramgrid[,1:2])
}

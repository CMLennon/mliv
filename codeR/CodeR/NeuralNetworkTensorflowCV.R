#Utility function to run a reproducible shallow or deep neural-network instrument

#Not CV'd yet. I think 

library(pacman)

p_load(tidymodels, glmnet, tidyverse, keras)

ff_net_cv = function(seed, datatable, folds, iter, wide = T, maxdepth = 18, maxwidth = 90){
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
  
  #weirdly, using the double pipe will break this command...
  datatable = datatable %>% arrange(id)
 
  zs = grep("z", names(datatable), value = "T")
  
  num_zs = length(zs)

  tfnnfs = function(foldnum, datatable, wide, zs, width, depth, drop, valid_split, epochs, mse = F){

    ## Training data and labels
    train_dt = datatable %>% filter(fold.id != foldnum) %>% data.table()
    train_dat = train_dt[,..zs]
    train_lab = train_dt[,x1]

    ## Testing data and labels
    test_dt = datatable %>% filter(fold.id == foldnum) %>% data.table()
    test_dat = test_dt[,..zs]
    test_lab = test_dt[,x1]

  if(wide){
    model_keras <- keras_model_sequential() %>%
    layer_dense(units = width,
      kernel_initializer = "uniform",
      activation = "relu",
      input_shape = ncol(train_dat)) %>%
    layer_dropout(rate = drop) %>%
    layer_dense(
      units = 1,
      activation = "linear") %>% compile(
      optimizer = 'adam',
      loss = 'MSE') %>% layer_batch_normalization()

      history = fit(
        model_keras,
        x = as.matrix(train_dat),
        y = train_lab,
        batch_size = 4,
        epochs = epochs,
        validation_split = valid_split)

      x1out = model_keras %>% predict(as.matrix(test_dat)) 
      if(mse){
      return(mean(x1out[,1]-test_lab))
     }
     else{
     test_dt$.preds = x1out[,1]
     return(test_dt)
   }
   
  } else{
    nlay = depth

      model_keras <- keras_model_sequential() %>% 
      layer_dense(units = 10,
      kernel_initializer = "uniform",
      activation = "relu",
      input_shape = ncol(train_dat)) %>% layer_batch_normalization()
      i = 0
      repeat{
    model_keras %>% layer_dense(units = 10,
      kernel_initializer = "uniform",
      activation = "linear") %>% layer_batch_normalization()
    i = i+1
    if(i > nlay){
      break
    }
  }
    model_keras %>% layer_dense(
      units = 1,
      activation = "linear") %>% compile(
      optimizer = 'sgd',
      loss = 'MSE')

      history = fit(
        model_keras,
        x = as.matrix(train_dat),
        y = train_lab,
        batch_size = 50,
        epochs = epochs,
        validation_split = valid_split,
        callbacks = list(
      callback_early_stopping(patience = 5),
      callback_reduce_lr_on_plateau(factor = 0.05)
      ))

     x1out = model_keras %>% predict(as.matrix(test_dat))

     if(mse){
      return(mean(x1out[,1]-test_lab))
     }
     else{
      test_dt$.preds = x1out[,1]
     return(test_dt)
   }
  }
}

nnetcv = function(datatable, gridcv, wide, valid_split, zs, fold_count){
  #gridcv is a grid of possible test values. column 1 corresponds to depth, 2 corresponds to width (width or depth will only be releavant one at a time) and drop is a dropout rate
  folds = 1:fold_count
  if(wide){
    widths = gridcv['width']
    dropouts = gridcv['dropout']
    depth = NaN
    mse = 10000000000
    for(width in widths){
      for(dropout in dropouts){
        mses = lapply(folds, tfnnfs, datatable, wide, zs, width, depth, drop, valid_split = .2, mse = T) %>% rbindlist()
        if(mse > mean(mses)){
          bestwidth = width
          bestdropout = dropout
        }
      }
    }
    datatable_done = lapply(folds, tfnnfs, datatable, wide, zs, bestwidth, drop = bestdropout, valid_split = .2, mse = F) %>% rbindlist() %>% arrange(id)
    return(datatable_done)
  } else{
    depths = gridcv['depth']
    dropouts = gridcv['dropout']
    width = NaN
    mse = 10000000000
    for(depth in depths){
      for(dropout in dropouts){
        mses = lapply(folds, tfnnfs, datatable, wide, zs, depth, width, drop, valid_split = .2, mse = T) %>% rbindlist()
        if(mse > mean(mses)){
          bestwidth = width
          bestdropout = dropout
        }
      }
    }
    datatable_done = lapply(folds, tfnnfs, datatable, wide, zs, bestwidth, drop = bestdropout, valid_split = .2, mse = F) %>% rbindlist() %>% arrange(id)
    return(datatable_done)

  }



}
  gridcv = expand.grid(depth = 3:maxdepth, width = 19:maxwidth, dropout = seq(from = 0, to = .3, by = .1))
  return(nnetcv(datatable, gridcv, wide, valud_split, zs, fold_count = 5))
}

}


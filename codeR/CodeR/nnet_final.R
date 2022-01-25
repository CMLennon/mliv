
library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor, tensorflow, magrittr)

#-----------------------------------------
#run final neural network once hyperparameters are learned. Must provide hyperparameter values.

#takes - datatable : all data
#		 width : width of hidden layers
#		 dropout: dropout rate for connnections
#		 depth: how many hidden layers
#		 epoch_num: how many epochs to train for
#		 batch_size: how large should minibatches be?
#		 activation: what function do we use for activation?
#-----------------------------------------	

final_nnet = function(datatable, folds, width, dropout, depth, epoch_num, batch_size, activationfcn, oos = FALSE, foldnum = 5, z_prob = FALSE){
  
  
  k_clear_session()
	n = nrow(datatable)
	


	foldmapper = function(i, folds){
    	fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    	data.table(fold)
    	fold[,fold.id := i]
    	return(fold)
 	 	}

 	datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
 	
 	if(z_prob == TRUE){
 	  #nonlinear excludability violation
 	  datatable = datatable %>% mutate(z_prob= ifelse(abs(sin(z_1)) > .995, round(sin(z_1)/10,2), 0))
 	  datatable %<>% mutate(x1 = x1 + z_prob, y = y + 2*z_prob)
 	}
 	
 	datatable %<>% arrange(id)
 	zs = grep("z", names(datatable), value = "T")

 	num_zs = length(zs)

 	if(oos == TRUE){
 		#when OOS is true, the model will be trained on the training portion of the 5 folds and use the final fold to produce predictions. All 5 folds stitched together will produce an out-of-sample prediction a-la split-sample IV.
 		out_dt = data.table()

 		#datatable$fold.id2 = sample(fold.id.vect, nrow(datatable))
 		for(fold in 1:foldnum){
 			train_dt = datatable %>% filter(fold.id != fold) %>% data.table() %>% arrange(id)
 			#input = input_fn(object = train_dt, features = zs, response = 'x1', batch_size = batch_size, epochs = epoch_num)
 			train_dat = train_dt[,..zs] %>% as.matrix()
    		train_lab = train_dt[,x1] %>% as.matrix()

    		#test and output data
    		test_dt = datatable %>% filter(fold.id == fold) %>% data.table() %>% arrange(id)
    		test_dat = test_dt[,..zs] %>% as.matrix()
    		test_lab = test_dt[,x1] %>% as.matrix()

			if(depth <= 1){
				model = keras_model_sequential() %>% layer_dense(units = 1, activation = 'linear', input_shape = c(num_zs))
			} else{
				model = keras_model_sequential() %>% layer_dense(units = width, input_shape = c(num_zs)) %>%
				layer_activation_leaky_relu()
				if(depth > 2){
					for(laycount in 1:(depth - 2)){
						model %>% layer_dense(units = width) %>%
						layer_activation_leaky_relu() %>%
						layer_dropout(rate = dropout)
							}
				
					}
				model %>% layer_dense(units = 1, activation = 'linear')}



		compile(model, optimizer = 'adam', loss = 'mse')


		history = model %>% fit(x= train_dat, y = train_lab, verbose = 0)

		#input = input_fn(test_dt, features = zs, response = 'x1')
		nnet_x1_hat = predict(model, x = test_dat)

		#nnet_x1_hat = predict(model, test_dat, test_lab, verbose = 0)

		test_dt$nnet_x1_hat = nnet_x1_hat

		out_dt = rbind(out_dt, test_dt)
 		}

 		out_dt = out_dt %>% arrange(id)
 		return(out_dt$nnet_x1_hat)
 	} else{
 		#when OOS is false, return in-sample predictions for model trained on all 5 folds with found optimal hyperparameters.
 		datatable = datatable %>% arrange(id)
 		train_dat = datatable[,..zs] %>% as.matrix()
 		train_lab = datatable[,x1] %>% as.matrix()
			if(depth <= 1){
				model = keras_model_sequential() %>% layer_dense(units = 1, activation = 'linear', input_shape = c(num_zs))
			} else{
				model = keras_model_sequential() %>% layer_dense(units = width, input_shape = c(num_zs), activation = 'relu') %>% 
				  layer_dropout(dropout = dropout)
				if(depth > 2){
					for(laycount in 1:(depth - 2)){
						model %>% layer_dense(units = width, activation = 'relu') %>%
						layer_dropout(rate = dropout)
							}
				
					}
				model %>% layer_dense(units = 1, activation = 'linear')}
 		
 		adamopt = keras::optimizer_adam(lr = .01, amsgrad = TRUE)
		compile(model, optimizer = adamopt, loss = 'mse')
		
		#early stopping to prevent model from burning useless CPU time. Plus, restore model to model of best fit
		erl_stop = keras::callback_early_stopping(monitor = 'val_loss', min_delta = .01, patience = 12, restore_best_weights = TRUE)
		reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5,
		                              patience=3, min_lr=0.0005)

		history = model %>% fit(x = train_dat, y =train_lab, epochs = epoch_num, batch_size = batch_size, verbose = 0, callbacks = list(erl_stop, reduce_lr))

		datatable$nnet_x1_hat = predict(model, x = train_dat)

		return(datatable$nnet_x1_hat)

 	}

}
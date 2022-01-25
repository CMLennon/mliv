
library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor)

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

final_nnet = function(datatable, z_prob = FALSE, folds, width, dropout, depth, epoch_num, batch_size, activationfcn, oos = FALSE, foldnum = 5){
	
	n = nrow(datatable)
	
	if(z_prob == TRUE){
	  datatable = datatable %>% mutate(z_prob = ifelse(abs(sin(datatable$z_3 + datatable$z_50 + datatable$z_100)) > .995, round(sin(datatable$z_1 + datatable$z_50 + datatable$z_100)/10,2), 0))
	  datatable %<>% mutate(x1 = x1 + z_prob, y = y + 2*z_prob)
	  }

	foldmapper = function(i, folds){
    	fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
    	data.table(fold)
    	fold[,fold.id := i]
    	return(fold)
 	 	}

 	datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
 	
 	
 	zs = grep("z", names(datatable), value = "T")

 	num_zs = length(zs)



 	if(oos == TRUE){
 		#when OOS is true, the model will be trained on the training portion of the 5 folds and use the final fold to produce predictions. All 5 folds stitched together will produce an out-of-sample prediction a-la split-sample IV.
 		out_dt = data.table()

 		#create new set of folds to build-model and then produce oos predictions
 		fold.id.vect = c(rep(1,nrow(datatable)/5 %>% floor()), 
 			rep(2,nrow(datatable)/5 %>% floor()),
 			rep(3,nrow(datatable)/5 %>% floor()),
 			rep(4,nrow(datatable)/5 %>% floor()),
 			rep(5, nrow(datatable)-(nrow(datatable)/5 %>% floor())*4))

 		datatable$fold.id2 = sample(fold.id.vect, nrow(datatable))
 		for(fold in 1:foldnum){
 			train_dt = datatable %>% filter(fold.id2 != foldnum) %>% data.table()
 			train_dat = train_dt[,..zs] %>% as.matrix()
    		train_lab = train_dt[,x1] %>% as.matrix()

    		#test and output data
    		test_dt = datatable %>% filter(fold.id2 == foldnum) %>% data.table()
    		test_dat = test_dt[,..zs] %>% as.matrix()
    		test_lab = test_dt[,x1] %>% as.matrix()

    		model = keras_model_sequential() %>% layer_dense(units = width1, activation = activationfcn, input_shape = c(num_zs))
 			#smallest model is 3 layers
			if(depth <= 1){
				model %>% layer_dense(units = width, activation = activationfcn) %>%
				layer_dropout(rate = dropout)
				} else{
					for(laycount in 1:depth){
					model %>% layer_dense(units = width, activation = activationfcn) %>%
					layer_dropout(rate = dropout)
					}
				}

		model %>% layer_dense(units = 1, activation = 'linear')

		compile(model, optimizer = 'adam', loss = 'mse')

		history = model %>% fit(x = train_dat, y =train_lab, epochs = epoch_num, batch_size = batch_size)

		nnet_x1_hat = evaluate(model, test_dat, test_lab, verbose = 0)

		test_dt$nnet_x1_hat = nnet_x1_hat

		out_dt = rbind(out_dt, test_dt)

 		}

 		out_dt = out_dt %>% arrange(id)
 		return(out_dt)
 	} else{
 		#when OOS is false, return in-sample predictions for model trained on all 5 folds with found optimal hyperparameters.
 		train_dat = datatable[,..zs] %>% as.matrix()
 		train_lab = datatable[,x1] %>% as.matrix()
 		model = keras_model_sequential() %>% layer_dense(units = width1, activation = activationfcn, input_shape = c(num_zs))
 		#smallest model is 3 layers
		if(depth <= 1){
			model %>% layer_dense(units = width, activation = activationfcn) %>%
			layer_dropout(rate = dropout)
			} else{
			for(laycount in 1:depth){
				model %>% layer_dense(units = width, activation = activationfcn) %>%
				layer_dropout(rate = dropout)
			}
			}

		model %>% layer_dense(units = 1, activation = 'linear')

		compile(model, optimizer = 'adam', loss = 'mse')

		history = model %>% fit(x = train_dat, y =train_lab, epochs = epoch_num, batch_size = batch_size)

		datatable$nnet_x1_hat = evaluate(model, train_dat, train_lab, verbose = 0)

		return(datatable)

 	}

}
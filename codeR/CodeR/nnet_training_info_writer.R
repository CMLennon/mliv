# build a file to figure out where the training happens,how the graphs diverge

## Key goals: find iteration-specifc hyper-parameters chosen

## Link to filepath-specific tuning run so we can see the training path.

library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor, magrittr, R.utils)

tsol = TRUE

training_info_nnet =function(iter, complex, shuffle, temp_solution = tsol, deep = TRUE, wide = TRUE){
	
	# First things first - we read in top-x object that matches the dataset and iteration given in arguments.

	# Then, we find out which one was chosen

	# Then, we read in loss-vector for said dataset and iteratively add those to that object

	if(temp_solution){

		if(complex){
			complex = 'T'
		} else{
			complex = 'F'
		}


		filename_topx = paste0('top.x_C', complex, '_s', shuffle, '_iter', iter, '.csv')
		read_file_topx = fread(here('topx_check', 'top_hyperparameters', filename_topx))

		filename_chosenx = paste0('nnet-finetune', complex, '-shuf', shuffle, '-', 'iter', iter, '.csv')
		read_file_chosen = fread(here('topx_check', 'hyperparameters_chosen', filename_chosenx))

	#create ranked hyperparameter combinations
		read_file_topx[,rank_topx := frank(mean_eval_loss)]
		read_file_chosen[,rank_chosen := frank(mean_eval_loss)]

		if(complex == TRUE | complex == 'T'){
			comp_lwr = 't'
		} else{
			comp_lwr = 'f'
		}
		res_file_is = fread(here('resources',paste0('sim-results-', complex, '-', shuffle, '-', 'nnetd', '.csv')))
		beta = res_file_is %>% filter(iter == iter) %>% select(estimate)

		jned_topx = merge(x = read_file_chosen, y = read_file_topx, by = c('flag_width1', 'flag_dropout', 'flag_depth', 'flag_batch_size'))

		jned_topx %<>% mutate(num_parameters = ifelse(flag_depth == 1, 101, 
			100*flag_width1 + ifelse(flag_depth-2 > 0,flag_width1*flag_width1*(flag_depth-2),0) + (flag_depth-1)*flag_width1 + 1))

		jned_topx$beta = rep(beta$estimate[iter],2)

		jned_topx %<>% mutate(simple_map = ifelse(flag_width1 == 1, 1, 0), rankpicked = frank(rank_chosen), rankorig = frank(rank_topx),
			picked_worse = ifelse(rankpicked != rankorig, 1, 0), iter = iter)
	}

	else{
		#pick nnet_type

		if(deep & wide){
			nnet_type = 'full'
			nnet_name = 'nnetf'
		} else if(deep){
			nnet_type = 'deep'
			nnet_name = 'nnetd'
		} else if(wide){
			nnet_type = 'wide'
			nnet_name= 'nnetw'
		}


		filename_topx = paste0('top.x_C', complex, '_s', shuffle, '_iter', iter, 'nnet_class_', nnet_type, '.csv')
		read_file_topx = fread(here('topx_check', 'top_hyperparameters', filename_topx))

		filename_chosenx = paste0('nnet-finetune', complex, '-shuf', shuffle, '-', 'iter', iter, '-', nnet_type, '.csv')
		read_file_chosen = fread(here('topx_check', 'hyperparameters_chosen', filename_chosenx))

	#create ranked hyperparameter combinations
		read_file_topx[,rank_topx := frank(mean_eval_loss)]
		read_file_chosen[,rank_chosen := frank(mean_eval_loss)]

		read_file_chosen %<>% mutate(iter = iter)

		if(complex == TRUE | complex == 'T'){
			comp_lwr = 't'
		} else{
			comp_lwr = 'f'
		}



		res_file_is = fread(here('resources',paste0('sim-results-', complex, '-', shuffle, '-', nnet_name, '.csv')))
		beta = res_file_is %>% filter(iter == iter) %>% select(estimate)

		jned_topx = merge(x = read_file_chosen, y = read_file_topx, by = c('flag_width1', 'flag_dropout', 'flag_depth'))

		jned_topx$beta = beta$estimate

		jned_topx %<>% mutate(simple_map = ifelse(flag_width1 == 1, 1, 0), rankpicked = frank(rank_chosen), rankorig = frank(rank_topx),
			picked_worse = ifelse(rankpicked != rankorig, 1, 0), iter = iter)

	}

	return(jned_topx)
	}

training_info_outputter = function(tsol = tsol, deep = TRUE, wide = TRUE, shuffle, complex){

	if(deep & wide){
			nnet_type = 'full'
			nnet_name = 'nnetf'
		} else if(deep){
			nnet_type = 'deep'
			nnet_name = 'nnetd'
		} else if(wide){
			nnet_type = 'wide'
			nnet_name= 'nnetw'
		}

	if(tsol){
		output = lapply(1:1000, training_info_nnet, wide = wide, deep = deep, shuffle = shuffle, complex = complex, temp_solution = tsol) %>% rbindlist()

		if(complex|complex == 'T'){
			comp = 'T'
		} else if(!complex| complex == 'F'){
			comp = 'F'
		}

		fwrite(output, here('traininfo',paste0('training_info_nnet_shuf_', shuffle, '_', 'complex', comp, '_', 'nnetd', '.csv')))


		output = data.table(output)

		print((setkey(output[, list(freq = .N), by =list(flag_width1, flag_depth, flag_dropout)], flag_width1, freq)[J(unique(flag_width1)), mult="last"]))
	
		print((setkey(output[, list(freq = .N), by =list(flag_width1, flag_depth, flag_dropout)], flag_depth, freq)[J(unique(flag_depth)), mult="last"]))

		print((output[,list(
			freq = .N),
		by = 'flag_depth,flag_width1,flag_dropout']))
	}

	else{
		output = lapply(1:1000, training_info_nnet, wide = wide, deep = deep, shuffle = shuffle, complex = TRUE, temp_solution = tsol,) %>% rbindlist()
	
		if(complex|complex == 'T'){
			comp = 'T'
		} else if(!complex| complex == 'F'){
			comp = 'F'
		}

		fwrite(output, here('traininfo',paste0('training_info_nnet_shuf_', shuffle, '_', 'complex', comp, '_', nnet_name, '.csv')))}

		output = data.table(output)

		output %<>% filter(rankpicked == 1)

		print('most common choice of depth given a width:')

		print((setkey(output[, list(freq = .N), by =list(flag_width1, flag_depth, flag_dropout)], flag_width1, freq)[J(unique(flag_width1)), mult="last"]))
		
		print('most common choice of width given a depth:')

		print((setkey(output[, list(freq = .N), by =list(flag_width1, flag_depth, flag_dropout)], flag_depth, freq)[J(unique(flag_depth)), mult="last"]))

		print('full frequency table:')

		print((output[,list(
			freq = .N),
		by = 'flag_depth,flag_width1,flag_dropout']))
	}

	tsol = TRUE

	training_info_outputter(tsol = tsol, deep = TRUE, wide = FALSE, shuffle = '2', complex =  TRUE)
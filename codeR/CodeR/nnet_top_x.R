#-----------------------------------------
#TOPLINE
# Extract top-x performing hyperparameter estimates from tuning directories for each dgp
#-----------------------------------------

#Types of choosing available:

# Top 5 can be measured by 1.) 'mse' (mean squared error, oos), 2.) 'a_rank' (average of inverse rank across all tuning runs) or 3.) probabilistically weighted by the OOS norm(z-score)

# 1.) MSE is calculated, where N = number of observations, y is some continuous prediction target, and y-hat is the regression NN predicted value of y, and i is observation, k
# is fold number, and j is model number

#1.)		MSE = Σ, i = 0 -> N, (y_i - y-hat_i)^2/N
#			returns: min_j(Σ, k = 1 ->5, (MSE_k)/5)

#2.) Uses MSE and ranks (as in the statisical analysis framework) the hyperparameter combinations over foldwise-MSE as defined above. For instance: if MSE in a given fold for combination 1 is 1.2, and 
# 	avg MSE for combination 2 is .79, then combination  1 has rank '2' and combination 2 has rank '1'.

# To choose, every hyperparameter set will be ranked in every fold, and then the best model, according to the model with the highest average inverse rank will be chosen.

# where model rank is r, and model number is j, ξ_j,k is the set of models that produce MSE less than or equal to model j in fold k, and k is fold number

# Rank Metric = max_j(Σ_k(1/||ξ_j,k||)/5)


#-----------------------------------------
# Other arguments:

# iter: iteration number

# num_folds: number of folds for k-fold cval

# shuffle: shuffle strategy

# complex: T/F. Belloni/Non-belloni

# top-x: number of combinations to return

#deep: cross-validate over all model depths or subset to models w/ depth 2 or less

#wide: cross validate over all model widths or subset to models w/ depth 64 or less.
#-----------------------------------------

#-----------------------------------------
#packages
#-----------------------------------------

library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor, magrittr)

#-----------------------------------------

#body

best_x_hyp = function(cv_runs = 25, numfolds = 5, shuffle, complex, top_x = 2, mech = 'randomized_choice', deep = TRUE, wide = TRUE, iter){
	#-----------------------------------------
	# Using iteration #, fold number, shuffle strategy and complex vs. simple dgp, find x-best (as determined by top_x) hyperparameter combinations for a given
	# dataset.
	#
	# 
	#-----------------------------------------

	extract_from_data = function(foldnum, shuffle, complex, iter){
		#utility function

			if(complex == FALSE | complex == 'F'){
				#make sure we can find our data
				shuffle = '1'
				complex = 'F'
			} else{
				complex = 'T'
			}
			
            runname = paste0('nnet-runon-complex', complex, '-shuf', shuffle)
            tunename = paste0('tuningdir-', iter, '-', foldnum)
            runstemp = ls_runs(runs_dir = here('tuningdir', '/',runname, '/', tunename))
            return(runstemp)
        }


    run_df = mapply(extract_from_data, iter = 1:cv_runs, foldnum = 1:numfolds, complex = complex, shuffle = shuffle, SIMPLIFY = FALSE) %>% rbindlist()

    if(deep){
    	max_depth = 6
    } else{max_depth = 2}

    if(wide){
    	max_width = 1024
    } else{max_width =64}

    if(tolower(mech) == 'mse'){
    	top = run_df %>% filter(flag_width1 <= max_width & flag_depth <= max_depth) %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss) %>% head(top_x)
    } else if (tolower(mech) == 'a_rank'){
    	top= run_df %>% filter(flag_width1 <= max_width & flag_depth <= max_depth) %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% arrange(flag_depth, flag_width1) %>% mutate(rank = min_rank(eval_)) %>% summarize(mean_inverse_rank = mean(1/rank)) %>% arrange(desc(mean_inverse_rank)) %>% head(top_x)
    } else if (tolower(mech) == 'randomized_choice'){
    	
    	#generate a z-score over MSE values, find pnorm, randomly select from parameters based on that.

        if(complex){comp = 'T'} else{
            comp = 'F'
        }

        
        if(deep & wide){
            nnet_type = 'full'
        } else if(deep){
            nnet_type = 'deep'
        } else if(wide){
            nnet_type = 'wide'
        }

        top = run_df %>% filter(flag_width1 <= max_width & flag_depth <= max_depth)%>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% summarize(mean_eval_loss = mean(eval_))  %>% ungroup() %>% mutate(z = (mean_eval_loss - mean(mean_eval_loss))/sd(mean_eval_loss))%>% 
        mutate(p = pnorm(z, 0 ,1, lower.tail = FALSE)) %>%
        mutate(p_n = p/sum(p)) %>% arrange(desc(p_n)) %>% mutate(nnet_type = nnet_type)
        
        newfile = paste0('top.x_c', comp, '_s',shuffle, '_iter', iter, 'nnet_class_', nnet_type, '.csv')
        dir.create(here('topx_check'), showWarnings = FALSE)
        dir.create(here('topx_check','top_hyperparameters'), showWarnings = FALSE)
        fwrite(x = top, file = here('topx_check', 'top_hyperparameters', newfile))

    	#sample top_x from the weighted distribution of hyperparameter candidates
    	top = top[sample.int(nrow(top),top_x, replace = FALSE, prob = top$p_n),]



    } else{
    	warning('You didn\'t enter one of the avaialble options: must choose \'mse\' \'a_rank\' or \'randomized_choice\': returning null.')
    	top = NULL
    }

    return(top)
}


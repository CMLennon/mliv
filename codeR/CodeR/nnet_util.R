
library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor, tfestimators, magrittr)


#-----------------------------------------
# Neural Network aux function for study
# Needs to take i_dt_corr (datatable) and a tidymodels folds object and
# crossvalidate to produce best oos predictor model given a set of restrictions.
#
# Choose from widths between 16, 32, 48, 64, ... etc... 256
# restricted versions only allow width to be 16 or 32

# Choose from depths between 1, 2, 3, 4, 5, 6
# restricted versions only allow depth of 1 or 2

#k-fold cv needs to be built in manually - which is challenging over so many models (it's challenging to work with in general, but tf makes it more manual than other procedures)
#-----------------------------------------
#
#options: 1) wide = T allows neural network to crossvalidate on full breadth of width of layers
#		  2) deep = T allows neural network to crossvalidate on a limited set of depths

#			When wide & deep = T, we call the resulting model 'fully xvalidated'

ff_net_cv = function(seed, datatable, folds, iter, wide = T, deep = T, numfolds, newrun = T, shuffle, complex){
    if(newrun){
	   n = nrow(datatable)
    
    	   foldmapper = function(i, folds){
            	fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
            	data.table(fold)
            	fold[,fold.id := i]
            	return(fold)
     	   	}
    
     	  if(wide){
     	  	endw = c(32, 64,256, 512)
     	  } else{endw=c(32, 64)}
    
     	  if(deep){
     	  	endd = 6
     	  } else{endd = 2}
    
    
     	  datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
    
     	  #grep names of instruments, find number
     	  zs = grep("z", names(datatable), value = "T")
     	  num_zs = length(zs)
    
     	  #function to run a set of tuning runs for a single fold in a datatable
     	  single_run = function(foldnum, datatable, iter, num_zs, complex, shuffle){
     	  	#separate training data and turn into matrix form so tf doesn't yell at us
     	  	train_dt = datatable %>% filter(fold.id != foldnum) %>% data.table()
            	train_dat = train_dt[,..zs] %>% as.matrix()
            	train_lab = train_dt[,x1] %>% as.matrix()
            	#separate testing data and """"
            	test_dt = datatable %>% filter(fold.id == foldnum) %>% data.table()
            	test_dat = test_dt[,..zs] %>% as.matrix()
            	test_lab = test_dt[,x1] %>% as.matrix()
    
            	#create temporary directory folder to hold tuning files: one per fold and one per iteration. These will be stored in super-folders that give complex/shuffle information as well

                if(complex){comp = 'T'} else{comp = 'F'}

                runname = paste0('nnet-runon-complex', comp, '-shuf', shuffle)
            	tunename = paste0('tuningdir-', iter, '-', foldnum)
            	tuning_run(here('codeR', 'nnutil.R'), runs_dir = here("tuningdir",runname, tunename), flags = list(
            		width1 = c(16, endw),
            		dropout = c(.1,.2),
            		depth = seq(from = 1, to = endd, by = 1),
            		epoch_num = c(40),
            		batch_size = c(10),
            		activationfcn = c('relu')
            		),
            	confirm = FALSE)
            	
            	runstemp = data.frame(ls_runs(runs_dir = here("tuningdir",runname, tunename)))
            	
            	#delete all run stats from disk storage
            	#unlink(paste0(here("tuningdir",tunename),"/*"), recursive = TRUE)
            	return(runstemp)
     	  }
    
     	  #use lapply to build a full list of run statisitcs
     	  run_df_list = lapply(1:numfolds, single_run, datatable = datatable, iter = iter, num_zs = num_zs, complex = complex, shuffle = shuffle)
     	  run_df = rbindlist(run_df_list)
            rm(run_df_list)
    
    
    } else{
        extract_from_data = function(iter, fold_num, shuffle, complex){
            if(complex){comp = 'T'} else{comp = 'F'}

            runname = paste0('nnet-runon-complex', comp, '-shuf', shuffle)
            tunename = paste0('tuningdir-', iter, '-', foldnum)
            runstemp = data.frame(ls_runs(runs_dir = here("tuningdir", runname, tunename)))
            return(runstemp)
        }

    run_df_list = lapply(1:numfolds, extract_from_data, iter = iter, complex = complex, shuffle = shuffle)
    run_df = rbindlist(run_df_list)
    rm(run_df_list)
}

    #find best runs (as defined as lowest MSE) under all three model conditions
    best_mod = run_df %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% 
    summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss) %>% head(1)

    best_modw = run_df %>% filter(flag_depth < 3) %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% 
    summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss) %>% head(1)

    best_modd = run_df %>% filter(flag_width1 < 65) %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% 
    summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss) %>% head(1)

    
    #create three named lists to hold hyperparameters for each subtype of model.

    nnfull = list(width = best_mod$flag_width1, dropout = best_mod$flag_dropout, depth = best_mod$flag_depth, epoch_num = best_mod$flag_epoch_num, 
        batch_size = best_mod$flag_batch_size, activationfcn = best_mod$flag_activationfcn)

    nnwide = list(width = best_modw$flag_width1, dropout = best_modw$flag_dropout, depth = best_modw$flag_depth, epoch_num = best_modw$flag_epoch_num, 
        batch_size = best_modw$flag_batch_size, activationfcn = best_modw$flag_activationfcn)

    nndeep = list(width = best_modd$flag_width1, dropout = best_modd$flag_dropout, depth = best_modd$flag_depth, epoch_num = best_modd$flag_epoch_num, 
        batch_size = best_modd$flag_batch_size, activationfcn = best_modd$flag_activationfcn)

    #return a named list of named lists to hold all hyperparameters for use in final model

    return(list(nnfull = nnfull, nnwide = nnwide, nndeep = nndeep))

}

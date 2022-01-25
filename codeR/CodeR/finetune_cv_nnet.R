library(pacman)

p_load(tidyverse, here, data.table, tidymodels, glmnet, tidyverse, keras, tfruns, tensor, magrittr, R.utils)


finetune_cv = function(top_x_obj,datatable, folds, iter, numfolds, shuffle, complex){

        p_load(R.utils)
    
	       n = nrow(datatable)
    
    	   foldmapper = function(i, folds){
            	fold = assessment(folds$splits[[i]]) %>% as.data.frame() %>% data.table()
            	data.table(fold)
            	fold[,fold.id := i]
            	return(fold)
     	   	}
    
     	  datatable = lapply(c(1:5), foldmapper, folds) %>% rbindlist() %>% data.table()
    
     	  #grep names of instruments, find number
     	  zs = grep("z", names(datatable), value = "T")
     	  num_zs = length(zs)
    
     	  #function to run a set of tuning runs for a single fold in a datatable
     	  single_run = function(top_x, foldnum, datatable, iter, num_zs, complex, shuffle){

                nnet_type = top_x$nnet_type[1]
     	  	   
               #separate training data and turn into matrix form so tf doesn't yell at us
     	  	   train_dt = datatable %>% filter(fold.id != fold) %>% data.table()
            	train_dat = train_dt[,..zs] %>% as.matrix()
            	train_lab = train_dt[,x1] %>% as.matrix()
            	#separate testing data and """"
            	test_dt = datatable %>% filter(fold.id == fold) %>% data.table()
            	test_dat = test_dt[,..zs] %>% as.matrix()
            	test_lab = test_dt[,x1] %>% as.matrix()
    
            	#create temporary directory folder to hold tuning files: one per fold and one per iteration. These will be stored in super-folders that give complex/shuffle information as well

                if(complex){comp = 'T'} else{comp = 'F'}

                hyp = top_x %>% select_if(grepl('flag_', names(.))) %>% unique()
                #print(hyp)

                #strip 'flag_' from names sonames(hyp) = gsub(names(hyp), pattern = 'flag_', replacement = '')
                names(hyp) = gsub(x=names(hyp), pattern='flag_')

                #programatically build new tuningrun names
                runname = paste0('nnet-finetune', comp, '-shuf', shuffle)
                tunename = paste0('tuningdir-finetune', iter, '-', foldnum, '-', nnet_type)


                map_hyp_to_tuning_run= function(hyp, hyp_row_num, runname, tunename){

                #takes list of parameters and evaluates the fine-tuning grid
                hyp = hyp[hyp_row_num,] %>% expand.grid()
                tuning_run(here('codeR', 'nnutil.R'), runs_dir = here("tuningdir",runname, tunename), flags = hyp,
                confirm = FALSE, echo = FALSE)

                
                runs = data.frame(ls_runs(runs_dir = here("tuningdir",runname, tunename)))
                
                #delete all run stats from disk storage
                storage_dir = here("tuningdir",paste0(runname, '_stored'), tunename)
                R.utils::copyDirectory(here("tuningdir",runname, tunename), storage_dir)

                unlink(paste0(here("tuningdir",runname, tunename),"/*"), recursive = TRUE)
                return(runs)
            }

                #utility fcn to map rows to individual runs
                
                runstemp = lapply(1:nrow(hyp), map_hyp_to_tuning_run, hyp = hyp, runname = runname, tunename = tunename) %>% rbindlist(fill = TRUE)

            	return(runstemp)
     	  }
    
     	  #call utility function fo find top-5 models by OOS MSE
     	  run_df = lapply(X  = 1:5, FUN = single_run, datatable = datatable, iter = iter, num_zs = num_zs, complex = complex, shuffle = shuffle, top_x = top_x_obj) %>% rbindlist(fill = TRUE)

    #find best runs (as defined as lowest MSE) under top_x model conditions

    nnet_type = top_x_obj$nnet_type[1]

    best_mod_ret = run_df %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% 
    summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss)

    if(complex){
        comp = 'T'
    } else{
        comp = 'F'
    }

    runname_out = paste0('nnet-finetune', comp, '-shuf', shuffle, '-iter', iter, '-', nnet_type, '.csv')

    dir.create(here('topx_check','hyperparameters_chosen'), showWarnings = FALSE)
    fwrite(file = here('topx_check','hyperparameters_chosen',runname_out), x = best_mod_ret)

    best_mod = run_df %>% group_by(flag_width1, flag_dropout, flag_depth, flag_epoch_num, 
        flag_batch_size, flag_activationfcn) %>% 
    summarize(mean_eval_loss = mean(eval_)) %>% 
    arrange(mean_eval_loss) %>% head(1)

    
    #create three named lists to hold hyperparameters for each subtype of model.

    nn_best = list(width = best_mod$flag_width1, dropout = best_mod$flag_dropout, depth = best_mod$flag_depth, epoch_num = best_mod$flag_epoch_num, 
        batch_size = best_mod$flag_batch_size, activationfcn = best_mod$flag_activationfcn)

    #return a named list of named lists to hold all hyperparameters for use in final model

    return(nn_best)

}


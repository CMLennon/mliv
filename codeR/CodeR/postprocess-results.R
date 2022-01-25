
# Function: Post-process simulation data -------------------------------------------------
# Arguments:
# 	- 'sim_dt': Output from the simulation (data.table)
# 	- 'is_belloni': Did sim_dt come from Belloni et al DGP? (T or F)
# 	- 'n_dgp': Which 'number' of the DGP options was used? (numeric)
	# The function
	post_process = function(sim_dt, is_belloni, n_dgp) {
		require(data.table, magrittr)
		# Define crosswalk between methods' full names and their abbreviations
		method_xwalk = rbindlist(list(
			data.table(abbr = "dgp", model_short = "oracle model"),
			data.table(abbr = "ovb", model_short = "naive ols"),
			data.table(abbr = "2sls", model_short = "2sls"),
			data.table(abbr = "lasso", model_short = "lasso selection"),
			data.table(abbr = "plasso", model_short = "post-lasso selection"),
			data.table(abbr = "rf", model_short = "random forest"),
			data.table(abbr = "rfcv", model_short = "random forest, cv"),
			data.table(abbr = "pca", model_short = "pca"),
			data.table(abbr = "gboost", model_short = "boosted trees"),
			data.table(abbr = "split", model_short = "ussiv"),
			data.table(abbr = "jive", model_short = 'jive'),
			data.table(abbr = "liml", model_short = 'liml (fuller)'),
			data.table(abbr = "nnetf", model_short = 'unrestricted neural net'),
			data.table(abbr = "nnetw", model_short = 'shallow neural net'),
			data.table(abbr = "nnetd", model_short = 'narrow neural net'),
			data.table(abbr = "nnetfos", model_short = 'unrestricted neural net oos'),
			data.table(abbr = "nnetwos", model_short = 'shallow neural net oos'),
			data.table(abbr = "nnetdos", model_short = 'narrow neural net oos')
		))
		# Generate shortened names for models
		sim_dt[, `:=`(
			model_short = model %>% stringr::str_to_lower() %>% str_remove_all("first stage: ")
		)]
		# Find the methods used to generate 'sim_dt'
		methods_used = sim_dt[,model_short]
		# Iterate over methods used, saving each method's results
		blah = lapply(
			X = methods_used,
			FUN = function(m) {
				# Create a filename using the DGP information and the method type
				m_filename = paste0(
					"sim-results-",
					ifelse(is_belloni, "t", "f"), "-",
					n_dgp, "-",
					method_xwalk[model_short == m, abbr],
					".csv"
				)
				# Save results from method 'm'
				fwrite(
					x = sim_dt[model_short == m, -"model_short"],
					file = here("Resources", m_filename)
				)
			}
		)
	}

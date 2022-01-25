

# Notes ----------------------------------------------------------------------------------
#   Goal: Produce figure on overfit and part a, B*cov(x-hat, e)


# Setup ----------------------------------------------------------------------------------
  # Packages
  library(pacman)
  p_load(
    MASS, tidyverse, tidymodels, 
    patchwork, latex2exp, extrafont,
    data.table, 
    parallel, future, future.apply,
    viridis,
    here, magrittr
  )


# Function: Random forest CV data --------------------------------------------------------
	produce_plot_grid = function(beta = 1) {
		# Set up multi-core processing
		plan(multicore, workers = detectCores() - 1)
    # Source the Belloni DGP function
    source(here("codeR", "Belloni_aux_function_shuffle.R"))
    # Generate data
    i_dt_corr = data_gen(2000,100,shuffle = '1')
    # Grab training data
    datatable = i_dt_corr[1:1000,]
    # Set up CV
    folds = vfold_cv(datatable, v = 5)
    # Grab test data
    datatable_test = i_dt_corr[1001:2000,]
    # Instruments
    zs = grep("z", names(datatable), value = "T")
    num_zs = length(zs)
    form = paste0("x1 ~ ", paste0(zs, collapse = '+')) %>% as.formula()
    # RF recipe
    forest_recipe = recipe(
      form, 
      data = datatable
    )
    # Set model, make the models tunable onnumber of trees, mtry, smallest leaf
    forest_model = rand_forest(others = list(inbag = 1000)) %>%
      set_engine("ranger") %>% 
      set_args(
        trees = tune(),
        mtry = tune(),
        min_n = tune()
      ) %>% set_mode("regression")
    # Forest workflow
		forest_wf = workflow() %>% add_recipe(forest_recipe) %>% add_model(forest_model)
		# Build forests
		build_forests = function(min_n) {
			hyparam = expand_grid(trees = 500, mtry = 100, min_n = min_n)
			# Emtpy output list
			model_out_data = list()
			# Temporary data
			dat_temp = datatable
			dattest_temp = datatable_test
			# Forest workflow
			forest_wf_fin = forest_wf %>% finalize_workflow(hyparam)
			# Fit forest
			forest_fit = parsnip::fit(forest_wf_fin, data = datatable)
			# In/out of sample predictions
			dat_temp$train_out = predict(forest_fit, new_data = dat_temp)
			dattest_temp$test_out = predict(forest_fit, new_data = dattest_temp)
			# Name for output
			outname = paste0("model_", min_n)
			# List to output
			model_out_data = list(dat_temp, dattest_temp, min_n)
			# Return output
			return(model_out_data)
		}
		# Run the function in parallel
		model_out_data = future_lapply(
			seq(from = 1, to = nrow(datatable)/2, by = 2),
			build_forests
		)
		# print(model_out_data[[1]])
		# Build bias dataset
		build_bias_data = function(model_data){
			is_dat = model_data[[1]]
			os_dat = model_data[[2]]
			min_n = model_data[[3]]
			# MSE
			mse_is = mean((is_dat$x1 - is_dat$train_out)^2)
			mse_os = mean((os_dat$x1 - os_dat$test_out)^2)
			# Define flexibility
			flexibility = 500 - min_n
			# Covariance, variance, and bias
			cov_e = cov(is_dat$train_out, is_dat$x1 - is_dat$train_out)
			cov_u = cov(is_dat$train_out, is_dat$en1)
			var_x = var(is_dat$train_out)
			bias = (beta*cov_e + cov_u)/(var_x)
			# Data table to output
			tmp = data.table(min_n = min_n,
				cov_e = cov_e,
				cov_u = cov_u,
				insample_mse = mse_is,
				outsample_mse = mse_os,
				variance_xhat = var_x,
				variance_x = var(is_dat$x1),
				calculated_bias = bias,
				model_flexibility = flexibility
			)
			return(tmp)
		}
		plot_data = future_lapply(model_out_data, build_bias_data) %>% rbindlist()
		return(plot_data)
	}


# Create and save data -------------------------------------------------------------------
	# Run the function
	outdata = produce_plot_grid()
	# Save as CSV
	write_csv(
		outdata,
		path = here("Resources", "example-data", "plot-cv2.csv")
	)


# Load data -------------------------------------------------------------------------------
	# Load data
	cv_dt = here("Resources", "example-data", "plot-cv.csv") %>% fread()


# Create figures -------------------------------------------------------------------------
	# Vector of colors
	col_v = magma(6, begin = 0, end = 0.9)
  # In-/out-of-sample MSE
	gg1 = ggplot(
		data = melt(
			data = cv_dt[, .(
				f = model_flexibility,
				mse_i = insample_mse,
				mse_o = outsample_mse
			)],
			id.vars = "f"
		),
		aes(x = f, y = value, color = factor(variable, levels = c("mse_o", "mse_i"), ordered = T))
	) +
	geom_line() +
  scale_x_continuous("Flexibility") +
  scale_y_continuous("") +
	scale_color_manual(
		"MSE",
		values = col_v[c(3, 5)],
		labels = c("Out of sample", "In sample")
	) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  )
	# Bias components
	gg2 = ggplot(
		data = melt(
			data = cv_dt[, .(
				f = model_flexibility,
				a = cov_e,
				b = cov_u,
				"1/c" = variance_xhat
			)],
			id.vars = "f"
		),
		aes(x = f, y = value, color = variable)
	) +
	geom_line() +
  scale_x_continuous("Flexibility") +
  scale_y_continuous("") +
	scale_color_manual(
		"Bias components",
		values = col_v[c(1, 4, 6)],
		labels = list(
      TeX("$a =\\beta\\mathrm{Cov}(\\hat{x},e)$"),
      TeX("$b=\\mathrm{Cov}(\\hat{x},u)$"),
      TeX("$1/c=\\mathrm{Var}(\\hat{x})$")
    )
	) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(legend.text = element_text(family = "Georgia"))
	# Bias
	gg3 = ggplot(
		data = melt(
			data = cv_dt[, .(
				f = model_flexibility,
				bias = calculated_bias
			)],
			id.vars = "f"
		),
		aes(x = f, y = value, color = variable)
	) +
	geom_line() +
  scale_x_continuous("Flexibility") +
  scale_y_continuous("") +
	scale_color_manual(
		TeX("Bias: $(a+b)c$"),
		values = col_v[2],
		labels = list(TeX("$\\beta - \\hat{\\beta}$"))
	) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(legend.text = element_text(family = "Georgia"))
  # Align the plots
  plots = align_plots(list(gg1, gg2, gg3))
  # Save
  ggsave(
    plot = plots[[1]],
    path = here("Figures"),
    filename = "cv-bias-rf-1.pdf",
    device = cairo_pdf,
    width = 6.5,
    height = 2.25
  )
  ggsave(
    plot = plots[[2]],
    path = here("Figures"),
    filename = "cv-bias-rf-2.pdf",
    device = cairo_pdf,
    width = 6.5,
    height = 2.25
  )
  ggsave(
    plot = plots[[3]],
    path = here("Figures"),
    filename = "cv-bias-rf-3.pdf",
    device = cairo_pdf,
    width = 6.5,
    height = 2.25
  )


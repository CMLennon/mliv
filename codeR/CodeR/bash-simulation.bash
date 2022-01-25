#!/bin/bash

# Notes

# 	- Post-script argument 1: complex: T/F. If true, this will run the belloni 100-instrument dataset. If 'F', runs 7 instrument dataset.
# 	- Post-script argument 2: Shuffle: options = (1,2,3). Only affects the complex dataset in any way. This is to isolate PCA performance changes due to the structure of the correlation in the instruments.
#			- "1": shuffle the beta-pattern to be random.
#			- "2": beta-pattern follows... instrument 1 = strongest instrument, 2 = 2nd strongest instrument, 3 = 3rd ...
#			- "3": beta-pattern is similar to above except instrument 50 = strongest instrument, instrument 51 = 2nd strongest, instrument 52 = 3rd... And, instrument 49 = weakest, 48 = 2nd weakest ...
#   - Post-script argument 3: number of observation per simulation (e.g., 1000).
#   - Post-script argument 4: number of simulations to run (e.g., 1000)
#   - Post-script argument 5: logical T/F. Whether or not we want to build the dataset.
#		- Beyond argument 5, follows a list of 'not-run' models. This will allow the user to disable certain models that will not be run. Available options here are:
# 		- "dgp" : run a model that uses x without shared error term
# 		- "ovb" : run a naive regression y ~ x1
# 		- "2sls": run two stage least squares
# 		- "lasso": run lasso on the first stage, with predicted results in 2nd stage
# 		- "plasso": as above, but run regression using all non-zero coefficients from a lasso model (post lasso method)
# 		- "rf": run a non-crossvalidated model in ranger (using ranger defaults) for first stage
# 		- "rfcv" : cross-validate first stage random forest.
# 		- "pca": build a set of synthetic instruments using priciple component analysis, cross-validated against x-hat/x mse.
# 		- "gboost": xgboost algorithm to estimate the first stage, cross-validated.
# 		- "split": split-sample 2SLS
#     -"JIVE": jackknife IV
#     -"LIML": limited-information maximum likelihood
#     - "nnetw": Neural net, wide
#     - "nnetd": Neural net, deep
#     - "nnetf": Neural net, full

# Note: If methods include any "nnet", then the script runs sequentially. Otherwise, the script runs in parallel (multicore) on all cores.

# Full list:
# dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML nnetf nnetw nnetd

# Run to generate data
# time Rscript --vanilla bash-friendly-sim.R F 1 1000 1000 T ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML nnetd nnetw nnetf

# Run for estimation

# Non-Belloni (type 1)
echo -e "F1\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla bash-friendly-sim.R F 1 1000 1000 F dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Belloni type 1
echo -e "T1\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla bash-friendly-sim.R T 1 1000 1000 F dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Belloni type 2
echo -e "T2\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla bash-friendly-sim.R T 2 1000 1000 F dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Belloni type 3
echo -e "T3\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla bash-friendly-sim.R T 3 1000 1000 F dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

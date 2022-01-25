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
#   - Post-script argument 6: integer 1-4 "data_modification." What manner does our nonlinearity arrive in?

#     - "1": sin-threshold-strategy. If sin of z1 exceeds threshold, it will add either .1 or -.1 to x and .5 or -.5 to y. This is learnable bad #variation.
#     - "2": Noise induced by nonlinear instrument. Exclusion restriction is violated by some power of an instrument(z1). Checks -4,4. Only whole #numbers are valid due to z1 being sometimes <0.
#     - "3": Noise induced by interaction of instruments. Checks interaction for up to 7 instruments.
#     - "4": As 1, but check the sum of z1-z4 instead of z1 alone.



#		- Beyond argument 6, follows a list of 'not-run' models. This will allow the user to disable certain models that will not be run. Available options here are:
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
#     -"jive": jackknife IV
#     -"liml": limited-information maximum likelihood
#     - "nnetw": Neural network


# Note: If methods include any "nnet", then the script runs sequentially. Otherwise, the script runs in parallel (multicore) on all cores.

# Full list:
# dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML nnetw

#time Rscript --vanilla bash-friendly-sim.R F 1 1000 1000 T 2 dgp ovb 2sls lasso plasso rf rfcv pca gboost split JIVE LIML nnetw

# Run for estimation

# # Non-Belloni (type 1 noise)
# echo -e "F1\n\nStarted: \n$(date +'%r\n%D')\n"
# time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 1 nnetw
# echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# # Non-Belloni (type 2 noise)
# echo -e "F2\n\nStarted: \n$(date +'%r\n%D')\n"
# time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 2 nnetw
# echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Non-Belloni (type 3 noise)
echo -e "F3\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 3 nnetw
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# # Non-Belloni (type 4 noise)
# echo -e "F4\n\nStarted: \n$(date +'%r\n%D')\n"
# time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 4 nnetw
# echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Belloni (type 1 noise)
echo -e "T1-1\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 1 nnetw
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# Belloni (type 3 noise)
echo -e "T1-1\n\nStarted: \n$(date +'%r\n%D')\n"
time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 3 nnetw
echo -e "Ended: \n$(date +'%r\n%D')\n\n"

# # Belloni (type 4 noise)
# echo -e "T3-1\n\nStarted: \n$(date +'%r\n%D')\n"
# time Rscript --vanilla codeR/nonlinear-results-new.R F 1 1000 1000 F 4 nnetw
# echo -e "Ended: \n$(date +'%r\n%D')\n\n"

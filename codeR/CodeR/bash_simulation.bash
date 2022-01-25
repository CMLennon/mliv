#Notes:
# Post-script argument 1: Complex: T/F. If true, this will run the belloni 100-instrument dataset (with 1000 observations.) If False, runs 7 instrument dataset (1000 observations)
#
# Post-script argument 2: Shuffle: options = (1,2,3). Only affects the complex dataset in any way. This is to isolate PCA performance changes due to the structure of the correlation in the instruments.
#
### 1: shuffle the beta-pattern to be random. 
#
### 2: beta-pattern follows... instrument 1 = strongest instrument, 2 = 2nd strongest instrument, 3 = 3rd ...
#
### 3: beta-pattern is similar to above except instrument 50 = strongest instrument, instrument 51 = 2nd strongest, instrument 52 = 3rd... And, instrument 49 = weakest, 48 = 2nd weakest ...

# Post-script argument 3: number of observation per simulation (normally 1000)

# Post-script argument 4: number of simulation iterations to run (normally 500)

  
#Post-script argument 5: data_transform - this can take on one of three values
    #1. - Slight learnable disturbance. This adds the noise term if(abs(sin(z1)) > .995), sign(sin(z1))*.1, else 0     #     twice over to y and once to x
    #2. - Quadratic form of bad-variation.
    #3. - Interaction form of bad-variation.

## Beyond argument 4, follows a list of 'NOT-RUN' models. This will allow the user to disable certain models that will not be run. Available options here are:

### "dgp" : run a model that uses x without shared error term
### "ovb" : run a naive regression y ~ x1
### "2sls": run two stage least squares
### "lasso": run lasso on the first stage, with predicted results in 2nd stage
### "plasso": as above, but run regression using all non-zero coefficients from a lasso model (post lasso method)
### "rf": run a non-crossvalidated model in ranger (using ranger defaults) for first stage
### "rfcv" : cross-validate first stage random forest.
### "pca": build a set of synthetic instruments using priciple component analysis, cross-validated against x-hat/x mse.
### "gboost": xgboost algorithm to estimate the first stage, cross-validated.
### "split": split-sample 2SLS
### "jive": jackknife IV
### "liml": limited-information maximum likelihood
### "nnetw": neural network

# simpler data
Rscript bash-friendly-sim.R F 1 1000 500 dgp ovb 2sls lasso plasso rf rfcv pca gboost split nnetw nnetd nnetf

# recreated belloni dataset
Rscript bash-friendly-sim.R T 1 1000 500 dgp ovb 2sls lasso plasso rf rfcv pca gboost split nnetw nnetd nnetf

Rscript bash-friendly-sim.R T 2 1000 500 dgp ovb 2sls lasso plasso rf rfcv pca gboost split nnetw nnetd nnetf

Rscript bash-friendly-sim.R T 3 1000 500 dgp ovb 2sls lasso plasso rf rfcv pca gboost split nnetw nnetd nnetf

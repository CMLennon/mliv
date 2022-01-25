generate_system = function(i_dt, n){
  #Utility function to generate our system of equations so 
  #that it can be done for multiple data sets
  β0 = 1
  β1 = 1
  β2 = 1
  num_zs = 7
  #corr = .2
  # z affects x1 (plus correlated noise)
  i_dt[, x1 := 1 + z_1 + z_2 + z_3 + z_4 + z_5 + z_6 + z_7 + e_common]
  # z^2 and affects x2 (plus correlated noise)
  # i_dt[, x2 := 1 + 10 * z^2 + y_err_full]
  i_dt[, x2 := 1 + e_common]
  # Calculate y
  i_dt[, y_err := rnorm(n)]
  i_dt[, y := β0 + β1 * x1 + β2 * x2 + y_err]
  i_dt[, y_err_full := y_err + e_common]
  i_dt[, true_x := x1]
  return(i_dt)
}

data_gen = function(n,k, shuffle = '1',norm= F,correlated= F){
  num_zs = k
  corr = .6
  #set up a covar mat with off-diagonals = corr^(|j-h|) where j is own-z (ie, for z2 = 2) and h is the correlated variable's z. Eg. z1 and z3 have corr = corr^2 and z1 has corr = corr^0=1.
  #cvmat = matrix(c(rep(c(1,rep(corr, num_zs)), num_zs-1)), nrow = num_zs, byrow = TRUE) {archived to have more variation in correlation levels across zs}
  
  ######USED FOR OLD SIMULATIONS -no correlation in new instruments. If -norm is specified the resulting instruments will have correlation implemented
  cvmat = matrix(0,num_zs,num_zs)
  for (j in (1:num_zs)) {
    for (h in (1:num_zs)) {
      cvmat[c(j),c(h)] = corr^abs(j-h)
    }
  }
  #######
  znames = sapply(1:num_zs, FUN = function(x) paste0('z_',as.character(x)))
  
  #pick normal correlated version or uniform normal
  if(correlated == T){
    i_dt_corr = data.table(
    mvrnorm(n = n, mu = rep(0, num_zs), Sigma = cvmat),
    rnorm(n))
  }
  else{
    i_dt_corr = data.table(rnorm(n))
    for(z in znames){
      i_dt_corr[,(z) := rnorm(n)]
    }
    i_dt_corr= i_dt_corr%>% relocate(V1, .after = last_col())
  }

  

  #replace names in multivariate matrix with z-names.
  i_dt_corr = setnames(i_dt_corr, c(znames, "e_common"))

  #replicate Ed's code using the new utility function (which is really just Ed's code)

  #generate system for correlated system, then split for cv process.
  i_dt_corr = generate_system(i_dt_corr, n)
  i_dt_corr[, id := c(1:n)]

  return(i_dt_corr)
}

data_gen = function(n, k, shuffle = '1', z_prob = FALSE){
  #n = number of observations to generate
  #k = number of instruments to generate
  #shuffle: takes values of '1', '2' ,'3'
  #z_prob: Make zs and part of correlated error binary?

    # '1' - shuffle beta pattern
    # '2' - set strongest beta at z_1 (Belloni reproduction)
    # '3' - set strongest beta at z_50, descends from there

  set_C <- function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
    #solving for a C value that changes behavior of lasso selectors.
    C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
    return(C)
  }
  #Generate two errors according to Bellloni (pg 26)

  #first need covariance and variance set, to produce covar matrix
  musq = 1250
  var_er = 1
  cor_er = .5
  var_inst = 1

  #produce Sigma Z. variance is 1, correlation (thus cov) is equal to .5^(|j-h|) for z_j and z_h
  Sig_z = matrix(0,nrow = k,ncol = k)
  for (j in (1:k)) {
    for (h in (1:k)) {
      Sig_z[c(j),c(h)] = .5^abs(j-h)
    }
  }
  z_n1 = mvrnorm(n = n, mu = rep(0,k), Sigma = Sig_z, empirical = TRUE)
  if(z_prob){
    z_n1 = t(apply(z_n1,1, function(x)(x-min(x))/(max(x) - min(x)))) %>% round()
  }
  

  #using exponential
  beta_pattern <-  unlist(lapply(c(0:(k-1)), function(x) (.7)^x), use.names = FALSE)


  Sig_z <- as.matrix(Sig_z)
  if(shuffle == '1'){
    beta_pattern <- as.vector(beta_pattern)
    beta_pattern = sample(beta_pattern)
  } else if(shuffle == '2'){
    beta_pattern <- as.vector(beta_pattern)
  } else{
    beta_pattern <- as.vector(beta_pattern)
    beta_pattern = c(tail(beta_pattern, -50), beta_pattern[1:50])
  }

  C = set_C(musq=musq, n = n, beta_pattern = beta_pattern,Sig_z = Sig_z)
  C = c(C)

  #now that C is known, set variance for error term on X

  sigma_v = abs(1 - C*crossprod(beta_pattern,Sig_z) %*% beta_pattern)

  #We can also set our full pi matrix to find our true instrument coefficients
  betas <- as.vector(C*beta_pattern)
  cov_er<-cor_er*sqrt(sigma_v)

  #use multivariate distribution given matrix above to get a matrix of values for sample
  #size of 'n1'
  covarmat <- matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)

  #we now need to get and separate our error terms
  evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
  
  if(z_prob){
    cov_er<-.999*sqrt(sigma_v)
    covarmat <- matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)
    evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
    evmat_n1 = t(apply(evmat_n1,1, function(x)(x-min(x))/(max(x) - min(x)))) %>% round()
  }
  
  errors = as.data.frame(evmat_n1)
  names(errors) <- c("e", "v")
  
  #separate errors
  en1 = evmat_n1[,1]
  vn1 = evmat_n1[,2]

  #let's construct x's and y's and add them to our z's to create a full data matrix

  X = z_n1 %*% betas + vn1
  #x created
  #x_true
  x_true = X - vn1
  Y = X + en1 + 1
  #y created
  znames = paste(rep("z", k), c(1:k), sep = "_")
  z = as.data.frame(z_n1)
  colnames(z) <- znames

  #X (or d, as in Belloni paper) name
  X = as.data.frame(X)

  colnames(X) <- "X"

  #Target value Y name
  Y = as.data.frame(Y, names = "Y")
  colnames(Y) <- "Y"

  #final matrix
  #true_x = X - vn1
  product<-cbind(z,X[,1],Y,errors, x_true)
  product = data.table(product)
  
  product = setnames(product, c(znames, "x1", "y", "en1", 'vn1', 'true_x'))
  product[,id := c(1:n)]
  return(product)
}

utility = function(x,y) {
  for(item in x){
    for(item in y){
      
    }
  }
}


#fig <- plot_ly() %>% add_surface(z = ~utilgrid)


#fig <- fig %>% layout(
#  scene = list(
#    camera=list(
#      eye = list(x=1.87, y=0.88, z=-0.64)
#    )
#  )
#)
#fig <- plot_ly() %>% add_surface(x = utilgrid[,1], y = utilgrid[,2], z = ~utilgrid, type = '3dmesh')
#plot_ly(z = ~utilgrid) %>% add_surface()

#p <- layout(fig, scene = list(xaxis = list(title = "good a", range = c(0,10)), yaxis = list(title = "good b", range = c(0,10)), zaxis = list(title = "utility")))



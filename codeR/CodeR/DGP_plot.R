## Build Correlation and Beta-strength diagrams


#modify data_gen to return betas and dataset

#-----------------------------------------
#simple case
seed = 12345

set.seed(seed)
k = 7
source(here("codeR", "simple_data_function.R"))
i_dt_corr = data_gen(n,k,shuffle)
betas = rep(1, 7)

z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]

betadata = data.frame(z_names = as.numeric(stringr::str_extract(z_names, "\\d+")), beta = betas)

i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]

corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()

corr_mat_dt$corrwith = rownames(corr_mat_dt)

corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))

corr_mat_dt %<>% mutate(value = ifelse(is.na(value), 0, value)) %>% mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))

corr_mat_dt = corr_mat_dt %>% mutate(corrwith =as.numeric(stringr::str_extract(corrwith, "\\d+")), variable =  as.numeric(stringr::str_extract(variable, "\\d+")), corr = ifelse(corrwith > variable, value, NA))

corrplot1 = ggplot(corr_mat_dt, aes(x = corrwith, y = variable, fill = corr)) + geom_raster() + scale_fill_viridis(option = 'magma', na.value = 'white', begin = 0, direction = -1) +
labs(title = "Correlation heatmap for instrumental variables", fill = 'correlation coefficient', subtitle = '\'simple case\' dgp') +
xlab('instrument number') +
ylab('instrument number') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none') +
coord_fixed()

betaplot1 = ggplot(betadata, aes(x = z_names, y = beta, fill = beta)) + geom_bar(stat = 'identity') + scale_fill_viridis(option = 'magma', na.value = 'white') + 
labs(title = "Magnitude of first-stage coefficients", fill = 'first stage coefficient on instrument x', subtitle = 'first stage coefficients, simple case') + xlab('instrument number') + ylab('coefficient of z on x') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none')

simplecase = plot_grid(corrplot1, betaplot1, nrow = 1)

#-----------------------------------------
#Complex case, shuffled betas
n = 1000
k = 100
shuffle = '1'
 #n = number of observations to generate
 #k = number of instruments to generate
 #shuffle: takes values of '1', '2' ,'3'

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
  

i_dt_corr = product

z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]

betadata = data.frame(z_names = as.numeric(stringr::str_extract(z_names, "\\d+")), beta = betas)

i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]

corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()

corr_mat_dt$corrwith = rownames(corr_mat_dt)

corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))

corr_mat_dt %<>% mutate(value = ifelse(is.na(value), 0, value)) %>% mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))

corr_mat_dt = corr_mat_dt %>% mutate(corrwith =as.numeric(stringr::str_extract(corrwith, "\\d+")), variable =  as.numeric(stringr::str_extract(variable, "\\d+")), corr = ifelse(corrwith > variable, value, NA))

corrplot2 = ggplot(corr_mat_dt, aes(x = corrwith, y = variable, fill = corr)) + geom_raster() + scale_fill_viridis(option = 'magma', na.value = 'white', begin = 0, direction = -1) + 
scale_x_continuous(limits = c(1,100), expand = c(0,0)) + 
scale_y_continuous(limits = c(1,100), expand = c(0,0)) +
labs(title = "Correlation heatmap for instrumental variables", fill = 'correlation coefficient', subtitle = 'all complex data generating processes') +
xlab('instrument number') +
ylab('instrument number') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none')+
coord_fixed()

betaplot2 = ggplot(betadata, aes(x = z_names, y = beta, fill = beta)) + geom_bar(stat = 'identity') + scale_fill_viridis(option = 'magma', na.value = 'white') + scale_x_continuous(limits = c(0,101), expand = c(0,0)) + 
labs(title = "Magnitude of first-stage coefficients", fill = 'first stage coefficient on instrument x', subtitle = 'randomized first stage coefficients, single draw, first complex dgp') + xlab('instrument number') + ylab('coefficient of z on x') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none')


corrline = ggplot(corr_mat_dt %>% filter(variable %in% c(1,50)), aes(x = corrwith, y = value, color = as.factor(variable))) + geom_line() + scale_color_viridis(option = 'magma', discrete = T, begin = .3, end = .9) + scale_x_continuous(limits = c(1,100), expand = c(0,0)) +
labs(title = "Correlation coefficients for instruments 1 and 50", color = 'instrument', subtitle = 'against all other instruments') + xlab('instrument number') + ylab('correlation coefficient') +
theme_minimal(
        base_size = 7.5
      )

#-----------------------------------------
#Complex case, baseline belloni (strongest instrument = z1, z2, z3 ...)
n = 1000
k = 100
shuffle = '2'
 #n = number of observations to generate
 #k = number of instruments to generate
 #shuffle: takes values of '1', '2' ,'3'

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

i_dt_corr = product

z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]

betadata = data.frame(z_names = as.numeric(stringr::str_extract(z_names, "\\d+")), beta = betas)

i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]

corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()

corr_mat_dt$corrwith = rownames(corr_mat_dt)

corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))

corr_mat_dt %<>% mutate(value = ifelse(is.na(value), 0, value)) %>% mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))

corr_mat_dt = corr_mat_dt %>% mutate(corrwith =as.numeric(stringr::str_extract(corrwith, "\\d+")), variable =  as.numeric(stringr::str_extract(variable, "\\d+")), corr = ifelse(corrwith > variable, value, NA))

betaplot3 = ggplot(betadata, aes(x = z_names, y = beta, fill = beta)) + geom_bar(stat = 'identity') + scale_fill_viridis(option = 'magma', na.value = 'white') + scale_x_continuous(limits = c(0,101), expand = c(0,0)) + 
labs(title = "Magnitude of first-stage coefficients", fill = 'first stage coefficient on instrument x', subtitle = 'strongest instrument is z1, second complex dgp') + xlab('instrument number') + ylab('coefficient of z on x') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none')
#-----------------------------------------
#Complex case, shifted belloni (strongest instrument = z50, z51, z52 ..., weakest instrument = z49, z48, ...)
n = 1000
k = 100
shuffle = '3'
 #n = number of observations to generate
 #k = number of instruments to generate
 #shuffle: takes values of '1', '2' ,'3'

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

i_dt_corr = product

z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]

betadata = data.frame(z_names = as.numeric(stringr::str_extract(z_names, "\\d+")), beta = betas)

i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]

corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()

corr_mat_dt$corrwith = rownames(corr_mat_dt)

corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))

corr_mat_dt %<>% mutate(value = ifelse(is.na(value), 0, value)) %>% mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))

corr_mat_dt = corr_mat_dt %>% mutate(corrwith =as.numeric(stringr::str_extract(corrwith, "\\d+")), variable =  as.numeric(stringr::str_extract(variable, "\\d+")), corr = ifelse(corrwith > variable, value, NA))

 betaplot4 = ggplot(betadata, aes(x = z_names, y = beta, fill = beta)) + geom_bar(stat = 'identity') + scale_fill_viridis(option = 'magma', na.value = 'white') + scale_x_continuous(limits = c(0,101), expand = c(0,0)) + 
labs(title = "Magnitude of first-stage coefficients", fill = 'first stage coefficient on instrument x', subtitle = 'strongest instrument is z50, third complex dgp') + xlab('instrument number') + ylab('coefficient of z on x') +
theme_minimal(
        base_size = 7.5
      ) +
theme(legend.position = 'none')

bcase1 = plot_grid(corrplot2, corrline, ncol = 1)
bcase2 = plot_grid(betaplot2, betaplot3, betaplot4, nrow = 3, ncol = 1)

bcase = plot_grid(bcase1, bcase2, ncol = 2, rel_widths = c(.65,1))
ggsave(here('Figures', 'dgpcomplex.png') ,plot = bcase, scale = 1.2, device = 'png')
ggsave(here('Figures', 'dgpsimple.png'), plot = simplecase, scale = 1, device = 'png')

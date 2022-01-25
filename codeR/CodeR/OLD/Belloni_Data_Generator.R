#let's load up our toolkit  
options(stringsAsFactors = F)
# Packages
library(pacman)
p_unload(plotly)  
p_load(
  ggplot2, tidymodels, e1071, ggthemes, latex2exp, Cairo,
  party, estimatr, lfe,
  tidyverse, haven, data.table, lubridate, magrittr, parallel, parsnip, glmnet, doParallel, MASS
)


BelloniData <- function(samplesize = c(100,250), muvec = c(30,180), typeof = "exponential", S = c(5, 50), k = 100) {
  
  # Function: 'BelloniData': Produces datasets with all possible combinations of samplesize, concentration parameters (mu) and in both 'styles' of beta.
  # Inputs: samplesize: vector of numeric dataset sample sizes
  
  #         muvec: vector of numeric density conditions: this will allow for multiple behaviors
  #               from an instrumental variable 
  
  #         typeof: string input, can be 'exponential' or 'cutoff'. See SPARSE MODELS AND METHODS 
  #                  FOR OPTIMAL INSTRUMENTS WITH AN APPLICATION TO EMINENT DOMAIN for more info
  
  #         s: vector of integer cutoff points for cutoff design (ie, 5 corresponds to 5 instruments with 
  #            coef = 1, 95 instruments where coef = 0)
  
  #         k: number of instruments. MUST BE LARGER THAN ALL VALUES IN VECTOR S.
  
  # Outputs: 1 csv file for each combination of samplesize, muevec, s (if cutoff instrument). Name of file
  #          follows 'Belloni_...' with mu, samplesize and s denoted in file name.

  #set simulation seed
  set.seed(42)
  true_beta = 1
  
  #function to find constant to produce empirical concentration, in order to produce
  #behavior.
  
  for (n in samplesize) {
    set_C <- function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
      C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
      return(C)
    }
    
    #using Belloni method, generate pi and beta values (according to page 26 of Belloni: SPARSE MODELS AND METHODS FOR OPTIMAL INSTRUMENTS WITH AN APPLICATION TO EMINENT DOMAIN
    if (typeof == "exponential") {
      for (musq in muvec) {
        print(paste0(musq, "concentrate"))
        
        
        #Generate two errors according to Bellloni (pg 26)
        
        #first need covariance and variamce set, to produce covar matrix
        var_er = 1
        cor_er = .6
        
        var_inst = 1
        
        #produce Sigma Z. variance is 1, correlation (thus cov) is equal to .5^(|j-h|) for z_j and z_h
        Sig_z = matrix(0,k,k)
        for (j in (1:k)) {
          for (h in (1:k)) {
            Sig_z[c(j),c(h)] = .5^abs(j-h)
          }
        }
        z_n1 = mvrnorm(n = n, mu = rep(0,k), Sigma = Sig_z, empirical = TRUE)
        
        #using exponential
        beta_pattern <-  unlist(lapply(c(0:(k-1)), function(x) (.7)^x), use.names = FALSE)
        #remove one here
        
        
        Sig_z <- as.matrix(Sig_z)
        beta_pattern <- as.vector(beta_pattern)
        
        print("first lin alg step")
        C = set_C(musq=musq, n = n, beta_pattern = beta_pattern,Sig_z = Sig_z)
        C = c(C)
        
        #now that C is known, set variance for error term on X
        
        print("lin alg step 2")
        sigma_v = abs(1 - C*crossprod(beta_pattern,Sig_z) %*% beta_pattern)
        
        #We can also set our full pi matrix to find our true instrument coefficients
        betas <- as.vector(C*beta_pattern)
        cov_er<-cor_er*sqrt(sigma_v)
        
        #use multivariate distribution given matrix above to get a matrix of values for sample
        #size of 'n1'
        covarmat <- matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)
        
        #we now need to get and separate our error terms
        evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
        
        #separate errors
        en1 = evmat_n1[,1]
        vn1 = evmat_n1[,2]
        
        #let's construct x's and y's and add them to our z's to create a full data matrix
        X = z_n1 %*% betas + vn1
        Y = X + en1
        
        #full matrix generation
        
        #instrument matrix names
        znames = paste(rep("z", k), c(1:k), sep = "_")
        z = as.data.frame(z_n1)
        colnames(z) <- znames
        
        #X (or d, as in Belloni paper) name
        X = as.data.frame(X, names = "X")
        colnames(X) <- "X"
        
        #Target value Y name
        Y = as.data.frame(Y, names = "Y")
        colnames(Y) <- "Y"
        
        #final matrix
        product<-cbind(z,X,Y)
        filename = paste("BelloniData","n",n,"typeof",typeof, 
                         "mu", musq, sep = "_")
        print(filename)
        write_csv(product, paste(".\\" , filename))
        print(paste0("Finsihed for sample size of ", n))
      }
    return(product)}
    #end of loop through vector of concentration parameters
    
    #end of 'exponential condition'
    
    else if (typeof == "cutoff") {
      for (s in S) {
        for (musq in muvec) {
          #Generate two errors according to Bellloni (pg 26)
          
          #first need covariance and variamce set, to produce covar matrix
          var_er = 1
          cor_er = .6
          
          var_inst = 1
          
          #produce Sigma Z. variance is 1, correlation (thus cov) is equal to .5^(|j-h|) for z_j and z_h
          Sig_z = matrix(0,k,k)
          for (j in (1:k)) {
            for (h in (1:k)) {
              Sig_z[c(j),c(h)] = .5^abs(j-h)
            }
          }
          
          z_n1 = mvrnorm(n = n, mu = rep(0,k), Sigma = Sig_z, empirical = TRUE)
          
          beta_pattern = c(rep(1,s), rep(0,k-s))
          
          C = set_C(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n)
          
          #now that C is known, set variance for error term on X
          sigma_v = abs(1 - C*(crossprod(beta_pattern,Sig_z) %*% beta_pattern))
          
          #We can also set our full pi matrix to find our true instrument coefficients
          betas <- C*beta_pattern
          cov_er<-cor_er*sqrt(sigma_v)
          
          #use multivariate distribution given matrix above to get a matrix of values for sample
          #size of 'n' pulled from vector of samples
          covarmat <- matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)
          
          #we now need to get and separate our error terms
          evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
          
          #separate errors
          en1 = evmat_n1[,1]
          vn1 = evmat_n1[,2]
          
          #let's construct x's and y's and add them to our z's to create a full data matrix
          X = z_n1 %*% betas + vn1
          Y = X + en1
          
          #full matrix generation
          
          #instrument matrix names
          znames = paste(rep("z", k), c(1:k), sep = "_")
          z = as.data.frame(z_n1)
          colnames(z) <- znames
          
          #X (or d, as in Belloni paper) name
          X = as.data.frame(X, names = "X")
          colnames(X) <- "X"
          
          #Target value Y name
          Y = as.data.frame(Y, names = "Y")
          colnames(Y) <- "Y"
          
          #final matrix
          product<-cbind(z,X,Y)
          filename = paste("BelloniData","n",n,"typeof",typeof, 
                           "mu", musq, "rlvvar", s, sep = "_")
          print(filename)
          write_csv(product, paste(".\\" , filename))
          print(paste0("Finsihed for sample size of ", n))
        }
        
        
        
        #end loop through different samplesizes
        
      }
    }
    
  }
return(product)}

BelloniData()

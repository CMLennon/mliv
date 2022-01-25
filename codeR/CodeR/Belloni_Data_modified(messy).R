#let's load up our toolkit  
options(stringsAsFactors = F)
# Packages
library(pacman)
p_unload(plotly)  
p_load(
  ggplot2, tidymodels, e1071, ggthemes, latex2exp, Cairo,
  party, estimatr, lfe,
  tidyverse, haven, data.table, lubridate, magrittr, parallel, parsnip, glmnet, doParallel, MASS, pracma
)


BelloniData <- function(samplesize = c(10000), muvec = c(180), typeof = "exponential", S = c(5, 50), k = 100, twst = "Nne", cov_mat_ret = FALSE) {
  
  # Function: 'BelloniData': Produces datasets with all possible combinations of samplesize, concentration parameters (mu) and in both 'styles' of beta.
  # Inputs: samplesize: vector of numeric dataset sample sizes
  
  #         muvec: vector of numeric density conditions: this will allow for multiple behaviors
  #               from an instrumental variable 
  
  #         typeof: string input, can be 'exponential' or 'cutoff'. See SPARSE MODELS AND METHODS 
  #                  FOR OPTIMAL INSTRUMENTS WITH AN APPLICATION TO EMINENT DOMAIN for more info
  
  #         S: vector of integer cutoff points for cutoff design (ie, 5 corresponds to 5 instruments with 
  #            coef = 1, 95 instruments where coef = 0)
  
  #         k: number of instruments. MUST BE LARGER THAN ALL VALUES IN VECTOR S.
  
  #         twst: takes Nne, quad, and noncont as arguments and makes f(z_i) for i in [1,3] either linear, quadratic and interacted, or a function of indicator variables
  
  #         cov_mat_ret: Returns the covariance matrix for generating the error structures. Necessary to produce true counterfactuals.
  
  # Outputs: 1 csv file for each combination of samplesize, muevec, s (if cutoff instrument). Name of file
  #          follows 'Belloni_...' with mu, samplesize and s denoted in file name.
  
  #set simulation seed
  set.seed(42)
  true_beta = 1
  
  #function to find constant to produce empirical concentration, in order to produce
  #behavior.
  
  for (n in samplesize) {
    set_C <- function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
      #solving for a C value that changes behavior of lasso selectors.
      C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
      return(C)
    }
    
    #using Belloni method, generate pi and beta values (according to page 26 of Belloni: SPARSE MODELS AND METHODS FOR OPTIMAL INSTRUMENTS WITH AN APPLICATION TO EMINENT DOMAIN
    if (typeof == "exponential") {
      for (musq in muvec) {
        print(paste0(musq, "concentrate"))
        
        for (n in samplesize) {
          set_C <- function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
            #solving for a C value that changes behavior of lasso selectors.
            C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
            return(C)
          }
        #Generate two errors according to Bellloni (pg 26)
        
        #first need covariance and variance set, to produce covar matrix
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
        errors = as.data.frame(evmat_n1)
        names(errors) <- c("e", "v")
        
        #separate errors
        en1 = evmat_n1[,1]
        vn1 = evmat_n1[,2]
        
        #let's construct x's and y's and add them to our z's to create a full data matrix
        
        if (twst == "Nne") {
          #this is really just a placeholder function. The dataframe could be returned here
          X = z_n1 %*% betas + vn1
          X_hat = X
          Y = X + en1
        }
        else if (twst == "quad") {
          #step 1, modify the betas vector to zero out for z1-z3
          betas = c(0,0,0, betas[ c(3:(length(betas)-1))])
          for (i in c(1:dim(z_n1)[1])) {
            print(dim(z_n1[1]))
            entry =  z_n1[i,] %*% betas + sum((z_n1[i,c(1:3)]) %*% t(z_n1[i,c(1:3)]))/3 + 
              sum(z_n1[i,1]*(z_n1[i,c(1:3)]) %*% t(z_n1[i,c(1:3)]))/3 + vn1[i]
            X[i] = entry
          }
          X_hat = X
          Y = X + en1
        }
        else if (twst == "noncont") {
          #Nullify X
          
          #create a logic matrix to zero out some values of z1-z3 based on some
          logic = cbind(matrix(as.integer(z_n1[,c(1:3)] >= 1),ncol = 3, nrow = n), ones(n=(dim(z_n1)[1]),m=(dim(z_n1)[2]-3)))
          
          z_n1_used = logic*z_n1
          X = z_n1 %*% betas + vn1
          X_hat = X
          Y = X + en1
        }
        
        
        #full matrix generation
        
        #instrument matrix names
        znames = paste(rep("z", k), c(1:k), sep = "_")
        z = as.data.frame(z_n1)
        colnames(z) <- znames
        
        #X (or d, as in Belloni paper) name
        X = as.data.frame(X, names = "x1")
        colnames(X) <- "x1"
        
        #Target value Y name
        Y = as.data.frame(Y, names = "Y")
        colnames(Y) <- "y"
        
        #final matrix
        true_x = x - vn1
        product<-cbind(z,X,Y,errors, true_x)
        twst = twst
        filename = paste("BelloniData","n",n,"typeof",typeof, 
                         "mu", musq, twst, "final",  sep = "_")
        print(filename)
        write_csv(product, paste0('/Users/connor/Desktop/',filename, '.csv'))
        print(paste0("Finsihed for sample size of ", n))
      }
    }
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
          errors = as.data.frame(evmat_n1)
          names(errors) <- c("e", "v")
          #separate errors
          en1 = evmat_n1[,1]
          vn1 = evmat_n1[,2]
          
          #let's construct x's and y's and add them to our z's to create a full data matrix
          if (twst == "Nne") {
            #Using the belloni betas and instruments, create a linear z->x relationship
            X = z_n1 %*% betas + vn1
            Y = X + en1
          }
          else if (twst == "quad") {
            #Keep the belloni data as is, but create new quadratic variables 
            #for the first three instruments. These include two-way interactions, three-way
            #interactions, squared terms, and cubed terms for z1. 
            
            #step 1, modify the betas vector to zero out for z1-z3
            X = NULL
            betas = betas
            for (i in c(1:dim(z_n1)[1])) {
              print(dim(z_n1[1]))
              entry =  z_n1[i,] %*% betas + sum((z_n1[i,c(1:3)]) %*% t(z_n1[i,c(1:3)]))+ 
                sum(z_n1[i,1]*(z_n1[i,c(1:3)]) %*% t(z_n1[i,c(1:3)])) + vn1[i]
              X[i] = entry
            }
            Y = X + en1
          }
          else if (twst == "noncont") {
            #noncont: create, for z1 through z3 a noncontinuous relationship between f(z) and X. In practice,
            #this is generated by an arbitrary cutoff point of z_i > 1 which is sufficiently placed to impact
            #the generated data. This is added on top of the normal beta effect.
            
            #create a logic matrix to zero out some values of z1-z3 based on some cutoff point
            logic = cbind(matrix(as.integer(z_n1[,c(1:3)] >= 1),ncol = 3, nrow = n), ones(n=(dim(z_n1)[1]),m=(dim(z_n1)[2]-3)))
            
            z_n1_used = logic*z_n1
            X = z_n1_used %*% betas + z_n1 %*% betas + vn1
            Y = X + en1
          }
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
          
          #betas so Connor doesn't have to reverse-engineer from the mu values in python
          bnames = paste(rep("beta_", k), c(1:k), sep = "_")
          bts = as.data.frame(ones(n,k))
          bts = transpose(transpose(bts) * betas)
          colnames(bts) <- bnames

          product<-cbind(z,X,Y,bts,errors)
          if(cov_mat_ret == TRUE){
            filename = paste("belcovmat","n",n,"typeof",typeof, 
                             "mu", musq, "rlvvar", s,twst, sep = "_") 
            print(filename)
            write_csv(as.data.frame(covarmat), paste0('/Users/connor/Desktop/',filename, '.csv'))
            print(paste0("finished covariance matrix for belloni covar ", filename))
          }
          else{
          filename = paste("BelloniData","n",n,"typeof",typeof, 
                           "mu", musq, "rlvvar", s,twst, "final", sep = "_")
          print(filename)
          write_csv(product, paste0('/Users/connor/Desktop/',filename, '.csv'))
          print(paste0("Finsihed for sample size of ", n))
          }
        }
        
        
        
        #end loop through different samplesizes
        
      }
    }
    
  }
  #return(betas)
  return(product)
}

BelloniData(samplesize = c(150000), mu = c(80), S =c(3,10,90), typeof = "exponential")



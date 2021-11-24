#' Personalized angle regression
#'
#' This function implements the algorithm to calculate the PAN estimate for a given tuning parameter path
#' specified by maximum value and step size. 
#' @param X.train An n x p matrix for the training set. The number of 
#' observations (rows) must be larger than the number of variables (columns).
#' @param X.test An n x p matrix for the test set. The number of 
#' observations (rows) can be one or more.
#' @param lambdaMax A scalar. The absolute value of the maximum tuning parameter. 
#' @param stepSize A scalar. The step size (either positive or negative) for the tuning parameter starting from 0.   
#' @export
#' @examples
#' cat_function()

PAN <- function(X.train,Y.train,X.test,lambdaMax=10,stepSize=1){
 
  if(!is.matrix(X.train)){
    stop('Training data is not a matrix')
  }
  
  if(!is.matrix(X.test)){
    stop('Test data is not a matrix')
  }
  
  if(dim(X.train)[2]>dim(X.train)[1]){
    stop('p > n: Number of rows in training data must be smaller than the number
          of columns')
  }
  
  norm <- function(x) sqrt(sum(x^2))
  
  # Find dimensions
  n <- dim(X.train)[1]
  p <- dim(X.train)[2]
  
  n.test <- dim(X.test)[1]
  test.tunings <- seq(0,lambdaMax,by=stepSize)
  
  #Find PAN train and test error 
  gammas_alg_test <- matrix(NA,p,length( test.tunings))
  rs_alg_test <- matrix(NA,1,length( test.tunings))
  betas_alg_test <- matrix(NA,p,length( test.tunings))
  
  PAN.pred.test <- matrix(NA,n.test,length(test.tunings))
  
  beta.ols <- solve(t(X.train)%*%X.train)%*%t(X.train)%*%Y.train
  gamma0 <- beta.ols/norm(beta.ols)
  
  for(k in 1:n.test){
    x0 <- c(X.test[k,])
    A <- x0 %*% t(x0)/norm(x0)^2
     
    XtX <- t(X.train) %*% X.train 
    B <- XtX %*% beta.ols %*% t(beta.ols) %*% XtX
    
    gamma_it <- gamma0
    for(i in 1:length(test.tunings)){
      matrix_it <- const(gamma_it,XtX)*B -  sign(lambdaMax)*test.tunings[i]*A - const(gamma_it,XtX)^2*c(t(gamma_it)%*%B%*%gamma_it)*XtX
      gammas_alg_test[,i] <- eigen(matrix_it,symmetric=TRUE)$vectors[,1]
      if(c(gammas_alg_test[,i]%*% XtX %*% beta.ols) < 0) gammas_alg_test[,i] <- -gammas_alg_test[,i]
      rs_alg_test[,i] <- c(gammas_alg_test[,i]%*% XtX %*% beta.ols)/c(gammas_alg_test[,i]%*% XtX %*% gammas_alg_test[,i])
      betas_alg_test[,i] <- rs_alg_test[,i]*gammas_alg_test[,i]
      
      #Update gamma
      gamma_it <- gammas_alg_test[,i] 
    }
    
    PAN.pred.test[k,] <- c(x0%*%betas_alg_test)
  }
  return(list(preditions = PAN.pred.test, lambdas = test.tunings))
}

#' Help function
#'
#' This function implements the algorithm to calculate the PAN estimate for a given tuning parameter path
#' specified by maximum value and step size. 
#' @param gamma A vector.   
#' @export

const <- function(gamma,XtX){
  return(1/c(t(gamma) %*% XtX %*% gamma))
}


#' Bootstrap
#'
#' This function implements the algorithm to calculate the PAN estimate for a given tuning parameter path
#' specified by maximum value and step size. 
#' @param X.train An n x p matrix for the training set. The number of 
#' observations (rows) must be larger than the number of variables (columns).
#' @param X.test An n x p matrix for the test set. The number of 
#' observations (rows) can be one or more.
#' @param lambdaMax A scalar. The absolute value of the maximum tuning parameter. 
#' @param stepSize A scalar. The step size (either positive or negative) for the tuning parameter starting from 0.   
#' @export
#' @examples
#' cat_function()

PAN.bootstrap <- function(X.train,Y.train,
                          lambdaMax=10,lambdaMaxNeg=NULL,stepSize=1,stepSizeNeg=NULL,
                          bootstrap.samples=2000){
  
  norm <- function(x) sqrt(sum(x^2))
  
  if(!is.matrix(X.train)){
    stop('Training data is not a matrix')
  }
  
  if(dim(X.train)[2]>dim(X.train)[1]){
    stop('p > n: Number of rows in training data must be smaller than the number
          of columns')
  }
  
  tunings <- seq(0,lambdaMax,by = stepSize)
  tunings_neg <- seq(0,-lambdaMax,by = -stepSize)[-1]
  if(!is.null(lambdaMaxNeg)){
    tunings_neg <- seq(0,-lambdaMaxNeg,by = -stepSizeNeg)
  }
  
  # Find dimensions
  n <- dim(X.train)[1]
  p <- dim(X.train)[2]
  
  ### Create outcomes by parametric bootstrap
  Y.boot <- matrix(,n,bootstrap.samples)
  beta.boot <- matrix(,p,bootstrap.samples)
  H <- solve(t(X.train)%*%X.train)%*%t(X.train)
  beta.ols <-  H%*%Y.train
  sigma_hat <- sqrt(1/(n-p)*t(Y.train-X.train%*%beta.ols)%*%(Y.train-X.train%*%beta.ols))
  
  for(r in 1:bootstrap.samples){
    Y.boot[,r] <- X.train%*%beta.ols + rnorm(n,0,sigma_hat)
    beta.boot[,r] <- H%*%Y.boot[,r]
  } 
  
  error.array <- array(,c(n,bootstrap.samples,length(tunings)+length(tunings_neg))) 
  bias.array <- array(,c(n,bootstrap.samples,length(tunings)+length(tunings_neg))) 
  
  for(r in 1:bootstrap.samples){
    if(r %% 50 == 0) print(paste0('Bootstrap samples: ',r))
    # Initialize error matrix over training data for each path
    gamma0 <- beta.boot[,r]/norm(beta.boot[,r])
    
    gammas_alg <- matrix(,p,length(tunings))
    rs_alg <- matrix(,1,length(tunings))
    betas_alg <- matrix(,p,length(tunings))
    
    gammas_alg_neg <- matrix(,p,length(tunings_neg))
    rs_alg_neg <- matrix(,1,length(tunings_neg))
    betas_alg_neg <- matrix(,p,length(tunings_neg))
    
    for(k in 1:n){
      
      # Define inputs
      x0 <- c(X.train[k,])
      A <- (x0 %*% t(x0))/norm(x0)^2
      XtX <- t(X.train) %*% X.train
      B <- XtX %*% beta.boot[,r] %*% t(beta.boot[,r]) %*% XtX
      
      #Initialize
      gamma_it <- gamma0
      
      for(i in 1:length(tunings_neg)){
        matrix_it <- const(gamma_it,XtX)*B - tunings_neg[i]*A - const(gamma_it,XtX)^2*c(t(gamma_it)%*%B%*%gamma_it)*XtX
        gammas_alg_neg[,i] <- eigen(matrix_it,symmetric=TRUE)$vectors[,1]
        if(c(gammas_alg_neg[,i]%*% XtX %*%  beta.boot[,r] ) < 0) gammas_alg_neg[,i] <- -gammas_alg_neg[,i]
        rs_alg_neg[,i] <- c(gammas_alg_neg[,i]%*% XtX %*% beta.boot[,r])/c(gammas_alg_neg[,i]%*% XtX %*% gammas_alg_neg[,i])
        betas_alg_neg[,i] <- rs_alg_neg[,i]*gammas_alg_neg[,i]
        
        #Update gamma
        gamma_it <- gammas_alg_neg[,i] 
      }
      
      #Initialize
      gamma_it <- gamma0
      
      for(i in 1:length(tunings)){
        matrix_it <- const(gamma_it,XtX)*B - tunings[i]*A - const(gamma_it,XtX)^2*c(t(gamma_it)%*%B%*%gamma_it)*XtX 
        gammas_alg[,i] <- eigen(matrix_it,symmetric=TRUE)$vectors[,1]
        if(c(gammas_alg[,i]%*% XtX %*%  beta.boot[,r] ) < 0) gammas_alg[,i] <- -gammas_alg[,i]
        rs_alg[,i] <- c(gammas_alg[,i]%*% XtX %*% beta.boot[,r] )/c(gammas_alg[,i]%*% XtX %*% gammas_alg[,i])
        betas_alg[,i] <- rs_alg[,i]*gammas_alg[,i]
        
        #Update gamma
        gamma_it <- gammas_alg[,i] 
      }
      
      #Record bootstrap error 
      error.array[k,r,] <- c(x0 %*% beta.ols) - c(rev(x0%*%betas_alg_neg),x0%*%betas_alg)
      bias.array[k,r,] <- c(x0 %*% beta.boot[,r]) - c(rev(x0%*%betas_alg_neg),x0%*%betas_alg)
      
    }
  }
  MSE_adj <- colMeans(apply(apply(error.array,c(1,3),function(x) mean(x)^2) - apply(bias.array,c(1,3),var),1:2,function(x) max(0,x)) + apply(error.array,c(1,3), function(x) var(x))*(bootstrap.samples-1)/(bootstrap.samples))
  lambdas <- c(rev(tunings_neg),tunings)
  return(list(MSE = MSE_adj, lambdas = lambdas, lambda.min = lambdas[which.min(MSE_adj)]))
}




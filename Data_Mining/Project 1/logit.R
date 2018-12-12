numer <- function(x, betas, k, K){
  return(exp(as.numeric(x) %*% betas))  
}

denom <- function(x, betas, K){
  esm <- 0
  for(k in 1:(K-1)){
    esm <- esm + exp(x %*% betas)
  }
  return(1+esm)
}

pi.fun <- function(x, betas, k, K, denom.v){
  if(k == K){
    return(1/denom.v)
  }
  else{
    return(numer(x, betas, k, K)/denom.v)
  }
}

class.ind <- function(cl) {
  n <- length(cl)   # number of observations
  x <- matrix(0, n, length(levels(cl))) # matrix of zeros nobs x nclasses
  x[(1L:n) + n * (as.integer(cl) - 1L)] <- 1 # Create indicator matrix
  dimnames(x) <- list(names(cl), levels(cl)) # add column names
  x # return indicator matrix
}

stopcrit <- function(beta.old, beta.new, iter, min.iter = 10, max.iter = 1000, ep = 0.01){
  beta.diff.norm <- sqrt(sum((beta.new - beta.old)^2))
  if((beta.diff.norm< ep | iter >= max.iter) & iter > min.iter){
    invisible(TRUE)
  }
  else{
    invisible(FALSE)
  }
}

logit.m <- function(X, y){
  N <- length(y)
  beta.old <- rep(0, dim(X)[2])
  beta.new <- rep(0, dim(X)[2])
  K <- 2
  iter <- 0
  while(stopcrit(beta.old, beta.new, iter)==FALSE){
    beta.old <- beta.new
    P <- matrix(NA, nrow =N, ncol = (K-1))
    # Denominators
    d.v <- apply(denom, X = X, MARGIN = 1, beta.old, K) 
    # Create Probability matrix
    for(i in 1:N){
      for(j in 1:(K-1)){
        P[i,j] <-pi.fun(X[i,],  betas = beta.old, k = j, K = K, denom.v = d.v[i])
      }
    }
    X.t <- t(t(X) %*% diag(as.vector(P)*as.vector(1-P), nrow = N))
    beta.new <- beta.old + solve(t(X)%*%X.t)%*%t(X)%*%(y - P)
    iter = iter +1
  }
  return(beta.new)
}


predict.logit <- function(betas, X){
  exp(X%*%betas)/(1 + exp(X%*%betas))
}

predict.logit.cl <- function(betas, X, thsh = 0.5){
  (predict.logit(betas, X) > thsh)
}


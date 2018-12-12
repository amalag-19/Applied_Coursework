sim.logistic <- function(n, n.cov, beta, Sigma = diag(.1,nrow = n.cov) ){
  # Make X multivariate normal
  X <- mvrnorm(n = n, mu = rep(0, n.cov), Sigma = Sigma )
  X <- cbind(rep(1, n), X)
  # Simulate response given X
  lpred <- X%*%beta
  nums <- exp(lpred)
  dens <- 1 + apply(X = nums, FUN = sum, MARGIN = 1)
  p <- cbind(nums/dens, 1- nums/dens)
  y <- apply(X = p, FUN = sample, MARGIN = 1, x = c(1,0), size = 1, replace = FALSE )
  data <- data.frame(cbind(y, X[,-1]))
  v.names <- rep(NA, n.cov)
  for(i in 1:n.cov){
    v.names[i] <- paste("X",i, sep  = "")
  }
  names(data)<- c("y", v.names)
  invisible(data)
}

## 
## Example
##
# 
# betas <- runif(n = (4 + 1))   # True betas (+ intercept beta)
# datas <- sim.logistic(n = 10, n.cov = 4, betas)

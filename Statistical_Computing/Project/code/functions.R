## loading required packages
{library(tmvtnorm)
library(ggplot2)
library(plyr)
library(mvtnorm)
library(reshape)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")}

#####################################################################################
## defining various distance measures
## defining the function to calculate euclidean distance between two vectors
euclid.dist<- function(x,y,dist.par=0){
  d<-norm(as.matrix(x-y),type = "f")
  return(d)
}


## defining the function to calculate mahalanobis distance between two vectors
maha.dist<-function(x,y,dist.par=0){
  s<-cov(x,y)
  x<-matrix(c(x),nrow = length(x),ncol = 1)
  y<-matrix(c(y),nrow = length(y),ncol = 1)
  value<-t(x-y)*(1/s)*(x-y)
  return(value)
}

## defining the function to calculate minkowski distance of degree p between two 
## vectors
mink.dist<-function(x,y,dist.par){
  p<-dist.par
  x<-matrix(c(x),nrow = length(x),ncol = 1)
  y<-matrix(c(y),nrow = length(y),ncol = 1)
  value<-(sum(abs(x-y)^p))^(1/p)
  return(value)
}

#####################################################################################
## defining the function for sufficient statistics for parameters in mixture
## normal distribution
suff.stat.MN<-function(data,suff.par=0){
  tp<-length(which(data[,1]==1))/nrow(data)
  tm1<-mean(data[which(data[,1]==1),2])
  tm2<-mean(data[which(data[,1]==2),2])
  return(c(tp,tm1,tm2))
}

#####################################################################################
## defining the uniform kernel density estimate that takes in the data frames x and y
## and returns indicator 1 or 0 depending on the distance between the
## sufficient statistics of the samples. defined on log scale
kernel.unif<-function(x,y,eps,suff.stat,suff.par,kernel.par=0,dist,dist.par){
  tx<-suff.stat(x,suff.par)
  ty<-suff.stat(y,suff.par)
  if(prod(!is.nan(tx))==1){
    if (dist(tx,ty,dist.par)<=eps){
      value<-0
    }
    else{value<-NA}
  }
  else{value<-NA}
  return(value)
}

## defining the gaussian kernel density estimate that takes in the data frames x and y
## and returns gaussian density estimate depending on the distance between the
## sufficient statistics of the samples. defined on log scale
kernel.gaussian<-function(x,y,eps,suff.stat,suff.par,kernel.par,dist,dist.par){
  sig<-kernel.par
  tx<-suff.stat(x,suff.par)
  ty<-suff.stat(y,suff.par)
  if(prod(!is.nan(tx))==1){
    if (dist(tx,ty,dist.par)<=eps){
      value<-dmvnorm(x = ty,mean = tx,sigma = (eps^2)*sig,log = TRUE)
    }
    else{value<-NA}
  }
  else{value<-NA}
  return(value)
}

## defining the epanechnikov kernel density estimate that takes in the data frames x
## and y and returns epanechnikov density estimate depending on the distance
## between the sufficient statistics of the samples. defined on log scale.
kernel.epan<-function(x,y,eps,suff.stat,suff.par,kernel.par,dist,dist.par){
  c<-kernel.par
  N<-sqrt((3/c)*(1-(1/(2*c^2))))
  sig<-kernel.par
  tx<-suff.stat(x,suff.par)
  ty<-suff.stat(y,suff.par)
  if(prod(!is.nan(tx))==1){
    if (dist(tx,ty,dist.par)<=min(eps,sqrt(c))){
      value<-log(N*(c-euclid.dist(tx,ty)^2))
    }
    else{value<-NA}
  }
  else{value<-NA}
  return(value)
}

#####################################################################################
## defining the function to simulate data from mixture normal density
data.sim<-function(m,p,m1,m2){
  comp<-sample(x = c(1,2),size = m,prob = c(p,(1-p)),replace = TRUE)
  id1<-which(comp==1)
  id2<-which(comp==2)
  samples<-rep(NA_real_,m)
  samples[id1]<-rnorm(n = length(id1),mean = m1,sd = 1)
  samples[id2]<-rnorm(n = length(id2),mean = m2,sd = 1)
  df<-data.frame(cbind(comp,samples))
  return(df)
}

#####################################################################################
## Augmented likelihood free Metropolis Hastings function function with inputs as the
## current state of the MC (theta.cs,df.cs), the variance of proposal (tuning parameter)
## & data and output as the next state (theta.ns,df.ns)
LF.MH<-function(theta.cs,df.cs,tune,data,eps,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par){
  ## sampling the next theta vector (p,m1,m2) from truncated multivariate proposal
  theta.star<-rtmvnorm(n=1, mean=theta.cs, sigma=diag(tune), lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf))
  
  ## generating augmented data from the model
  df.star<-data.sim(nrow(data),theta.star[1],theta.star[2],theta.star[3])
  kern.value<-kernel.func(df.star,data,as.numeric(eps),suff.stat,suff.par,kernel.par,dist,dist.par)
  
  if (!is.na(kern.value)){
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior is uniform and this
    ## is likelihood free.
    num<-kern.value+dtmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(tune),lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf),log = TRUE)
    denom<-kernel.func(df.cs,data,as.numeric(eps),suff.stat,suff.par,kernel.par,dist,dist.par)+dtmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(tune),lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf),log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      df.ns<-df.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      df.ns<-df.cs
    }
  }
  else {
    theta.ns<-theta.cs
    df.ns<-df.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,df.ns,flag))
}

## Augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug<-function(n,start,tune,data,eps,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par){
  theta<-matrix(NA_real_,3,n)
  ## Defining the initial value for the chain
  theta[,1]<-start
  ## initializing the augmented variable
  df.temp<-data
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH(theta.cs = theta[,i],df.cs = df.temp,tune = tune,data = data,eps = eps,suff.stat = suff.stat,suff.par = suff.par,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par)
    theta[,i+1]<-temp[[1]]
    df.temp<-temp[[2]]
    accept<-accept+temp[[3]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(theta)[,1],"X2"=t(theta)[,2],"X3"=t(theta)[,3])
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,3,ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,3,ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    MCMC.means[,j]<-rowSums(theta[,1:m[j]])/m[j]
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(theta[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(theta[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(theta[3,1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

###################################################################################
## Error augmented likelihood free Metropolis Hastings function function with inputs
## as the current state of the MC (theta.cs,df.cs), the variance of proposal (tuning
## parameter) & data and output as the next state (theta.ns,df.ns)
LF.MH.eps<-function(theta.cs,df.cs,eps.cs,theta.tune,eps.tune,data,kernel.func,kernel.par,dist,dist.par){
  ## sampling the next theta vector (p,m1,m2) from truncated multivariate proposal
  theta.star<-rtmvnorm(n=1, mean=theta.cs, sigma=diag(theta.tune), lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf))
  eps.star<-rtmvnorm(n=1, mean=eps.cs, sigma=diag(eps.tune), lower=0.05, upper = 10)
  ## generating augmented data from the model
  df.star<-data.sim(nrow(data),theta.star[1],theta.star[2],theta.star[3])
  kern.value<-kernel.func(df.star,data,as.numeric(eps.star),kernel.par,dist,dist.par)
  
  if (!is.na(kern.value)){
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior for theta is uniform 
    ## and for epsilon is gamma. Also this is likelihood free.
    num<-kern.value+dgamma(x = eps.star,shape = 10,scale = 1,log = TRUE)+
      dtmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(theta.tune),lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf),log = TRUE)+
      dtmvnorm(x = eps.cs,mean = as.vector(eps.star),sigma = eps.tune,lower=0.05, upper = 10,log = TRUE)
    denom<-kernel.func(df.cs,data,as.numeric(eps.cs),kernel.par,dist,dist.par)+dgamma(x = eps.cs,shape = 10,scale = 1,log = TRUE)+
      dtmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(theta.tune),lower=c(0.01,-Inf,-Inf), upper = c(0.99,Inf,Inf),log = TRUE)+
      dtmvnorm(x = eps.star,mean = as.vector(eps.cs),sigma = eps.tune,lower=0.05, upper = 10,log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      df.ns<-df.star
      eps.ns<-eps.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      df.ns<-df.cs
      eps.ns<-eps.cs
    }
  }
  else {
    theta.ns<-theta.cs
    df.ns<-df.cs
    eps.ns<-eps.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,df.ns,eps.ns,flag))
}

## Error augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug.eps<-function(n,start,theta.tune,eps.tune,data,kernel.func,kernel.par,dist,dist.par){
  theta<-matrix(NA_real_,3,n)
  eps<-rep(NA_real_,n)
  ## Defining the initial value for the chain
  theta[,1]<-start[1:3]
  ## initializing the augmented variable
  df.temp<-data
  ## Defining the initial value for the epsilon chain
  eps[1]<-start[4]
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH.eps(theta.cs = theta[,i],df.cs = df.temp,eps.cs = eps[i], theta.tune = theta.tune,eps.tune=eps.tune,data = data,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par)
    theta[,i+1]<-temp[[1]]
    df.temp<-temp[[2]]
    eps[i+1]<-temp[[3]]
    accept<-accept+temp[[4]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(theta)[,1],"X2"=t(theta)[,2],"X3"=t(theta)[,3],"X4"=eps)
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,4,ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,4,ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    MCMC.means[,j]<-as.vector(c(rowSums(theta[,1:m[j]])/m[j],mean(eps[1:m[j]])))
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(theta[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(theta[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(theta[3,1:m[j]])$se
    MCMC.se[4,j]<-bm(eps[1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]),ess(samples[,5]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

#####################################################################################
## Marginal space likelihood free Metropolis Hastings function function with inputs as the
## current state of the MC (theta.cs,df.cs), the variance of proposal (tuning parameter)
## & data and output as the next state (theta.ns,df.ns)
LF.MH.MS<-function(theta.cs,kern.cs,tune,data,eps,S,kernel.func,kernel.par,dist,dist.par){
  ## sampling the next theta vector (p,m1,m2) from truncated multivariate proposal
  theta.star<-rtmvnorm(n=1, mean=theta.cs, sigma=diag(tune), lower=c(0.05,-Inf,-Inf), upper = c(0.95,Inf,Inf))
  
  ## generating augmented data 1:S from the model
  df.len<-rep(NA_real_,S)
  df.star<-lapply(X = df.len,FUN = function(x){
    value<-data.sim(nrow(data),theta.star[1],theta.star[2],theta.star[3])
    invisible(value)
  })
  
  kern.value<-lapply(X = df.star,FUN = function(x){
    value<-kernel.func(x,data,as.numeric(eps),kernel.par,dist,dist.par)
    invisible(value)
  })
  kern.star<-as.vector(do.call(rbind,kern.value))
  
  if (prod(!is.na(kern.star))==1){
    kern.star<-log(mean(exp(kern.star)))
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior is uniform and this
    ## is likelihood free.
    num<-kern.star+dtmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(tune),lower=c(0.05,-Inf,-Inf), upper = c(0.95,Inf,Inf),log = TRUE)
    denom<-kern.cs+dtmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(tune),lower=c(0.05,-Inf,-Inf), upper = c(0.95,Inf,Inf),log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      kern.ns<-kern.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      kern.ns<-kern.cs
    }
  }
  else {
    theta.ns<-theta.cs
    kern.ns<-kern.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,kern.ns,flag))
}

## Augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug.MS<-function(n,start,tune,data,eps,S,kernel.func,kernel.par,dist,dist.par){
  theta<-matrix(NA_real_,3,n)
  ## Defining the initial value for the chain
  theta[,1]<-start
  ## initializing the augmented variable
  kern.temp<-0
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH.MS(theta.cs = theta[,i],kern.cs = kern.temp,tune = tune,data = data,eps = eps,S = S,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par)
    theta[,i+1]<-temp[[1]]
    kern.temp<-temp[[2]]
    accept<-accept+temp[[3]]
    print(i)
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(theta)[,1],"X2"=t(theta)[,2],"X3"=t(theta)[,3])
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,3,ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,3,ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    MCMC.means[,j]<-rowSums(theta[,1:m[j]])/m[j]
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(theta[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(theta[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(theta[3,1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

##############################################################################################
## defining the function for sufficient statistics for multivariate normal using gaussian copula
suff.stat.G<-function(data,suff.par){
  ## suff.par contains the id's of the components of data to be returned
  return(data[suff.par])
}

log.prior<-function(theta,b,p){
  if (p>2){
    value<--((theta[1]^2)/200)-(((theta[2]-b*(theta[1]^2)+100*b)^2)/2)-sum(theta[3:p]^2)
  }
  else{value<--((theta[1]^2)/200)-(((theta[2]-b*(theta[1]^2)+100*b)^2)/2)}
  return(value)
}

## Gaussian copula for multivariate normal: PAIRS: Augmented likelihood free Metropolis Hastings 
## function function with inputs as the current state of the MC (theta.cs,y.cs), the variance of
## proposal (tuning parameter) & data and output as the next state (theta.ns,y.ns)
LF.MH.G<-function(theta.cs,y.cs,tune,data,eps,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  ## defining the dimensionality of the multivariate normal
  p<-length(data)
  ## sampling the next theta vector from multivariate proposal
  theta.star<-rmvnorm(n=1, mean=theta.cs, sigma=diag(tune))
  
  ## generating augmented data from the likelihood model
  s0<-1
  y.star<-rmvnorm(n = 1,mean = theta.star,sigma = diag(rep(s0,p)))
  kern.value<-kernel.func(y.star,data,as.numeric(eps),suff.stat,suff.par,kernel.par,dist,dist.par)
  
  if (!is.na(kern.value)){
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior is uniform and this
    ## is likelihood free.
    num<-kern.value+log.prior(theta = theta.star,b = b,p = p)+dmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(tune),log = TRUE)
    denom<-kernel.func(y.cs,data,as.numeric(eps),suff.stat,suff.par,kernel.par,dist,dist.par)+log.prior(theta = theta.cs,b = b,p = p)+dmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(tune),log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      y.ns<-y.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      y.ns<-y.cs
    }
  }
  else {
    theta.ns<-theta.cs
    y.ns<-y.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,y.ns,flag))
}

## Augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug.G<-function(n,start,tune,data,eps,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  if (nrow(suff.par)>1){
    if (suff.par[which(suff.par!=2)]>2){
      suff.par.m<-c(1,2,suff.par[which(suff.par!=2)])
    }
    else{suff.par.m<-suff.par}
  }
  else if(suff.par==2){
    suff.par.m<-c(1,2)
  }
  else{suff.par.m<-suff.par}
  kernel.par<-diag(1,length(suff.par.m))
  ## matrix for storing the parameter values
  theta<-matrix(NA_real_,length(data),n)
  ## Defining the initial value for the chain
  theta[,1]<-start
  ## initializing the augmented variable
  y.temp<-data
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH.G(theta.cs = theta[,i],y.cs = y.temp,tune = tune,data = data,eps = eps,suff.stat = suff.stat,suff.par = suff.par.m,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par,b=b)
    theta[,i+1]<-temp[[1]]
    y.temp<-temp[[2]]
    accept<-accept+temp[[3]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame(cbind("iterations"=1:n,t(theta)))
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,nrow(suff.par),ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,nrow(suff.par),ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    for (i in 1:nrow(suff.par)){
      MCMC.means[i,j]<-sum(theta[suff.par[i],1:m[j]])/m[j]
      ## Using the batchmeans function to calculate MCMC standard errors
      MCMC.se[i,j]<-bm(theta[i,1:m[j]])$se
    } 
  }
  ESS.10K<-ess(samples[(1:10000),2])
  ESS<-c(ess(samples[,2]))
  if (nrow(suff.par)>1){
    for (j in 2:nrow(suff.par)){
      ESS.10K<-c(ESS.10K,ess(samples[(1:10000),j]))
      ESS<-c(ESS,ess(samples[,j]))
    }
  }
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate,"time"=sys.time[3]))
}

################################################################################################
## Error augmented for Gaussian copula for multivariate normal: PAIRS: Likelihood free Metropolis Hastings 
## function function with inputs as the current state of the MC (theta.cs,y.cs,eps.cs), the variance of
## proposal (tuning parameter) & data and output as the next state (theta.ns,y.ns,eps.ns)
LF.MH.Geps<-function(theta.cs,y.cs,eps.cs,theta.tune,eps.tune,data,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  ## defining the dimensionality of the multivariate normal
  p<-length(data)
  ## sampling the next theta vector from multivariate proposal
  theta.star<-rmvnorm(n=1, mean=theta.cs, sigma=diag(theta.tune))
  ## sampling the next epsilon from truncated normal proposal
  eps.star<-rtmvnorm(n=1, mean=eps.cs, sigma=eps.tune, lower=0.05, upper = 20)
  
  ## generating augmented data from the likelihood model
  s0<-1
  y.star<-rmvnorm(n = 1,mean = theta.star,sigma = diag(rep(s0,p)))
  kern.value<-kernel.func(y.star,data,as.numeric(eps.star),suff.stat,suff.par,kernel.par,dist,dist.par)
  
  if (!is.na(kern.value)){
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior is uniform and this
    ## is likelihood free.
    num<-kern.value+dgamma(x = eps.star,shape = 10,scale = 1,log = TRUE)+
      log.prior(theta = theta.star,b = b,p = p)+dmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(theta.tune),log = TRUE)+
      dtmvnorm(x = eps.cs,mean = as.vector(eps.star),sigma = eps.tune,lower=0.05, upper = 20,log = TRUE)
    denom<-kernel.func(y.cs,data,as.numeric(eps.cs),suff.stat,suff.par,kernel.par,dist,dist.par)+dgamma(x = eps.cs,shape = 10,scale = 1,log = TRUE)+
      log.prior(theta = theta.cs,b = b,p = p)+dmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(theta.tune),log = TRUE)+
      dtmvnorm(x = eps.star,mean = as.vector(eps.cs),sigma = eps.tune,lower=0.05, upper = 10,log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      y.ns<-y.star
      eps.ns<-eps.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      y.ns<-y.cs
      eps.ns<-eps.cs
    }
  }
  else {
    theta.ns<-theta.cs
    y.ns<-y.cs
    eps.ns<-eps.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,y.ns,eps.ns,flag))
}

## Augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug.Geps<-function(n,start,theta.tune,eps.tune,data,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  if (nrow(suff.par)>1){
    if (suff.par[which(suff.par!=2)]>2){
      suff.par.m<-c(1,2,suff.par[which(suff.par!=2)])
    }
    else{suff.par.m<-suff.par}
  }
  else if(suff.par==2){
    suff.par.m<-c(1,2)
  }
  else{suff.par.m<-suff.par}
  ## matrix for storing the parameter values
  theta<-matrix(NA_real_,length(data),n)
  ## vector for epsilon
  eps<-rep(NA_real_,n)
  ## Defining the initial value of theta for the chain
  theta[,1]<-start[[1]]
  ## Defining the initial value for the epsilon chain
  eps[1]<-start[[2]]
  ## initializing the augmented variable
  y.temp<-data
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH.Geps(theta.cs = theta[,i],y.cs = y.temp,eps.cs = eps[i],theta.tune = theta.tune,eps.tune = eps.tune,data = data,suff.stat = suff.stat,suff.par = suff.par.m,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par,b=b)
    theta[,i+1]<-temp[[1]]
    y.temp<-temp[[2]]
    eps[i+1]<-temp[[3]]
    accept<-accept+temp[[4]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame(cbind("iterations"=1:n,t(theta),eps))
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,nrow(suff.par)+1,ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,nrow(suff.par)+1,ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    for (i in 1:nrow(suff.par)){
      MCMC.means[i,j]<-sum(theta[suff.par[i],1:m[j]])/m[j]
      ## Using the batchmeans function to calculate MCMC standard errors
      MCMC.se[i,j]<-bm(theta[i,1:m[j]])$se
    }
    MCMC.means[nrow(MCMC.means),j]<-sum(eps[1:m[j]])/m[j]
    MCMC.se[nrow(MCMC.means),j]<-bm(eps[1:m[j]])$se
  }
  ESS.10K<-ess(samples[(1:10000),2])
  ESS<-c(ess(samples[,2]))
  if (nrow(suff.par)>1){
    for (j in 2:nrow(suff.par)){
      ESS.10K<-c(ESS.10K,ess(samples[(1:10000),j]))
      ESS<-c(ESS,ess(samples[,j]))
    }
  }
  ESS.10K<-c(ESS.10K,ess(samples[(1:10000),ncol(samples)]))
  ESS<-c(ESS,ess(samples[,ncol(samples)]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate,"time"=sys.time[3]))
}

###############################################################################################
## Multiple augmented LF MCMC via Gaussian copula for multivariate normal: PAIRS: Augmented likelihood free Metropolis Hastings 
## function function with inputs as the current state of the MC (theta.cs,y.cs), the variance of
## proposal (tuning parameter) & data and output as the next state (theta.ns,y.ns)
LF.MH.Gmult<-function(theta.cs,kern.cs,tune,data,eps,S,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  ## defining the dimensionality of the multivariate normal
  p<-length(data)
  
  ## sampling the next theta vector from multivariate proposal
  theta.star<-rmvnorm(n=1, mean=theta.cs, sigma=diag(tune))
  
  ## generating S augmented data from the likelihood model
  s0<-1
  y.star<-rmvnorm(n = S,mean = theta.star,sigma = diag(rep(s0,p)))
  kern.star<-rep(NA_real_,S)
  for (i in 1:S){
    kern.star[i]<-kernel.func(as.vector(y.star[i,]),data,as.numeric(eps),suff.stat,suff.par,kernel.par,dist,dist.par)
  }
  
  if (prod(!is.na(kern.star))==1){
    kern.star<-log(mean(exp(kern.star)))
    ## defining the acceptance probability for sampled (theta.star,x.star) on log
    ## scale. note that here there is only proposal term. prior is uniform and this
    ## is likelihood free.
    num<-kern.star+log.prior(theta = theta.star,b = b,p = p)+dmvnorm(x = theta.cs,mean = as.vector(theta.star),sigma = diag(tune),log = TRUE)
    denom<-kern.cs+log.prior(theta = theta.cs,b = b,p = p)+dmvnorm(x = as.vector(theta.star),mean = theta.cs,sigma = diag(tune),log = TRUE)
    accept.probab<-num-denom
    ## sampling u from uniform(0,1) to check for acceptance
    u<-runif(1, min=0, max=1)
    ## initializing the indicator flag=0 to check if the sampled x.star
    ## will be accepted
    flag<-0
    ## if-else to define the next state of the chain based on acceptance probability
    if(log(u)<=accept.probab){
      theta.ns<-theta.star
      kern.ns<-kern.star
      flag<-1
    }
    else {
      theta.ns<-theta.cs
      kern.ns<-kern.cs
    }
  }
  else {
    theta.ns<-theta.cs
    kern.ns<-kern.cs
    flag<-0
  }
  ## returning the next state and indicator if the sampled value was accepted
  return(list(theta.ns,kern.ns,flag))
}

## Augmented likelihood free sampler function to generate the chains, calculate 
## the expectation i.e. Monte Carlo estimates and MCMC standard errors
LF.aug.Gmult<-function(n,start,tune,data,eps,S,suff.stat,suff.par,kernel.func,kernel.par,dist,dist.par,b){
  if (nrow(suff.par)>1){
    if (suff.par[which(suff.par!=2)]>2){
      suff.par.m<-c(1,2,suff.par[which(suff.par!=2)])
    }
    else{suff.par.m<-suff.par}
  }
  else if(suff.par==2){
    suff.par.m<-c(1,2)
  }
  else{suff.par.m<-suff.par}
  ## matrix for storing the parameter values
  theta<-matrix(NA_real_,length(data),n)
  ## Defining the initial value for the chain
  theta[,1]<-start
  ## initializing the augmented variable
  kern.temp<-0
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for LF.MH updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-LF.MH.Gmult(theta.cs = theta[,i],kern.cs = kern.temp,tune = tune,data = data,eps = eps,S=S,suff.stat = suff.stat,suff.par = suff.par.m,kernel.func = kernel.func,kernel.par = kernel.par,dist=dist,dist.par = dist.par,b=b)
    theta[,i+1]<-temp[[1]]
    kern.temp<-temp[[2]]
    accept<-accept+temp[[3]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame(cbind("iterations"=1:n,t(theta)))
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,nrow(suff.par),ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,nrow(suff.par),ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    for (i in 1:nrow(suff.par)){
      MCMC.means[i,j]<-sum(theta[suff.par[i],1:m[j]])/m[j]
      ## Using the batchmeans function to calculate MCMC standard errors
      MCMC.se[i,j]<-bm(theta[i,1:m[j]])$se
    } 
  }
  ESS.10K<-ess(samples[(1:10000),2])
  ESS<-c(ess(samples[,2]))
  if (nrow(suff.par)>1){
    for (j in 2:nrow(suff.par)){
      ESS.10K<-c(ESS.10K,ess(samples[(1:10000),j]))
      ESS<-c(ESS,ess(samples[,j]))
    }
  }
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate,"time"=sys.time[3]))
}
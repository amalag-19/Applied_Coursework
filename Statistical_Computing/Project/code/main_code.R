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
}## loading required packages
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

###########################################
## defining the parameters
p<-0.5
m1<--5
m2<-5

## simulating the data
m<-500 ## sample size
data<-data.sim(m = m,p = p,m1 = m1,m2 = m2)
plot(density(data[,2]))

###########################################
## defining the length of the chains
n<-10000

## defining the initial values of the chains
start<-c(0.4,-5.5,5.5)
start
eps<-2
var.tune.unif<-c(0.06,1,1)

## running 1 chain
undebug(LF.aug)
undebug(LF.MH)
source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
M<-LF.aug(n = n,start = start,tune=var.tune.unif,data = data,eps = eps, suff.stat = suff.stat.MN, kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3)

## checking sample variances for tuning
var(M[[1]][,2])
var(M[[1]][,3])
var(M[[1]][,4])

## MCMC estimates for 3 chains
M[[3]][,nrow(M[[3]])]

## acceptance rates
M[[5]]

#undebug(LF.MH)
# LF.MH(start,data,var.tune,data,eps)

################################################################################################
#############################################
## defining the parameters
p<-0.5
m1<--5
m2<-5

## simulating the data
m<-500 ## sample size
data<-data.sim(m = m,p = p,m1 = m1,m2 = m2)
plot(density(data[,2]))

##############################################
## running multiple chains
n<-10000
start<-matrix(c(0.6,-2,7,0.4,-3,3,0.45,-4,4),3,3)
##start<-c(0,3,4)
start
eps<-2
sig<-diag(3)
var.tune.unif<-c(0.06,1.02,1)
var.tune.gaussian<-c(0.06,1,1)
var.tune.epan<-c(0.06,0.98,0.99)
M<-vector("list",ncol(start))

## running multiple chains
cl <- makeCluster(3)
registerDoParallel(cl)
M<-vector("list",ncol(start))
M<-foreach (j = 1:ncol(start))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug(n = n,start = start[,j],tune=var.tune.gaussian,data = data,eps = eps,suff.stat.MN,suff.par=0,kernel.func = kernel.gaussian,kernel.par = sig,dist = mink.dist,dist.par = 3)
}
stopCluster(cl)

## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])

## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")

## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:ncol(start)){
  f[[j]]<-data.frame(M[[j]][[2]],t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","mean.X3","MCMCse.X1","MCMCse.X2","MCMCse.X3","start.label")

g<-do.call(rbind,g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X1 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X2 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X3 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
multiplot(p1, p2, p3,cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X1",colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X2",colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X3",colour="Label for starting values")+thema
multiplot(p4, p5, p6,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X1",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X1",title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X2",title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X2",title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X3",title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X3",title="After N")+thema

multiplot(q1,q3,q5,q2,q4,q6, cols=2)

## estimated density for different starting values
q7<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X1",colour="Label of starting values")+thema
q8<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X2",colour="Label of starting values")+thema
q9<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X3",colour="Label of starting values")+thema
multiplot(q7,q8,q9, cols=1)

###########################################
## running the chain and plotting the acceptance rates for different absolute values

## defining the length of the chains
n<-10000

## defining the initial values of the chains
start<-c(0.4,-5.5,5.5)

eps<-seq(0.1,5,by = 0.05)
var.tune.unif<-c(0.06,1.02,1)
var.tune.gaussian<-c(0.06,1,1)
var.tune.epan<-c(0.06,0.98,0.99)

## to learn the dependence of acceptance rate with scale parameter in the kernel
## density running the various chains in parallel for uniform kernel 
cl <- makeCluster(6)
registerDoParallel(cl)
M.unif<-vector("list",length(eps))
M.unif<-foreach (j = 1:length(eps)) %dopar% {
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  LF.aug(n = n,start = start,tune=var.tune.unif,data = data,eps = eps[j],kernel.func = kernel.unif,kernel.par = 0,dist = euclid.dist,dist.par = 0)
}
stopCluster(cl)

## to learn the dependence of acceptance rate with scale parameter in the kernel
## density running the various chains in parallel for Gaussian kernel 
cl <- makeCluster(6)
registerDoParallel(cl)
M.gaussian<-vector("list",length(eps))
M.gaussian<-foreach (j = 1:length(eps)) %dopar% {
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  LF.aug(n = n,start = start,tune=var.tune.gaussian,data = data,eps = eps[j],kernel.func = kernel.gaussian,kernel.par = sig)
}
stopCluster(cl)

## to learn the dependence of acceptance rate with scale parameter in the kernel
## density running the various chains in parallel for Epanechnikov kernel 
cl <- makeCluster(6)
registerDoParallel(cl)
M.epan<-vector("list",length(eps))
M.epan<-foreach (j = 1:length(eps)) %dopar% {
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  LF.aug(n = n,start = start,tune=var.tune.epan,data = data,eps = eps[j],kernel.func = kernel.epan,kernel.par = 100)
}
stopCluster(cl)

## creating dataframes for plotting
df.accept<-list()
df.ESS10K<-list()
df.ESSs<-list()
for(j in 1:length(eps)){
  df.accept[[j]]<-data.frame(cbind(M.unif[[j]][[5]],M.gaussian[[j]][[5]],M.epan[[j]][[5]]))
  df.ESS10K[[j]]<-data.frame(cbind(M.unif[[j]][[6]],M.gaussian[[j]][[6]],M.epan[[j]][[6]]))
  df.ESSs[[j]]<-data.frame(cbind(M.unif[[j]][[7]],M.gaussian[[j]][[7]],M.epan[[j]][[7]]))
}

## plotting the acceptance rates vs. epsilon for different kernels
df.temp<-do.call(rbind,df.accept)
df.accept<-data.frame(cbind(eps,as.matrix(df.temp)))
names(df.accept)<-c("epsilon","uniform","gaussian","epanechnikov")
g<-melt(data = df.accept,id=c("epsilon"))
p<-ggplot(data = g)
p+geom_line(mapping = aes(x=epsilon,y = value,colour=variable))+labs(x=expression(epsilon),y="Acceptance rate",colour="Kernel")

## plotting ESS per 10,000 vs. epsilon for different kernels
df.temp<-do.call(rbind,df.ESS10K)
df.ESS10K<-data.frame(cbind(eps,as.matrix(df.temp)))
names(df.ESS10K)<-c("epsilon","uniform","gaussian","epanechnikov")
g<-melt(data = df.ESS10K,id=c("epsilon"))
p<-ggplot(data = g)
p+geom_line(mapping = aes(x=epsilon,y = value,colour=variable))+labs(x=expression(epsilon),y="ESS per 10,000",colour="Kernel")

## plotting the ESS per second vs. epsilon for different kernels
df.temp<-do.call(rbind,df.ESSs)
df.ESSs<-data.frame(cbind(eps,as.matrix(df.temp)))
names(df.ESSs)<-c("epsilon","uniform","gaussian","epanechnikov")
g<-melt(data = df.ESSs,id=c("epsilon"))
p<-ggplot(data = g)
p+geom_line(mapping = aes(x=epsilon,y = value,colour=variable))+labs(x=expression(epsilon),y="ESS per second",colour="Kernel")

################################################################################################
## Testing the error augmented likelihood free MCMC sampler
## defining the parameters
p<-0.5
m1<--5
m2<-5

## simulating the data
m<-500 ## sample size
data<-data.sim(m = m,p = p,m1 = m1,m2 = m2)
plot(density(data[,2]))

## defining the chain size
n<-10000
start<-matrix(c(0.6,-2,7,2,0.4,-3,3,3,0.45,-4,4,8),4,3)
sig<-diag(3)
var.tune.unif<-c(0.06,1.02,1)
var.tune.gaussian<-c(0.06,16.4,15.7)
var.tune.epan<-c(0.06,0.98,0.99)
eps.tune<-1.79

# debug(LF.aug.eps)
# debug(LF.MH.eps)
# undebug(kernel.gaussian)
# 
# undebug(LF.aug.eps)
# undebug(LF.MH.eps)
# M<-LF.aug.eps(n,start[,1],theta.tune=var.tune.gaussian,eps.tune=eps.tune,data=data,kernel.func = kernel.gaussian,kernel.par = sig,dist = euclid.dist,dist.par = 0)
# 

M<-vector("list",ncol(start))

## running multiple chains for error augmented likelihood free MCMC sampler
cl <- makeCluster(3)
registerDoParallel(cl)
M<-vector("list",ncol(start))
M<-foreach (j = 1:ncol(start))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.eps(n,start[,j],theta.tune=var.tune.gaussian,eps.tune=eps.tune,data=data,kernel.func = kernel.gaussian,kernel.par = sig,dist = euclid.dist)
}
stopCluster(cl)

## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])
var(M[[1]][[1]][,5])

## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")

## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:ncol(start)){
  f[[j]]<-data.frame((M[[j]][[2]]),t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","mean.X3","mean.X4","MCMCse.X1","MCMCse.X2","MCMCse.X3","MCMCse.X4","start.label")

g<-do.call(rbind,g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","samples.X4","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X1 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X2 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X3 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p4<-p+ geom_line(mapping = aes(x=iterations,y=mean.X4,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X4 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X4,ymin=mean.X4-MCMCse.X4,ymax=mean.X4+MCMCse.X4,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
multiplot(p1, p2, p3, p4, cols=2)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X1",colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X2",colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X3",colour="Label for starting values")+thema
p7<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X4),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X4",colour="Label for starting values")+thema

multiplot(p4, p5, p6,p7,cols=2)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X1",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X1",title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X2",title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X2",title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X3",title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X3",title="After N")+thema

q7<-q+ geom_density(mapping = aes(x=samples.X4),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X4",title="After N/2")+thema
q8<-q+geom_density(mapping = aes(x=samples.X4),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X4",title="After N")+thema
multiplot(q1,q3,q5,q7,q2,q4,q6,q8, cols=2)

## estimated density for different starting values
q9<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X1",colour="Label of starting values")+thema
q10<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X2",colour="Label of starting values")+thema
q11<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X3",colour="Label of starting values")+thema
q12<-q+ geom_density(mapping = aes(x=samples.X4,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X4",colour="Label of starting values")+thema
multiplot(q9,q10,q11,q12, cols=2)

################################################################################################
## Testing the marginal space likelihood free MCMC sampler
## defining the parameters
p<-0.5
m1<--5
m2<-5

## simulating the data
m<-500 ## sample size
data<-data.sim(m = m,p = p,m1 = m1,m2 = m2)
plot(density(data[,2]))

## defining the chain size
n<-10000
start<-matrix(c(0.6,-2,7,0.3,-3,3,0.45,-4,4),3,3)
sig<-diag(3)
var.tune.unif<-c(0.06,3.8,3.9)
var.tune.gaussian<-c(0.06,6.4,4.95)
var.tune.epan<-c(0.06,0.98,0.99)
eps<-4
S<-50

# undebug(LF.aug.MS)
# undebug(LF.MH.MS)
# undebug(kernel.unif)
# 
# undebug(LF.aug.eps)
# undebug(LF.MH.eps)
# start<-c(0.6,-3,3)
# source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
# M<-LF.aug.MS(n = n,start=start,tune=var.tune.unif,data,eps,S=S,kernel.func = kernel.unif,kernel.par = 0,dist = euclid.dist,dist.par = 0)


## running multiple chains for marginal space likelihood free MCMC sampler
cl <- makeCluster(3)
registerDoParallel(cl)
M<-vector("list",ncol(start))
M<-foreach(j = 1:ncol(start))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.MS(n = n,start = start[,j],tune=var.tune.unif,data = data,eps = eps,S=S,kernel.func = kernel.unif,kernel.par = 0,dist = euclid.dist,dist.par = 0)
}
stopCluster(cl)

## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])

## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")

## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:ncol(start)){
  f[[j]]<-data.frame(M[[j]][[2]],t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","mean.X3","MCMCse.X1","MCMCse.X2","MCMCse.X3","start.label")

g<-do.call(rbind,g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X1 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X2 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X3 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
multiplot(p1, p2, p3,cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X1",colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X2",colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X3",colour="Label for starting values")+thema
multiplot(p4, p5, p6,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X1",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X1",title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X2",title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X2",title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X3",title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X3",title="After N")+thema

multiplot(q1,q3,q5,q2,q4,q6, cols=2)

## estimated density for different starting values
q7<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X1",colour="Label of starting values")+thema
q8<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X2",colour="Label of starting values")+thema
q9<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X3",colour="Label of starting values")+thema
multiplot(q7,q8,q9, cols=1)


################################################################################################
## Part 2 
library(mvtnorm)
## sampling from matrix normal distributions

## defining the multinormal dimensionality
p<-4
## observed data
data<-c(10,rep(0,(p-1)))

## defining the parameters
b<-0.1
s0<-1
var.tune.unif<-rep(1,p)
eps<-15

## defining the hyperparameters
A<-diag(c(100,rep(1,(p-1))))

## defining the chain size for pairs
n.pair<-10000

# start<-rep(0,p)
# source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
# undebug(LF.aug.pair)
# undebug(LF.MH.pair)
# undebug(kernel.unif)
# suff.par<-c(2,3)
# 
# 
# M<-LF.aug.pair(n = n.pair,start = start,tune=var.tune.unif,data = data,eps = eps, suff.stat = suff.stat.G, suff.par=suff.par, kernel.func = kernel.unif,kernel.par = 100,dist = euclid.dist,dist.par = 3,b=b)


# ## defining the starting value of the chain
# start<-matrix(0,p,3)
# 
# suff.par<-c(2,3)
# 
# M<-vector("list",3)
# 
# ## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
# cl <- makeCluster(3)
# registerDoParallel(cl)
# M<-vector("list",ncol(start))
# M<-foreach(j=1:ncol(start))%dopar%{
#   source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
#   M[[j]]<-LF.aug.G(n = n.pair,start = start[,j],tune=var.tune.unif,data = data,eps = eps, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
# }
# stopCluster(cl)


## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
comb<-combn(p,2)
start<-rep(0,p)
b<-0.1
sig<-diag(nrow(comb))
cl <- makeCluster(6)
registerDoParallel(cl)
M<-vector("list",ncol(comb))
M<-foreach(j = 1:ncol(comb))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.G(n = n.pair,start = start,tune=var.tune.unif,data = data,eps = eps, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.gaussian,kernel.par = sig,dist = mink.dist,dist.par = 3,b=b)
}
stopCluster(cl)

var.tune.unif<-matrix(NA_real_,p,ncol(comb))
start.unif<-matrix(NA_real_,p,ncol(comb))
for (j in 1:ncol(comb)){
  var.temp<-var(M[[j]][[1]][,2])
  for (i in 3:(p+1)){
    var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
  }
  var.tune.unif[,j]<-var.temp
  start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
}

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
i<-1
while(i<5){
  M<-vector("list",ncol(comb))
  cl <- makeCluster(6)
  registerDoParallel(cl)
  M<-foreach(j = 1:ncol(comb))%dopar%{
    source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
    M[[j]]<-LF.aug.G(n = n.pair,start = start.unif[,j],tune=var.tune.unif[,j],data = data,eps = eps, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.gaussian,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
  }
  stopCluster(cl)
  for (j in 1:ncol(comb)){
    var.temp<-var(M[[j]][[1]][,2])
    for (i in 3:(p+1)){
      var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
    }
    var.tune.unif[,j]<-var.temp
    start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
  }
  i<-i+1
}

## extracting the ESS.10K, ESS/s and accept rate for pairs
ESS.10K<-rep(NA_real_,ncol(comb))
ESS.rate<-rep(NA_real_,ncol(comb))
accept.rate<-rep(NA_real_,ncol(comb))
sys.time<-0
for (j in 1:ncol(comb)){
  ESS.10K[j]<-min(M[[j]][[6]])
  ESS.rate[j]<-min(M[[j]][[7]])
  accept.rate[j]<-min(M[[j]][[5]])
  sys.time<-sys.time+M[[j]][[8]]
}
sys.time

## calculating the mean ESS.10K, ESS/s and accept rate for pairs
ESS.10K.mean<-mean(ESS.10K)
ESS.10K.mean
ESS.rate.mean<-mean(ESS.rate)
ESS.rate.mean
accept.rate.mean<-mean(accept.rate)
accept.rate.mean

## estimating the correlation matrix
Lambda<-diag(1,p)
for (j in 1:ncol(comb)){
  r<-rank(M[[j]][[1]][,comb[1,j]+1])
  q<-rank(M[[j]][[1]][,comb[2,j]+1])
  eta<-matrix(NA_real_,nrow = n.pair,ncol = 2)
  eta[,1]<-qnorm(r/(n.pair+1))
  eta[,2]<-qnorm(q/(n.pair+1))
  Lambda[comb[1,j],comb[2,j]]<-cor(eta)[1,2]
  Lambda[comb[2,j],comb[1,j]]<-Lambda[comb[1,j],comb[2,j]]
}
Lambda

##################################################
## chains for single theta for univariate density estimations

## defining the multinormal dimensionality
p<-6
## observed data
data<-c(10,rep(1,(p-1)))

## defining the parameters
b<-0.1
s0<-1
var.tune.unif<-rep(1,p)
eps<-15

## defining the hyperparameters
A<-diag(c(100,1,1,1))

## defining the chain size for single theta
n.sing<-10000

## defining the starting value of the chain
start<-rep(0,p)
comb<-matrix(c(1:p),nrow = 1,ncol = p)

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
cl <- makeCluster(6)
registerDoParallel(cl)
M<-vector("list",ncol(comb))
M<-foreach(j = 1:ncol(comb))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.G(n = n.sing,start = start,tune=var.tune.unif,data = data,eps = eps, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
}
stopCluster(cl)

var.tune.unif<-matrix(NA_real_,p,ncol(comb))
start.unif<-matrix(NA_real_,p,ncol(comb))
for (j in 1:ncol(comb)){
  var.temp<-var(M[[j]][[1]][,2])
  for (i in 3:(p+1)){
    var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
  }
  var.tune.unif[,j]<-var.temp
  start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
}

## running multiple chains for theta individually in multinormal using likelihood free MCMC sampler
i<-1
while(i<10){
  M<-vector("list",ncol(comb))
  cl <- makeCluster(6)
  registerDoParallel(cl)
  M<-foreach(j = 1:ncol(comb))%dopar%{
    source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
    M[[j]]<-LF.aug.G(n = n.sing,start = start.unif[,j],tune=var.tune.unif[,j],data = data,eps = eps, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
  }
  stopCluster(cl)
  for (j in 1:ncol(comb)){
    var.temp<-var(M[[j]][[1]][,2])
    for (i in 3:(p+1)){
      var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
    }
    var.tune.unif[,j]<-var.temp
    start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
  }
  i<-i+1
}

## density estimation
require(ks)
dens<-list()
for (j in 1:ncol(comb)){
  dens[[j]]<-kde(x = as.matrix(M[[j]][[1]][j+1]),gridsize = 10000)
}

GC.sampler<-function(n,Lambda,dens){
  theta.prop<-matrix(NA_real_,p,n)
  weight<-rep(NA_real_,n)
  for(i in 1:n){
    theta.prop[,i]<-rmvnorm(n = 1, mean = rep(0,p),sigma = diag(1,p))
    eta<-rep(NA_real_,p)
    g_i<-rep(NA_real_,p)
    for (j in 1:p){
      eta[j]<-qnorm(pkde(theta.prop[j,i],dens[[j]]),mean = 0,sd = 1)
      g_i[j]<-predict(object = dens[[j]],x = theta.prop[j,i])
    }
    eta<-as.matrix(eta)
    num<-(1/sqrt(det(Lambda)))*as.numeric(exp(0.5*t(eta)%*%(diag(1,p)-solve(Lambda))%*%eta))*prod(g_i)
    denom<-dmvnorm(x = theta.prop,mean = rep(0,p),sigma = diag(1,p))
    weight[i]<-num/denom
  }
  theta.ind<-sample(x = c(1:n),prob = weight,replace = TRUE)
  theta<-theta.prop[,theta.ind]
  return(theta)
}

n<-10
theta<-GC.sampler(n=n,Lambda=Lambda,dens=dens)
head(t(theta))

f<-kde(x = rnorm(n = 100000,mean = 0,sd = 1))
predict(object = f,x = 0.1)
dnorm(x = 0.1,mean = 0,sd = 1)



###################################################
## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])
var(M[[1]][[1]][,5])


## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")


## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:ncol(start)){
  f[[j]]<-data.frame(M[[j]][[2]],t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","MCMCse.X1","MCMCse.X2","start.label")

g<-do.call(rbind,g)
head(g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","samples.X4","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X1 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y="Estimate of the expectation of X2 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
multiplot(p1, p2, cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p3<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X1",colour="Label for starting values")+thema
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X2",colour="Label for starting values")+thema
multiplot(p3, p4,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X1",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X1",title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X2",title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X2",title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X3",title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X3",title="After N")+thema

q7<-q+ geom_density(mapping = aes(x=samples.X4),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X4",title="After N/2")+thema
q8<-q+geom_density(mapping = aes(x=samples.X4),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X4",title="After N")+thema

multiplot(q1,q3,q5,q7,q2,q4,q6,q8, cols=2)

## estimated density for different starting values
q9<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X1",colour="Label of starting values")+thema
q10<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X2",colour="Label of starting values")+thema
q11<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X3",colour="Label of starting values")+thema
q12<-q+ geom_density(mapping = aes(x=samples.X4,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X4",colour="Label of starting values")+thema
multiplot(q9,q10,q11,q12, cols=2)

###############################################################################################
## Implementing error augmented sampler for Gaussian copula
h=0
Nt<-5
grand.ESS<-rep(NA_real_,Nt)
MSE
while (h<Nt){
library(mvtnorm)
## sampling from matrix normal distributions

## defining the multinormal dimensionality
p<-4
## observed data
data<-c(10,rep(0,(p-1)))

## defining the parameters
b<-0.1
s0<-1
var.tune.unif<-rep(1,p)
eps.tune<-1

## defining the hyperparameters
A<-diag(c(100,rep(1,(p-1))))

## defining the chain size for pairs
n.pair<-10000

start<-list()
start[[1]]<-rep(0,p)
start[[2]]<-10
source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
undebug(LF.aug.Geps)
undebug(LF.MH.Geps)
undebug(kernel.unif)
comb<-matrix(c(1,2),2,1)


M<-LF.aug.Geps(n = n.pair,start = start,theta.tune=var.tune.unif,eps.tune=eps.tune,data = data, suff.stat = suff.stat.G, suff.par=as.matrix(comb), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)


## defining the starting value of the chain
start<-list()
for (j in 1:3){
  start[[j]]<-list()
  start[[j]][[1]]<-rep(0,p)
  start[[j]][[2]]<-10
}
start

var.tune.unif<-c(84,59,723,530)
eps.tune<-12

comb<-matrix(c(1,2),2,1)

M<-vector("list",3)

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
cl <- makeCluster(3)
registerDoParallel(cl)
M<-vector("list",3)
M<-foreach(j=1:3)%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.Geps(n = n.pair,start = start[[j]],theta.tune=var.tune.unif,eps.tune=eps.tune,data = data, suff.stat = suff.stat.G, suff.par=as.matrix(comb), kernel.func = kernel.unif,kernel.par = 100,dist = euclid.dist,dist.par = 3,b=b)
}
stopCluster(cl)


## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
comb<-combn(p,2)
## defining the starting value of the chain
start.unif<-list()
for (j in 1:ncol(comb)){
  start.unif[[j]]<-list()
  start.unif[[j]][[1]]<-rep(0,p)
  start.unif[[j]][[2]]<-10
}

library(foreach)
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
M<-vector("list",ncol(comb))
M<-foreach(j = 1:ncol(comb))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.Geps(n = n.pair,start = start.unif[[j]],theta.tune=var.tune.unif,eps.tune=eps.tune,data = data, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = euclid.dist,dist.par = 3,b=b)
}
stopCluster(cl)

var.tune.unif<-matrix(NA_real_,p,ncol(comb))
for (j in 1:ncol(comb)){
  var.temp<-var(M[[j]][[1]][,2])
  for (i in 3:(p+1)){
    var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
  }
  var.tune.unif[,j]<-var.temp
  start.unif[[j]][[1]]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-c(1,ncol(M[[j]][[1]]))])
  start.unif[[j]][[2]]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),ncol(M[[j]][[1]])])
}
eps.tune<-var(M[[j]][[1]][,p+2])

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
i<-1
while(i<5){
  M<-vector("list",ncol(comb))
  cl <- makeCluster(6)
  registerDoParallel(cl)
  M<-foreach(j = 1:ncol(comb))%dopar%{
    source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
    M[[j]]<-LF.aug.Geps(n = n.pair,start = start.unif[[j]],theta.tune=var.tune.unif[,j],eps.tune=eps.tune,data = data, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
  }
  stopCluster(cl)
  for (j in 1:ncol(comb)){
    var.temp<-var(M[[j]][[1]][,2])
    for (i in 3:(p+1)){
      var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
    }
    var.tune.unif[,j]<-var.temp
    start.unif[[j]][[1]]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-c(1,ncol(M[[j]][[1]]))])
    start.unif[[j]][[2]]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),ncol(M[[j]][[1]])])
  }
  i<-i+1
}

## extracting the ESS.10K, ESS/s and accept rate for pairs
ESS.10K<-rep(NA_real_,ncol(comb))
ESS.rate<-rep(NA_real_,ncol(comb))
accept.rate<-rep(NA_real_,ncol(comb))
sys.time<-0
for (j in 1:ncol(comb)){
  ESS.10K[j]<-min(M[[j]][[6]])
  ESS.rate[j]<-min(M[[j]][[7]])
  accept.rate[j]<-min(M[[j]][[5]])
  sys.time<-sys.time+M[[j]][[8]]
}
sys.time

## calculating the mean ESS.10K, ESS/s and accept rate for pairs
ESS.10K.mean<-mean(ESS.10K)
ESS.10K.mean
ESS.rate.mean<-mean(ESS.rate)
grand.ESS[h]<-ESS.rate.mean
accept.rate.mean<-mean(accept.rate)
accept.rate.mean
h=h+1
}
grand.ESS

## estimating the correlation matrix
Lambda<-diag(1,p)
for (j in 1:ncol(comb)){
  r<-rank(M[[j]][[1]][,comb[1,j]+1])
  q<-rank(M[[j]][[1]][,comb[2,j]+1])
  eta<-matrix(NA_real_,nrow = n.pair,ncol = 2)
  eta[,1]<-qnorm(r/(n.pair+1))
  eta[,2]<-qnorm(q/(n.pair+1))
  Lambda[comb[1,j],comb[2,j]]<-cor(eta)[1,2]
  Lambda[comb[2,j],comb[1,j]]<-Lambda[comb[1,j],comb[2,j]]
}
Lambda

###################################################
## checking sample variances for tuning
n<-n.pair
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])
var(M[[1]][[1]][,5])
var(M[[1]][[1]][,6])


## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,c(2,3,6)],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")


## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:3){
  f[[j]]<-data.frame(M[[j]][[2]],t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","mean.X3","MCMCse.X1","MCMCse.X2","MCMCse.X3","start.label")

g<-do.call(rbind,g)
head(g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","samples.X4","samples.X5","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y=expression(paste("parameter i")), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y=expression(paste("parameter j")), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))+labs(x="Number of samples",y=expression(paste(epsilon)), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%100)==0))
multiplot(p1, p2, p3,cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for i")),colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for j")),colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for ",epsilon)),colour="Label for starting values")+thema
multiplot(p4, p5, p6,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For i")),title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For i")),title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For j")),title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For j")),title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X5),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For ", epsilon)),title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X5),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For ", epsilon)),title="After N")+thema

multiplot(q1,q3,q5,q2,q4,q6, cols=2)

## estimated density for different starting values
q7<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For i")),colour="Label of starting values")+thema
q8<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For j")),colour="Label of starting values")+thema
q9<-q+ geom_density(mapping = aes(x=samples.X5,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For ", epsilon)),colour="Label of starting values")+thema
multiplot(q7,q8,q9, cols=1)
#############################################################################################
## Implementing the multiple augmented LF MCMC sampler via gaussian copula
library(mvtnorm)
## sampling from matrix normal distributions

## defining the multinormal dimensionality
p<-16
## observed data
data<-c(10,rep(1,(p-1)))

## defining the parameters
b<-0.1
s0<-1
var.tune.unif<-rep(1,p)
eps<-15
S=10

## defining the hyperparameters
A<-diag(c(100,rep(1,(p-1))))

## defining the chain size for pairs
n.pair<-10000

start<-rep(0,p)
source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
undebug(LF.aug.pair)
undebug(LF.MH.pair)
undebug(kernel.unif)
comb<-as.matrix(c(1,2),2,1)


M<-LF.aug.Gmult(n = n.pair,start = start,tune=var.tune.unif,data = data,eps = eps, S=S, suff.stat = suff.stat.G, suff.par=as.matrix(comb), kernel.func = kernel.unif,kernel.par = 100,dist = euclid.dist,dist.par = 3,b=b)


# defining the starting value of the chain
start<-matrix(0,p,3)

comb<-as.matrix(c(1,2),2,1)

M<-vector("list",3)

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
cl <- makeCluster(3)
registerDoParallel(cl)
M<-vector("list",ncol(start))
M<-foreach(j=1:ncol(start))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.Gmult(n = n.pair,start = start[,j],tune=var.tune.unif,data = data,eps = eps, S=S,suff.stat = suff.stat.G, suff.par=as.matrix(comb), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
}
stopCluster(cl)


## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
comb<-combn(p,2)
start<-rep(0,p)
b<-0.1
cl <- makeCluster(6)
registerDoParallel(cl)
M<-vector("list",ncol(comb))
M<-foreach(j = 1:ncol(comb))%dopar%{
  source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
  M[[j]]<-LF.aug.Gmult(n = n.pair,start = start,tune=var.tune.unif,data = data,eps = eps,S=S, suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
}
stopCluster(cl)

var.tune.unif<-matrix(NA_real_,p,ncol(comb))
start.unif<-matrix(NA_real_,p,ncol(comb))
for (j in 1:ncol(comb)){
  var.temp<-var(M[[j]][[1]][,2])
  for (i in 3:(p+1)){
    var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
  }
  var.tune.unif[,j]<-var.temp
  start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
}

## running multiple chains for theta pairs in multinormal using likelihood free MCMC sampler
i<-1
while(i<5){
  M<-vector("list",ncol(comb))
  cl <- makeCluster(6)
  registerDoParallel(cl)
  M<-foreach(j = 1:ncol(comb))%dopar%{
    source(file = "Box Sync/PSU/Fall 2015/Statistical Computing (STAT 540)/Project/code/functions.R")
    M[[j]]<-LF.aug.Gmult(n = n.pair,start = start.unif[,j],tune=var.tune.unif[,j],data = data,eps = eps, S=S,suff.stat = suff.stat.G, suff.par=as.matrix(comb[,j]), kernel.func = kernel.unif,kernel.par = 100,dist = mink.dist,dist.par = 3,b=b)
  }
  stopCluster(cl)
  for (j in 1:ncol(comb)){
    var.temp<-var(M[[j]][[1]][,2])
    for (i in 3:(p+1)){
      var.temp<-c(var.temp,var(M[[j]][[1]][,i]))
    }
    var.tune.unif[,j]<-var.temp
    start.unif[,j]<-as.numeric(M[[j]][[1]][nrow(M[[j]][[1]]),-1])
  }
  i<-i+1
}

## extracting the ESS.10K, ESS/s and accept rate for pairs
ESS.10K<-rep(NA_real_,ncol(comb))
ESS.rate<-rep(NA_real_,ncol(comb))
accept.rate<-rep(NA_real_,ncol(comb))
sys.time<-0
for (j in 1:ncol(comb)){
  ESS.10K[j]<-min(M[[j]][[6]])
  ESS.rate[j]<-min(M[[j]][[7]])
  accept.rate[j]<-min(M[[j]][[5]])
  sys.time<-sys.time+M[[j]][[8]]
}
sys.time

## calculating the mean ESS.10K, ESS/s and accept rate for pairs
ESS.10K.mean<-mean(ESS.10K)
ESS.10K.mean
ESS.rate.mean<-mean(ESS.rate)
ESS.rate.mean
accept.rate.mean<-mean(accept.rate)
accept.rate.mean

## estimating the correlation matrix
Lambda<-diag(1,p)
for (j in 1:ncol(comb)){
  r<-rank(M[[j]][[1]][,comb[1,j]+1])
  q<-rank(M[[j]][[1]][,comb[2,j]+1])
  eta<-matrix(NA_real_,nrow = n.pair,ncol = 2)
  eta[,1]<-qnorm(r/(n.pair+1))
  eta[,2]<-qnorm(q/(n.pair+1))
  Lambda[comb[1,j],comb[2,j]]<-cor(eta)[1,2]
  Lambda[comb[2,j],comb[1,j]]<-Lambda[comb[1,j],comb[2,j]]
}
Lambda

## creating the plots
## defining the matrix sizes vector
x<-c(4,6,9,12,16)
## ESS rates
Basic<-c(20.05,20.51,17.2,16.05,14.15)
Error.aug.<-c(21.21,25.35,22.44,23.12,18.98)
Multiple.aug.<-c(15.50,15.58,12.18,11.56,10.64)
df<-data.frame(cbind(x,Basic,Error.aug.,Multiple.aug.))
library(reshape2)
df<-melt(df,id=c("x"))
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 16, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 16, colour = "black"),
        legend.text  = element_text(size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "black", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 16, colour = "black"),
        title =  element_text(size = 8, face = "bold"))
p<-ggplot(data=df)
p+geom_line(mapping = aes(x=x, y=value,colour=factor(variable)))+thema+labs(x="Matrix sizes (n*p)",y="ESS/s", colour="LF MCMC Alg.")+scale_x_discrete(breaks = c(4,6,9,12,16), labels=c("2*2","3*2","3*3","4*3","4*4"),expand=c(0.05,-2.5))#xlim(c(3.8,16.2))

## TMSE
Basic<-c(2.3*10^(-4),8.5*10^(-3),8.3*10^(-2),8.7*10^(-2),9.5*10^(-2))
Error.aug.<-c(2.5*10^(-4),9.1*10^(-3),8.9*10^(-2),1.3*10^(-1),1.5*10^(-1))
Multiple.aug.<-c(2.3*10^(-4),8.5*10^(-3),7.3*10^(-2),8.5*10^(-2),9.8*10^(-2))
df<-data.frame(cbind(x,Basic,Error.aug.,Multiple.aug.))
library(reshape2)
df<-melt(df,id=c("x"))
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 16, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 16, colour = "black"),
        legend.text  = element_text(size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "black", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 16, colour = "black"),
        title =  element_text(size = 8, face = "bold"))
p<-ggplot(data=df)
p+geom_line(mapping = aes(x=x, y=value,colour=factor(variable)))+thema+labs(x="Matrix sizes (n*p)",y="TMSE for 50 simulations", colour="LF MCMC Alg.")+scale_x_discrete(breaks = c(4,6,9,12,16), labels=c("2*2","3*2","3*3","4*3","4*4"),expand=c(0.05,-2.5))#xlim(c(3.8,16.2))

## Total Run time
Basic<-c(46.5,112.1,302.2,576.8,1185.7)
Error.aug.<-c(106.5,252.6,798.4,1361.3,2455.9)
Multiple.aug.<-c(156.5,396.9,1023.5,1832.5,3006.4)
df<-data.frame(cbind(x,Basic,Error.aug.,Multiple.aug.))
library(reshape2)
df<-melt(df,id=c("x"))
thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 16, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 16, colour = "black"),
        legend.text  = element_text(size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "black", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 16, colour = "black"),
        title =  element_text(size = 8, face = "bold"))
p<-ggplot(data=df)
p+geom_line(mapping = aes(x=x, y=value,colour=factor(variable)))+thema+labs(x="Matrix sizes (n*p)",y="Mean total run time per simulation", colour="LF MCMC Alg.")+scale_x_discrete(breaks = c(4,6,9,12,16), labels=c("2*2","3*2","3*3","4*3","4*4"),expand=c(0.05,-2.5))#xlim(c(3.8,16.2))


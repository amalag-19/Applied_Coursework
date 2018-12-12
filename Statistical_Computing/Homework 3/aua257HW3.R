## Please note that each question is seperated by a full line of '#'
## and subparts are seperated by half lines of '#'. Also run each 
## question from the begining as a whole to get the correct output
## since variable names/function names may overlap for different
## questions.

## Please install the necessary packages mentioned in the libraries beforehand
#####################################################################
## Q1 (a)
## libraries
{library(ggplot2)
library(pracma)
library(reshape)
library(fExtremes)
library(plyr)
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
}}

## Reading the data
data<-read.table("http://sites.stat.psu.edu/~mharan/540/hwdir/AtlanticGEV.dat",header = FALSE)
head(data)

## defining the negative of log likelihood for calculating MLE from given data
log.gev.dens<-function(x,mu,zeta,sigma){
  t4<--log(sigma)
  t1<-(1+zeta*((x-mu)/sigma))
  t2<--(t1^(-1/zeta))
  t3<-(-1-(1/zeta))*log(t1)
  value<-t2+t3+t4
  invisible(value)
}
log.lik<-function(param){
  mu<-param[1]
  zeta<-param[2]
  sigma<-param[3]
  if((sigma>0&&(zeta>0&&(prod(as.matrix(data)>(mu-(sigma/zeta)))==1))||(zeta<0&&(prod(as.matrix(data)<(mu-(sigma/zeta)))==1))||(zeta==0))){
    t4<--log(sigma)
    value<-sum(log.gev.dens(data[,1],mu = mu,zeta = zeta, sigma = sigma))
  }
  else {value<--(10^5)}
  invisible(-value)
}

## estimating initial parameter values by evaluating the log likelihood over a grid
## and choosing the maximum
mu.grid<-seq(1,5,length.out = 10)
zeta.grid<-seq(0.1,5,length.out = 10)
sigma.grid<-seq(0.05,5.5,length.out = 10)
par.grid<-expand.grid(mu.grid,zeta.grid,sigma.grid)
head(par.grid)
initial.index<-which.min(apply(X = par.grid,MARGIN = 1,FUN = log.lik))
initial.value<-par.grid[initial.index,]

## optimization using BFGS over above intial value
BFGS<-optim(par = initial.value, fn = log.lik, method = "BFGS",hessian = TRUE)
Info.inv<-solve(BFGS$hessian)

## extracting the estimates of mu, zeta and sigma for BFGS
mu.hat.BFGS<-as.numeric(BFGS$par[1])
mu.hat.BFGS
zeta.hat.BFGS<-as.numeric(BFGS$par[2])
zeta.hat.BFGS
sigma.hat.BFGS<-as.numeric(BFGS$par[3])
sigma.hat.BFGS

## caclulating the standard errors of mu, zeta and sigma for BFGS
mu.se.BFGS<-sqrt(Info.inv[1,1])
mu.se.BFGS
zeta.se.BFGS<-sqrt(Info.inv[2,2])
zeta.se.BFGS
sigma.se.BFGS<-sqrt(Info.inv[3,3])
sigma.se.BFGS

## caclulating the 95% confidence bounds of mu, zeta and sigma for BFGS
mu.lowc.BFGS<-mu.hat.BFGS-1.96*mu.se.BFGS
mu.upc.BFGS<-mu.hat.BFGS+1.96*mu.se.BFGS
zeta.lowc.BFGS<-zeta.hat.BFGS-1.96*zeta.se.BFGS
zeta.upc.BFGS<-zeta.hat.BFGS+1.96*zeta.se.BFGS
sigma.lowc.BFGS<-sigma.hat.BFGS-1.96*sigma.se.BFGS
sigma.upc.BFGS<-sigma.hat.BFGS+1.96*sigma.se.BFGS
  

## optimization using Nelder-Mead over above intial value
NM<-optim(par=c(1,3,1), fn = log.lik, method = "Nelder-Mead", hessian = TRUE)
Info.inv<-solve(NM$hessian)

## extracting the estimates of mu, zeta and sigma for Nelder-Mead
mu.hat.NM<-as.numeric(NM$par[1])
mu.hat.NM
zeta.hat.NM<-as.numeric(NM$par[2])
zeta.hat.NM
sigma.hat.NM<-as.numeric(NM$par[3])
sigma.hat.NM

## caclulating the standard errors of mu, zeta and sigma for Nelder-Mead
mu.se.NM<-sqrt(Info.inv[1,1])
mu.se.NM
zeta.se.NM<-sqrt(Info.inv[2,2])
zeta.se.NM
sigma.se.NM<-sqrt(Info.inv[3,3])
sigma.se.NM

## caclulating the 95% confidence bounds of mu, zeta and sigma for Nelder-Mead
mu.lowc.NM<-mu.hat.NM-1.96*mu.se.NM
mu.upc.NM<-mu.hat.NM+1.96*mu.se.NM
zeta.lowc.NM<-zeta.hat.NM-1.96*zeta.se.NM
zeta.upc.NM<-zeta.hat.NM+1.96*zeta.se.NM
sigma.lowc.NM<-sigma.hat.NM-1.96*sigma.se.NM
sigma.upc.NM<-sigma.hat.NM+1.96*sigma.se.NM

## optimization using Newton-Raphson over intial value obtained from Nelder Mead
## defining log likelihood for newton raphson
log.lik.n<-function(param){
  mu<-param[1]
  zeta<-param[2]
  sigma<-param[3]
  t4<--log(sigma)
  value<-sum(apply(X = as.matrix(data), MARGIN = 1, FUN = function(y){
    t1<-(1+zeta*((y-mu)/sigma))
    t2<--(t1^(-1/zeta))
    t3<-(-1-(1/zeta))*log(t1)
    t5<-t2+t3+t4
    invisible(t5)
    }))
  invisible(value)
}

## number of iterations for newton raphson
n<-100
## initializing the parameter estimates matrix 
param.NR<-matrix(NA_real_,3,n)
#param.NR[,1]<-c(jitter(NM$par,amount = 0))
param.NR[,1]<-c(1.64,0.18,0.09)

## initializing list for hessian matrix at each iterate
H<-list()
for (i in 1:(n-1)){
  H[[i]]<-hessian(log.lik.n,param.NR[,i],h=10^(-7))
  param.NR[,i+1]<-param.NR[,i]-(solve(H[[i]]))%*%grad(log.lik.n,param.NR[,i],heps=10^(-7))
  print(param.NR[,i+1])
}

Info.inv<--solve(hessian(log.lik.n,param.NR[,n],h=10^(-7)))

## extracting the estimates of mu, zeta and sigma for Newton Raphson
mu.hat.NR<-param.NR[1,n]
mu.hat.NR
zeta.hat.NR<-param.NR[2,n]
zeta.hat.NR
sigma.hat.NR<-param.NR[3,n]
sigma.hat.NR

## caclulating the standard errors of mu, zeta and sigma for Newton Raphson
mu.se.NR<-sqrt(Info.inv[1,1])
mu.se.NR
zeta.se.NR<-sqrt(Info.inv[2,2])
zeta.se.NR
sigma.se.NR<-sqrt(Info.inv[3,3])
sigma.se.NR

## caclulating the 95% confidence bounds of mu, zeta and sigma for Newton Raphson
mu.lowc.NR<-mu.hat.NR-1.96*mu.se.NR
mu.upc.NR<-mu.hat.NR+1.96*mu.se.NR
zeta.lowc.NR<-zeta.hat.NR-1.96*zeta.se.NR
zeta.upc.NR<-zeta.hat.NR+1.96*zeta.se.NR
sigma.lowc.NR<-sigma.hat.NR-1.96*sigma.se.NR
sigma.upc.NR<-sigma.hat.NR+1.96*sigma.se.NR

#########################################
## Q1 (c)
## plotting the histogram of samples and comparing with true GEV
x.vec.BFGS<-seq(mu.hat.BFGS-(sigma.hat.BFGS/zeta.hat.BFGS)+0.001,3,length.out = nrow(data))
y.vec.BFGS<-dgev(x.vec.BFGS, xi = zeta.hat.BFGS, mu = mu.hat.BFGS, beta = sigma.hat.BFGS)

x.vec.NM<-seq(mu.hat.NM-(sigma.hat.NM/zeta.hat.NM)+0.001,3,length.out = nrow(data))
y.vec.NM<-dgev(x.vec.NM, xi = zeta.hat.NM, mu = mu.hat.NM, beta = sigma.hat.NM)

x.vec.NR<-seq(mu.hat.NR-(sigma.hat.NR/zeta.hat.NR)+0.001,3,length.out = nrow(data))
y.vec.NR<-dgev(x.vec.NR, xi = zeta.hat.NR, mu = mu.hat.NR, beta = sigma.hat.NR)

df1<-data.frame(data,"BFGS"=x.vec.BFGS,"NM"=x.vec.NM,"NR"=x.vec.NR)
df2<-data.frame(data,"BFGS"=y.vec.BFGS,"NM"=y.vec.NM,"NR"=y.vec.NR)
df3<-melt(df1, id=c("V1"))
names(df3)<-c("data","method","x.vec")
df4<-data.frame(cbind(df3, "y.vec"=melt(df2, id=c("V1"))$value))
head(df4)

p<-ggplot(data = df4)
p+geom_histogram(mapping = aes(x=data,y=..density..),fill="steelblue",binwidth=0.04)+labs(x="x",y="Density")+geom_line(mapping = aes(x = x.vec,y = y.vec,colour=method))

#####################################################################################
## Q2
## First M-H algo:
## defining the function Expectation1 for Variable at a time Metropolis Hastings
## algortithm and to calculate Monte Carlo estimates and standard errors
## for x
Expectation1<-function(n,start,var.tune){
  ## libraries
  {library(ggplot2)
    library(tmvtnorm)
    library(fExtremes)
    library(plyr)
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
    }}
  
  ## Reading the data
  data<-read.table("http://sites.stat.psu.edu/~mharan/540/hwdir/AtlanticGEV.dat",header = FALSE)
  head(data)
  
  ## sourcing batchmeans function used to calculate MCMC standard
  ## errors later
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  log.gev.dens<-function(x,mu,zeta,sigma){
    t4<--log(sigma)
    t1<-(1+zeta*((x-mu)/sigma))
    t2<--(t1^(-1/zeta))
    t3<-(-1-(1/zeta))*log(t1)
    value<-t2+t3+t4
    invisible(value)
  }
  
  ## defining log likelihood function as a part of target
  log.lik<-function(param){
    mu<-param[1]
    zeta<-param[2]
    sigma<-param[3]
    if((sigma>0&&(zeta>0&&(prod(as.matrix(data)>(mu-(sigma/zeta)))==1))||(zeta<0&&(prod(as.matrix(data)<(mu-(sigma/zeta)))==1))||(zeta==0))){
      t4<--log(sigma)
      value<-sum(log.gev.dens(data,mu = mu,zeta = zeta, sigma = sigma))
    }
    else {value<-NA}
    invisible(value)
  }
  
  ## Variable at a time MCMC
  ## Metropolis Hastings function for mu with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.mu<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rnorm(n=1, mean=current.state[1], sd = sqrt(var.tune))
    if (!is.na(log.lik(c(x.star,current.state[2],current.state[3])))){
      num<-log.lik(c(x.star,current.state[2],current.state[3]))+dnorm(current.state[1], mean=x.star, sd = sqrt(var.tune),log = TRUE)
      denom<-log.lik(current.state)+dnorm(x.star, mean=current.state[1], sd=sqrt(var.tune),log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[1]}
    }
    else {
      next.state<-current.state[1]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  ## Metropolis Hastings function for zeta with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.zeta<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rnorm(n=1, mean=current.state[2], sd = sqrt(var.tune))
    if(!is.na(log.lik(c(current.state[1],x.star,current.state[3])))){
      num<-log.lik(c(current.state[1],x.star,current.state[3]))+dnorm(current.state[2], mean=x.star, sd = sqrt(var.tune),log = TRUE)
      denom<-log.lik(current.state)+dnorm(x.star, mean=current.state[2], sd=sqrt(var.tune),log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[2]}
    }
    else {
      next.state<-current.state[2]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  ## Metropolis Hastings function for sigma with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.sigma<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rtmvnorm(n=1, mean=current.state[3], sigma = var.tune, lower = 0)
    
    if(!is.na(log.lik(c(current.state[1],current.state[2],x.star)))){
      num<-log.lik(c(current.state[1],current.state[2],x.star))-log(x.star)+dtmvnorm(current.state[3], mean=as.vector(x.star), sigma = var.tune, lower = 0,log = TRUE)
      denom<-log.lik(current.state)-log(current.state[3])+dtmvnorm(x.star, mean=as.vector(current.state[3]), sigma = var.tune, lower = 0, log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[3]}
    }
    else {
      next.state<-current.state[3]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  
  ## Initializing the Markov chain for beta_1
  x<-matrix(NA_real_,3,n)
  ## Defining the initial value for the chain
  x[,1]<-start
  ## Initializing the accept count for mu, zeta, sigma used to calculate acceptance
  ## rate of x
  accept.mu<-0
  accept.zeta<-0
  accept.sigma<-0
  ## loop for RWM updates
  sys.time<-system.time(for(i in 1:(n-1)){
    curr.state<-x[,i]
    
    temp<-MH.mu(curr.state,var.tune[1])
    curr.state[1]<-temp[1]
    accept.mu<-accept.mu+temp[2]
    
    temp<-MH.zeta(curr.state,var.tune[2])
    curr.state[2]<-temp[1]
    accept.zeta<-accept.zeta+temp[2]
    
    temp<-MH.sigma(curr.state,var.tune[3])
    curr.state[3]<-temp[1]
    accept.sigma<-accept.sigma+temp[2]
    
    x[,i+1]<-curr.state
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(x)[,1],"X2"=t(x)[,2],"X3"=t(x)[,3])
  ## calculating the acceptance rates for mu, zeta and sigma
  acceptance.rate<-c((accept.mu/n),(accept.zeta/n),(accept.sigma/n))
  
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
    MCMC.means[,j]<-rowSums(x[,1:m[j]])/m[j]
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(x[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(x[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(x[3,1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

## Second M-H algo:
## defining the function Expectation2 for Variable at a time Metropolis Hastings
## algortithm and to calculate Monte Carlo estimates and standard errors
## for x
Expectation2<-function(n,start,var.tune){
  ## libraries
  {library(ggplot2)
    library(tmvtnorm)
    library(fExtremes)
    library(plyr)
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
    }}
  
  ## Reading the data
  data<-read.table("http://sites.stat.psu.edu/~mharan/540/hwdir/AtlanticGEV.dat",header = FALSE)
  head(data)
  
  ## sourcing batchmeans function used to calculate MCMC standard
  ## errors later
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
  
  log.gev.dens<-function(x,mu,zeta,sigma){
    t4<--log(sigma)
    t1<-(1+zeta*((x-mu)/sigma))
    t2<--(t1^(-1/zeta))
    t3<-(-1-(1/zeta))*log(t1)
    value<-t2+t3+t4
    invisible(value)
  }
  
  ## defining log likelihood function as a part of target
  log.lik<-function(param){
    mu<-param[1]
    zeta<-param[2]
    sigma<-param[3]
    if((sigma>0&&(zeta>0&&(prod(as.matrix(data)>(mu-(sigma/zeta)))==1))||(zeta<0&&(prod(as.matrix(data)<(mu-(sigma/zeta)))==1))||(zeta==0))){
      t4<--log(sigma)
      value<-sum(log.gev.dens(data,mu = mu,zeta = zeta, sigma = sigma))
    }
    else {value<-NA}
    invisible(value)
  }
  
  ## Variable at a time MCMC
  ## Metropolis Hastings function for mu with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.mu<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rnorm(n=1, mean=current.state[1], sd = sqrt(var.tune))
    if (!is.na(log.lik(c(x.star,current.state[2],current.state[3])))){
      num<-log.lik(c(x.star,current.state[2],current.state[3]))+dnorm(current.state[1], mean=x.star, sd = sqrt(var.tune),log = TRUE)
      denom<-log.lik(current.state)+dnorm(x.star, mean=current.state[1], sd=sqrt(var.tune),log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[1]}
    }
    else {
      next.state<-current.state[1]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  ## Metropolis Hastings function for zeta with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.zeta<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rnorm(n=1, mean=current.state[2], sd = sqrt(var.tune))
    if(!is.na(log.lik(c(current.state[1],x.star,current.state[3])))){
      num<-log.lik(c(current.state[1],x.star,current.state[3]))+dnorm(current.state[2], mean=x.star, sd = sqrt(var.tune),log = TRUE)
      denom<-log.lik(current.state)+dnorm(x.star, mean=current.state[2], sd=sqrt(var.tune),log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[2]}
    }
    else {
      next.state<-current.state[2]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  ## Metropolis Hastings function for sigma with inputs as the variance of proposal
  ## (tuning parameter) and current state of the MC
  MH.sigma<-function(current.state,var.tune){
    ## sampling x.star from normal proposal with mean as current state
    ## of x and specified tuning parameter
    x.star<-rgamma(n=1, shape=(current.state[3]^2)/var.tune, rate = current.state[3]/var.tune)
    
    if(!is.na(log.lik(c(current.state[1],current.state[2],x.star)))){
      num<-log.lik(c(current.state[1],current.state[2],x.star))-log(x.star)+dgamma(current.state[3], shape = (x.star^2)/var.tune, rate = x.star/var.tune,log = TRUE)
      denom<-log.lik(current.state)-log(current.state[3])+dgamma(x.star, shape = (current.state[3]^2)/var.tune, rate = current.state[3]/var.tune, log = TRUE)
      ## defining the acceptance probability for sampled x.star on log
      ## scale
      accept.probab<-num-denom
      ## sampling u from uniform(0,1) to check for acceptance
      u<-runif(1, min=0, max=1)
      ## initializing the indicator flag=0 to check if the sampled x.star
      ## will be accepted
      flag<-0
      ## if-else to define the next state of the chain based on acceptance
      ## probability
      if(log(u)<=accept.probab){
        flag<-1
        next.state<-x.star
      }
      else {next.state<-current.state[3]}
    }
    else {
      next.state<-current.state[3]
      flag<-0
    }
    ## returning the next state and indicator if the sampled value was
    ## accepted
    return(c(next.state,flag))
  }
  
  ## Initializing the Markov chain for beta_1
  x<-matrix(NA_real_,3,n)
  ## Defining the initial value for the chain
  x[,1]<-start
  ## Initializing the accept count for mu, zeta, sigma used to calculate acceptance
  ## rate of x
  accept.mu<-0
  accept.zeta<-0
  accept.sigma<-0
  ## loop for RWM updates
  sys.time<-system.time(for(i in 1:(n-1)){
    curr.state<-x[,i]
    
    temp<-MH.mu(curr.state,var.tune[1])
    curr.state[1]<-temp[1]
    accept.mu<-accept.mu+temp[2]
    
    temp<-MH.zeta(curr.state,var.tune[2])
    curr.state[2]<-temp[1]
    accept.zeta<-accept.zeta+temp[2]
    
    temp<-MH.sigma(curr.state,var.tune[3])
    curr.state[3]<-temp[1]
    accept.sigma<-accept.sigma+temp[2]
    
    x[,i+1]<-curr.state
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(x)[,1],"X2"=t(x)[,2],"X3"=t(x)[,3])
  ## calculating the acceptance rates for mu, zeta and sigma
  acceptance.rate<-c((accept.mu/n),(accept.zeta/n),(accept.sigma/n))
  
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
    MCMC.means[,j]<-rowSums(x[,1:m[j]])/m[j]
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(x[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(x[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(x[3,1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

## defining the chain size
n<-10000
start<-matrix(NA_real_,3,3)
start[,1]<-c(1.65,0.172,0.098)
for (j in 2:3){
  start[,j]<-jitter(start[,1],amount = 0)
}

## Tuning parameters chosen after many pilot runs
var.tune1<-c(0.00000436,0.000420,2.88*(10^(-6)))
var.tune2<-c(0.00000446,0.000433,2.94*(10^(-6)))

M<-vector("list",ncol(start))

library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
## For running second M-H algo, comment the Expectation1 line in following for each
## loop and uncomment the Expectation2 line
M1<-foreach (j = 1:3) %dopar% {
  Expectation1(n,start[,j],var.tune1)
  #Expectation2(n,start[,j],var.tune2)
}
stopCluster(cl)
M<-M1

## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])

## calculating the sample quantiles for credible intervals
mu.lowc<-quantile(M[[1]][[1]][,2],0.025)
mu.upc<-quantile(M[[1]][[1]][,2],0.975)
zeta.lowc<-quantile(M[[1]][[1]][,3],0.025)
zeta.upc<-quantile(M[[1]][[1]][,3],0.975)
sig.lowc<-quantile(M[[1]][[1]][,4],0.025)
sig.upc<-quantile(M[[1]][[1]][,4],0.975)

## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

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

## acceptance rate rate for 3 chains
M[[1]][[5]]
M[[2]][[5]]
M[[3]][[5]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## Effective samples for every 10,000 samples generated
min(M[[1]][[6]])

## Number of effective samples generated per second
min(M[[1]][[7]])

###########################################
## Q2 (b)
## plotting
library(ggplot2)
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
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y=expression(paste("Estimate of the E(", mu ,") with MCMCse")), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y=expression(paste("Estimate of the E(", xi ,") with MCMCse")), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y=expression(paste("Estimate of the E(", sigma ,") with MCMCse")), colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
multiplot(p1, p2, p3,cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for", mu)),colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for", xi)),colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y=expression(paste("MCMCse for", sigma)),colour="Label for starting values")+thema
multiplot(p4, p5, p6,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For", mu)),title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For", mu)),title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For", xi)),title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For", xi)),title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y=expression(paste("For", sigma)),title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y=expression(paste("For", sigma)),title="After N")+thema

multiplot(q1,q3,q5,q2,q4,q6, cols=2)

## estimated density for different starting values
q7<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For", mu)),colour="Label of starting values")+thema
q8<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For", xi)),colour="Label of starting values")+thema
q9<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y=expression(paste("For", sigma)),colour="Label of starting values")+thema
multiplot(q7,q8,q9, cols=1)

####################################################################################
## Q3 
library(tmvtnorm)
library(fExtremes)
## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

## defining the parameters
{mu<-3
zeta<-0.4
sigma<-0.8}

## log likelihood for calculating MLE from simulated data
log.gev.dens<-function(x,mu,zeta,sigma){
  t4<--log(sigma)
  t1<-(1+zeta*((x-mu)/sigma))
  t2<--(t1^(-1/zeta))
  t3<-(-1-(1/zeta))*log(t1)
  value<-t2+t3+t4
  invisible(value)
}

log.lik<-function(param){
  mu<-param[1]
  zeta<-param[2]
  sigma<-param[3]
  if((sigma>0&&(zeta>0&&(prod(data>(mu-(sigma/zeta)))==1))||(zeta<0&&(prod(data<(mu-(sigma/zeta)))==1))||(zeta==0))){
    t4<--log(sigma)
    value<-sum(log.gev.dens(data,mu = mu,zeta = zeta,sigma = sigma))
  }
  else {value<--(Inf)}
  return(-value)
}

## Monte Carlo sample size chosen based on 5% tolerance for maximum of se for
## mu, zeta and sigma
n<-10000

## initializing the matrix of estimates
p.hat<-matrix(NA_real_,3,n)

## initializing the estimate vectors for all three parameters
mu.hat<-rep(NA_real_,n)
zeta.hat<-rep(NA_real_,n)
sigma.hat<-rep(NA_real_,n)

## initializing the standard errors, MCse and coverage vectors for all three parameters
mu.se<-rep(NA_real_,n)
zeta.se<-rep(NA_real_,n)
sigma.se<-rep(NA_real_,n)
mu.MCse<-rep(NA_real_,n)
zeta.MCse<-rep(NA_real_,n)
sigma.MCse<-rep(NA_real_,n)

mu.cov<-rep(NA_real_,n)
zeta.cov<-rep(NA_real_,n)
sigma.cov<-rep(NA_real_,n)

i<-1
## initializing the grid endpoints to choose intial value for 1st simulation
mu.lower<--10
mu.upper<-10
zeta.lower<--10
zeta.upper<-10
sigma.lower<-0.1
sigma.upper<-10.1


while(i<=n){
  ## simulating 1000 samples from given parameters 
  data<-rgev(n = 1000, xi = zeta, mu = mu, beta = sigma)
  ## evaluating the log likelihood on a grid to choose intial value
  mu.grid<-seq(mu.lower,mu.upper,length.out = 5)
  zeta.grid<-seq(zeta.lower,zeta.upper,length.out = 5)
  sigma.grid<-seq(sigma.lower,sigma.upper,length.out = 5)
  par.grid<-expand.grid(mu.grid,zeta.grid,sigma.grid)
  head(par.grid)
  initial.index<-which.min(apply(X = as.matrix(par.grid),MARGIN = 1,FUN = log.lik))
  initial.value<-par.grid[initial.index,]

  ## estimates
  param.hat<-optim(par = initial.value, fn = log.lik, method = "Nelder-Mead",hessian = TRUE)
  p.hat[,i]<-param.hat$par
  mu.hat[i]<-mean(p.hat[1,1:i])
  zeta.hat[i]<-mean(p.hat[2,1:i])
  sigma.hat[i]<-mean(p.hat[3,1:i])
  
  Info.mat<-param.hat$hessian
  Info.inv<-solve(Info.mat)
  mu.se[i]<-sqrt(Info.inv[1,1])
  mu.MCse[i]<-mu.se[i]/sqrt(i)
  mu.lower<-p.hat[1,i]-(mu.se[i]*qnorm(0.975,mean = 0, sd = 1))
  mu.upper<-p.hat[1,i]+(mu.se[i]*qnorm(0.975,mean = 0, sd = 1))
  mu.cov[i]<-(mu.lower<=mu)&(mu<=mu.upper)
  
  zeta.se[i]<-sqrt(Info.inv[2,2])
  zeta.MCse[i]<-zeta.se[i]/sqrt(i)
  zeta.lower<-p.hat[2,i]-(zeta.se[i]*qnorm(0.975,mean = 0, sd = 1))
  zeta.upper<-p.hat[2,i]+(zeta.se[i]*qnorm(0.975,mean = 0, sd = 1))
  zeta.cov[i]<-(zeta.lower<=zeta)&(zeta<=zeta.upper)
  
  sigma.se[i]<-sqrt(Info.inv[3,3])
  sigma.MCse[i]<-sigma.se[i]/sqrt(i)
  sigma.lower<-p.hat[3,i]-(sigma.se[i]*qnorm(0.975,mean = 0, sd = 1))
  sigma.upper<-p.hat[3,i]+(sigma.se[i]*qnorm(0.975,mean = 0, sd = 1))
  sigma.cov[i]<-(sigma.lower<=sigma)&(sigma<=sigma.upper)
  
  if (max(mu.MCse[i]/mu.hat[i],zeta.MCse[i]/zeta.hat[i],sigma.MCse[i]/sigma.hat[i])<=0.005){
    n.star<-i
    break
  }
  print(max(mu.se[i]/mu.hat[i],zeta.se[i]/zeta.hat[i],sigma.se[i]/sigma.hat[i]))
  i<-i+1
}

############################################
## Q3 (a)

n.star
mu.hat[n.star]
zeta.hat[n.star]
sigma.hat[n.star]

MSE.mu<-mean((p.hat[1,1:n.star]-mu)^2)
MSE.mu
MSE.mu.MCse<-sqrt(var((p.hat[1,1:n.star]-mu)^2)/n.star)
MSE.mu.MCse

MSE.zeta<-mean((p.hat[2,1:n.star]-zeta)^2)
MSE.zeta
MSE.zeta.MCse<-sqrt(var((p.hat[2,1:n.star]-zeta)^2)/n.star)
MSE.zeta.MCse

MSE.sigma<-mean((p.hat[3,1:n.star]-sigma)^2)
MSE.sigma
MSE.sigma.MCse<-sqrt(var((p.hat[3,1:n.star]-sigma)^2)/n.star)
MSE.sigma.MCse

############################################
## Q3 (b)

mu.covp<-sum(mu.cov[1:n.star])/n.star
mu.covp
mu.covp.MCse<-sqrt(var(mu.cov[1:n.star])/n.star)
mu.covp.MCse

zeta.covp<-sum(zeta.cov[1:n.star])/n.star
zeta.covp
zeta.covp.MCse<-sqrt(var(zeta.cov[1:n.star])/n.star)
zeta.covp.MCse

sigma.covp<-sum(sigma.cov[1:n.star])/n.star
sigma.covp
sigma.covp.MCse<-sqrt(var(sigma.cov[1:n.star])/n.star)
sigma.covp.MCse

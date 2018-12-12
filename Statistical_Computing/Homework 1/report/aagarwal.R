## Please note that each question is seperated by a full line of '#'
## and subparts are seperated by half lines of '#'. Also run each 
## question from the begining as a whole to get the correct output
## since variable names/function names may overlap for different
## questions.

#####################################################################
## Q1
## Reading the data
ys<-read.table("http://sites.stat.psu.edu/~mharan/515/hwdir/EMG1.dat")
data<-data.frame(ys)
names(data)<-c("X","Y")
n.rows<-nrow(data)

## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

## defining fixed parameter values
lambda<-0.4
beta.0<-5
sigma.i<-1

## defining the erfc function
erfc<-function(x){
  value<-(2/sqrt(pi))*(1-pnorm(x,mean=0,sd=(1/sqrt(2))))
  return(value)
}

## defining the posterior kernel on log scale
log.pr<-function(b){
  z<-(beta.0+(data$X)*b+(lambda*sigma.i^2)-data$Y)/(sqrt(2)*sigma.i)
  value<--(b^2/200)+(lambda*sum(data$X)*b)+sum(log(erfc(z)))
  return(value)
}

## Random Walk Metropolis function with inputs as the variance of proposal
## (tuning parameter) and current state of the MC
RWM<-function(variance,current.state){
  ## sampling beta1.star from normal proposal with mean as current state
  ## of beta_1 and specified tuning parameter
  beta1.star<-rnorm(1,mean=current.state,sd=sqrt(variance))
  ## defining the acceptance probability for sampled beta1.star on log
  ## scale
  accept.probab<-log.pr(beta1.star)-log.pr(current.state)
  ## sampling u from uniform(0,1) to check for acceptance
  u<-runif(1, min=0, max=1)
  ## initializing the indicator flag=0 to check if the sampled beta1.star
  ## will be accepted
  flag<-0
  ## if-else to define the next state of the chain based on acceptance
  ## probability
  if(log(u)<=accept.probab){
    flag<-1
    next.state<-beta1.star
  }
  else {next.state<-current.state}
  ## returning the next state and indicator if the sampled value was
  ## accepted
  return(c(next.state,flag))
}

## defining the function Expectation for running Random Walk Metropolis
## algortithm and to calculate Monte Carlo estimates and standard errors
## for beta_1
Expectation<-function(n,start,variance){
  ## Initializing the Markov chain for beta_1
  beta.1<-rep(NA_real_,n)
  ## Defining the initial value for the chain
  beta.1[1]<-start
  ## Initializing the accept count used to calculate acceptance rate
  ## of beta_1
  accept<-0
  ## loop for RWM updates
  for(i in 1:(n-1)){
    temp<-RWM(variance,beta.1[i])
    beta.1[i+1]<-temp[1]
    accept<-accept+temp[2]
  }
  ## samples obtained from the running the chain for given n
  samples<-beta.1
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix R for storing Monte Carlo estimates and MCMC
  ## standard errors for beta_1
  R<-matrix(NA_real_,ml,2)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    R[j,1]<-mean(beta.1[1:m[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,2]<-bm(beta.1[1:m[j]])$se
  }
  ## returning the sampled values of the chain, estimates & standard
  ## errors of beta_1 at varying number of iterations and acceptance rate
  ## of the chain
  return(list("samples"=samples,"R"=R,"acceptance.rate"=acceptance.rate))
}

## defining the function plotter to plot (1) the estimates and MCMC standard
## errors vs. number of iterations for beta_1 and (2) histograms and 
## smoothed posterior densities using MCMC samples for beta_1 along with
## autocorrelation of samples from the chain
plotter<-function(lists){
  par(mfrow=c(2,2), mar=c(2.8, 2.5, 2.8, 2.5))
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  lists.len<-length(lists)
  if (lists.len==1){
    plot(m,lists[[1]][[2]][,1],type="l",main="(a)",ylab=expression(paste("E(",beta[1],")")),ylim=c(0,7.5))
  }
  else if(lists.len>=2){
    ylim.upper<-0
    for (i in 1:lists.len){
      ylim.upper<-max(ylim.upper,lists[[i]][[2]][,1])
    }
    plot(m,lists[[1]][[2]][,1],type="l",main="(a)",ylab=expression(paste("E(",beta[1],")")),ylim=c(7.25,ylim.upper))
    for (i in 2:lists.len){
      E.beta1<-lists[[i]][[2]][,1]
      lines(m,E.beta1)
    }
  }
  E.beta1.se<-lists[[1]][[2]][,2]
  plot(m,E.beta1.se,type="l",main="(b)",ylab=expression(paste("MCMCse(E(",beta[1],"))")))
  samples.beta1<-lists[[1]][[1]]
  plot(density(samples.beta1), col="RED",main="(c)",lwd = 2,xlim=c(6,8.5),xlab="(c)")
  lines(density(samples.beta1[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(7.6,1.3, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  acf(samples.beta1,ylab=expression(paste("Autocorrelation")),main="(d)")
}

## defining the sample size
n<-30000

## defining the tuning parameter
var<-2

## Running the RWM algorithm for 3 times with different starting values
start<-c(7,7.2,6.6)
## creating a vector of empty lists
result<-vector("list", 3) 
for (j in 1:3){
  result[[j]]<-Expectation(n,start[j],var)
}

##################################
## Q1 (b)
E.beta1<-result[[1]][[2]][nrow(result[[1]][[2]]),1]
E.beta1.se<-result[[1]][[2]][nrow(result[[1]][[2]]),2]

##################################
## Q1 (c)
CI<-quantile(result[[1]][[1]],c(0.025,0.975))

##################################
## Q1 (d) and (e) 
## plotting estimates, MCMCse, smoothed density and autocorrelation
plotter(result)

## Effective sample size
ESS<-ess(result[[1]][[1]])

## Acceptance rate of the chain
Accept.rate<-result[[1]][[3]]

########################################################################
## Q2
## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

## defining the fixed parameter sigma_i
sigma.i<-1

## defining the dexpgauss density function on a log scale using the given
## link
dexpgauss<-function(x, mu=0, sigma=1, lambda=1, log=TRUE){
  l<-max(length(x), length(mu), length(sigma), length(lambda))
  x<-rep(x, times=ceiling(l/length(x)), length.out=l)
  mu<-rep(mu, times=ceiling(l/length(mu)), length.out=l)
  sigma<-rep(sigma, times=ceiling(l/length(sigma)), length.out=l)
  lambda<-rep(lambda, times=ceiling(l/length(lambda)), length.out=l)
  if (min(sigma)<=0){
    stop("Sigma must be greater than zero")
  }
  if (min(lambda)<=0){
    stop("Lambda must be greater than zero")
  }
  erfc<-pnorm((mu+lambda*sigma*sigma-x)/sigma, lower.tail=FALSE, log.p=log)
  if (log){
    result<-lambda/2*(2*mu+lambda*sigma*sigma-2*x)+Re(erfc)+log(lambda)
  }
  else{
    result<-exp(lambda/2*(2*mu+lambda*sigma*sigma-2*x))*Re(erfc)*lambda
  }
  result[is.nan(result)]<-0
  return(result)
}

## defining the erfc function
erfc<-function(x){
  value<-(2/sqrt(pi))*(1-pnorm(x,mean=0,sd=(1/sqrt(2))))
  return(value)
}

## defining the log of beta_0 & beta_1 joint posterior density
## conditional on everything else
log.b0b1<-function(b0,b1,lamb,data){
  z<-(b0+(data$X)*b1+(lamb*sigma.i^2)-data$Y)/(sqrt(2)*sigma.i)
  n.rows<-nrow(data)
  value<-lamb*b0*n.rows+lamb*b1*sum(data$X)+sum(log(erfc(z)))
  return(value)
}

## Random Walk Metropolis function for beta_0 and beta_1 with inputs 
## as the covariance matrix of (tuning parameter) bivariate normal 
## proposal with mean as vector of current states of beta_0 and beta_1
## Markov chain
library(MASS)
RWM.b0b1<-function(cov,b0.cs,b1.cs,lamb.cs,data){
  ## sampling b0.star from normal proposal with mean as current state of 
  ## beta_0 and specified tuning parameter
  b0b1<-mvrnorm(1,c(b0.cs,b1.cs),cov)
  b0.star<-b0b1[1]
  b1.star<-b0b1[2]
  ## defining the acceptance probability for sampled b0.star on log
  ## scale
  accept.probab<-log.b0b1(b0.star,b1.star,lamb.cs,data)-log.b0b1(b0.cs,b1.cs,lamb.cs,data)
  ## sampling u0 from uniform(0,1) to check for acceptance
  u0<-runif(1, min=0, max=1)
  ## initializing the indicator flag=0 to check if the sampled 
  ## b0.star will be accepted
  flag<-0
  ## if-else to define the next state of the chain based on acceptance
  ## probability
  if(log(u0)<=accept.probab){
    flag<-1
    b0.ns<-b0.star
    b1.ns<-b1.star
  }
  else {
    b0.ns<-b0.cs
    b1.ns<-b1.cs
  }
  ## returning the next state and indicator if the sampled value was
  ## accepted
  return(c(b0.ns,b1.ns,flag))
}

## Metropolis-Hastings update for lambda using a gamma proposal with
## parametrization as mean being current state and variance being 
## tuning parameter
RWM.lamb<-function(tune.lamb,b0.cs,b1.cs,lamb.cs,data){
  lamb.star<-rgamma(1,shape=(lamb.cs^2)/tune.lamb,scale=tune.lamb/lamb.cs)
  ul<-runif(1, min=0, max=1)
  n.rows<-nrow(data)
  log.lik<-rep(NA_real_,n.rows)
  for (j in 1:n.rows){
    num<-dexpgauss(data$Y[j], mu=(b0.cs+b1.cs*data$X[j]), sigma=1, lambda=lamb.star, log=TRUE)
    denom<-dexpgauss(data$Y[j], mu=(b0.cs+b1.cs*data$X[j]), sigma=1, lambda=lamb.cs, log=TRUE)
    log.lik[j]<-(num-denom)
  }
  num.p<--0.99*(log(lamb.star))-(lamb.star/100)
  denom.p<--0.99*(log(lamb.cs))-(lamb.cs/100)
  num.q<-dgamma(lamb.cs,shape=(lamb.star^2)/tune.lamb,scale=tune.lamb/lamb.star,log=TRUE)
  denom.q<-dgamma(lamb.star,shape=(lamb.cs^2)/tune.lamb,scale=tune.lamb/lamb.cs,log=TRUE)
  accept.probab<-sum(log.lik)+(num.p-denom.p)+(num.q-denom.q)
  ##accept.probab<-((lamb.star/lamb.cs)^(n.rows-0.99))*(exp(log.lamb(b0.cs,b1.cs,lamb.star)-log.lamb(b0.cs,b1.cs,lamb.cs)))
  flag<-0
  if(log(ul)<=accept.probab){
    flag<-1
    lamb.ns<-lamb.star
  }
  else {lamb.ns<-lamb.cs}
  return(c(lamb.ns,flag))
}

## defining the function Expectations for running variable one at a time
## Metropolis Hastings algortithm and to calculate Monte Carlo estimates
## and standard errors for beta_0, beta_1 and lambda
Expectations<-function(n,start,cov,scale.lamb,data){
  ## Initializing the 3 chains for beta_0, beta_1 and lambda
  beta.0<-rep(NA_real_,n)
  beta.1<-rep(NA_real_,n)
  lambda<-rep(NA_real_,n)
  ## Defining the initial values for the 3 chains
  beta.0[1]<-start[1]
  beta.1[1]<-start[2]
  lambda[1]<-start[3]
  ## Initializing the accept counts used to calculate acceptance rate
  ## of beta_0, beta_1 and lambda
  accept.b0b1<-0
  accept.lamb<-0
  ## loop for M-H updates
  for(i in 1:(n-1)){
    v.b0b1<-RWM.b0b1(cov,beta.0[i],beta.1[i],lambda[i],data)
    v.lamb<-RWM.lamb(scale.lamb,v.b0b1[1],v.b0b1[2],lambda[i],data)
    beta.0[i+1]<-v.b0b1[1]
    beta.1[i+1]<-v.b0b1[2]
    lambda[i+1]<-v.lamb[1]
    accept.b0b1<-accept.b0b1+v.b0b1[3]
    accept.lamb<-accept.lamb+v.lamb[2]
  }
  ## samples obtained from the running the chains for given n
  samples<-cbind(beta.0,beta.1,lambda)
  ## calculating the acceptance rates for each of the chains
  accept.rate.b0b1<-accept.b0b1/n
  accept.rate.lamb<-accept.lamb/n
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix R for storing Monte Carlo estimates and MCMC
  ## standard errors for all parameters 
  R<-matrix(NA_real_,ml,6)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    R[j,1]<-mean(beta.0[1:m[j]])
    R[j,2]<-mean(beta.1[1:m[j]])
    R[j,3]<-mean(lambda[1:m[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,4]<-bm(beta.0[1:m[j]])$se
    R[j,5]<-bm(beta.1[1:m[j]])$se
    R[j,6]<-bm(lambda[1:m[j]])$se
  }
  ## returning the sampled values of the 3 chains, estimates & standard
  ## errors of beta_0, beta_1 and lambda at varying number of iterations
  ## and acceptance rates of the 3 chains
  return(list("samples"=samples,"R"=R,"acceptance.rates"=c(accept.rate.b0b1,accept.rate.lamb)))
}

## defining the function plotter1 to plot the estimates and MCMC standard
## errors vs. number of iterations for beta_0, beta_1 and lambda for 
## 1 chain
plotter1<-function(list1){
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  E.beta0<-list1[[2]][,1]
  E.beta1<-list1[[2]][,2]
  E.lambda<-list1[[2]][,3]
  par(mfcol=c(3,2), mar=c(1, 4, 1, 1))
  plot(m,E.beta0,type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")))
  plot(m,E.beta1,type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[1],")")))
  plot(m,E.lambda,type="l",xlab="Number of iterations",ylab=expression(paste("E(",lambda,")")))
  E.beta0.se<-list1[[2]][,4]
  E.beta1.se<-list1[[2]][,5]
  E.lambda.se<-list1[[2]][,6]
  plot(m,E.beta0.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[0],"))")))
  plot(m,E.beta1.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[1],"))")))
  plot(m,E.lambda.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",lambda,"))")))
}

## defining the function plot.dens.acf1 for plotting histograms and 
## smoothed posterior densities using MCMC samples for beta_0, beta_1 and
## lambda along with  autocorrelation of samples for 1 chain
plot.dens.acf1<-function(list1){
  samples.beta0<-list1[[1]][,1]
  samples.beta1<-list1[[1]][,2]
  samples.lambda<-list1[[1]][,3]
  par(mfcol=c(3,2), mar=c(2.7, 2.5, 2.7, 2.5))
  hist(list1[[1]][,1], freq=FALSE, breaks=50,xlim=c(1.8,2.8),ylim=c(0,8),main=expression(paste("(a) Smoothed density plot for ",beta[0])))
  lines(density(list1[[1]][,1]), col="RED",  lwd = 2)
  hist(list1[[1]][,2], freq=FALSE, breaks=50,xlim=c(2.8,4.1),ylim=c(0,4),main=expression(paste("(b) Smoothed density plot for ",beta[1])))
  lines(density(list1[[1]][,2]), col="RED",  lwd = 2)
  hist(list1[[1]][,3], freq=FALSE, breaks=50,xlim=c(0.6,1),ylim=c(0,15),main=expression(paste("(c) Smoothed density plot for ",lambda)))
  lines(density(list1[[1]][,3]), col="RED",  lwd = 2)
  acf(samples.beta0,ylab=expression(paste("Autocorrelation")),main=expression(paste("(d) Autocorrelation of ",beta[0]," samples")),xlab="(d)")
  acf(samples.beta1,ylab=expression(paste("Autocorrelation")),main=expression(paste("(e) Autocorrelation of ",beta[1]," samples")),xlab="(e)")
  acf(samples.lambda,ylab=expression(paste("Autocorrelation")),,main=expression(paste("(f) Autocorrelation of ",lambda," samples")),xlab="(f)")
}

## For diagnostics: defining the function plotter to plot the estimates
## and MCMC standard errors vs. number of iterations for beta_0, beta_1
## and lambda
plotter<-function(lists){
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  par(mfcol=c(3,2), mar=c(1.6, 2, 1.6, 2))
  lists.len<-length(lists)
  for (j in 1:3){
    if (lists.len==1){
      plot(m,lists[[1]][[2]][,j],type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")))
    }
    else if(lists.len>=2){
      ylim.lower<-100
      ylim.upper<-0
      for (i in 1:lists.len){
        ylim.lower<-min(ylim.lower,lists[[i]][[2]][,j])
        ylim.upper<-max(ylim.upper,lists[[i]][[2]][,j])
      }
      plot(m,lists[[1]][[2]][,j],type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")),ylim=c(ylim.lower,ylim.upper+0.03))
      for (i in 2:lists.len){
        E.beta0<-lists[[i]][[2]][,j]
        lines(m,E.beta0)
      }
    }
  }
  E.beta0.se<-lists[[1]][[2]][,4]
  E.beta1.se<-lists[[1]][[2]][,5]
  E.lambda.se<-lists[[1]][[2]][,6]
  plot(m,E.beta0.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[0],"))")))
  plot(m,E.beta1.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[1],"))")))
  plot(m,E.lambda.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",lambda,"))")))
}

## For diagnostics: defining the function plot.dens.acf by plotting
## histograms and smoothed posterior densities using MCMC samples for
## beta_0, beta_1 and lambda after sample size n/2 and n 
## along with  autocorrelation of samples 
plot.dens.acf<-function(list1){
  samples.beta0<-list1[[1]][,1]
  samples.beta1<-list1[[1]][,2]
  samples.lambda<-list1[[1]][,3]
  par(mfcol=c(3,2), mar=c(2.5, 2.3, 2.5, 2.3))
  plot(density(samples.beta0), col="RED", lwd = 2,xlim=c(1.8,2.8),ylim=c(0,4),main=expression(paste("(a) Smoothed density plot for ",beta[0])))
  lines(density(samples.beta0[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(2.55,3.5, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  plot(density(samples.beta1), col="RED", lwd = 2,xlim=c(2.8,4.1),ylim=c(0,2.5),main=expression(paste("(a) Smoothed density plot for ",beta[1])))
  lines(density(samples.beta1[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(3.75,2.2, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  plot(density(samples.lambda), col="RED", lwd = 2,xlim=c(0.6,1),ylim=c(0,10),main=expression(paste("(a) Smoothed density plot for ",lambda)))
  lines(density(samples.lambda[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(0.89,8.5, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  acf(samples.beta0,ylab=expression(paste("Autocorrelation")),main=expression(paste("(d) Autocorrelation of ",beta[0]," samples")),xlab="(d)")
  acf(samples.beta1,ylab=expression(paste("Autocorrelation")),main=expression(paste("(e) Autocorrelation of ",beta[1]," samples")),xlab="(e)")
  acf(samples.lambda,ylab=expression(paste("Autocorrelation")),,main=expression(paste("(f) Autocorrelation of ",lambda," samples")),xlab="(f)")
}

## Reading the data
ys=read.table("http://sites.stat.psu.edu/~mharan/515/hwdir/EMG2.dat")
data2<-data.frame(ys)
names(data2)<-c("X","Y")

####################################
## checking if variable one at a time M-H algorithm works for arbitrary
## starting values and arbitrary tuning parameters
## defining the sample size
n<-1000

## defining initial values of the 3 chains
start<-c(1,1,1)

## defining the tuning parameters (now)
s11<-0.064
s22<-0.112
s12<-0.023
cov.b0b1<-matrix(c(s11,s12,s12,s22),2,2)
var.lambda<-0.0079

## Running variable one at a time M-H algorithm
list1<-Expectations(n,start,cov.b0b1,var.lambda,data2)

## plots of estimates and MCMC standard errors of beta_0, beta_1 and
## lambda vs. number of iterations
plotter1(list1)

## plotting histograms with smoothed posterior densities using MCMC 
## samples for beta_0, beta_1 and lambda along with the autocrrelation
## plots
plot.dens.acf1(list1)

## Acceptance rates of the 3 chains
list1[[3]]

## Effective sample size
ess(list1[[1]])

##################################
## NOTE THAT THIS IS AN INTERMEDIATE TIME CONSUMING STEP AND IS NOT 
## REQUIRED LATER ONCE BEST SET OF TUNING PARAMETERS IS KNOWN 

## Running the MH-algorithm over a grid of tuning parameters to choose
## the best tuning parameters. 

## defining the sample size
n<-5000
## defining a sequence for number of iterations
m<-seq(100,n,by=100)
## defining initial values of the 3 chains
start<-c(0.1,0.1,0.1)
## defining the tuning parameters
s11<-seq(1,5,by=0.2)
s22<-seq(1,5,by=0.2)
s12<-seq(0.1,0.36,by=0.02)
var.lambda<-0.3
grid<-expand.grid(s11,s22,s12,var.lambda)
mhl<-nrow(grid)

## creating a vector of empty lists
result<-vector("list", mhl) 
## loop for running variable one at a time M-H algorithm for different
## combinations of tuning parameters
for (i in 1:mhl){
  cov.b0b1<-matrix(c(grid[i,1],grid[i,3],grid[i,3],grid[i,2]),2,2)
  result[[i]]<-Expectations(n,start,cov.b0b1,grid[i,4],data2)
}

## checking plots for different sets of tuning parameters and identifying
## the best one with highest ESS and lowest autocorrelation
plotter(result)
plot.dens.acf(result[[12]])

## Effective sample size
temp<-rep(NA_real_,mhl)
for (i in 1:mhl){
  temp[i]<-ess(result[[i]][[1]][,2])
}
s12[which.max(temp)]
## This does not give conclusive results. Tuning done in part (d)

##################################
## From Q2 part (d) (later) it was found that the best set of tuning
## parameters are s11=0.0187,s22=0.0447,s12=-0.0226, var.lamb=0.0036

## Running the M-H algorithm for n=50000, by selecting the best set of
## tuning parameters
## defining the sample size
n<-50000

## defining initial values of the final chain
start<-c(2.3,3.5,0.8)

## defining the tuning parameters
s11<-0.0187
s22<-0.0447
s12<--0.0226
cov.b0b1<-matrix(c(s11,s12,s12,s22),2,2)
var.lambda<-0.0036

## Running variable one at a time M-H algorithm
list.final<-Expectations(n,start,cov.b0b1,scale.lambda,data2)

## plots of estimates and MCMC standard errors of beta_0, beta_1 and
## lambda vs. number of iterations
plotter1(list.final)

## plotting histograms with smoothed posterior densities using MCMC 
## samples for beta_0, beta_1 and lambda along with the autocrrelation
## plots
plot.dens.acf1(list.final)

## Acceptance rates of the 3 chains
list.final[[3]]

## Effective sample size
ess(list.final[[1]])

##################################
## Q2 (b)
## Posterior means and standard errors
E.beta0<-list.final[[2]][nrow(list.final[[2]]),1]
E.beta0
E.beta1<-list.final[[2]][nrow(list.final[[2]]),2]
E.beta1
E.lambda<-list.final[[2]][nrow(list.final[[2]]),3]
E.lambda
E.beta0.se<-list.final[[2]][nrow(list.final[[2]]),4]
E.beta0.se
E.beta1.se<-list.final[[2]][nrow(list.final[[2]]),5]
E.beta1.se
E.lambda.se<-list.final[[2]][nrow(list.final[[2]]),6]
E.lambda.se

## 95% Credible intervals
CI.beta0<-quantile(list.final[[1]][,1],c(0.025,0.975))
CI.beta0
CI.beta1<-quantile(list.final[[1]][,2],c(0.025,0.975))
CI.beta1
CI.lambda<-quantile(list.final[[1]][,3],c(0.025,0.975))
CI.lambda

##################################
## Q2 (c)
corr.b0b1<-cor(list.final[[1]][,1],list.final[[1]][,2])
corr.b0b1

##################################
## Q2 (d) 
## Diagnostics plots of estimates and MCMC standard errors of beta_0, beta_1 and
## lambda vs. number of iterations

## defining the sample size
n<-10000

## defining the number of chains
mhl<-3

## defining initial values of the 3 chains
start<-matrix(c(2.2,3.2,0.7,2.3,3.5,0.8,2.6,3.8,0.9),3,mhl)

## defining the tuning parameters
s11<-0.0187
s22<-0.0447
s12<--0.0226
cov.b0b1<-matrix(c(s11,s12,s12,s22),2,2)
var.lambda<-0.0036

## creating a vector of empty lists
result.final<-vector("list", mhl) 
## loop for running variable one at a time M-H algorithm for different
## combinations of tuning parameters
for (i in 1:mhl){
  result.final[[i]]<-Expectations(n,start[,i],cov.b0b1,var.lambda,data2)
}
plotter(result.final)

## plotting histograms with smoothed posterior densities using MCMC 
## samples for beta_0, beta_1 and lambda along with the autocrrelation
## plots
plot.dens.acf(result.final[[1]])

## checking the ESS for first chain
ess(result.final[[1]][[1]][,1])
ess(result.final[[1]][[1]][,2])
ess(result.final[[1]][[1]][,3])

## Tuning the parameters based on sample variances of beta_0,beta_1 and
## lambda and covariance of beta_0 and beta_1
mean(var(result.final[[1]][[1]][,1]),var(result.final[[2]][[1]][,1]),var(result.final[[3]][[1]][,1]))
mean(var(result.final[[1]][[1]][,2]),var(result.final[[2]][[1]][,2]),var(result.final[[3]][[1]][,2]))
mean(var(result.final[[1]][[1]][,3]),var(result.final[[2]][[1]][,3]),var(result.final[[3]][[1]][,3]))

c1<-cov(result.final[[1]][[1]][,1],result.final[[1]][[1]][,2])
c2<-cov(result.final[[2]][[1]][,1],result.final[[2]][[1]][,2])
c3<-cov(result.final[[3]][[1]][,1],result.final[[3]][[1]][,2])
mean(c1,c2,c3)

#######################################################################
## Q3
## Please run all the functions defined in Q2

## Defining the function plotter1 again (to accomodate for different
## x and y limits) to plot the estimates and MCMC standard errors vs.
## number of iterations for beta_0, beta_1 and lambda for 1 chain
plotter1<-function(list1){
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  E.beta0<-list1[[2]][,1]
  E.beta1<-list1[[2]][,2]
  E.lambda<-list1[[2]][,3]
  par(mfcol=c(3,2), mar=c(1, 4, 1, 1))
  plot(m,E.beta0,type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")))
  plot(m,E.beta1,type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[1],")")))
  plot(m,E.lambda,type="l",xlab="Number of iterations",ylab=expression(paste("E(",lambda,")")))
  E.beta0.se<-list1[[2]][,4]
  E.beta1.se<-list1[[2]][,5]
  E.lambda.se<-list1[[2]][,6]
  plot(m,E.beta0.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[0],"))")))
  plot(m,E.beta1.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[1],"))")))
  plot(m,E.lambda.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",lambda,"))")))
}

## defining the function plot.dens.acf1 again (to accomodate for different
## x and y limits) for plotting histograms and smoothed posterior
## densities using MCMC samples for beta_0, beta_1 and lambda along with
## autocorrelation of samples for 1 chain
plot.dens.acf1<-function(list1){
  samples.beta0<-list1[[1]][,1]
  samples.beta1<-list1[[1]][,2]
  samples.lambda<-list1[[1]][,3]
  par(mfcol=c(3,2), mar=c(2.7, 2.5, 2.7, 2.5))
  hist(list1[[1]][,1], freq=FALSE, breaks=50,xlim=c(-0.5,1),ylim=c(0,3),main=expression(paste("(a) Smoothed density plot for ",beta[0])))
  lines(density(list1[[1]][,1]), col="RED",  lwd = 2)
  hist(list1[[1]][,2], freq=FALSE, breaks=50,xlim=c(1.5,3.5),ylim=c(0,2),main=expression(paste("(b) Smoothed density plot for ",beta[1])))
  lines(density(list1[[1]][,2]), col="RED",  lwd = 2)
  hist(list1[[1]][,3], freq=FALSE, breaks=50,xlim=c(0.12,0.2),ylim=c(0,90),main=expression(paste("(c) Smoothed density plot for ",lambda)))
  lines(density(list1[[1]][,3]), col="RED",  lwd = 2)
  acf(samples.beta0,ylab=expression(paste("Autocorrelation")),main=expression(paste("(d) Autocorrelation of ",beta[0]," samples")),xlab="(d)")
  acf(samples.beta1,ylab=expression(paste("Autocorrelation")),main=expression(paste("(e) Autocorrelation of ",beta[1]," samples")),xlab="(e)")
  acf(samples.lambda,ylab=expression(paste("Autocorrelation")),,main=expression(paste("(f) Autocorrelation of ",lambda," samples")),xlab="(f)")
}

## For diagnostics: defining the function plotter again (to accomodate for different
## x and y limits) to plot the estimates and MCMC standard errors vs. 
## number of iterations for beta_0, beta_1 and lambda
plotter<-function(lists){
  ## defining a sequence for number of iterations
  m<-seq(100,n,by=100)
  par(mfcol=c(3,2), mar=c(1.6, 2, 1.6, 2))
  lists.len<-length(lists)
  for (j in 1:3){
    if (lists.len==1){
      plot(m,lists[[1]][[2]][,j],type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")))
    }
    else if(lists.len>=2){
      ylim.lower<-100
      ylim.upper<-0
      for (i in 1:lists.len){
        ylim.lower<-min(ylim.lower,lists[[i]][[2]][,j])
        ylim.upper<-max(ylim.upper,lists[[i]][[2]][,j])
      }
      plot(m,lists[[1]][[2]][,j],type="l",xlab="Number of iterations",ylab=expression(paste("E(",beta[0],")")),ylim=c(ylim.lower,ylim.upper+0.03))
      for (i in 2:lists.len){
        E.beta0<-lists[[i]][[2]][,j]
        lines(m,E.beta0)
      }
    }
  }
  E.beta0.se<-lists[[1]][[2]][,4]
  E.beta1.se<-lists[[1]][[2]][,5]
  E.lambda.se<-lists[[1]][[2]][,6]
  plot(m,E.beta0.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[0],"))")))
  plot(m,E.beta1.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",beta[1],"))")))
  plot(m,E.lambda.se,type="l",xlab="Number of iterations",ylab=expression(paste("MCMCse(E(",lambda,"))")))
}

## For diagnostics: defining the function plot.dens.acf again (to 
## accomodate for different x and y limits)by plotting histograms and
## smoothed posterior densities using MCMC samples for beta_0, beta_1
## and lambda after sample size n/2 and n along with  autocorrelation
## of samples 
plot.dens.acf<-function(list1){
  samples.beta0<-list1[[1]][,1]
  samples.beta1<-list1[[1]][,2]
  samples.lambda<-list1[[1]][,3]
  par(mfcol=c(3,2), mar=c(2.5, 2.3, 2.5, 2.3))
  plot(density(samples.beta0), col="RED", lwd = 2,xlim=c(-0.5,1),ylim=c(0,3),main=expression(paste("(a) Smoothed density plot for ",beta[0])))
  lines(density(samples.beta0[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(0.6,2.7, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  plot(density(samples.beta1), col="RED", lwd = 2,xlim=c(1.5,3.5),ylim=c(0,2),main=expression(paste("(a) Smoothed density plot for ",beta[1])))
  lines(density(samples.beta1[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(2.9,1.8, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  plot(density(samples.lambda), col="RED", lwd = 2,xlim=c(0.12,0.2),ylim=c(0,90),main=expression(paste("(a) Smoothed density plot for ",lambda)))
  lines(density(samples.lambda[1:(n/2)]),col="BLUE",  lwd = 2)
  legend(0.175,85, bty="n",c("After n/2","After n"), cex=0.65,lty=c(1,1), lwd=c(1,1),col=c("blue","red"))
  
  acf(samples.beta0,ylab=expression(paste("Autocorrelation")),main=expression(paste("(d) Autocorrelation of ",beta[0]," samples")),xlab="(d)")
  acf(samples.beta1,ylab=expression(paste("Autocorrelation")),main=expression(paste("(e) Autocorrelation of ",beta[1]," samples")),xlab="(e)")
  acf(samples.lambda,ylab=expression(paste("Autocorrelation")),,main=expression(paste("(f) Autocorrelation of ",lambda," samples")),xlab="(f)")
}

## Reading the data
ys=read.table("http://sites.stat.psu.edu/~mharan/515/hwdir/EMG3.dat")
data3<-data.frame(ys)
names(data3)<-c("X","Y")

#######################################
## Tuning the parameters using heuristics
## defining the sample size
n<-10000

## defining initial values of the 3 chains
start<-matrix(c(0,2,0.15,0.2,2.5,0.17,0.5,3,0.16),3,mhl)

## defining the tuning parameters
s11<-0.0257
s22<-0.0765
s12<--0.0371
cov.b0b1<-matrix(c(s11,s12,s12,s22),2,2)
var.lambda<-0.0000299

## Number of chains
mhl<-3

## creating a vector of empty lists
result.final<-vector("list", mhl) 

## loop for running variable one at a time M-H algorithm for different
## combinations of tuning parameters
for (i in 1:mhl){
  result.final[[i]]<-Expectations(n,start[,i],cov.b0b1,var.lambda,data3)
}
plotter(result.final)

## plotting histograms with smoothed posterior densities using MCMC 
## samples for beta_0, beta_1 and lambda along with the autocrrelation
## plots
plot.dens.acf(result.final[[2]])

## checking the ESS for first chain
ess(result.final[[1]][[1]][,1])
ess(result.final[[1]][[1]][,2])
ess(result.final[[1]][[1]][,3])

## Tuning the parameters based on sample variances of beta_0,beta_1 and
## lambda and covariance of beta_0 and beta_1
mean(var(result.final[[1]][[1]][,1]),var(result.final[[2]][[1]][,1]),var(result.final[[3]][[1]][,1]))
mean(var(result.final[[1]][[1]][,2]),var(result.final[[2]][[1]][,2]),var(result.final[[3]][[1]][,2]))
mean(var(result.final[[1]][[1]][,3]),var(result.final[[2]][[1]][,3]),var(result.final[[3]][[1]][,3]))

c1<-cov(result.final[[1]][[1]][,1],result.final[[1]][[1]][,2])
c2<-cov(result.final[[2]][[1]][,1],result.final[[2]][[1]][,2])
c3<-cov(result.final[[3]][[1]][,1],result.final[[3]][[1]][,2])
mean(c1,c2,c3)

##################################
## From above it was found that the best set of tuning parameters are 
## s11=0.0257,s22=0.0765,s12=-0.0371, var.lamb=0.0000299

## Running the M-H algorithm for n=50000, by selecting the best set of
## tuning parameters
## defining the sample size
n<-10000

## defining initial values of the final chain
start<-c(0.2,0.15,0.2)

## defining the tuning parameters
s11<-0.0257
s22<-0.0765
s12<--0.0371
cov.b0b1<-matrix(c(s11,s12,s12,s22),2,2)
var.lambda<-0.0000299

## Running variable one at a time M-H algorithm
list.final<-Expectations(n,start,cov.b0b1,var.lambda,data3)

##################################
## Q3 (a)
## Posterior means and standard errors
E.beta0<-list.final[[2]][nrow(list.final[[2]]),1]
E.beta0
E.beta1<-list.final[[2]][nrow(list.final[[2]]),2]
E.beta1
E.lambda<-list.final[[2]][nrow(list.final[[2]]),3]
E.lambda

E.beta0.se<-list.final[[2]][nrow(list.final[[2]]),4]
E.beta0.se
E.beta1.se<-list.final[[2]][nrow(list.final[[2]]),5]
E.beta1.se
E.lambda.se<-list.final[[2]][nrow(list.final[[2]]),6]
E.lambda.se

## 95% credible intervals
beta0.CI<-quantile(list.final[[1]][,1],c(0.025,0.975))
beta0.CI
beta1.CI<-quantile(list.final[[1]][,2],c(0.025,0.975))
beta1.CI
lambda.CI<-quantile(list.final[[1]][,3],c(0.025,0.975))
lambda.CI

###################################
## Q3 (b)
## plots of estimates and MCMC standard errors of beta_0, beta_1 and
## lambda vs. number of iterations
plotter1(list.final)

## plotting histograms with smoothed posterior densities using MCMC 
## samples for beta_0, beta_1 and lambda along with the autocrrelation
## plots
plot.dens.acf1(list.final)

## Acceptance rates of the 3 chains
list.final[[3]]

## Effective sample sizes
ess(list.final[[1]][,1])
ess(list.final[[1]][,2])
ess(list.final[[1]][,3])
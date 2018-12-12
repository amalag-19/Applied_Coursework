## Please note that each question is seperated by a full line of '#'
## and subparts are seperated by half lines of '#'. Also run each 
## question from the begining as a whole to get the correct output
## since variable names/function names may overlap for different
## questions.

## Please install the necessary packages mentioned in the libraries beforehand
#####################################################################
## Q1

## loading the required libraries
library(ggplot2)

## reading the data
df<-read.table("http://www.stat.psu.edu/~mharan/540/hwdir/bulbsA.dat",header = FALSE,sep = "")
df<-t(df)

## observed data
Y<-as.vector(df[,1])

## defining the function to calculate EM estimate
EM<-function(start,tol,data){
  ## initializing the parameter vector
  theta<-vector()
  ## defining the initial value
  theta[1]<-start
  m<-length(which(data==100))
  s<-sum(data)
  diff<-Inf
  i<-1
  while (diff>tol){
    theta[i+1]<-(s+(m*theta[i]))/length(data)
    diff<-abs(theta[i+1]-theta[i])
    i<-i+1
  }
  return(theta)
}

## defining the EM number of iterations
n<-100
## defining the starting value
start<-1
## defining the tolerance
tol<-10^(-6)

theta<-EM(start = start,tol = tol,data = Y)

f<-data.frame("iterations"=1:length(theta),theta)
p<-ggplot(data=f)
p+geom_line(mapping = aes(x=iterations,y=theta))+labs(x="iterations",y=expression(theta))

#####################################################################################
## Q2

## defining the EM number of iterations
n<-50
## defining the starting value
start<-50
## defining the tolerance for EM
tol<-10^(-6)

theta<-EM(start = start,tol = tol,data = Y)

f<-data.frame("iterations"=1:length(theta),theta)
p<-ggplot(data=f)
p+geom_line(mapping = aes(x=iterations,y=theta))+labs(x="iterations",y=expression(theta))

## defining the truncated exponential function to calculated the MLE for uncensored
## data (<100) used as an initial value
trunc.exp<-function(mu){
  temp<-dexp(x = Y[which(Y<100)],rate  = 1/mu,log = TRUE)
  value<-sum(temp)-length(which(Y<100))*log(pexp(q = 100,rate = 1/mu))
  return(-value)
}
## calculating the MLE to be used as the starting value in the EM algorithm
start<-optim(par = 400,fn = trunc.exp,method = "L-BFGS-B",lower = 0, upper = Inf)$par
## defining the EM number of iterations
n<-50

## defining the tolerance for EM
tol<-10^(-6)

## running the EM algorithm again for this new starting value
theta.st<-EM(start = start,tol = tol,data = Y)

## plotting
f<-data.frame("iterations"=1:length(theta.st),theta.st)
p<-ggplot(data=f)
p+geom_line(mapping = aes(x=iterations,y=theta.st))+labs(x="iterations",y=expression(theta))

#####################################################################################
## Q3
theta[length(theta)]

#####################################################################################
## Q4

## initializing the vector of parameter estimates for bootstrap samples
theta.boot<-vector()

## defining the tolerance for standard error estimate
tol<-10^(-6)

## initializing the vector of standard error estimate for the parameter for bootstrap
## samples
boot.se.theta<-vector()
boot.se.theta[1:99]<--Inf
boot.se.theta[100]<-Inf

## initializing the vector of quantiles for calculating bootstrap t confidence
## intervals
r<-vector()

## bootstrap loop
i<-1

while (abs(boot.se.theta[length(boot.se.theta)]-boot.se.theta[length(boot.se.theta)-99])>tol){
  Y.boot<-sample(x = Y,size = length(Y),replace = TRUE)
  ## defining the EM number of iterations
  n<-30
  ## defining the starting value
  start<-451
  ## defining the tolerance for EM
  tol<-10^(-6)
  ## caculating parameter estimates for bootstrap samples via EM
  temp<-EM(start = start,tol = tol,data = Y.boot)
  theta.boot[i]<-temp[length(temp)]
  if (i>1){
    boot.se.theta[i-1]<-sqrt(var(theta.boot[1:i]))
  }
  ## bootstrap sample size
  B<-length(theta.boot)
  ## quantile value for the bootstrap sample
  r[i]<-(theta.boot[i]-theta[n])/(boot.se.theta[length(boot.se.theta)])
  i<-i+1
}
boot.se<-boot.se.theta[length(boot.se.theta)]
B

## plotting the estimate of standard error with bootstrap sample size
g<-data.frame("iterations"=1:(B-1),boot.se.theta)
q<-ggplot(data=g)
q+geom_line(mapping = aes(x=iterations,y=boot.se.theta))+labs(x="Bootstrap sample size",y=expression(paste("Bootstrap estimate for se(",theta,")")))

#####################################################################################
## Q5
## calculating the t bootstrap CI

rl<-quantile(x = r,probs = 0.025)
ru<-quantile(x = r,probs = 0.975)

lower<-theta[n]-abs(ru)*(boot.se)
upper<-theta[n]+abs(rl)*(boot.se)

CI<-c(lower,upper)

##Poisson Regression
n<-50 ## sample size
nv<-seq(from=100,to=1000,by=100) ## vector of different sample sizes
n.iter<-rep(NA,length(nv)) ## Placeholder vector for number of 
## iterations required for different sample sizes
for (k in 1:length(nv)){
  ## Generating the data: Design matrix from normal distribution 
  ## for a particular sample size
  n<-nv[k]
  x1<-rnorm(n,mean=1,sd=0.5)
  x2<-rnorm(n,mean=1,sd=0.5)
  X<-cbind(1,x1,x2) ## Design matrix
  p<-3 ## Number of predictors
  beta<-c(1,2,3) ## true parameter vector
  eta<-X%*%beta ## Linear predictor
  lambda<-exp(eta) ## Poisson parameter vector
  Y<-rep(NA,n) ## placeholder vector for response data
  
  ## response sampled from Poisson according to lambda vector
  for (i in 1:n){
    Y[i]<-rpois(1,lambda[i])
  }
  
  ## Fisher Scoring Algorithm
  m<-100 ## size of the iteration matrix
  
  ## score function
  s <- function(b){
    score<-t(X)%*%(Y-exp(X%*%b))
    return(score)
  }
  
  ## Fisher matrix
  f<-function(b){
    w<-exp(X%*%b)
    W<-diag(c(w))
    fisher<-t(X)%*%W%*%X
    return(fisher)
  }
  
  ## placeholder matrix for different iterates of parameter vector 
  beta.iterates<-matrix(0,nrow=p, ncol=m+1)
  ## Initial beta vector randomly chosen
  beta.iterates[,1]<-c(3,2,1)
  ## threshold for stopping the iterations
  epsilon<-10^(-2)
  ## placeholder vector for the difference of iterated beta vector 
  ## and true beta vector for checking convergence criteria
  d<-rep(NA,p)
  ## indicator for stopping the algorithm when convergence criteria is met
  stop<-0
  ## iterations
  for (i in 1:m){
    if (stop==0){
      beta.iterates[,i+1]<-beta.iterates[,i]+solve(f(beta.iterates[,i]))%*%s(beta.iterates[,i])
      for (j in 1:p){
        d[j]<-abs((beta.iterates[j,i+1]-beta[j])) 
      }
      if (prod(d<=epsilon)==TRUE){
        stop<-1
        stop.idx<-i+1
      }
    }
  }
  n.iter[k]<-stop.idx
}
n.iter
plot(nv,n.iter,type="l") ## shaky plot nothing concrete can be inferred


## Bernoulli Regression

n<-10000
x1=rnorm(n,mean=1,sd=0.5)
x2=rnorm(n,mean=1,sd=0.5)
X=cbind(1,x1,x2)
pn<-3
beta<-c(1,-2,3)
eta<-X%*%beta
p<-exp(eta)/(1+exp(eta))
## placeholder vector for response data
Y<-matrix(NA,nrow=n,ncol=1)

## response sampled from Poisson according to lambda vector
for (i in 1:n){
  Y[i]<-rbinom(1,1,p[i])
}

## Fisher Scoring Algorithm
m=100
s <- function(b){
  score<-t(X)%*%(Y-c(((exp(X%*%b))/(1+exp(X%*%b)))))
  return(score)
}
s(beta)

f<-function(b){
  w<-c(exp(X%*%b))/c(((1+exp(X%*%b))^2))
  W<-diag(c(w))
  fisher<-t(X)%*%W%*%X
  return(fisher)
}

beta.iterates<-matrix(0,nrow=pn, ncol=m+1)
beta.iterates[,1]<-c(1,-1,2)
sum<-beta.iterates[,1]
## threshold
epsilon<-10^(-2)
d<-rep(NA,pn)
stop<-0
stop.idx<-0
for (i in 1:m){
  if (stop==0){
    beta.iterates[,i+1]<-beta.iterates[,i]+solve(f(beta.iterates[,i]))%*%s(beta.iterates[,i])
    sum<-sum+beta.iterates[,i+1]
    for (j in 1:pn){
      d[j]<-abs((beta.iterates[j,i+1]-beta.iterates[j,i])) 
    }
    if (prod(d<=epsilon)==TRUE){
      stop<-1
      stop.idx<-i+1
    }
  }
}
stop.idx
beta.iterates[,stop.idx]
beta.iterates
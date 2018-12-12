##
## 
## MVNexample.r
##
## Introduction to the Multivariate Normal distribution
##


## 1. Preliminaries (install R package if needed)
## 2. Bivariate Normal Distribution
## 3. Groups of correlated observations
## 4. Correlation based on proximity (Moving Average in Time Series)
## 5. Correlation based on proximity (Exponential Decay in Covariance Function)


#######################################################
##
## 1. Preliminaries (download and install R package)
##
#######################################################

## this only needs to be run once on a computer
install.packages("mvtnorm")




#######################################################
##
## 2. Bivariate Normal distribution
##    (x1,x2)' ~ N( (mu1,mu2)', Sigma )
##
#######################################################

library(mvtnorm)

## vector mean
mu=c(1,2)

Sigma=matrix(NA,nrow=2,ncol=2)
## marginal variance of x1
s2.1=1
## marginal variance of x2
s2.2=3
## correlation
rho=.8

Sigma[1,1]=s2.1
Sigma[2,2]=s2.2
Sigma[1,2]=rho*sqrt(s2.1*s2.2)
Sigma[2,1]=rho*sqrt(s2.1*s2.2)
Sigma

X=rmvnorm(10,mean=mu,sigma=Sigma)
X
plot(X)

X=rmvnorm(1000,mean=mu,sigma=Sigma)
plot(X)

##
## sample moments
##
x1=X[,1]
x2=X[,2]

par(mfrow=c(1,2))
hist(x1)
hist(x2)

mean(x1) ## should be mu.1
mu[1]
mean(x2) ## should be mu.2
mu[2]
var(x1)  ## should be s2.1
s2.1
var(x2)  ## should be s2.2
s2.2
cor(X)   ## off diagonal elements should be rho




#######################################################
##
## 3. Groups of correlated observations
##    x1,...,x5 come from Group A, x6,...,x10 come from Group B
##    Group A and Group B are uncorrelated
##    Within groups, RVs are correlated
##
#######################################################

## mean vector
mu=rep(0,10)
## group correlation
rho=.95


## initialize covariance matrix
Sigma=diag(10) ## diagonal covariance matrix
## correlations for group A
Sigma[1:5,1:5] <- rho
## correlations for group B
Sigma[6:10,6:10] <- rho
## fix diagonal
diag(Sigma) <- 1
Sigma

## plot 9 realizations
par(mfrow=c(3,3))
X=rmvnorm(9,mean=mu,sigma=Sigma)
for(i in 1:9){
    plot(X[i,],pch=20,cex=0.1,mar=c(5,4,4,2))
}




#######################################################
##
## 4. Moving Average Correlation (time series)
##
#######################################################


n.obs=10 ## try for 10 and for 20

## mean vector
mu=rep(0,n.obs-1)

## variance of epsilon
sigma=1
## moving average parameter
phi=.9

## simulate epsilon ~ N(0,sigma^2*I)
epsilon=rnorm(n.obs,mean=mu,sd=sigma)
epsilon
## make "A" matrix
A=matrix(0,nrow=n.obs-1,ncol=n.obs)
for(i in 1:(n.obs-1)){
    A[i,i]=1
    A[i,i+1]=phi
}
A


## get x=mu+A*epsilon
x=mu+1/(1+phi^2)*A%*%epsilon 
x
par(mfrow=c(1,1))
plot(x,type="b")

## plot epsilon and x side by side
par(mfrow=c(4,2))
for(i in 1:4){
    epsilon=rnorm(n.obs,mean=mu,sd=sigma)
    x=mu+1/(1+phi^2)*A%*%epsilon
    plot(epsilon,type="b",main="epsilon",ylim=c(-2,2))
    plot(x,type="b",main="x",ylim=c(-2,2))
}




##
## Alternate formulation  x~N(mu,sigma^2*AA')
##

Sigma=sigma^2*A%*%t(A)
Sigma
x=rmvnorm(n=1,mean=mu,Sigma)
par(mfrow=c(1,1))
plot(x[1,],type="b")



#######################################################
##
## 5. Correlations based on proximity (time series or spatial)
##    cor(x_i,x_j)=exp(-d_{ij}/phi)
##
#######################################################


n.obs=100 ## try for 10 and for 100

## mean vector
mu=rep(0,n.obs)
## range parameter (large phi means higher correlation)
##                 (small phi means little correlation) 
phi=50  ## try for phi=.001 (nearly independent data) and for phi=n.obs/2

## specify locations:
locs=1:n.obs
locs
?dist
dist(locs)
## get distance matrix
D=as.matrix(dist(locs))
D

## get covariance matrix
Sigma=exp(-D/phi)
Sigma


## plot 9 realizations
par(mfrow=c(3,3))
X=rmvnorm(9,mean=mu,sigma=Sigma)
for(i in 1:9){
    plot(X[i,],type="b")
}





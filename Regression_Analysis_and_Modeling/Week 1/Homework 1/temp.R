y=rbinom(10000,50,3/5)
##make row of two plots
par(mfrow=c(1,2))
##plot the RVs
plot(y)
##plot a histogram of the RVs
hist(y,col="yellow")
z<-asin(y/50)
qqnorm(z)
qqline(z)
hist(z,col="blue")

load("interarrival.Rdata")
interarrival
y<-interarrival
par(mfrow=c(1,1))
hist(y,col="yellow")
##sample mean 
sm <- sum(y)/length(y)
##sample variance
sv <- (sum((y-sm)^2))/(length(y)-1)
min(y)
max(y)

dist=cars$dist
dist



## STAT 514:

m=5000 ## number of samples
n=225 ## sample size

A=matrix(0,nrow=m,ncol=n)

MLE.est1=matrix(0,nrow=m,ncol=1)

Med=matrix(0,nrow=m,ncol=1)
Med.est2=matrix(0,nrow=m,ncol=1)

xx<-matrix(0,nrow=m,ncol=n)
yy<-(1:n)/n
cdf.est3<-matrix(0,nrow=m,ncol=1)

pdf.est4.1<-matrix(0,nrow=m,ncol=1)
pdf.est4.2<-matrix(0,nrow=m,ncol=1)
pdf.est4.3<-matrix(0,nrow=m,ncol=1)
pdf.est4.4<-matrix(0,nrow=m,ncol=1)

for (i in 1:5000) {
  A[i,]<-rexp(225,1/10)
  
  MLE.est1[i,]<-mean(A[i,])
  MSE.est1<-mean((MLE.est1-10)^2)
  
  Med[i,]<-median(A[i,])
  Med.est2<-1.44*Med
  MSE.est2<-mean((Med.est2-10)^2)
  
  xx[i,]<-sort(A[i,])
  t1<-lm(log(1-yy+1/n) ~ xx[i,])
  cdf.est3[i,]<- -1/t1$coef[2]
  MSE.est3<-mean((cdf.est3-10)^2)
  
  h<-hist(A[i,])
  counts<-h$counts
  dens<-h$density[counts>0]
  midpts<-h$mids[counts>0]
  m1<-lm(log(dens) ~ midpts)
  m2<-lm(log(dens) ~ midpts, weights=counts[counts>0])
  pdf.est4.1[i,]<-exp(-m1$coef[1])
  MSE.est4.1<-mean((pdf.est4.1-10)^2)
  pdf.est4.2[i,]<--1/m1$coef[2]
  MSE.est4.2<-mean((pdf.est4.2-10)^2)
  pdf.est4.3[i,]<-exp(-m2$coef[1])
  MSE.est4.3<-mean((pdf.est4.3-10)^2)
  pdf.est4.4[i,]<--1/m2$coef[2]
  MSE.est4.4<-mean((pdf.est4.4-10)^2)
}

## STAT 511 Homework 2:
## Q1

n<-c(29,53,61,62,51,62,53,49,71,26)
y<-c(18,31,34,33,27,33,28,23,33,12)
p.MLE=sum(y)/sum(n)
p<-(c(1:1000))/1000
l<-matrix(0,nrow=length(p),ncol=1)
like <- function (x) x^(sum(y))*(1-x)^(sum(n)-sum(y))
lmax<-0
for (i in 1:length(p)){
  l[i]=like(p[i])
  if(l[i]>lmax){
    p.MLE.graph<-p[i]
    lmax<-l[i]
  }
}

plot(p,l, type="o", col="black", lwd=1)
title(main = "Plot of likelihood function vs. parameter p" , sub = NULL, xlab = NULL, ylab = "likelihood")
abline(v=p.MLE.graph,col=10,lty=1,lwd=4, xlab="p.MLE")
p.MLE
p.MLE.graph

## Q4(f)
n=20
mu=rep(0,n)
sig=1
tau=1
R=matrix(0,nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:n){
    R[i,j]=exp(-abs(i-j)/10)
  }
}
I=diag(n)
Sigma=((sig^2)*I)+((tau^2)*R)
library(mvtnorm)
Y=rmvnorm(9,mean=mu,sigma=Sigma)
Y

## Q4(g)
n=20
mu=rep(0,n)
sig=10
tau=1
R=matrix(0,nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:n){
    R[i,j]=exp(-abs(i-j)/10)
  }
}
I=diag(n)
Sigma=((sig^2)*I)+((tau^2)*R)
library(mvtnorm)
Y=rmvnorm(9,mean=mu,sigma=Sigma)
Y

## Q4(h)
n=20
mu=rep(0,n)
sig=1
tau=10
R=matrix(0,nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:n){
    R[i,j]=exp(-abs(i-j)/10)
  }
}
I=diag(n)
Sigma=((sig^2)*I)+((tau^2)*R)
library(mvtnorm)
Y=rmvnorm(9,mean=mu,sigma=Sigma)
Y

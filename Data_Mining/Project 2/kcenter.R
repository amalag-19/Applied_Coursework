library(MASS)
library(ggplot2)

euclid.dist<- function(x,y){
  d<-norm(as.matrix(x-y),type = "f")
  invisible(d)
}

## y is the matrix with columns as cluster objects h1, h2, ... hk
cluster.dist<- function(x,y){
  d<-min(apply(X = y, MARGIN = 2, FUN = function(t) euclid.dist(x,t)))
  invisible(d)
}

## simulating data
setwd(dir = "/Users/Amal/Box Sync/PSU/Fall 2015/Data Mining (STAT 557)/Project 2")
source("simulation_function.R")

## defining parameters
n <- 10000
n.cov <- 2
n.class <- 3
p.vec <- c(0.3, 0.2, 0.5)
mu <- cbind(c(1,1), c(2,2), c(3,3))
data<-sim(n,n.cov,n.class,p.vec,mu)

df<-data[,-1]
K<-n.class
head(df)
n<-nrow(df)
initial<-sample(1:n)[1]

H<-matrix(NA_real_,n.cov,K)
H[,1]<-as.matrix(df[initial,])
H.index<-rep(NA_integer_,K)
H.index[1]<-initial
H
H.index

dist<-apply(X = df,MARGIN = 1, FUN = function(x) euclid.dist(x,H[,1]))
clust<-rep(1,n)
for (i in 2:K){
  D<-max(dist[-which(is.na(H.index)==FALSE)])
  H.index[i]<-which.max(dist[-which(is.na(H.index)==FALSE)])
  H[,i]<-as.matrix(df[H.index[i],])
  for(j in 1:n){
    dist.new<-euclid.dist(df[j,],H[,i])
    if(dist.new<=cluster.dist(df[j,],H[,1:i])){
      dist[j]<-dist.new
      clust[j]<-i
    }
  }
}
df$clust<-clust
head(df)
p<-ggplot(data=df)
p+geom_point(mapping = aes(x = X1,y = X2,group = factor(clust), colour = factor(clust)))


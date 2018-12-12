library(MASS)
library(ggplot2)

## simulating a dataset
samp.size <- 10000
n.cov <- 3
n.class <- 3
p.vec <- c(0.3, 0.2, 0.5)
y <- sample(1:n.class, size = samp.size, prob = p.vec, replace = TRUE)
y
y.list <- split(x = y,f = y)
y.list[[1]]
mu <- cbind(c(1,1,1), c(2,2,2), c(3,3,3))
x.list <- list()
for(i in 1:n.class){
  x.list[[i]]<- mvrnorm(n = length(y.list[[i]]), mu = mu[,i], Sigma = diag(rep(1, n.cov)) )
}

for(i in 1:n.class){
  x.list[[i]]<- mvrnorm(n = length(y.list[[i]]), mu = mu[,i], Sigma = diag(.1,nrow = n.cov) )
}

temp2<-lapply(X = 1:n.class,FUN = function(x){
  temp<-data.frame(y.list[[x]],x.list[[x]])
  names(temp)<-c("Y","X1","X2","X3")
  return(temp)
})

datas<-do.call(what = rbind,args = temp2)
head(datas)
datas <- data.frame(datas)

p<-ggplot(data = datas)
p+geom_point(mapping = aes(x = X1,y = X2,group = factor(Y), colour = factor(Y)))

## Example multinomial logistitic model
library(nnet)

QDA.est<-function(dframe){
  n.class<-nlevels(as.factor(dframe[,1]))
  n.cov<-ncol(dframe)-1
  N<-nrow(dframe)
  prior<-rep(NA,n.class)
  mk<-matrix(NA,n.cov,n.class)
  Sk<-list()
  for(i in 1:n.class){
    index<-which(as.factor(dframe[,1])==(i-1))
    prior[i]<-length(index)/N
    mk[,i]<-colMeans(dframe[index,2:(n.cov+1)])
    df.centered<-as.matrix(dframe[index,2:(n.cov+1)]-matrix(mk[,i],length(index),n.cov,byrow=TRUE))
    Sk[[i]]<-(t(df.centered)%*%df.centered)/(N-n.class)
  }
  invisible(list("classes"=n.class,"prior"=prior,"mean"=mk,"variance"=Sk))
}

predict.QDA.1<-function(x,n.class,prior,mk,Sk){
  del<-rep(NA,n.class)
  for (i in 1:n.class){
    eig<-eigen(Sk[[i]])
    temp<-diag(eig$values^(-1/2),length(eig$values))%*%t(eig$vectors)
    temp<-temp%*%cbind(c((x-mk[,i])))
    quad.term<-t(temp)%*%temp
    del[i]<-((-1/2)*(sum(log(eig$values))+quad.term))+log(prior[i])
  }
  class<-which.max(del)
  invisible(class)
}

QDA.1<-function(dframe,x){
  temp1<-QDA.est(dframe)
  QDA.class<-predict.QDA.1(x,temp1[[1]],temp1[[2]],temp1[[3]],temp1[[4]])
  return(QDA.class-1)
}

predict.QDA.2<-function(x,n.class,prior,mk,Sk){
  x<-as.matrix(unlist(x))
  del<-rep(NA,n.class)
  for (i in 1:n.class){
    temp<-ginv(as.matrix(Sk[[i]]))
    del[i]<-t(x)%*%temp%*%mk[,i]-((1/2)*t(mk[,i])%*%temp%*%mk[,i])+log(prior[i])
    }
  QDA.class<-which.max(del)
  invisible(QDA.class)
}

QDA.2<-function(dframe,x){
  temp1<-QDA.est(dframe)
  QDA.class<-predict.QDA.2(x,temp1[[1]],temp1[[2]],temp1[[3]],temp1[[4]])
  return(QDA.class-1)
}

## QDA on simulated dataset
QDA.est(datas)
QDA(datas,as.matrix(datas[9999,-1]))
datas[1,-1]

## QDA.1 applied on aerial images dataset as a whole
load("cv_data.RData")
head(cv.data)
head(cv.data[[1]][[1]][,ncol(cv.data[[1]][[1]])])
df<-do.call(rbind,cv.data[[1]])
df[,1]<-df$Class2
df<-df[,-ncol(df)]
head(df)

c.pred<-rep(NA,nrow(df))

for(i in 1:nrow(df)){
  c.pred[i]<-QDA.1(df,as.matrix(df[i,-1]))
}

MSE<-sum((c.pred-df[,1])^2)/nrow(df)

## defining cross validation function
{# Function: trsfrm.gen 
#
# Purpose: Generate transformation functions. Normalization during 
# cross validation needs to be based on training data only. So, since
# sample means and sds from the training data will change depending
# on which partition of the data is used for training, different
# transformations will be used for each iteration.
#
# Arguments: 
# train.data - training data matrix with column names (should all
#             be numeric). 
# Returns:
# trsfrm - function for log transformation and normalization based
# on sample means and sds from trainin data. trsfrm() takes a matrix
# with the same column names passed to trsfrm.gen (can be both 
# training and test data)
#
# Example:
# tr.fun <- trsfrm.gen(train.data)
# tr.fun(test.data)


trsfrm.gen <- function(train.data){
  ncols <- dim(train.data)[2]
  for(i in 1:ncols){
    if(names(train.data)[i] %in% c("BrdIndx", "Area", 
                                   "Compact", "ShpIndx","SD_NIR", "LW")){
      train.data[,i] <- log(train.data[,i])
    }
  }
  train.data <- as.matrix(train.data)
  mns <- apply(X = train.data, MARGIN = 2, FUN = mean)
  sds <- apply(X = train.data, MARGIN = 2, FUN = mean)
  trsfrm <- function(x, mn = mns, sd = sds){
    ncols <- dim(x)[2]
    for(i in 1:ncols){
      if(names(x)[i] %in% c("BrdIndx", "Area", 
                            "Compact", "ShpIndx","SD_NIR", "LW")){
        x[,i] <- log(x[,i])
      }
    }
    x <- as.matrix(x)
    nrml.data <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    xcols <- dim(x)[2]
    for(i in 1:xcols){
      nrml.data[,i] <- (x[,i] - mn[i])/sd[i]
    }
    return(nrml.data)
  }
  return(trsfrm)
}}

trsfrm.gen()
p<-ggplot(data = df)
p+geom_point(mapping = aes(x = df[,2],y = df[,3],group = factor(df[,1]), colour = factor(df[,1])))
ins

library(MASS)
library(ggplot2)

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

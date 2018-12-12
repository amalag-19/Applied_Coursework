library(MASS)
library(ggplot2)
require(reshape2)
require(Matrix)

samp.size <- 10000
n.cov <- 3
n.class <- 3
p.vec <- c(0.3, 0.2, 0.5)
y <- sample(1:n.class, size = samp.size, prob = p.vec, replace = TRUE)
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
tail(datas)

names(datas)<-c("Y","X1","X2","X3")
nm.df<-names(datas)
datas$Intercept<-1

datas<-datas[,c("Intercept",nm.df)]
head(datas)

mm<-model.matrix(~factor(Y),data=datas)
head(mm)
datas.n<-data.frame(mm[,-1],datas[,-2])
datas.n$ID<-1:nrow(datas.n)

names(datas.n)<-c("Y2","Y3","Intercept","X1","X2","X3","ID")
head(datas.n)
datas.n1<-melt(datas.n,id.vars = c("ID","Intercept","X1","X2","X3"))
head(datas.n1)
nrow(datas.n1)
X<-datas.n[,-c(1,2,7)]
X<-as.matrix(X)
head(X)
s<-(paste0("bdiag(X",rep(",X",length(unique(datas$Y))-2),")"))
s<-parse(text=s)
X.til<-eval(s)
X.til<-as.matrix(X.til)

head(X.til)
beta.X<-matrix(0,nrow=ncol(X.til),ncol=1)
l.pred<-X.til%*%beta.X
head(l.pred)

datas.n2<-data.frame(datas.n1$ID,datas.n1$value,X.til,l.pred)
names(datas.n2)[c(1,2)]<-c("ID","Y")
names(datas.n2)
head(datas.n2)
datas.n2$e.pred<-exp(datas.n2$l.pred)
aggdata<-aggregate(datas.n2$e.pred,by=list(datas.n2$ID),FUN=function(x)1+sum(x))
head(aggdata)
pos<-match(datas.n2$ID,aggdata$Group.1)
datas.n2$s.pred<-aggdata[pos,"x"]
head(datas.n2)
datas.n2$p<-datas.n2$e.pred/datas.n2$s.pred
temp<-as.matrix(datas.n2$p)%*%t(as.matrix(datas.n2$p))


head(temp[,1:6])

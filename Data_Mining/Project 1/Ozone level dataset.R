## ozone level dataset
data<-read.csv("onehr.data",na.strings = "?")
data<-data[which(complete.cases(data)==TRUE),]
head(data)
head(data[,1])
data$X0.<-as.factor(data$X0.)
levels(data$X0.)
data[,1]<-data[,ncol(data)]
data<-data[,-ncol(data)]

head(data[,1:6])

undebug(QDA.est)
y<-QDA.est(data)
str(y)
str(y[[4]])


debug(predict.QDA.2)
predict.QDA.2(data[1,-1],y[[1]],y[[2]],y[[3]],y[[4]])

undebug(QDA.2)

QDA.predictions<-apply(X = as.matrix(1:nrow(data)), MARGIN = 1, FUN = function(x) QDA.2(data,data[x,-1]))
MSE<-mean((QDA.predictions-as.numeric(data[,1]))^2)

df<-data.frame(QDA.predictions,Y=data[,1])
df<-split(df,df$Y)
sens<-sum(df[[2]]$QDA.predictions==df[[2]]$Y)/nrow(df[[2]])
spec<-sum(df[[1]]$QDA.predictions==df[[1]]$Y)/nrow(df[[1]])
spec
df$Y

which(QDA.predictions==0)

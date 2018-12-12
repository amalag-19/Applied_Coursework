oly=read.csv("Olympics.csv")
str(oly)

fit=lm(goldtime~year*gender,data=oly)
summary(fit)
yhat=fit$fitted
res=fit$res

plot(oly$year,res)
plot(yhat,res)




idx.M=which(oly$gender=="M")
col=rep("red",nrow(oly))
col[idx.M] <- "black"

plot(oly$year,oly$goldtime,col=col,pch=20)


##
## regression line for men / women
##

oly.m=oly[idx.M,]
oly.w=oly[-idx.M,]

fit.m=lm(goldtime~year,data=oly.m)
fit.w=lm(goldtime~year,data=oly.w)

xvals=1900:2250
xvals.df=data.frame(year=xvals)
yvals.m=predict(fit.m,newdata=xvals.df)
yvals.w=predict(fit.w,newdata=xvals.df)


plot(xvals,yvals.m,type="l",col="black",lwd=3,ylim=c(7,14))
points(xvals,yvals.w,type="l",col="red",lwd=3,ylim=c(7,14))
points(oly$year,oly$goldtime,col=col,pch=20)


##
## 95% CI of mean (pointwise confidence band)
##

X.m=cbind(1,oly.m$year)
X.w=cbind(1,oly.w$year)

XtX.inv.m=solve(t(X.m)%*%X.m)
XtX.inv.w=solve(t(X.w)%*%X.w)


sigma.hat.m=summary(fit.m)$sigma
sigma.hat.w=summary(fit.w)$sigma

X.unobs=cbind(1,xvals)
se.pred.m=rep(0,length(xvals))
for(i in 1:length(se.pred.m)){
    se.pred.m[i]=sigma.hat.m*sqrt( t(X.unobs[i,])%*%XtX.inv.m%*%X.unobs[i,])
}
se.pred.w=rep(0,length(xvals))
for(i in 1:length(se.pred.w)){
    se.pred.w[i]=sigma.hat.w*sqrt( t(X.unobs[i,])%*%XtX.inv.w%*%X.unobs[i,])
}

n.m=nrow(X.m)
n.w=nrow(X.w)

upper.m=yvals.m+qt(1-.05/2,df=n.m-2)*se.pred.m
lower.m=yvals.m-qt(1-.05/2,df=n.m-2)*se.pred.m

upper.w=yvals.w+qt(1-.05/2,df=n.w-2)*se.pred.w
lower.w=yvals.w-qt(1-.05/2,df=n.w-2)*se.pred.w


##
## Now 95% prediction interval
##

X.unobs=cbind(1,xvals)
se.pred.m=rep(0,length(xvals))
for(i in 1:length(se.pred.m)){
    se.pred.m[i]=sigma.hat.m*sqrt(1+ t(X.unobs[i,])%*%XtX.inv.m%*%X.unobs[i,])
}
se.pred.w=rep(0,length(xvals))
for(i in 1:length(se.pred.w)){
    se.pred.w[i]=sigma.hat.w*sqrt(1+ t(X.unobs[i,])%*%XtX.inv.w%*%X.unobs[i,])
}

n.m=nrow(X.m)
n.w=nrow(X.w)

upper.pred.m=yvals.m+qt(1-.05/2,df=n.m-2)*se.pred.m
lower.pred.m=yvals.m-qt(1-.05/2,df=n.m-2)*se.pred.m

upper.pred.w=yvals.w+qt(1-.05/2,df=n.w-2)*se.pred.w
lower.pred.w=yvals.w-qt(1-.05/2,df=n.w-2)*se.pred.w

##
## Plot both
##

par(mfrow=c(1,2))
plot(xvals,yvals.m,type="l",col="black",lwd=3,xlim=c(1900,2050),ylim=c(9,13),main="CI on MEAN")
points(xvals,yvals.w,type="l",col="red",lwd=3)
points(oly$year,oly$goldtime,col=col,pch=20)
points(xvals,upper.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,lower.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,upper.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))
points(xvals,lower.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))

plot(xvals,yvals.m,type="l",col="black",lwd=3,xlim=c(1900,2050),ylim=c(9,13),main="Prediction CI")
points(xvals,yvals.w,type="l",col="red",lwd=3)
points(oly$year,oly$goldtime,col=col,pch=20)
points(xvals,upper.pred.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,lower.pred.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,upper.pred.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))
points(xvals,lower.pred.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))

##
## Predicting into the future
##

par(mfrow=c(1,2))
plot(xvals,yvals.m,type="l",col="black",lwd=3,ylim=c(7,14),main="CI on MEAN")
points(xvals,yvals.w,type="l",col="red",lwd=3)
points(oly$year,oly$goldtime,col=col,pch=20)
points(xvals,upper.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,lower.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,upper.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))
points(xvals,lower.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))

plot(xvals,yvals.m,type="l",col="black",lwd=3,ylim=c(7,14),main="Prediction CI")
points(xvals,yvals.w,type="l",col="red",lwd=3)
points(oly$year,oly$goldtime,col=col,pch=20)
points(xvals,upper.pred.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,lower.pred.m,type="l",lty=2,col="black",lwd=1,ylim=c(7,14))
points(xvals,upper.pred.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))
points(xvals,lower.pred.w,type="l",lty=2,col="red",lwd=1,ylim=c(7,14))





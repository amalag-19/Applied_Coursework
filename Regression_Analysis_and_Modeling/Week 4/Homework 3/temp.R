data(trees)
?trees
str(trees)
x<-trees$Height
y<-trees$Girth

## EDA
plot(x,y)
plot(x,sqrt(y))
plot(x^2,y)
plot(x,log(y))
plot(log(x),y)

fit=lm(y~x)
plot(x,y)
abline(fit)

fit=lm(sqrt(y)~x)
plot(x,sqrt(y))
abline(fit)

fit=lm(y~I(x^2))
plot(x^2,y)
abline(fit)

res=fit$resid
hist(res)
yhat=fit$fitted.values
plot(yhat,y)
summary(fit)$sigma
## the table given in summary(fit)
summary(fit)$coef

Munich=read.csv("rent99.raw",sep=" ")
Munich3<-Munich[,c(1,3,4)]
str(Munich3)
pairs(Munich3)
area<-Munich3[,2]
yearc<-Munich3[,3]

fit=lm(rent~area+yearc,data=Munich)
summary (fit)

residuals=fit$resid
yhat=fit$fitted.values
plot(yhat,residuals)
beta=fit$coef
residuals.area=residuals+beta[2]*Munich[,2]

x1<-rep(1,length(area))
X=matrix(0,nrow=length(area), ncol=3)
X[,1]<-x1
X[,2]<-area
X[,3]<-yearc
X
str(X)
names(X)

## plot data

plot(Nile,main="Annual flow of the river Nile at Ashwan")


## date when Ashwan dam was constructed

abline(v=1902,col="blue")

## create covariate
y=as.numeric(Nile)
y
x=rep(0,length(y))
1902-1871
x[32:length(x)]=1

plot(x)
pairs(cbind(y,x))

## fit linear model

fit=lm(y~x)
summary(fit)

ypred=predict(fit)
ypred
res=resid(fit)

## examine model fit graphically
plot(Nile,main="Annual flow of the river Nile at Ashwan - Data with Fitted Model")
points(1871:1970,ypred,type="l",col="red")

## examine residuals for normality

par(mfrow=c(1,3))
plot(1871:1970,res,pch=20,type="b",main="Model Residuals Over Time")
hist(res,main="Histogram of Model Residuals")
qqnorm(res)
qqline(res,col="red")

## examine residuals for serial correlation

r=res[-1]
r
r.prev=res[-length(res)]
r.prev

par(mfrow=c(1,2))
plot(r,r.prev)
summary(lm(r~r.prev))

acf(res)


## Specify and fit a Linear Model with AR(1) time series correlated random effect

library(nlme)

group=rep(1,length(y))
t=1871:1970

fit.gls=gls(y~x,correlation=corAR1(form= ~t))
summary(fit.gls)
intervals(fit.gls)


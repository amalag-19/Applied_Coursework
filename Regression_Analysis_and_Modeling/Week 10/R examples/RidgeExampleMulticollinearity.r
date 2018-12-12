####
#### RidgeExample.r
####


## 
## load diabetes data
##
## source:      Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle
##              Regression" (with discussion) _Annals of Statistics_
##

load("diab.Rdata")

str(diab)

## note how correlated some of the predictors are
pairs(diab)
cor(diab)


## fit some small models - note the effect of correlated variables
##   NOTE: compare betahat for "tc" in all four models below!
##         (including "full" model)
##
fit.lm=lm(y~ldl,data=diab)
summary(fit.lm)
fit.lm=lm(y~tc,data=diab)
summary(fit.lm)
fit.lm=lm(y~tc+ldl,data=diab)
summary(fit.lm)


## fit a full model with many predictors
fit.lm=lm(y~.,data=diab)
summary(fit.lm)

## variance inflation factors
library(car)
vif(fit.lm)

##
## ridge regression
##
library(glmnet)

## glmnet takes an "X" matrix and a "y" vector separately
head(diab)
X=as.matrix(diab[,-1])
y=diab$y

rr=glmnet(x=X,y,alpha=0,lambda=exp(seq(-5,10,by=.1)))

##
## choose tuning parameter by cross-validation
##
fit=cv.glmnet(x=X,y=y,alpha=0,lambda=exp(seq(-5,10,by=.1)))
plot(fit)


## get lambda and best rr fit
lambda.rr=fit$lambda.1se
lambda.rr

## some plots
##   NOTE: multicollinearity does NOT adversely affect prediction!
##           
par(mfrow=c(1,2))
plot(fit)
abline(v=log(lambda.rr))
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter",label=T)
abline(v=log(lambda.rr))


## betas from best ridge regression
coef(fit)
## betas from lm
coef(fit.lm)

cbind(as.numeric(coef(fit)),coef(fit.lm))


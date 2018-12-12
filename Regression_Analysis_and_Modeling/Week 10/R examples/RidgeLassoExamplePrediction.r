####
#### RidgeLassoExamplePrediction.r
####
#### Outline:
####  1. Simulating data and splitting into train / test sets
####  2. AIC stepwise selection
####  3. Ridge regression
####  4. Lasso regression
####  5. Comparison of CVMSPE on training set and MSPE on test set
####  6. Ridge / Lasso for GLMs (logistic regression example)


#######################################################
##
## (1) Simulating 50 predictors and one response
##     v1-v10 are important
##     v11-v50 are not
##
#######################################################

## generate a dummy dataset with 30 predictors (10 useful & 20 useless) 
## note that each predictor has mean zero and sd=1
N=150
P=50

X=matrix(NA,nrow=N,ncol=P)


##
## do this for correlated predictors
##
covmat=matrix(rnorm(P^2,sd=2),nrow=P)
covmat

covmat=covmat+t(covmat)
U=eigen(covmat)$vectors
D=diag(rexp(P,rate=10))
covmat=U%*%D%*%t(U)

##
## Do this for uncorrelated predictors
##
covmat=diag(P)

install.packages("mvtnorm")
library(mvtnorm)
for(i in 1:N){
    X[i,]=rmvnorm(1,mean=rep(0,P),sigma=covmat)
}
X=data.frame(X)
head(X)

## true betas
betas.true=c(1,2,3,4,5,-1,-2,-3,-4,-5,rep(0,P-10))

## simulating "y"
sigma=15.7
X=as.matrix(X)
y=X%*%betas.true+rnorm(N,mean=0,sd=sigma)


## splitting into "train" and "test" data

alldata=data.frame(cbind(y,X))
names(alldata)[1] <- "y"
head(alldata)
train=alldata[1:100,]
test=alldata[101:150,]


## linear regression
fit=lm(y~.,data=train)
summary(fit)

betas.lm=coef(fit)

yhat.lm=predict(fit,newdata=test)
mspe.lm=mean((test$y-yhat.lm)^2)
mspe.lm

#######################################################
##
## (2) AIC Stepwise Selection
##
#######################################################
fit.AIC=step(fit)
summary(fit.AIC)

yhat.AIC=predict(fit.AIC,newdata=test)
mspe.AIC=mean((test$y-yhat.AIC)^2)
mspe.AIC

#######################################################
##
## (3) ridge regression
##
#######################################################

## first check VIFs
library(car)
vif(fit)
library(glmnet)
## alpha=0 gives ridge regression
## alpha=1 gives lasso regression

## fit ridge (trying 100 different lambda values)
rr=glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=0,nlambda=100)
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")

## use 10-fold crossvalidation to find the best lambda
cv.rr=cv.glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=0,nfolds=10,nlambda=100)

## getting cvmspe from best value of lambda
cvmspe.rr=min(cv.rr$cvm)

## get lambda and best rr fit
lambda.rr=cv.rr$lambda.min
lambda.rr

## some plots
par(mfrow=c(1,2))
plot(cv.rr)
abline(v=log(lambda.rr))
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")
abline(v=log(lambda.rr))

## beta estimates for best lambda
betas.rr=coef(cv.rr,s="lambda.min")

plot(betas.rr,betas.lm,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

yhat.rr=predict(cv.rr,s="lambda.min",newx=as.matrix(test[,-1]))
mspe.rr=mean((test$y-yhat.rr)^2)
mspe.rr



#######################################################
##
## (4) lasso regression
##
#######################################################

## alpha=0 gives ridge regression
## alpha=1 gives lasso regression

## fit lasso (trying 100 different lambda values)
lasso=glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=1,nlambda=100)
plot(lasso,xvar="lambda",main="Lasso Regression Betas for Different Values of the Tuning Parameter")
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")

## use 10-fold crossvalidation to find the best lambda
cv.lasso=cv.glmnet(x=as.matrix(train[,-1]),y=as.numeric(train[,1]),alpha=1,nfolds=10)

## get lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.min
lambda.lasso

## getting cvmspe from best value of lambda
cvmspe.lasso=min(cv.lasso$cvm)

## some plots
par(mfrow=c(1,2))
plot(cv.lasso)
abline(v=log(lambda.lasso))
plot(lasso,xvar="lambda")
abline(v=log(lambda.lasso))

## beta estimates for best lambda
betas.lasso=coef(cv.lasso,s="lambda.min")
betas.lasso

yhat.lasso=predict(cv.lasso,newx=as.matrix(test[,-1]),s="lambda.min")
mspe.lasso=mean((test$y-yhat.lasso)^2)
mspe.lasso

#############################################################
##
## (5) comparison of CVMSPE on training set and MSPE on test set
##
#############################################################

##
## V-fold Cross Validation for best AIC model 
##

formula(fit.AIC)

## divide training data into V equal groups
n=nrow(train)
V=10
idx.group=rep(1:V,floor(n/V+1))
idx.group=idx.group[1:n]
idx.group=sample(idx.group)
idx.group

## cross-validation
cvmspe=rep(NA,V)
for(v in 1:V){
    idx.holdout=which(idx.group==v)
    fit=lm(formula(fit.AIC),data=train[-idx.holdout,])
    yhat=predict(fit,newdata=train[idx.holdout,])
    cvmspe[v]=mean((train$y[idx.holdout]-yhat)^2)
}
cvmspe

cvmspe.AIC=mean(cvmspe)

##
## compare CVMSPE
##

cvmspe.AIC
cvmspe.rr
cvmspe.lasso


##
## Comparison of MSPE on test set
##

mspe.AIC
mspe.rr
mspe.lasso


## table
T=matrix(c(cvmspe.AIC,cvmspe.rr,cvmspe.lasso,mspe.AIC,mspe.rr,mspe.lasso),ncol=2)
rownames(T)=c("AIC","RR","Lasso")
colnames(T)=c("CVMSPE-Train","MSPE-Test")
T






#######################################################
##
## (6) GLM with ridge / lasso 
##
#######################################################


## simulate logistic regression

N=150
P=50

X=matrix(rnorm(N*P),nrow=N,ncol=P)
## true betas
betas.true=c(1,2,3,4,5,-1,-2,-3,-4,-5,rep(0,P-10))

## simulating "y"
X=as.matrix(X)
eta.true=X%*%betas.true
pi.true=1/(1+exp(eta.true))
y=rbinom(N,size=1,prob=pi.true)
y

## lasso
lasso=glmnet(x=X,y=y,family="binomial",alpha=1,nlambda=100)
## use 10-fold crossvalidation to find the best lambda
cv.lasso=cv.glmnet(x=X,y=y,alpha=1,nfolds=10)

## get lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.1se
lambda.lasso

## some plots
par(mfrow=c(1,2))
plot(cv.lasso)
abline(v=log(lambda.lasso))
plot(lasso,xvar="lambda")
abline(v=log(lambda.lasso))

## beta estimates for best lambda
betas.lasso=coef(cv.lasso)
betas.lasso

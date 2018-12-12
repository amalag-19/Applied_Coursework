# Note all the comments in this file precede the code they explain.
# Notable results from code chunks are highlighted in the report.

# Loading the dataset
load("fire.Rdata")
head(fire)

# creating data frame X consisting the response and the subset of predictors
# described as "characteristics near the home" that homeowners have control over
X<-data.frame(fire$burnt,fire$planted,fire$buildings,fire$perc.woody,fire$perc.cleared,fire$distance2bush,fire$distance2tree)

# Exploratory data analysis for "characteristics near a home" predictors
# Pairwise scatter plots in X
pairs(X)

# Since the response is binary, it seems best to use logistic regression.
# Model 1: Fitting a logistic model with all the predictors
fit=glm(burnt~.,data=fire,family="binomial")
summary(fit)
for.m1<-formula(fit)

# Logistic regression with ONE predictor at a time for all the 
# "characteristics near the home" to see their relationship with the 
# probability that a home is burned in a wildfire

# Fitting with planted
fit1=glm(burnt~planted,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=as.factor(c('r','p'))
newdata=data.frame(planted=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
par(mfrow=c(2,3))
plot(xvals,eta,main="Linear Predictor",xlab="planted",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="planted",ylab=expression(mu),ylim=c(0,1))

# Fitting with buildings
fit1=glm(burnt~buildings,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=seq(1,10)
newdata=data.frame(buildings=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
par(mfrow=c(1,1))
plot(xvals,eta,main="Linear Predictor",xlab="buildings",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="buildings",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(fire$buildings),fire$burnt)

## Now fit with perc.woody
fit1=glm(burnt~perc.woody,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=seq(1,100)
newdata=data.frame(perc.woody=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
plot(xvals,eta,main="Linear Predictor",xlab="perc.woody",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="perc.woody",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(fire$perc.woody),fire$burnt)

## Now fit with perc.cleared
fit1=glm(burnt~perc.cleared,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=seq(1,100)
newdata=data.frame(perc.cleared=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
plot(xvals,eta,main="Linear Predictor",xlab="perc.cleared",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="perc.cleared",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(fire$perc.cleared),fire$burnt)

## Now fit with distance2bush
fit1=glm(burnt~distance2bush,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=seq(1,100)
newdata=data.frame(distance2bush=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
plot(xvals,eta,main="Linear Predictor",xlab="distance2bush",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="distance2bush",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(fire$distance2bush),fire$burnt)

## Now fit with distance2tree
fit1=glm(burnt~distance2tree,family="binomial",data=fire)
summary(fit1)
##get estimated linear predictor
xvals=seq(1,100)
newdata=data.frame(distance2tree=xvals)
eta=predict(fit1,newdata=newdata,type="link")

## get estimated mean
plot(xvals,eta,main="Linear Predictor",xlab="distance2tree",ylab=expression(eta),type="l")
mu=predict(fit1,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="distance2tree",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(fire$distance2tree),fire$burnt)

# Seperating the observations based on planted status "r" or "p"
n<-length(fire$planted)
p<-ncol(fire)
idx<-rep(0,n)
j<-1
for (i in 1:n){
  if (fire$planted[i]=="r"){
    idx[j]<-i
    j<-j+1
  }
}
idx<-idx[which(idx!=0)]

# dataframe for observations with planted status "r"
fire.r<-fire[idx,]
fire.r$planted
# Total number of observations with planted status "r"
N1<-nrow(fire.r)

# dataframe for observations with planted status "p"
fire.p<-fire[-idx,]
fire.p$planted
# Total number of observations with planted status "p"
N2<-n-N1

# Creating short names for columns (just for convenience)
Y<-fire$burnt
pw<-fire$perc.woody
d2t<-fire$distance2tree

# Placeholder vector for the difference in response (burnt value) for
# the two houses when they have similar nearby characteristics
# but different planted status
diff<-rep(0,N1*N2)
k<-1
for (i in 1:N1){
  for (j in 1:N2){
    if ((abs(d2t[i]-d2t[j])<=2)&&(abs(pw[i]-pw[j])<=5)){
      diff[k]<-fire.r$burnt[i]-fire.p$burnt[j]
      k<-k+1
    }
  }
}
mean(diff)

# Model Diagnostics: Checking for multicollinearity
vif(fit)
# No problems of multicollinearity since all the the values are much
# less than 10

# Outliers
r.star<-rstudent(fit)
plot(r.star)
cutoff=qt(.975,df=n-p-1)
abline(h=cutoff)
abline(h=-cutoff)

# Influential points (Leverage criteria)
h=hatvalues(fit)
cutoff.h=2*p/n
cutoff.h
plot(h)
abline(h=cutoff)

# Influential Points (Cook's Distance criteria)
cd=cooks.distance(fit)
plot(cd)

# Checking AIC
AIC(fit)

# Residual Plots for the Model 1
par(mfrow=c(2,2))
plot(fit)

# Partial Residual Plots
library(car)
crPlots(fit)
# These plots look pretty good. We could leave the model as it is
# (based on partial residual plots)

# Model Diagnostics: Simulating the data from model 1 and matching
# the residuals and QQ plots:
p.hat=predict(fit,type="response")
Y.sim=rbinom(length(p.hat),size=1,prob=p.hat)
fire.sim<-fire
fire.sim$burnt<-Y.sim
fit.sim=glm(burnt~.,data=fire,family="binomial")
summary(fit.sim)
par(mfrow=c(2,2))
plot(fit.sim)
# Residuals and QQ Plots look the same. Verifies our model 1 to some extent

#############################################################################
# Model Selection: Trying different approaches of regularization to estimate
# parameters and select the best model: 

# Interaction terms

# Modified Model 1
fit=glm(burnt ~ ffdi + slope + aspect + topo + perc.cleared + amt.burntless5yrs + 
          perc.burntless5yrs + amt.not.burnt5to10yrs + perc.burnt5to10yrs + 
          amt.unlogged + perc.logged + amt.not.NP + amt.not.SF + adj.for.type + 
          edge + distance2tree + planted + buildings + perc.woody + 
          distance2bush+perc.woody*planted,data=fire,family="binomial")
summary(fit)

# Dividing the data into training and test sets in the ratio 2:1 (approx.)
train=fire[1:325,]
test=fire[326:n,]
fit.t=glm(burnt ~ ffdi + slope + aspect + topo + perc.cleared + amt.burntless5yrs + 
                perc.burntless5yrs + amt.not.burnt5to10yrs + perc.burnt5to10yrs + 
                amt.unlogged + perc.logged + amt.not.NP + amt.not.SF + adj.for.type + 
                edge + distance2tree + planted + buildings + perc.woody + 
                distance2bush+perc.woody*planted,data=train,family="binomial")
betas.fit<-coef(fit.t)

# MSPE for Model 1
mu=predict(fit,newdata=test,type="response")
yhat<-rep(0,nrow(test))
yhat[which(c(mu)>=0.5)]<-1
mspe.model1=mean((test$burnt-yhat)^2)
mspe.model1

# Model 2: Stepwise AIC model selection
fit.AIC=step(fit.t)
summary(fit.AIC)
f.AIC<-formula(fit.AIC)
coef(fit.AIC)

# MSPE for Model 2 (AIC)
Yhat.AIC=predict(fit.AIC,newdata=test[,-1])
yhat<-rep(0,nrow(test))
yhat[which(c(Yhat.AIC)>=0.5)]<-1
mspe.AIC=mean((test$burnt-yhat)^2)
mspe.AIC

##########################################
# Model 3: Ridge regression
# Checking VIF's
vif(fit.t)
# vif of edge is greater than 14.45>10. Can use ridge to alleviate this 
# multicollinearity

# Modifying the data structure to make it suitable for glmnet
xfactors<-model.matrix(fire$burnt~fire$adj.for.type+fire$planted+fire$planted*fire$perc.woody)[1:325,-1]
x<-as.matrix(data.frame(fire[1:325,c(-1,-15,-18)], xfactors))

# Fitting ridge (trying 100 different lambda values)
library(glmnet)
rr=glmnet(x,y=as.factor(train[,1]),alpha=0,family='binomial')
par(mfrow=c(1,1))
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")

# Using 10-fold crossvalidation to find the best lambda
cv.rr=cv.glmnet(x,y=as.numeric(train[,1]),alpha=0,nfolds=10,nlambda=100)

# Getting cvmspe from best value of lambda
cvmspe.rr=min(cv.rr$cvm)

## get lambda and best rr fit
lambda.rr=cv.rr$lambda.min
lambda.rr

## Some relevant plots
par(mfrow=c(1,1))
plot(cv.rr,main="Plot of MSE vs. log of Tuning Parameter")
abline(v=log(lambda.rr))
par(mfrow=c(1,1))
plot(rr,xvar="lambda",main="Betas for Different Values of Tuning Parameter")
abline(v=log(lambda.rr))

## Beta estimates for best lambda
betas.rr=coef(cv.rr,s="lambda.min")
betas.rr

# Prediction at test set
xfactors<-model.matrix(fire$burnt~fire$adj.for.type+fire$planted+fire$planted*fire$perc.woody)[326:n,-1]
x<-as.matrix(data.frame(fire[326:n,c(-1,-15,-18)], xfactors))
yhat.rr=predict(cv.rr,s="lambda.min",newx=x)


yhat<-rep(0,nrow(test))
yhat[which(c(yhat.rr)>=0.5)]<-1
mspe.rr=mean((test$burnt-yhat)^2)
mspe.rr

##########################################
# Model 4: Lasso regression
# Modifying the data structure to make it suitable for glmnet
xfactors<-model.matrix(fire$burnt~fire$adj.for.type+fire$planted+fire$planted*fire$perc.woody)[1:325,-1]
x<-as.matrix(data.frame(fire[1:325,c(-1,-15,-18)], xfactors))

# Fitting lasso (trying 100 different lambda values)
lasso=glmnet(x,y=as.factor(train[,1]),alpha=1,family='binomial')
par(mfrow=c(1,1))
plot(lasso,xvar="lambda",main="Lasso Regression Betas for Different Values of the Tuning Parameter")

# Comparing lasso and ridge plots side by side
par(mfrow=c(1,2))
plot(rr,xvar="lambda",main="Ridge Regression Betas for Different Values of the Tuning Parameter")
plot(lasso,xvar="lambda",main="Lasso Regression Betas for Different Values of the Tuning Parameter")

# Using 10-fold crossvalidation to find the best lambda
cv.lasso=cv.glmnet(x,y=as.numeric(train[,1]),alpha=1,nfolds=10)

# Getting lambda and best lasso fit
lambda.lasso=cv.lasso$lambda.min
lambda.lasso

# Getting cvmspe from best value of lambda
cvmspe.lasso=min(cv.lasso$cvm)

# Some relevant plots
par(mfrow=c(1,1))
plot(cv.lasso,main="Plot of MSE vs. log of Tuning Parameter")
abline(v=log(lambda.lasso))
par(mfrow=c(1,1))
plot(lasso,xvar="lambda",main="Betas for different Values of Tuning Parameter")
abline(v=log(lambda.lasso))

# Beta estimates for best lambda
betas.lasso=coef(cv.lasso,s="lambda.min")
betas.lasso

# Prediction at test set
xfactors<-model.matrix(fire$burnt~fire$adj.for.type+fire$planted+fire$planted*fire$perc.woody)[326:n,-1]
x<-as.matrix(data.frame(fire[326:n,c(-1,-15,-18)], xfactors))
yhat.lasso=predict(cv.lasso,newx=x,s="lambda.min")
yhat<-rep(0,nrow(test))
yhat[which(c(yhat.lasso)>=0.5)]<-1
mspe.lasso=mean((test$burnt-yhat)^2)
mspe.lasso

#################################################
# V-fold Cross Validation for best AIC model 

# Dividing training data into V equal groups
nt=nrow(train)
V=10
idx.group=rep(1:V,floor(nt/V+1))
idx.group=idx.group[1:nt]
idx.group=sample(idx.group)
idx.group

# Cross-validation
# placeholder vector for cvmspe for AIC model
cvmspe1=rep(NA,V)
# placeholder vector for cvmspe for model 1
cvmspe2=rep(NA,V)
for(v in 1:V){
  idx.holdout=which(idx.group==v)
  fit1<-glm(formula(fit.AIC),data=train[-idx.holdout,],family='binomial')
  fit2<-glm(formula(fit),data=train[-idx.holdout,],family='binomial')
  yhat1=predict(fit1,newdata=train[idx.holdout,])
  yhat2=predict(fit2,newdata=train[idx.holdout,])
  cvmspe1[v]=mean((train$burnt[idx.holdout]-yhat1)^2)
  cvmspe2[v]=mean((train$burnt[idx.holdout]-yhat2)^2)
}
cvmspe.AIC=mean(cvmspe1)
cvmspe.model1<-mean(cvmspe2)

# comparing CVMSPE of all the models
cvmspe.model1
cvmspe.AIC
cvmspe.rr
cvmspe.lasso

# Comparing MSPE on test set for all models
mspe.model1
mspe.AIC
mspe.rr
mspe.lasso

# Comparing multiple approaches for estimating regression parameters in a table
T=matrix(c(cvmspe.model1,cvmspe.AIC,cvmspe.rr,cvmspe.lasso, mspe.model1, mspe.AIC,mspe.rr,mspe.lasso),ncol=2)
rownames(T)=c("Model 1","AIC","RR","Lasso")
colnames(T)=c("CVMSPE-Train","MSPE-Test")
T

# Conclusion: Best predictive model obtained from AIC according to MSPE on test set 
# i.e. Model 2.

###########################################
# Fitting the model from AIC penalization over all dataset

fitfull.AIC<-glm(f.AIC,data=fire,family='binomial')

# Beta estimates
betas.AIC=coef(fitfull.AIC)
betas.AIC
p<-length(betas.AIC)

# Final Model Diagnostics: Checking for multicollinearity
vif(fitfull.AIC)
# No problems of multicollinearity since all the the values are much
# less than 10

# Outliers
par(mfrow=c(1,1))
r.star<-rstudent(fitfull.AIC)
plot(r.star)
cutoff=qt(.975,df=n-p-1)
abline(h=cutoff)
abline(h=-cutoff)

# Influential points (Leverage criteria)
h=hatvalues(fit)
cutoff.h=2*p/n
cutoff.h
plot(h)
abline(h=cutoff)

# Influential Points (Cook's Distance criteria)
cd=cooks.distance(fit)
plot(cd)

# Residual Plots for the final model
par(mfrow=c(2,2))
plot(fit)

# Model Diagnostics: Simulating the data from final model and matching
# the residuals and QQ plots:
p.hat=predict(fitfull.AIC,type="response")
Y.sim=rbinom(length(p.hat),size=1,prob=p.hat)
fire.sim<-fire
fire.sim$burnt<-Y.sim
fit.sim=glm(burnt~.,data=fire,family="binomial")
summary(fit.sim)
par(mfrow=c(2,2))
plot(fit.sim)
# Residuals and QQ Plots look the same. Verifies our final model to a good extent

# Prediction by final model
yhat.AIC=predict(fitfull.AIC,newdata=fire,type="response")
yhat<-rep(0,n)
yhat[which(yhat.AIC>=0.5)]<-1
yhat
# Difference in the predicted and the original response
diff<-abs(yhat-fire$burnt)

# Percentage of time the prediction is correct by the best predictive model
(1-mean(diff))*100

#############################################################################
# Q3 (a)

# Forming the subset of observations where there are any trees within 10 meters
# of the house
fire.t<-fire[which(d2t<=10),]

# Removing trees within 10 meters of all houses in this subset 
fire.t$distance2tree<-10
m<-nrow(fire.t)

# Prediction of the response (houses burnt) with this managament practice 
# under Model 2
mu<-predict(fitfull.AIC,newdata=fire.t,type="response")
Y.hat<-rep(0,m)
Y.hat[which(mu>=0.5)]<-1

# Difference in the predicted and the original response
diff<-Y.hat-Y[which(d2t<=10)]
#diff<-diff[which(abs(diff)==1)]
mean(diff)
# This shows that the probability of getting a house burnt decreases by 0.0479
# with this management practice under model 2

# Q3 (b)

# Forming the subset of observations where there is greater than 50%
# woody vegetation within 40 meters of the house
fire.w<-fire[which(fire$perc.woody>=50),]

# Removing woody vegetation until all the houses only have 50% woody vegetation
# in this subset 
fire.w$perc.woody<-50
m<-nrow(fire.w)

# Prediction of the response (houses burnt) with this managament practice 
# under Model 4
mu<-predict(fitfull.AIC,newdata=fire.w,type="response")
Y.hat<-rep(0,m)
Y.hat[which(mu>=0.5)]<-1

# Difference in the predicted and the original response
diff<-Y.hat-Y[which(fire$perc.woody>=50)]
#diff<-diff[which(abs(diff)==1)]
mean(diff)

# This shows that the probability of getting a house burnt decreases by 0.0487
# with this management practice under model 2

# Q3 (c)

# Forming the subset of observations where there is remnant vegetation
# within 40 meters of the house

fire.p<-fire[which(fire$planted=="r"),]
m<-nrow(fire.p)

# Removing remnant vegetation until all the houses only have planted vegetation
# in this subset 
for (i in 1:m){
  fire.p$planted[i]<-as.factor('p')
}

# Prediction of the response (houses burnt) with this managament practice 
# under Model 2
mu<-predict(fitfull.AIC,newdata=fire.p,type="response")
Y.hat<-rep(0,m)
Y.hat[which(mu>=0.5)]<-1

# Difference in the predicted and the original response
diff<-Y.hat-Y[which(fire$planted=="r")]
#diff<-diff[which(abs(diff)==1)]
mean(diff)

# This shows that the probability of getting a house burnt decreases by 0.68
# with this management practice under model 2

# Total No. of houses burnt in 2009 fires 
N<-sum(Y)
# Total no. of houses that would have been burned in the 2009 fires if 
# this management practice had been implemented before 2009
N.m<-sum(Y[which(fire$planted=="p")])+sum(Y.hat)

# Semiparametric Modeling: gam
library(mgcv)
fit.gam=gam(burnt~s(perc.woody+perc.cleared)+s(distance2bush+distance2tree)+buildings+planted,data=fire,family=binomial)
summary(fit.gam)
plot(fit.gam)

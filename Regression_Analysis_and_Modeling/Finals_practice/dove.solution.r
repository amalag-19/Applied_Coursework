## read in data

dove=read.csv("dove.csv")
head(dove)

## create P.dev = P-P.clim
##  and   lT.dev
##  and   eT.dev

dove$P.dev=dove$P-dove$P.clim
dove$eT.dev=dove$eT-dove$eT.clim
dove$lT.dev=dove$lT-dove$lT.clim


##
## Exploratory data analysis
##

pairs(dove)

#########################################################
##
##
## Model 1 : global adaptation
## y_i ~ Binom(N_i,p_i)
## logit(p_i) = eta_i
##
#########################################################

fit=glm(cbind(y,N-y)~P+eT+lT+Lat+Lon,data=dove,family=binomial)
summary(fit)
par(mfrow=c(2,3))
plot(fit)

##
## potential outlier / influential point
##

dove[35,]

## highly suspicious.  Both "y" and "N" are "999" which seems very odd
hist(dove$y,breaks=100)
hist(dove$N,breaks=100)

## Assume: 999=missing value.  Remove row 35 from analysis.

which(dove$y==999)

dove=dove[-35,]

## now re-fit

fit.global=glm(cbind(y,N-y)~P+eT+lT+Lat+Lon,data=dove,family=binomial)
summary(fit.global)
par(mfrow=c(2,3))
plot(fit.global)

##
## look for nonlinear relationships
##

library(car)
crPlots(fit.global)
## looks pretty good.  Could leave as is (based on partial residual plots),
##  or could add in polynomials:
## Only a quadratic in "P" is significant:

fit.global=glm(cbind(y,N-y)~poly(P,2)+eT+lT+Lat+Lon,data=dove,family=binomial)
summary(fit.global)
par(mfrow=c(2,3))
plot(fit.global)

crPlots(fit.global)

##
## Check for multicollinearity
##

vif(fit.global)
## borderline.  Could use ridge regression.  At least need to acknowledge multicollinearity


##
## Check residuals for appropriateness:
##
par(mfrow=c(2,2))
plot(fit.global)

## no problems here.  QQ plot has long tails.  Could try simulating data from model to
## see if this pattern persists:
p.hat=predict(fit.global,type="response")
y.sim=rbinom(length(p.hat),size=dove$N,prob=p.hat)
dove.sim=dove
dove.sim$y=y.sim
fit.global.sim=glm(cbind(y,N-y)~poly(P,2)+eT+lT+Lat+Lon,data=dove.sim,family=binomial)
summary(fit.global.sim)
par(mfrow=c(2,2))
plot(fit.global.sim)
## pattern persists - so is not a problem at all.



#########################################################
##
##
## Model 2 : local adaptation
## y_i ~ Binom(N_i,p_i)
## logit(p_i) = eta_i
##
#########################################################


fit.local=glm(cbind(y,N-y)~P.dev+eT.dev+lT.dev+Lat+Lon,data=dove,family=binomial)
summary(fit.local)
par(mfrow=c(2,3))
plot(fit.local)

##
## look for nonlinear relationships
##

library(car)
crPlots(fit.local)
## looks pretty good.  Could leave as is (based on partial residual plots),
##  or could add in polynomials, but none are significant:

fit.local=glm(cbind(y,N-y)~P.dev+eT.dev+lT.dev+Lat+Lon,data=dove,family=binomial)
summary(fit.local)
par(mfrow=c(2,3))
plot(fit.local)

crPlots(fit.local)

##
## Check for multicollinearity
##

vif(fit.local)
## just fine.  No problems here


##
## Check residuals for appropriateness:
##
par(mfrow=c(2,2))
plot(fit.local)

## no problems here.  QQ plot has long tails.  Could try simulating data from model to
## see if this pattern persists:
p.hat=predict(fit.local,type="response")
y.sim=rbinom(length(p.hat),size=dove$N,prob=p.hat)
dove.sim=dove
dove.sim$y=y.sim
fit.local.sim=glm(cbind(y,N-y)~P.dev+eT.dev+lT.dev+Lat+Lon,data=dove.sim,family=binomial)
summary(fit.local.sim)
par(mfrow=c(2,2))
plot(fit.local.sim)
## pattern persists - so is not a problem at all.





##
## Differential Adaptation Model
##
#########################################################
##

##
## Model 3 : Differential  adaptation
## y_i ~ Binom(N_i,p_i)
## logit(p_i) = eta_i
##
#########################################################


## try a model with interactions
fit.diff=glm(cbind(y,N-y)~P.clim*P.dev+lT.clim*lT.dev+eT.clim*eT.dev+Lat+Lon,data=dove,family=binomial)
summary(fit.diff)
## lT interaction is not significant - remove
fit.diff=glm(cbind(y,N-y)~P.clim*P.dev+lT.clim+lT.dev+eT.clim*eT.dev+Lat+Lon,data=dove,family=binomial)
summary(fit.diff)

## stepwise AIC or p-value model selection will remove both lT covariates
fit.diff=step(fit.diff)
summary(fit.diff)


##
## look for nonlinear relationships (nothing significant)
##

fit.diff.nonsig=glm(cbind(y,N-y)~P.clim*P.dev+lT.clim+lT.dev+eT.clim*poly(eT.dev,2)+Lat+Lon,data=dove,family=binomial)
summary(fit.diff.nonsig)

##
## Check for multicollinearity
##

vif(fit.diff)
## strong multicollinearity, partly due to the interactions.
## try centering / scaling covariates
X=dove[,-c(1,2)]
X.sc=scale(X)
summary(X)
summary(X.sc)

dove.sc=cbind(dove[,1:2],X.sc)

fit.diff=glm(cbind(y,N-y)~P.clim*P.dev+eT.clim*eT.dev+Lat+Lon,data=dove.sc,family=binomial)
summary(fit.diff)

vif(fit.diff)
## much better.  Alternative is to use ridge regression!

##
## Check residuals for appropriateness:
##
par(mfrow=c(2,2))
plot(fit.diff)

## no problems here.  QQ plot has long tails.  Could try simulating data from model to
## see if this pattern persists:
p.hat=predict(fit.diff,type="response")
y.sim=rbinom(length(p.hat),size=dove$N,prob=p.hat)
dove.sim=dove
dove.sim$y=y.sim
fit.diff.sim=glm(cbind(y,N-y)~P.clim*P.dev+eT.clim*eT.dev+Lat+Lon,data=dove.sc,family=binomial)
summary(fit.diff.sim)
par(mfrow=c(2,2))
plot(fit.diff.sim)
## pattern persists - so is not a problem at all.



## plot of 2D mean function
P.clim.vals=0:80
P.dev.vals=-40:40
xy=expand.grid(P.clim.vals,P.dev.vals)
names(xy) <- c("P.clim","P.dev")
xy$eT.dev=0
xy$eT.clim=10
xy$Lat=0
xy$Lon=0
z=predict(fit.diff,newdata=xy)
zmat=matrix(z,nrow=length(P.clim.vals),ncol=length(P.dev.vals))

filled.contour(P.clim.vals,P.dev.vals,zmat)
par(mfrow=c(1,1))
image(P.clim.vals,P.dev.vals,zmat)
points()
## interpretation: deviation from mean is very important at high precip areas
##  but not very important at low precip areas.


## a similar plot for eT would be needed.




##
## alternative using gam
##

library(mgcv)
fit.diff.gam=gam(cbind(y,N-y)~s(P.clim,P.dev)+s(eT.clim,eT.dev)+Lat+Lon,data=dove,family=binomial)
summary(fit.diff.gam)
plot(fit.diff.gam)


##
## compare models
##

AIC(fit.global)
AIC(fit.local)
AIC(fit.diff)

##
## Prediction
##

## Do on your own! 


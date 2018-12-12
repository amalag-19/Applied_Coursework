#################################################
##
## SemiparExamples.r
##
## Contents:
##   (1) Semiparametric modeling in one dimension
##       y = f(x) + e
##   (2) Semiparametric modeling in two dimensions
##       y = f(x,z) + e
##   (3) Varying Coefficient models
##       y = x*f(z) + e
##
##
#################################################
library(mgcv)

Munich=read.csv("rent99.raw",sep=" ")

## Exploratory Data Analysis

str(Munich)       ## gives general structure
summary(Munich)   ## computes summary statistics of each column
head(Munich)      ## shows the first 6 rows of the data frame

## Plotting relevant variables

pairs(Munich[,2:4])


#####################################################
##
##   (1) Semiparametric modeling in one dimension
##       y = f(x) + e
##
#####################################################

##
## (1a) rentsqm = f(area)+e
##

fit.area=gam(rentsqm~s(area),data=Munich)
summary(fit.area)

## plot of mean function (and CI band)
plot(fit.area,resid=T,pch=2,lwd=3,cex=.015)

## plot of residuals
res=resid(fit.area)
plot(Munich$area,res)
abline(h=0)

##
## (1b) adding in other predictor variables (linear predictor)
## rentsqm = beta0 + beta*yearc + f(area)+e
##

fit=gam(rentsqm~yearc+s(area),data=Munich)
summary(fit)
## plot of mean function (and CI band)
plot(fit,resid=T,pch=2,lwd=3,cex=.015)
## plot of residuals
res=resid(fit)
plot(Munich$area,res)
abline(h=0)


##
## (1c) Now just yearc
## rentsqm = f(yearc)+e
##

fit.yearc=gam(rentsqm~s(yearc),data=Munich)
summary(fit.yearc)
## plot of mean function (and CI band)
plot(fit.yearc,resid=T,pch=2,lwd=3,cex=.015)
## plot of residuals
res=resid(fit.yearc)
plot(Munich$yearc,res)
abline(h=0)


##
## (1d) Now both (plus some covariates)
## rentsqm = Xbeta + f(area) + g(yearc)+e
##

fit.both=gam(rentsqm~s(area)+s(yearc),data=Munich)

## plot of mean function (and CI band)
par(mfrow=c(1,2))
plot(fit.both,resid=T,pch=2,lwd=3,cex=.005)
## plot of residuals
res=resid(fit.both)
par(mfrow=c(1,1))
plot(Munich$yearc,res)
abline(h=0)


#####################################################
##
##   (2) Semiparametric modeling in two dimensions
##       y = f(x,z) + e
##
#####################################################

##
## (2a) rentsqm = f(area,yearc)+e
##

fit.ay=gam(rentsqm~s(area,yearc),data=Munich)
summary(fit.ay)
## plot of 2D mean function
plot(fit.ay)
## plot of residuals
res=resid(fit.ay)
plot(Munich$area,res)
abline(h=0)
plot(Munich$yearc,res)
abline(h=0)
qqnorm(res)
qqline(res)


##
## (2b) Interaction between continuous variables
##    rentsqm = B0 + B1*area + B2*yearc + B3*area*yearc + e
##

fit.int=lm(rentsqm~area*yearc,data=Munich)
summary(fit.int)

## plot of 2D mean function
area.vals=20:160
yearc.vals=1920:2000
xy=expand.grid(area.vals,yearc.vals)
xy
names(xy) <- c("area","yearc")
z=predict(fit.int,newdata=xy)
zmat=matrix(z,nrow=length(area.vals),ncol=length(yearc.vals))
filled.contour(area.vals,yearc.vals,zmat)
image(area.vals,yearc.vals,zmat)

## plot of residuals
plot(fit.int)
res=resid(fit.int)
plot(Munich$area,res)
abline(h=0)
plot(Munich$yearc,res)
abline(h=0)
qqnorm(res)
qqline(res)


#####################################################
##
##   (3) Varying Coefficient Model
##       y = f(x)*z + e
##
#####################################################

##
## (3a) rentsqm = B0+area*f(yearc)+e
##

fit.vc=gam(rentsqm~s(yearc,by=area),data=Munich)
summary(fit.vc)
plot(fit.vc)




##
## comparison of all models by AIC
##
AIC(fit.area)
AIC(fit.yearc)
AIC(fit.both)
AIC(fit.ay)
AIC(fit.int)
AIC(fit.vc)


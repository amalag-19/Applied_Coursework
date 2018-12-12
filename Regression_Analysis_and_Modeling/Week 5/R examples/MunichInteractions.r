#################################################
##
##
## MunichInteractions.r
##
##
#################################################


## read in rent99.raw data
Munich=read.csv("rent99.raw",sep=" ")

str(Munich)       ## gives general structure
summary(Munich)   ## computes summary statistics of each column
head(Munich)      ## shows the first 6 rows of the data frame

fit=lm(rent~area,data=Munich)
summary(fit)

##
## Diagnostic plots
##

## plot fitted values vs residuals (check for heteroscedasticity)
par(mfrow=c(1,1))
plot(yhat,res)
## plot residuals vs each covariate (check for heteroscedasticity and nonlinearity)
plot(Munich$area,res)
abline(0,0,col="red")
## partial residual plots (requires the "car" package)
library(car)
crPlots(fit)
## QQ-plot of residuals (check normality)
qqnorm(res)
qqline(res)



##
## adding in a categorical covariate
##

fit=lm(rent~area+kitchen,data=Munich)
summary(fit)

##
## Diagnostic plots
##

## plot fitted values vs residuals (check for heteroscedasticity)
plot(yhat,res)
## plot residuals vs each covariate (check for heteroscedasticity and nonlinearity)
plot(Munich$area,res)
abline(0,0,col="red")
## partial residual plots (requires the "car" package)
library(car)
crPlots(fit)
## QQ-plot of residuals (check normality)
qqnorm(res)
qqline(res)


##
## adding in an interaction term
##

fit=lm(rent~area+kitchen+area*kitchen,data=Munich)
summary(fit)

##
## Diagnostic plots
##

## plot fitted values vs residuals (check for heteroscedasticity)
plot(yhat,res)
## plot residuals vs each covariate (check for heteroscedasticity and nonlinearity)
plot(Munich$area,res)
abline(0,0,col="red")
## partial residual plots (requires the "car" package)
library(car)
crPlots(fit)
## QQ-plot of residuals (check normality)
qqnorm(res)
qqline(res)


beta.hat=fit$coeff
beta.hat
##
## plotting interaction effect
##

plot(Munich$area,Munich$rent,pch=20,col=Munich$kitchen+1)
abline(beta.hat[1],beta.hat[2],lwd=3)
abline(beta.hat[1]+beta.hat[3],beta.hat[2]+beta.hat[4],lwd=3,col=2)





##
## Now an interaction between two categorical covariates
##

## categorical dummy coding for "location" 3=best, 2=mid-range location, 1=worst

Munich$loc2=rep(0,nrow(Munich))
Munich$loc2[Munich$location==2]=1

Munich$loc3=rep(0,nrow(Munich))
Munich$loc3[Munich$location==3]=1

fit=lm(rent~area+loc2*kitchen+loc3*kitchen,data=Munich)

summary(fit)

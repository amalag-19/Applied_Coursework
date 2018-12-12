#################################################
##
##
## RegressionPlots.r
##
##
#################################################


## These are examples of descriptive plots to show 
## the effect of covariates, especially when the  
## effect is nonlinear
##
## I will use the "predict" command in R
## (these plots can be done by hand as well)

#################################################
##
## Outline
## (1) Plots for one predictor variable
## (2) Partial residual plots for 2+ predictor variables
## (3) Plots for interactions between categorical and continuous covariates
##
#################################################

## read in rent99.raw data
Munich=read.csv("rent99.raw",sep=" ")


#################################################
##
##
## (1) Plots for one predictor variable
##
## Idea: a) plot scatterplot
##       b) create "x" values
##       c) use "predict" to get "f(x)" values
##       d) plot the ordered pairs (x,f(x))
#################################################

## Use polynomial basis functions to model nonlinear effect 

fit=lm(rentsqm~poly(area,3),data=Munich)
summary(fit)


## (a) scatterplot
plot(Munich$area,Munich$rentsqm)

## (b) create predictor variables to use
x.vals=seq(from=20,160,by=1)
x.vals
## turn this into a data.frame object so that R can use the "predict" function
## NOTE: call your variable the same thing as your predictor variable is called
##       in your regression data.frame
area.vals=data.frame(area=x.vals)
str(area.vals)

## (c) get "f(area)" values (regression line evaluated at "x.vals"
f.vals=predict(fit,newdata=area.vals)

## (d) plot (x,f(x))
points(x.vals,f.vals,type="l",lwd=3,col="red")









#################################################
##
##
## (2) Partial residual plots for 2+ predictor variable
##
## Idea: a) get eps.k  =  k-th partial residuals 
##       b) fit regression   eps.k ~ f(x.k)
##       c) plot scatterplot and regression line as in part (1)
#################################################

## Use inverse functions and polynomial basis functions to model nonlinear effect 

fit=lm(rentsqm~I(1/area)+poly(yearc,2),data=Munich)
summary(fit)

##
##
## Partial residuals for area
##
##

## fit using all variables except for area
fit.yearc=lm(rentsqm~poly(yearc,2),data=Munich)
summary(fit.yearc)

## get k-th partial residuals
eps.k=fit.yearc$resid

## fit regression with eps.k as response and f(area) as predictor
part.res.fit=lm(eps.k~I(1/area),data=Munich)
summary(part.res.fit)

## scatterplot of area vs partial residuals
plot(Munich$area,eps.k)

## create predictor variables to use
x.vals=seq(from=20,160,by=1)
area.vals=data.frame(area=x.vals)

## (c) get "f(area)" values (regression line evaluated at "x.vals"
f.vals=predict(part.res.fit,newdata=area.vals)

## (d) plot (x,f(x))
points(x.vals,f.vals,type="l",lwd=3,col="red")



##
##
## Partial residuals for yearc
##
##

## fit using all variables except for yearc
fit.area=lm(rentsqm~I(1/area),data=Munich)
summary(fit.area)

## get k-th partial residuals
eps.k=fit.area$resid

## fit regression with eps.k as response and f(area) as predictor
part.res.fit=lm(eps.k~poly(yearc,2),data=Munich)
summary(part.res.fit)

## scatterplot of yearc vs partial residuals
plot(Munich$yearc,eps.k)

## create predictor variables to use
x.vals=seq(from=1910,2000,by=1)
yearc.vals=data.frame(yearc=x.vals)

## (c) get "f(yearc)" values (regression line evaluated at "x.vals"
f.vals=predict(part.res.fit,newdata=yearc.vals)

## (d) plot (x,f(x))
points(x.vals,f.vals,type="l",lwd=3,col="red")




#################################################
##
##  (3) Plots for interactions between categorical and continuous covariates
##
## Idea: a) fit regression  
##       b) plot scatterplot of response vs continuous covariate
##       c) plot one regression line for each category
#################################################




##
## (a) fitting regression with interaction
##

fit=lm(rentsqm~I(1/area)*kitchen,data=Munich)
summary(fit)
#par(mfrow=c(2,2))
plot(fit)

##
## (b) scatterplot
##
par(mfrow=c(1,1))
plot(Munich$area,Munich$rentsqm)

##
## (c) plot one regression line for each category
##

x.vals=seq(from=20,to=160,by=1)

## create data.frame with kitchen=0
## (for regression line when kitchen=0)
kitchen.vals=rep(0,length(x.vals))
df.0=data.frame(area=x.vals,kitchen=kitchen.vals)


## create data.frame with kitchen=1
## (for regression line when kitchen=1)
kitchen.vals=rep(1,length(x.vals))
df.1=data.frame(area=x.vals,kitchen=kitchen.vals)

## look at the two data.frames
head(df.0)
head(df.1)

## get f(area) when kitchen=0
f.vals.0=predict(fit,newdata=df.0)

## get f(area) when kitchen=1
f.vals.1=predict(fit,newdata=df.1)

## plot regression line for kitchen=0
points(x.vals,f.vals.0,type="l",col="red",lwd=3)

## plot regression line for kitchen=1
points(x.vals,f.vals.1,type="l",col="blue",lwd=3)

##
## add legend to plot
##
legend("topright",legend=c("Kitchen=0","Kitchen=1"),col=c("red","blue"),lwd=3)

legend("topright",col=c("red","blue"),legend=c("Kitchen=0","Kitchen=1"),lwd=10)

##
## Interpretation:
##  
## When kitchen=0, we have a strong decreasing trend in rentsqm as area increases,
##  so houses without a designer kitchen typically have lower rentsqm as area
##  increases
##
## When kitchen=1, this trend is almost non-existant.  That is, houses WITH a 
##  designer kitchen show almost no effect between area and rentsqm.
##

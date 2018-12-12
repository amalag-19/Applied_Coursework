golf=read.table("golffull.raw",header=TRUE)
head(golf)
pairs(golf)

## price     - in euros
## age       - in months
## kilometer - kilometer reading in 1000 kilometers
## TIA       - number of months until next required inspection
## extras1   - 1=ABS brake, 0=no ABS brake
## extras2   - 1=sunroof, 0=no sunroof


fit.all=lm(price~age+kilometer+TIA+extras1+extras2,data=golf)
summary(fit.all)

plot(fit.all)
library(car)
crPlots(fit.all)


##############################################################
##
## Approach #1: Adding polynomial terms one at a time
##              Accept them if the regression parameter is
##               significantly different from zero
##
##############################################################

## for "age", age^2 is significant but age^3 is not

fit=lm(price~poly(age,2)+kilometer+TIA+extras1+extras2,data=golf)
summary(fit)

fit=lm(price~poly(age,3)+kilometer+TIA+extras1+extras2,data=golf)
summary(fit)

## for "kilometer", km^2 is significant but km^3 is not

fit=lm(price~poly(age,2)+poly(kilometer,2)+TIA+extras1+extras2,data=golf)
summary(fit)

fit=lm(price~poly(age,2)+poly(kilometer,3)+TIA+extras1+extras2,data=golf)
summary(fit)


##############################################################
##
## Approach #2: AIC model selection 
##              (stepwise model selection)
##
##############################################################

##
## Starting from a null model
##

fit.0=lm(price~1,data=golf)
summary(fit.0)

## add one
add1(fit.0,scope=~ageop1+ageop2+kilometerop1+kilometerop2+extras1+extras2+TIA)

fit.1=lm(price~ageop1,data=golf)
summary(fit.1)

## add another one
add1(fit.1,scope=~ageop1+ageop2+kilometerop1+kilometerop2+extras1+extras2+TIA)

##
## could continue this, or could also consider "drop1" in which one term is dropped
##
## However, "step" does this automatically.  It considers add1, drop1 repeatedly
##  and searches to find the best model via AIC
##
## Simple approach: fit a BIG model (too many terms)
##                  then do stepwise selection via the "step" command:
##

fit.all=lm(price~.-kilometer-age,data=golf)
summary(fit.all)

best.AIC=step(fit.all)
summary(best.AIC)






##############################################################
##
## Comparing Models via AIC  
##             - Must have SAME RESPONSE VARIABLE (can't transform)
##             - Otherwise models can be very different
##
##############################################################

## load polio data
load("polio.Rdata")
str(polio)

plot(polio$y,type="b")

## fit LM and get AIC
fit.lm <- lm(y ~ .,data=polio)
summary(fit.lm)
AIC(fit.lm)

## fit Poisson GLM and get AIC

fit.glm <- glm(y ~ .,data=polio,family=poisson)
summary(fit.glm)
AIC(fit.glm)

## AIC for GLM (558) is better (lower) than AIC for LM (680)
## So GLM fits the data better and is preferred.

## stepwise selection using AIC
summary(fit.glm)
best.fit=step(fit.glm)
summary(best.fit)

## note that dropping any variable leads to worse AIC:
drop1(fit.glm)

## also note that AIC model selection often leaves in variables
## that are NOT significantly different from zero
##
## AIC is a global measure of model fit
## p-values are a local measure of model fit

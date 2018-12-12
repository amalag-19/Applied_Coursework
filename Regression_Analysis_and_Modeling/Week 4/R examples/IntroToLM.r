#################################################
##
##
## IntroToLM.r
##
##
#################################################

## Outline:
## (1) Getting data into R using "read.csv"
## (2) Manipulating data frames
## (3) Basics of "lm"
## (4) Residual plots
## (5) Manipulating variables within "lm"


#################################################
##
##
## (1) Reading data into R
##
##
#################################################

## General procedure:
##  (a) get data into ".csv" or similar tabular format using Excel / etc
##  (b) read into R using "read.csv"
##      Make sure that you are in the right directory
##  (c) result is a data frame in R


?read.csv
## important read.csv options:
##    sep=      specify what separates entries.  Could be commas (default) or blank space (shown below)
##    header=   either "TRUE" or "FALSE".  Set to true if the first row contains variable names.

## read in rent99.raw data
Munich=read.csv("rent99.raw",sep=" ")


#################################################
##
##
## (2) Manipulating Data Frames in R
##
##
#################################################

##
## Important summary functions
##

str(Munich)       ## gives general structure
summary(Munich)   ## computes summary statistics of each column
head(Munich)      ## shows the first 6 rows of the data frame
pairs(Munich)     ## plots a matrix of pairwise scatterplots
names(Munich)     ## gives the names of the columns in the data frame
dim(Munich)       ## gives the number of rows and columns in the data frame

##
## Subsetting data frames
##

## selecting only first 10 rows
m10row=Munich[1:10,]
m10row

## selecting only first 3 columns
m3col=Munich[,1:3]
str(m3col)

## selecting the 1st, 2nd, and 4th columns
Msmall=Munich[,c(1,2,4)]

## calling columns by name
y=Munich$rent
x=Munich$area
par(mfrow=c(1,1))
plot(x,y)

## adding columns to data frame
RV=rnorm(n=nrow(Munich),mean=0,sd=1)
Munich$newcolumn=RV
str(Munich)

## adding a column that is a transformation of an existing column
Munich$area.squared=Munich$area^2
str(Munich)

## using an index to keep only rows with a kitchen=1
idx.kitchen=which(Munich$kitchen==1)
idx.kitchen
Mk=Munich[idx.kitchen,]
str(Mk)



#################################################
##
##
## (3) Basics of the "lm" command in R
##
##
#################################################

?lm
## Important inputs to lm:
##
##   formula:   in the form of   response ~ linear predictor 
##       For example, to fit the model:
##       y_i ~ N( mu + beta_1*x_i1 + beta_2*x_i2 , sigma^2 ) 
##       use the formula:
##       y~x1+x2
##
##   data:      data frame containing response and predictor variables
##

##
## ex: fit rent_i ~ N(beta0 + beta1*area_i , sigma^2)
##
fit=lm(rent~area,data=Munich)


##
## ex: fit rent_i ~ N(beta0 + beta1*area_i +beta2*yearc_i , sigma^2)
##
fit=lm(rent~area+yearc,data=Munich)




## Most of the important information about the fit is given by the "summary" command

summary(fit)

##
## Extracting residuals, regression parameter estimates, and predicted values from an lm object
##

## estimates of regression parameters
fit$coef

## residuals
res=fit$resid
hist(res)

## fitted values (y.hat = E(y|X,beta.hat))
yhat=fit$fitted.values
plot(yhat,y)


##
## Extracting estimates of variance and test statistics for common hypothesis tests
##

## estimate of sigma (not sigma^2)
summary(fit)$sigma

## the table given in summary(fit)
summary(fit)$coef

## t-statistic values for H0: beta_k=0 vs Ha: beta_k nonzero
summary(fit)$coef[,3]

## p values for H0: beta_k=0 vs Ha: beta_k nonzero
summary(fit)$coef[,4]



#################################################
##
##
## (3) Residual diagnostic plots
##
##
#################################################

##
## Diagnostic plots
##

## plot fitted values vs residuals (check for heteroscedasticity)
plot(yhat,res)
## plot residuals vs each covariate (check for heteroscedasticity and nonlinearity)
plot(Munich$area,res)
abline(h=0,col="red")
plot(Munich$yearc,res)
abline(0,1,col="red")
## partial residual plots (requires the "car" package)
library(car)
crPlots(fit)
## QQ-plot of residuals (check normality)
qqnorm(res)
qqline(res)





#################################################
##
##
## (5) Manipulating variables within lm
##
##
#################################################

##
## lm allows you to manipulate variables within the "formula" statement
##

## EX1: log-transform response
fit=lm(log(rent)~area+yearc,data=Munich)
summary(fit)

## EX2: log-transform "area"
fit=lm(rent~log(area)+yearc,data=Munich)
summary(fit)


##
## The "I" function is used to transform variables when the
##   transformation is not written as a function in R
##   EX: polynomial powers, inverse, ...
##

## EX3: inverse of "area" DOESNT WORK WITHOUT "I"
fit=lm(rent~1/area+yearc,data=Munich)
summary(fit)

## EX4: inverse of "area"    USING "I"
fit=lm(rent~I(1/area)+yearc,data=Munich)
summary(fit)

## EX5  polynomial powers of area
fit=lm(rent~area+I(area^2)+I(area^3)+yearc,data=Munich)
summary(fit)


##
## Instead of using the polynomials in EX4 and EX5 above,
##  it is often preferable to use ORTHOGONAL polynomial 
##  basis functions.  These can be obtained using "poly"
##

w=1:100
X=poly(w,degree=4)
matplot(X)
## the resulting polynomial basis functions are all orthogonal:
t(X)%*%X
## and zero-centered
summary(X)


## EX6 polynomial powers of area using "poly"
fit=lm(rent~poly(area,3)+yearc,data=Munich)
summary(fit)


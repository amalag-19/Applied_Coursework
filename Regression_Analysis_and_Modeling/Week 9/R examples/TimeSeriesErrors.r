####################################################################################
## Polio Time Series
## Description:
##      Time series of Polio incidences in U.S.A. from 1970 to 1983.
## Format:
##      A data frame with the 168 monthly observations (from January 1970
##      to December 1983) with the following variables
##        ‘y’                 time series of polio incidences.
##        ‘t*10^( -3 )’       linear trend multiplied by factor 10^{(-3)}.
##        ‘cos( 2*pi*t/12 )’  cosine annual seasonal component.
##        ‘sin( 2*pi*t/12 )’  sine annual seasonal component.
##        ‘cos( 2*pi*t/6 )’   cosine semi-annual seasonal component.
##        ‘sin( 2*pi*t/6 )’   sine semi-annual seasonal component.
## Source:
##      Zeger, S.L. (1988). A regression model for time series of counts.
##      _Biometrika_ *75*, 822-835.
# More info.: http://artax.karlin.mff.cuni.cz/r-help/library/gcmr/html/polio.html


# Load data

load("polio.Rdata")
str(polio)

# Look at data
cases <- polio$y
plot.ts(cases)

# Many zeros in data set.  I'm not sure the data is zero-inflated though.



# Fit the linear model

fit <- lm(y ~ .,data=polio)
summary(fit)


# Look at residuals
res=resid(fit)
plot(polio$t,res,type="l")
abline(h=0)

qqnorm(res)
qqline(res)

## Test for autocorrelation
acf(res)
pacf(res)
## lag-1 autocorrelation test is significant





##
## Fitting regression with AR1 error covariance via GLS in "nlme"
##

library(nlme)
X=as.matrix(polio[,-1])
X=cbind(1,X)
head(X)

fit.ar1=gls(y~0+X,data=polio,correlation=corAR1(),method="REML")
summary(fit.ar1)
summary(fit)

res.ar1=resid(fit.ar1)
plot(polio$t,res.ar1,type="l")
abline(h=0)

plot(res,res.ar1)
abline(0,1)

###########################################
##
## Now fitting via GLS by hand
##
###########################################

lag.one.res=cbind(res[-length(res)],res[-1])
head(lag.one.res)

cor(lag.one.res)

summary(lm(lag.one.res[,2]~lag.one.res[,1]))

rho.hat=.23548


## make correlation matrix

C.ar1=corAR1(rho.hat)
C.ar1
C.ar1=Initialize(C.ar1,data=polio)
C.ar1
R=corMatrix(C.ar1)
R[1:10,1:10]


## precision matrix
W=solve(R)
W[1:10,1:10]


## take spectral decomposition
E=eigen(W)
P=E$vectors
D=diag(E$values)

## get square root of eigenvectors
D.sqrt=sqrt(D)
D[1:10,1:10]
D.sqrt[1:10,1:10]

W.sqrt=P%*%D.sqrt%*%t(P)

## checking propoerties of W.sqrt
max(abs(W.sqrt%*%W.sqrt-W))
I.tilde=W.sqrt%*%R%*%W.sqrt
I.tilde[1:10,1:10]
diag(I.tilde)

##
## creating y.star, X.star for GLS
##

y.star=W.sqrt%*%polio$y
X.star=W.sqrt%*%X

fit.gls=lm(y.star~0+X.star)
summary(fit.gls)

## compare to OLS fit
summary(fit)

## Different approaches to estimating variance parameters can
##  lead to different results.
## There is considerable current research in estination of
##  variance parameters in different types of models









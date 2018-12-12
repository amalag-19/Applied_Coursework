##
##
##
##   Distribuions.r
##
##     Exploration of distributions through sampling and transformation
##





##
##   A. Exponential Distribution
##

?rexp

## simulate exponential RVs with mean = 3
y=rexp(1000,1/3)
y
mean(y)
var(y)

##
## Plot
##

## make row of three plots
par(mfrow=c(1,3))
## plot the RVs
plot(y)
## plot a histogram of the RVs
hist(y,col="yellow")
## add a line at the mean
abline(v=3,lwd=4)

## try various transformations to make the data approximately normally-distributed
qqnorm(sqrt(y))
qqline(sqrt(y))
hist(sqrt(y))

##
##   B. Gamma Distribution
##

?rgamma

## simulate gamma RVs with mean = 3 and var=1
y=rgamma(1000,9,3)
y
mean(y)
var(y)

##
## Plot
##

## make two rows of two plots
par(mfrow=c(2,2))
## plot the RVs
plot(y)
## plot a histogram of the RVs
hist(y,col="yellow")
## add a line at the mean
abline(v=3,lwd=4)

## try various transformations to make the data approximately normally-distributed
qqnorm(sqrt(y))
qqline(sqrt(y))
hist(sqrt(y))





##
##   C. Geometric Distribution
##

?rgeom

## simulate gamma RVs with mean = 3 
y=rgeom(1000,1/4)
y
mean(y)
var(y)

##
## Plot
##

## make two rows of two plots
par(mfrow=c(2,2))
## plot the RVs
plot(y)
## plot a histogram of the RVs
hist(y,col="yellow")
## add a line at the mean
abline(v=3,lwd=4)

## try various transformations to make the data approximately normally-distributed
qqnorm((y)^(1/3))
qqline((y)^(1/3))
hist((y)^(1/3))





##
##   D. Simulating from the Laplace (double exponential) distribution
##

rlaplace <- function(n,mu=0,b=1){
    ## n = number of replicates
    ## mu = location parameter
    ## b = scale parameter
    v=rexp(n)
    z=rnorm(n)
    mu+b*sqrt(2*v)*z
}

par(mfrow=c(1,1))
y=rlaplace(10000,mu=3,b=2)
hist(y,breaks=100)


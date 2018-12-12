##
##
## QQ plot examples
##
## (interpreting QQ-plots)
##


## Examples:
##
## 1. Normal distribution, large sample size (matches exactly)
## 2. Normal distribution, small sample size (see random variation)
## 3. Truncated Normal (see data with shorter tails)
## 4. T distribution (see data with longer tails)
## 5. Exponential distribution (right/positive skewed)
## 6. Negative Exponential distribution (left/negative skew)

par(mfrow=c(1,2))

###############################################
##
## 1. Normal (large sample size)
##
###############################################

y=rnorm(1000)
hist(y)
qqnorm(y)
qqline(y)



###############################################
##
## 2. Normal (small sample size)
##
###############################################

y=rnorm(10)
hist(y)
qqnorm(y)
qqline(y)


## look at nine examples to see typical variation
par(mfrow=c(3,3))
for(i in 1:9){
    y=rnorm(20)
    qqnorm(y)
    qqline(y)
}


###############################################
##
## 3. Truncated Normal (shorter tails than Normal Distribution)
##
###############################################

## function to simulate from truncated normal

rtnorm <- function(n,mu,sigma,a,b){
    ## n = number of samples
    ## mu= location
    ## sigma = scale
    ## a = lower bound
    ## b = upper bound
    x=rep(NA,n)
    for(i in 1:n){
        accept=0
        while(accept==0){
            ## step 1: draw y~N(mu,sigma^2)
            y=rnorm(1,mu,sigma)
            ## step 2: accept if y is in right region
            if(y>a & y<b){
                x[i]=y
                accept=1
            }
        }
    }
    x
}


##
## 3.1 Truncated only on left
##

par(mfrow=c(1,2))

y=rtnorm(1000,0,1,-2,Inf)
hist(y)
qqnorm(y)
qqline(y)

##
## 3.2 Truncated only on left
##

y=rtnorm(1000,0,1,-Inf,2)
hist(y)
qqnorm(y)
qqline(y)

##
## 3.3 Truncated on both right and left
##

y=rtnorm(1000,0,1,-2,2)
hist(y)
qqnorm(y)
qqline(y)






###############################################
##
## 4. T-distribution (Fatter Tails than Gaussian)
##
###############################################

y=rt(1000,df=3)
hist(y)
qqnorm(y)
qqline(y)

##
## 4.1 try to pull in tails with a transformation (not very successful!)
##

f <- function(x){
    log(abs(x))*sign(x)
}

hist(f(y))
qqnorm(f(y))
qqline(f(y))

g <- function(x){
    sign(x)*abs(x)^-.5
}

hist(g(y))
qqnorm(g(y))
qqline(g(y))


###############################################
##
## 5. Exponential distribution (right skewed)
##
###############################################

y=rexp(1000,1)
hist(y)
qqnorm(y)
qqline(y)

## transform to approximately Gaussian with cube root transform

hist((y)^(1/3))
qqnorm((y)^(1/3))
qqline((y)^(1/3))




###############################################
##
## 6. Negative Exponential distribution (left skewed)
##
###############################################

y=-rexp(1000,1)
hist(y)
qqnorm(y)
qqline(y)





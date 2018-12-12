
##
## (1) Read in data
##
##  Precipitation data for State College PA for the months of April and May 2013
##  Source: http://www.ncdc.noaa.gov/IPS/coop/coop.html
##
##  1 = rain
##  0 = no rain

rain=c(1,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,0)
rain

## plot data over time
plot(rain,type="b",lwd=2)

## plot a histogram of the data
hist(rain,col="yellow")

## calculate the sample mean and variance
mean(rain)
var(rain)

##
##
## (2) Fit Bernoulli model to data
##
##

## (2.1) Find the likelihood of the data given "p"= Prob(y=1)


p=.3

## find the likelihood (pmf) of the first data point under the model:
##    x~binom(1,p)
x=rain[1]
dbinom(x,1,p)

## look at the "dbinom" function
?dbinom

## find the likelihood (pmf) of all data points
dbinom(rain,1,p)
prod(dbinom(rain,1,p))

##
## (2.2) Now we'll try this for a range of values of "p"
##

## first create a sequence of candidate "p" values
p.vals=seq(.01,.99,.01)
n.p=length(p.vals)
p.vals
n.p

## create a placeholder vector to store the likelihood values 
lik.vals=rep(NA,n.p)

## calculate the likelihood function for each "p" value in p.vals
for(i in 1:n.p){
  lik.vals[i] <- prod(dbinom(rain,1,p.vals[i]))
}

## find the "p" which maximizes the likelihood function
max.p.idx=which.max(lik.vals)
max.p=p.vals[max.p.idx]
max.p

## Plot
plot(p.vals,lik.vals)
abline(v=max.p,col="red",lwd=2)


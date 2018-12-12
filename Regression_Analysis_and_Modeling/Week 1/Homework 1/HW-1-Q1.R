y=rbinom(1000,3,2/3)
## make row of three plots
par(mfrow=c(1,2))
## plot the RVs
plot(y)
## plot a histogram of the RVs
hist(y,col="yellow")

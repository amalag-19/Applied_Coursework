#####################################################
##
##   BsplineEx.r
##
##   Examples of building B-spline basis functions
##    and fitting using penalized methods
##
##   Outline:
##   (1) Semiparametric modeling in one dimension using mgcv
##       y = f(x) + e
##   (2) Creating a set of  B-spline basis functions
##        f(x) = sum_i gamma_i * B_i(x)
##   (3) fit unpenalized LM:
##        y = f(x) + e = Z %*% gamma + e
##   (4) fit penalized LM:
##        pen(gamma) = lambda * gamma' K gamma
##
#####################################################




## load packages and data
install.packages("fda")
library(fda)
library(mgcv)

Munich=read.csv("rent99.raw",sep=" ")

#####################################################
##
##   (1) Semiparametric modeling in one dimension using gam
##       y = f(x) + e
##
#####################################################

##
## (1a) rentsqm = f(area)+e
##

fit.area=gam(rentsqm~s(area),data=Munich)
summary(fit.area)
plot(fit.area)

#####################################################
##
##   (2) Creating a set of  B-spline basis functions
##        f(x) = sum_i gamma_i * B_i(x)
##
#####################################################

mn=min(Munich$area)
mx=max(Munich$area)

## note: the "order" of a B-spline is 1 more than the degree
##       => cubic B splines need norder=4

## try the following code for nbasis=7, nbasis=30, nbasis=80
## and for norder=1 through norder=4

B=create.bspline.basis(rangeval=c(mn,mx),nbasis=10,norder=1)
plot(B)

##
## make "Z" matrix
##

Z=eval.basis(Munich$area,B)
str(Z)
image(Z)

##############################################################
##
## (3) fit unpenalized LM:
##      rentsqm = f(area) + e
##              = Z %*% gamma + e
##
##############################################################


fit.lm=lm(Munich$rentsqm~0+Z)
summary(fit.lm)

gamma.hat=fit.lm$coef
pred.vals=seq(mn,mx,by=.1)
Zpred=eval.basis(pred.vals,B)
yhat=Zpred%*%gamma.hat


## plot data
plot(Munich$area,Munich$rentsqm,pch=20,cex=.3)
## plot basis functions
matplot(pred.vals,Zpred,type="l",col="blue",add=T,lwd=2,lty=1)
## plot weighted basis functions
Zpred.weighted=Zpred
for(i in 1:ncol(Zpred.weighted)){
    Zpred.weighted[,i]=Zpred.weighted[,i]*gamma.hat[i]
}
matplot(pred.vals,Zpred.weighted,type="l",col="red",add=T,lwd=2,lty=1)
## add up all weighted basis functions to get fitted line
points(pred.vals,yhat,col="black",type="l",lwd=4)


##############################################################
##
## (4) fit penalized LM:
##      rentsqm = f(area) + e
##              = Z %*% gamma + e
##
##     penalty function:
##        pen(gamma) = lambda * gamma' K gamma
##        where K = D%*%D  and "D" is a matrix of 2nd differences (FKLM p437)  
##############################################################

## matrix of second differences

p=ncol(Z)
D=matrix(0,nrow=p-2,ncol=p)
for(i in 1:(p-2)){
    D[i,i]=1
    D[i,i+1]=-2
    D[i,i+2]=1
}
D


## K
K=t(D)%*%D
dim(K)
K

##
## choose lambda and estimate gamma
##

lambda=10

y=Munich$rentsqm
gamma.hat=solve(t(Z)%*%Z+lambda*K)%*%t(Z)%*%y
gamma.hat

## plot estimated function

f.hat=Zpred%*%gamma.hat
plot(Munich$area,Munich$rentsqm,pch=20,cex=.3)
points(pred.vals,f.hat,col="black",type="l",lwd=4)


##
## Plotting fits for different lambdas
##

par(mfrow=c(2,3))
for(lambda in exp(c(-2,2,4,7,10,20))){
    gamma.hat=solve(t(Z)%*%Z+lambda*K)%*%t(Z)%*%y
    f.hat=Zpred%*%gamma.hat
    plot(Munich$area,Munich$rentsqm,pch=20,cex=.3,main=paste("log(lambda)=",log(lambda)))
    points(pred.vals,f.hat,col="black",type="l",lwd=4)
}



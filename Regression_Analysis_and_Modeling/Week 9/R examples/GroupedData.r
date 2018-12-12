##
##
## Grouped Data example
##
##

##Create four groups

x1=c(1,-2,3)
x2=c(1,0,1)
x3=c(1,4,-4)
x4=c(1,0,0)

n1=3
n2=3
n3=3
n4=5

X=rbind(x1,x2,x3,x4)
X=X[c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4)),]
X

n=nrow(X)

## fix parameters
beta=c(1,-2,3)
sigma=.3

## simulate
y=rnorm(n,mean=X%*%beta,sd=sigma)
y

## OLS regression
summary(lm(y~0+X))

##
## weighted regression
##

## get group means
ybar1=mean(y[1:3])
ybar2=mean(y[4:6])
ybar3=mean(y[7:9])
ybar4=mean(y[10:14])

ybarvec=c(ybar1,ybar2,ybar3,ybar4)
Xbar=rbind(x1,x2,x3,x4)

w=c(n1,n2,n3,n4)

ybarvec
Xbar

## weighted regression
summary(lm(ybarvec~0+Xbar,weights=w))

## compare to full regression
summary(lm(y~0+X))

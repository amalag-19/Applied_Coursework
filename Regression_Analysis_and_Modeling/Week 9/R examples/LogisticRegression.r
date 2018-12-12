#################################################
##
## Logistic Regression
##
#################################################

##
## Read in data on grad school admissions
##
## admit = admission status (1=admitted, 0=not)
## gre = GRE score
## gpa = undergraduate GPA
## rank = "rank" of undergraduate institution
##         (four categories 1=best, 4=worst)

admissions <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
## view the first few rows of the data
head(admissions)
pairs(admissions)

##
## Let's try a linear regression
##

fit.lm=lm(admit~gre+gpa+rank,data=admissions)
summary(fit.lm)

## residual checks
res=resid(fit.lm)
yhat=fit.lm$fitted
plot(yhat,res)

qqnorm(res)
qqline(res)

##
## Now try logistic regression
##

fit=glm(admit~.,family="binomial",data=admissions)
summary(fit)

##
## Let's do logistic regression with ONE predictor to see what is happening
##

fit2=glm(admit~gre,family="binomial",data=admissions)
summary(fit2)

## get estimated linear predictor
xvals=seq(220,800)
newdata=data.frame(gre=xvals)
eta=predict(fit2,newdata=newdata,type="link")

## get estimated mean
par(mfrow=c(1,2))
plot(xvals,eta,main="Linear Predictor",xlab="gre",ylab=expression(eta),type="l")
mu=predict(fit2,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="gre",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(admissions$gre),admissions$admit)


fit2=glm(admit~gpa,family="binomial",data=admissions)
summary(fit2)

## get estimated linear predictor
xvals=seq(2,4,by=.01)
newdata=data.frame(gpa=xvals)
eta=predict(fit2,newdata=newdata,type="link")

## get estimated mean
par(mfrow=c(1,2))
plot(xvals,eta,main="Linear Predictor",xlab="gre",ylab=expression(eta),type="l")
mu=predict(fit2,newdata=newdata,type="response")
plot(xvals,mu,main="Mean Response as a Function of the Predictor",xlab="gre",ylab=expression(mu),ylim=c(0,1)) 
points(jitter(admissions$gpa),admissions$admit)



library(car)
crPlots(fit)

fit=glm(admit~.,family="binomial",data=admissions)
summary(fit)

res=resid(fit)
plot(res)

idx=which(res>0)
admissions$minority=0
admissions$minority[idx] <- 1


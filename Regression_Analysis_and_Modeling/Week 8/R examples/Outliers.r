##################################################################
##
## (1) Teacher Pay
##    "PAY" = Average public school teacher annual salary ($)
##    "SPEND" = Spending on public schools per pupil
##    "AREA" = "1" = Northeast/North Central, "2"=South, "2"=West
##
##################################################################

pay=read.csv("TeacherPay.csv")
pairs(pay)

fit=lm(PAY~SPEND,data=pay)
summary(fit)

plot(pay$SPEND,pay$PAY,xlab="State Spending Per Student",ylab="Avg Teacher Salary")

plot(fit)

r.star=rstudent(fit)

r.star
n=nrow(pay)
p=2


## plotting outlier test:

max(r.star)
plot(r.star)
cutoff=qt(.975,df=n-p-1)
abline(h=cutoff)
abline(h=-cutoff)

## hat matrix diagonals

h=hatvalues(fit)

cutoff.h=2*p/n
cutoff.h
plot(h)
abline(h=cutoff)


## Cook's Distance

cd=cooks.distance(fit)
cd

plot(cd)


##
## Conclusion: Outliers are not very influential.  Probably fine to leave them in.  
## For comparison, below is the model without outliers
##

outlier.idx=which(abs(r.star)>cutoff)
outlier.idx

## model with outliers
fit=lm(PAY~SPEND,data=pay)
summary(fit)
## model without outliers
fit.no.outliers=lm(PAY~SPEND,data=pay[-outlier.idx,])
summary(fit.no.outliers)


##
##
## (2) Ad spending
##
##

ad.spending=read.csv("AdSpending.csv")
pairs(ad.spending)
summary(ad.spending)
str(ad.spending)

fit=lm(ad.spending$MILIMP~(ad.spending$SPEND))
summary(fit)

plot(ad.spending$SPEND,ad.spending$MILIMP)
abline(fit)


r.star=rstudent(fit)
yhat=fit$fitted

## slightly longer tails than nornmal
qqnorm(r.star)
qqline(r.star)

## residuals vs fitted valuese
plot(yhat,r.star)
abline(h=0)


n=nrow(ad.spending)
p=2


## plotting outlier test:

max(r.star)
plot(r.star)
cutoff=qt(.975,df=n-p-1)
abline(h=cutoff)
abline(h=-cutoff)

## hat matrix diagonals

h=hatvalues(fit)

cutoff.h=2*p/n
cutoff.h
plot(h)
abline(h=cutoff.h)


## Cook's Distance

cd=cooks.distance(fit)
cd

plot(cd)


## find the potentially-influential observations
which(cd>.5)
which(abs(h)>cutoff.h)

ad.spending[c(7,10,13),]


## plotting potential outlier
plot(ad.spending$SPEND,ad.spending$MILIMP)
points(ad.spending$SPEND[10],ad.spending$MILIMP[10],col="red",pch=8,cex=2,lwd=2)

## fitting model without this point
fit.no.FORD=lm(MILIMP~SPEND,data=ad.spending[-10,])
summary(fit)
summary(fit.no.FORD)

plot(ad.spending$SPEND,ad.spending$MILIMP)
points(ad.spending$SPEND[10],ad.spending$MILIMP[10],col="red",pch=3,cex=3,lwd=2)
abline(fit)
abline(fit.no.FORD,col="red")

## residual analysis
res=resid(fit.no.FORD)
yhat=fit.no.FORD$fitted

plot(yhat,res)
abline(h=0)

qqnorm(res)
qqline(res)


## checking to see how bad this qq-plot actually is
par(mfrow=c(3,3))
for(i in 1:9){
    sim=rnorm(length(res))
    qqnorm(sim)
    qqline(sim)
}

##
## Conclusion: outlier (FORD) is influential.  
##  However: (a) qualitative results don't change (positive correlation between
##               spend and milimp)
##           (b) residuals are probably ok
##           (c) might want to consider quadratic trend (not significant - see below)
##           (d) report slopes from both model with and model without FORD


fit.quad=lm(MILIMP~SPEND+I(SPEND^2),data=ad.spending[-10,])
summary(fit.quad)
## quadratic term is not significant.


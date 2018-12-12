fluid=read.csv("fluid.csv",sep=",")
#str(fluid)
t<-fluid$t
x<-fluid$x
x5=0
n5=0
x10=0
n10=0
x15=0
n15=0
x20=0
n20=0
j=1
y=rep(NA,200)
for (i in 1:length(t)) {
  if (t[i]==5) {
    y[j]<-x[i]
    x5<-x5+(x[i])^2
    n5<-n5+1
    j<-j+1
  }
  if (t[i]==10) {
    x10<-x10+(x[i])^2
    n10<-n10+1
  }
  if (t[i]==15) {
    x15<-x15+(x[i])^2
    n15<-n15+1
  }
  if (t[i]==20) {
    x20<-x20+(x[i])^2
    n20<-n20+1
  }
}
x5.mean<-x5/n5
x10.mean<-x10/n10
x15.mean<-x15/n15
x20.mean<-x20/n20
xsq.mean=c(x5.mean,x10.mean,x15.mean,x20.mean)
tm=c(5,10,15,20)
xsq.mean


# answer 4
oly=read.csv("Olympics.csv")
str(oly)
#fit<-lm(goldtime~year+gender,data=oly)
#fit<-lm(goldtime~year+gender+year*gender,data=oly)
#fit<-lm(goldtime~I(1/year)+gender,data=oly)
fit<-lm(goldtime~I(1/year)+(year*gender),data=oly)
#fit<-lm(goldtime~I(1/year)+gender+(year*gender),data=oly)
#fit<-lm(goldtime~I(sqrt(year))+gender,data=oly)
#fit<-lm(goldtime~I(sqrt(year))+(year*gender),data=oly)
#fit<-lm(goldtime~I(sqrt(year))+gender+(year*gender),data=oly)
#fit<-lm(goldtime~I((year)^(1/3))+gender,data=oly)
#fit<-lm(goldtime~I((year)^(1/3))+(year*gender),data=oly)
#fit<-lm(goldtime~I((year)^(1/3))+gender+(year*gender),data=oly)
#fit<-lm(goldtime~poly(year,3)+I(1/year)+gender+year*gender,data=oly)
#fit<-lm(goldtime~I(1/year)+(poly(year,2)*gender),data=oly)

summary(fit)
summary(fit)$sigma

par(mfrow=c(1,1))
plot(oly$year,oly$goldtime,pch=as.integer(oly$gender))

x.vals=seq(from=1900,to=2010,by=1)

gender.vals=rep("W",length(x.vals))
df.0=data.frame(year=x.vals,gender=gender.vals)

gender.vals=rep("M",length(x.vals))
df.1=data.frame(year=x.vals,gender=gender.vals)

## look at the two data.frames
head(df.0)
head(df.1)

f.vals.0=predict(fit,newdata=df.0)

f.vals.1=predict(fit,newdata=df.1)

points(x.vals,f.vals.0,type="l",col="red",lwd=3)

points(x.vals,f.vals.1,type="l",col="blue",lwd=3)

par(mfrow=c(2,2))
plot(fit)

library(car)
crPlots(fit)

for (i in 1:length(x.vals)) {
  if (x.vals[i]==1944) {
    Women.1944.predicted<-f.vals.0[i]
    Men.1944.predicted<-f.vals.1[i]
  }
}
Women.1944.predicted
Men.1944.predicted


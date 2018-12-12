library(mgcv)
##Q1
mist=read.csv("mistletoe.csv",sep=",")
head(mist)
## fit using gam
fit.gm=gam(infected.mndnr~mortal+phys+si+usize+height+dbh+s(x,y),data=mist,family=binomial)
summary(fit.gm)
etahat.gm=predict(fit.gm)
muhat.gm=predict(fit.gm,type="response")
muhat.gm
plot(fit.gm,scheme=1)
plot(fit.gm,scheme=2)
## fit using linear predictor
fit.lp=glm(infected.mndnr~mortal+phys+si+usize+height+dbh+x+y,data=mist,family=binomial)
summary(fit.lp)
etahat.lp=predict(fit.lp)
muhat.lp=predict(fit.lp,type="response")
muhat.lp
hrisk=rep(0,length(mist$infected.mndnr))
count=0
for (i in 1:length(mist$infected.mndnr)){
  if (muhat.gm[i]>muhat.lp[i]){
    hrisk.idx[i]=1
    count<-count+1
  }
}
count

##Q2
admissions <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
head(admissions)
fit1=gam(admit~rank+s(gre)+s(gpa),data=admissions,family=binomial)
summary(fit1)
fit2=gam(admit~rank+gre+s(gpa),data=admissions,family=binomial)
summary(fit2)
AIC(fit1)
AIC(fit2)
plot(fit1)

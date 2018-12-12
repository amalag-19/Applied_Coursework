ice=read.csv("icecream.csv",sep=",")
ice<-ice[-length(ice$date),]
ice
pairs(ice)
plot.ts(ice$IC)
library(nlme)
X=cbind(1,ice$temp)
X
n=29
p=2
rho.hat=0.7557313
C.ar1=corAR1(rho.hat)
C.ar1
C.ar1=Initialize(C.ar1,data=ice)
C.ar1
R=corMatrix(C.ar1)
R[1:3,1:3]

W=solve(R[1:29,1:29]/(1-(rho.hat)^2))
Y<-t(X)%*%W%*%X
Y
F<-((0.002603578^2)*Y[2,2])/(0.001708249)

sqrt(F)
2*(1-pt(sqrt(F),27))
2*(pt(-4.382617,27))

ice=read.csv("icecream.csv",sep=",")
ice
newdata=ice[30,]
newdata
predict(fit3, newdata)

H=X%*%solve(Y)%*%t(X)%*%W
h29=H[29,29]
X30=matrix(c(1,ice$temp[30]),nrow=1,ncol=2)
cf=(1/(1-rho.hat^2))+(rho.hat^2*(1-h29))+(X30%*%solve(Y)%*%t(X30))
s=sqrt(fit3$sigma^2*cf)
U=Y30_cap+s*qt(0.975,(n-p))
L=Y30_cap-s*qt(0.975,(n-p))




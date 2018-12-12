data(anscombe)
x1<-anscombe$x1
x2<-anscombe$x2
x3<-anscombe$x3
x4<-anscombe$x4
y1<-anscombe$y1
y2<-anscombe$y2
y3<-anscombe$y3
y4<-anscombe$y4
n<-length(x1)
# fitting the simple linear model for the pairs (x_i,y_i), i=1,2,3,4
fit1=lm(y1~x1)
fit2=lm(y2~x2)
fit3=lm(y3~x3)
fit4=lm(y4~x4)
fit1$coeff
fit2$coeff
fit3$coeff
fit4$coeff
p<-length(fit1$coeff)

## (x1,y1)
plot(x1,y1)
abline(fit1)

# level of significance
alpha<-0.05
# observation numbers of the outliers
out1<-rep(0,length(x1))
rst1<-rstudent(fit1)
j=1
for (i in 1:length(x1)){
  if (rst1[i]>qt((1-(alpha/(2*length(x1)))),(n-p-1))){
    out1[j]<-i
    j<-j+1
  }
}
influence.measures(fit1)


hii1<-hatvalues(fit1)
cd1<-cooks.distance(fit1)
inf1.1<-rep(0,length(x1))
inf1.2<-rep(0,length(x1))
j=1
for (i in 1:length(x1)){
  if (hii1[i]>(2*p/n)){
    inf1.1[j]<-i
    j<-j+1
  }
}
j=1
for (i in 1:length(x1)){
  if (cd1[i]>(4/n)){
    inf1.2[j]<-i
    j<-j+1
  }
}
inf1.1
inf1.2

## (x2,y2)
plot(x2,y2)
plot(log(x2),y2)
fit2.2=lm(y2~x2+I(x2^2))
plot(x2,y2)
vals<-seq(from=3,15,by=0.1)
vals
f.vals<-fit2.2$coeff[1]+fit2.2$coeff[2]*vals+fit2.2$coeff[3]*(vals^2)
f.vals
lines(vals,f.vals)
length(x2)

rst2<-rstudent(fit2.2)
out2<-rep(0,length(x2))
rst2<-rstudent(fit2.2)
j=1
for (i in 1:length(x2)){
  if (rst2[i]>qt((1-(alpha/(2*length(x2)))),(n-p-1))){
    out2[j]<-i
    j<-j+1
  }
}
out2

rst4<-rstudent(fit4)
j=1
for (i in 1:length(x4)){
  if (rst4[i]>qt((1-(alpha/(2*length(x4)))),(n-p-1))){
    out4[j]<-i
    j<-j+1
  }
}
out4

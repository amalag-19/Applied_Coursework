fuel=read.csv("FuelEfficiency.csv")
str(fuel)
head(fuel)
pairs(fuel)

fit=lm(MPG~Weight+Horsepower,data=fuel)
summary(fit)
plot(fit)
r.star=rstudent(fit)
r.star

n=nrow(fuel)
p=3

max(r.star)
plot(r.star)
cutoff=qt(.975,df=n-p-1)
abline(h=cutoff)
abline(h=-cutoff)


library(car)

fit=lm(MPG~Weight+Horsepower+Drive_Ratio+Displacement+Cylinders,data=fuel)
vif(fit)
summary(fit)

## Q1 (a)
## To check the code, please install the following packages
library(ggplot2)
library(plyr)
library(reshape)

## length of the vector to be sorted
N<-c(20,200,2000,20000,200000,2000000,20000000)
## Number of iterations
M<-1000
## tolerance
tol<-0.3

## creating an empty list
sys<-list()

## Generating samples and storing up the sample times for different iteratons.
for (m in 1:M){
  samples<-runif(N[length(N)])
  sys[[m]]<-list()
  for (n in 1:length(N)){
    sys.q<-system.time(sort(samples[1:N[n]],method="quick"))
    sys.s<-system.time(sort(samples[1:N[n]],method="shell"))
    temp<-rbind(sys.q,sys.s)
    temp<-data.frame(temp,"method"=c("quick","shell"),"vector.length"=N[n],"iteration"=m)
    sys[[m]][[n]]<-temp
  }
  if (m%%100==0){
    save(sys,file="data1.Rdata")
    print(m)
  }
}

## Appending the list with cumulative means and variances
a<-do.call(rbind,do.call(rbind,sys))
a<-a[-which(a$iteration==469),] ## removing one bad observation since my computer slept at this iteration
a$iteration[which(a$iteration>468)]<-a$iteration[which(a$iteration>468)]-1

data<-split(a,a$method)
MCse<-lapply(data,function(x){split(x,x$vector.length)})
s.mean.var<-lapply(MCse,function(x){
  lapply(x,function(y){
    y$cumean<-cumsum(y$elapsed)/y$iteration
    y$cuvar<-(cumsum((y$elapsed)^2)-(y$iteration*(y$cumean)^2))/(y$iteration-1)
    y$M.star<-((2*1.96*sqrt(y$cuvar))/(tol*y$cumean))^2
    return(y)
    })
  })
f<-do.call(rbind,do.call(rbind,s.mean.var))
save(f,file="data2.Rdata")

M.star.q<-f$M.star[which(f$iteration==999&f$method=="quick")]
M.star.s<-f$M.star[which(f$iteration==999&f$method=="shell")]
M.star.q
M.star.s

## (b)
f$cumean[which(f$iteration==999&f$method=="quick")]
f$cuvar[which(f$iteration==999&f$method=="quick")]
f$cumean[which(f$iteration==999&f$method=="shell")]
f$cuvar[which(f$iteration==999&f$method=="shell")]

## plotting
p<-ggplot(data=f)

thema<-theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 8, colour = "black"), 
        axis.text.x  = element_text(angle = 0, size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"), 
        axis.text.y  = element_text(angle = 0, size = 8, colour = "black"),
        legend.text  = element_text(size = 8, colour = "black"), 
        legend.title = element_text(size = 8, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white", linetype = NULL),
        panel.grid.minor = element_line(colour = "white", linetype = NULL),
        text = element_text(size = 8, colour = "black"),
        title =  element_text(size = 8, face = "bold"))

## plotting mean vs. no. of iterations for different n and methods
p+geom_line(mapping = aes(x = iteration,y = log(cumean),colour = factor(vector.length),linetype = method))+labs(x="Number of iterations",y="Log of cumulative sample mean (MC estimates of time complexity)",colour="Vector length",linetype="Sort method")+thema

## plotting variance vs. no. of iterations for different n and methods
p+geom_line(mapping = aes(x = iteration,y = log(cuvar),colour = factor(vector.length),linetype = method))+
  labs(x="Number of iterations",y="Log of cumulative sample variance",colour="Vector length",linetype="Sort method")+thema

## plotting no. of iterations required vs. different vector lengths for both methods
df<-data.frame(N,M.star.q,M.star.s)
names(df)<-c("N","quick","shell")
df<-melt.data.frame(data = df,id.vars = "N")

q<-ggplot(data=df)
q+geom_line(mapping = aes(x = factor(N),y = value,group=variable,colour=variable))+
  geom_point(mapping = aes(x = factor(N),y = value,group=variable,colour=variable))+
  labs(x="Vector Length",y="Number of optimal iterations required",colour="method")+thema

## (c) 
## plotting mean vs. n with error bars for both methods
p+geom_errorbar(mapping = aes(x=factor(vector.length),y=cumean,ymin=cumean-1.96*cuvar/iteration,ymax=cumean+1.96*cuvar/iteration,group=method,colour=method),subset=.(iteration==999))+
  geom_line(mapping = aes(x=factor(vector.length),y=cumean,group=method,colour=method),subset=.(iteration==999))+labs(x="Length of the vector",y="Computational cost",colour="Sort method")+thema

## plotting log(mean) vs. n with error bars for both methods
p+geom_errorbar(mapping = aes(x=factor(vector.length),y=log(cumean),ymin=log(cumean-1.96*cuvar/iteration),ymax=log(cumean+1.96*cuvar/iteration),group=method,colour=method),subset=.(iteration==999))+
  geom_line(mapping = aes(x=factor(vector.length),y=log(cumean),group=method,colour=method),subset=.(iteration==999))+
  labs(x="Length of the vector",y="Computational cost",colour="Sort method")+thema

## plotting mean vs. no. of iterations for different n and methods
p+geom_line(mapping = aes(x = iteration,y = cumean,colour = method,group=method))+facet_grid(vector.length~.,scales="free")+
  labs(x="Number of iterations",y="Computational cost",colour="Sort method")+thema

## (e)
ratio<-function(df){
  rat<-rep(NA,nrow(df))
  for (i in 1:nrow(df)){
    if (df$method[i]=="quick"){
      rat[i]<-df$cumean[i]/(df$vector.length[i]*log(df$vector.length[i]))
    }
    else if (df$method[i]=="shell"){
      rat[i]<-df$cumean[i]/(df$vector.length[i]*(log(df$vector.length[i]))^2)
    }
  }
  df[, "rat"] <- rat
  return(df)
}
f.r<-ratio(f)
r<-ggplot(data=f.r)
r+geom_point(mapping = aes(x=factor(vector.length),y=rat,group=method,colour=method),subset=.(iteration==999))+
  geom_line(mapping = aes(x=factor(vector.length),y=rat,group=method,colour=method),subset=.(iteration==999))+labs(x="Length of the vector (n)",y=expression("Computational cost divided by n log (n)"),colour="Sort method")+thema

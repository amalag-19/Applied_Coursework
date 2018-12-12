## Please note that each question is seperated by a full line of '#'
## and subparts are seperated by half lines of '#'. Also run each 
## question from the begining as a whole to get the correct output
## since variable names/function names may overlap for different
## questions.

## Please install the necessary packages mentioned in the libraries beforehand
#####################################################################

## libraries
{library(ggplot2)
  library(tmvtnorm)
  library(astsa)
  library(plyr)
  library(mvtnorm)
  library(Matrix)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }}

## Q1 (b)
## defining the parameters
{mu<-3
sigma<-0.8
zeta<-0.4}

## defining square root of target kernel function
sqrt.pr<-function(x){
  t1<-(1+zeta*((x-mu)/sigma))
  t2<-exp(-(t1^(-1/zeta)))
  t3<-t1^(-1-(1/zeta))
  value<-(t2*t3)^(1/2)
  return(-value)
}

## defining x times square root of target kernel function
x.sqrt.pr<-function(x){
  t1<-(1+zeta*((x-mu)/sigma))
  t2<-exp(-(t1^(-1/zeta)))
  t3<-t1^(-1-(1/zeta))
  value<-x*((t2*t3)^(1/2))
  return(value)
}

## defining negative of x times square root of target kernel function
x.sqrt.pr2<-function(x){
  t1<-(1+zeta*((x-mu)/sigma))
  t2<-exp(-(t1^(-1/zeta)))
  t3<-t1^(-1-(1/zeta))
  value<-x*((t2*t3)^(1/2))
  return(-value)
}

## defining log of target kernel function
log.pr<-function(x){
  if (x>=1){
    t1<-(1+zeta*((x-mu)/sigma))
    t2<--(t1^(-1/zeta))
    t3<-(-1-(1/zeta))*log(t1)
    value<-t2+t3
  }
  else {value<--Inf}
  return(value)
}

## calculating the bounds b,c,d
b<--optim(par = 1.5,fn = sqrt.pr,method = "L-BFGS-B",lower = 1, upper = Inf)$value
c<-optim(par = 1.5,fn = x.sqrt.pr,method = "L-BFGS-B",lower = 1, upper = Inf)$value
d<--optim(par = 1.5,fn = x.sqrt.pr2,method = "L-BFGS-B",lower = 1, upper = Inf)$value

## defining the ratio of uniforms function for GEV
ROU.GEV<-function(n){
  u.f<-rep(NA_real_,n)
  v.f<-rep(NA_real_,n)
  i<-1
  j<-0
  while(i<=n){
    u<-runif(1, min=c, max=d)
    v<-runif(1,min=0,max=b)
    if (log(v)<(0.5*log.pr(u/v))){
      u.f[i]<-u
      v.f[i]<-v
      i<-i+1
    }
    else {j<-j+1}
  }
  samples<-u.f/v.f
  return(list(u.f,v.f,samples))
}

## defining the function expectation.trial to calculate the required Monte Carlo sample 
## size
Expectation.trial<-function(n){
  s<-ROU.GEV(n)
  iterations<-1:n
  ## Initializing a matrix R for storing Monte Carlo estimates and MC
  ## standard errors for samples
  R<-matrix(NA_real_,n,2)
  ## loop for storing Monte Carlo estimates and MC standard errors in R
  for (j in 1:n){
    R[j,1]<-mean(s[[3]][1:iterations[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,2]<-sqrt(var(s[[3]][1:iterations[j]])/iterations[[j]])
  }
  temp<-(R[,2]/R[,1])<0.05
  n.star<-max(which(temp==FALSE))+1
  return(n.star)
}

## choosing the optimum sample size according to 5% tolerance
n<-10000
n.star<-rep(NA_integer_,100)

for (i in 1:100){
  n.star[i]<-Expectation.trial(n)
}
n.opt<-mean(n.star)

## defining the function expectation to calculate Monte Carlo estimates of expectation
## and MonteCarlo standard errors
Expectation<-function(n){
  t<-system.time(
    s<-ROU.GEV(n)
  ) 
  sampling.rate<-n/as.numeric(t[3])
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix R for storing Monte Carlo estimates and MC
  ## standard errors for x
  R<-matrix(NA_real_,ml,2)
  ## loop for storing Monte Carlo estimates and MC standard errors in R
  for (j in 1:ml){
    R[j,1]<-mean(s[[3]][1:m[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,2]<-sqrt(var(s[[3]][1:m[j]])/m[[j]])
  }
  return(list(s,R,sampling.rate))
}

## defining the sample size as the optimum chosen as above
n<-100000

## For optimal sample size, generating samples from Ratio of uniforms, estimates and MCse
s<-Expectation(n)

############################################
## Q1 (c)
## creating vectors to plot true target kernel
x.vec<-seq(1,15,length.out = n)
y.vec<-rep(NA_real_,length(x.vec))
for (i in 1:length(x.vec)){
  y.vec[i]<-(1/sigma)*exp(log.pr(x.vec[i]))
}

## defining preliminaries for plotting using ggplot
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

f<-data.frame("samples"=s[[1]][[3]],"u"=s[[1]][[1]],"v"=s[[1]][[2]],"x.vec"=x.vec,"y.vec"=y.vec)
p<-ggplot(data=f)

## plotting the histogram of samples and comparing with true target kernel
p+ geom_histogram(mapping = aes(x=samples,y=..density..),fill="steelblue",subset=.(samples<15))+labs(x="x",y="Density")+thema+geom_line(mapping = aes(x = x.vec,y = y.vec))

## plotting the ratio of unifroms region
p+geom_point(mapping = aes(x=u,y=v),col="steelblue")

############################################
## Q1(d)
## Approximate expected value
s[[2]][nrow(s[[2]]),1]

## plotting the estimates vs. sample size
g<-data.frame("sample.size"=seq(100,n,by=100),"estimates"=s[[2]][,1],"MCse"=s[[2]][,2])
q<-ggplot(data=g)
q+geom_line(mapping = aes(x=sample.size,y=estimates),col="steelblue")

## plotting the Monte carlo standard error vs. sample size
q+geom_line(mapping = aes(x=sample.size,y=MCse),col="steelblue")

############################################
## Q1 (e)
sampling.rate<-s[[3]]

########################################################################################
## Q2 (a) and (b)

## libraries
{library(ggplot2)
  library(tmvtnorm)
  library(astsa)
  library(plyr)
  library(mvtnorm)
  library(Matrix)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }}

## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

## defining the parameters
{mu<-3
sigma<-0.8
zeta<-0.4
a<-mu-(sigma/zeta)}

## defining the posterior kernel on log scale
log.pr<-function(x){
  t1<-(1+zeta*((x-mu)/sigma))
  t2<--(t1^(-1/zeta))
  t3<-(-1-(1/zeta))*log(t1)
  value<-t2+t3
  return(value)
}

## defining the proposal kernel on log scale
log.prop<-function(x,m,s,a){
  t1<-dnorm(x, mean=m, sd=s,log = TRUE)
  t2<-log(1-pnorm(a, mean=m, sd=s))
  value<-t1-t2
  return(value)
}

## Metropolis Hastings function with inputs as the variance of proposal
## (tuning parameter) and current state of the MC
MH<-function(variance,current.state){
  ## sampling x.star from normal proposal with mean as current state
  ## of x and specified tuning parameter
  x.star<-rtmvnorm(n=1, mean=current.state, sigma=variance, lower=mu-(sigma/zeta))
  ## defining the acceptance probability for sampled x.star on log
  ## scale
  accept.probab<-(log.pr(x.star)+log.prop(current.state,x.star,sqrt(variance),a))-(log.pr(current.state)+log.prop(x.star,current.state,sqrt(variance),a))
  ## sampling u from uniform(0,1) to check for acceptance
  u<-runif(1, min=0, max=1)
  ## initializing the indicator flag=0 to check if the sampled x.star
  ## will be accepted
  flag<-0
  ## if-else to define the next state of the chain based on acceptance
  ## probability
  if(log(u)<=accept.probab){
    flag<-1
    next.state<-x.star
  }
  else {next.state<-current.state}
  ## returning the next state and indicator if the sampled value was
  ## accepted
  return(c(next.state,flag))
}

## defining the function Expectation.trial for running Metropolis-Hastings
## algortithm and to calculate optimum sample size for 5% relative error
Expectation.trial<-function(n,start,variance){
  ## Initializing the Markov chain
  x<-rep(NA_real_,n)
  ## Defining the initial value for the chain
  x[1]<-start
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for RWM updates
  for(i in 1:(n-1)){
    temp<-MH(variance,x[i])
    x[i+1]<-temp[1]
    accept<-accept+temp[2]
  }
  iter<-1001:n
  ## Initializing a matrix R for storing Monte Carlo estimates and MCMC
  ## standard errors for x
  m<-n-1000
  R<-matrix(NA_real_,m,2)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:m){
    R[j,1]<-mean(x[1:iter[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,2]<-bm(x[1:iter[j]])$se
  }
  temp<-(R[,2]/R[,1])<0.05
  n.star<-min(which(temp==TRUE))+1000
  return(n.star)
}
## choosing the optimum sample size according to 5% tolerance

n<-10000
n.star<-rep(NA_integer_,10)
start<-3.5
variance<-9.4

for (i in 1:10){
  n.star[i]<-Expectation.trial(n,start,variance)
}
n.opt<-mean(n.star)
n.opt

## defining the function Expectation for running Metropolis-Hastings
## algortithm and to calculate Monte Carlo estimates and standard errors
## for x
Expectation<-function(n,start,variance){
  ## Initializing the Markov chain
  x<-rep(NA_real_,n)
  ## Defining the initial value for the chain
  x[1]<-start
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for RWM updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-MH(variance,x[i])
    x[i+1]<-temp[1]
    accept<-accept+temp[2]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"samples"=x)
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix R for storing Monte Carlo estimates and MCMC
  ## standard errors for x
  R<-matrix(NA_real_,ml,2)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    R[j,1]<-mean(x[1:m[j]])
    ## Using the batchmeans function to calculate MCMC standard errors
    R[j,2]<-bm(x[1:m[j]])$se
  }
  n.star<-m[which.min((R[,2]/R[,1])<0.5)]
  ESS<-ess(samples[,2])
  ESS.rate<-ESS/as.numeric(sys.time[3])
  ## returning the sampled values of the chain, estimates & standard
  ## errors of x at varying number of iterations and acceptance rate
  ## of the chain
  return(list("samples"=samples,"iterations"=m,"R"=R,"acceptance.rate"=acceptance.rate,"n.star"=n.star,"ESS"=ESS,"ESS.rate"=ESS.rate))
}

## defining the sample size
n<-100000

## defining the tuning parameter (chosen after many pilot runs)
var.tune<-9.4

## Running the MH algorithm for 3 times with different starting values
start<-c(3.5,3.7,3.9)
## creating a vector of empty lists
result<-vector("list", length(start)) 
t<-system.time(
for (j in 1:length(start)){
  result[[j]]<-Expectation(n,start[j],var.tune)
})

## estimates
result[[1]][[3]][50,1]
result[[2]][[3]][50,1]
result[[3]][[3]][50,1]

## verifying the tuning parameter using acf
acf(result[[1]][[1]][,2],lag.max = 50,main="Plot of ACF of samples")
acf(result[[2]][[1]][,2],lag.max = 50,main="Plot of ACF of samples")
acf(result[[3]][[1]][,2],lag.max = 50,main="Plot of ACF of samples")

## verifying ESS
result[[1]][[6]]/n
result[[2]][[6]]/n
result[[3]][[6]]/n

## verifying ESS rate
result[[1]][[7]]
result[[2]][[7]]
result[[3]][[7]]

## new starting values
result[[1]][[1]][n,2]
result[[2]][[1]][n,2]
result[[3]][[1]][n,2]

##########################################
## Q2 (c)
## plotting
f<-list()
g<-list()
for(j in 1:length(start)){
  f[[j]]<-data.frame(result[[j]][[2]],result[[j]][[3]],start[j])
  g[[j]]<-data.frame(result[[j]][[1]],start[j])
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean","MCMCse","start")

g<-do.call(rbind,g)
names(g)<-c("iterations","samples","start")

p<-ggplot(data=f)
q<-ggplot(data=g)
head(g)

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

## plotting mean vs. sample size for different starting values
p+ geom_line(mapping = aes(x=iterations,y=mean,group=factor(start),colour=factor(start)))+labs(x="Number of samples",y="Estimate of the expectation with MCMC standard errors", colour="Starting values")+thema+
geom_errorbar(mapping = aes(x=iterations,y=mean,ymin=mean-MCMCse,ymax=mean+MCMCse,group=factor(start),colour=factor(start)))

## plotting MCMCse vs. sample size for different starting values
p+ geom_line(mapping = aes(x=iterations,y=MCMCse,group=factor(start),colour=factor(start)))+labs(x="Number of samples",y="MCMC standard errors",colour="Starting values")+thema

## trace plot
q+geom_line(mapping = aes(x=iterations,y=samples,group=factor(start),colour=factor(start)),subset=.(start==g[1,3]))+labs(x="",y="MCMC standard errors",colour="Starting values")+thema

## estimated density after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples),fill="steelblue",subset=.(start==g[1,3]&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="Estimated density from the samples",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples),fill="tomato",subset=.(start==g[1,3]))+
  labs(x="",y="Estimated density from the samples",title="After N")+thema

multiplot(q1, q2, cols=2)

## estimated density for different starting values
q+ geom_density(mapping = aes(x=samples,group=factor(start),colour=factor(start)))+
  labs(x="",y="Estimated density from the samples",colour="Starting values")+thema

##########################################
## Q2 (d)
result[[1]][[7]]

#######################################################################################
## Q3 (b)

## libraries
{library(ggplot2)
  library(tmvtnorm)
  library(astsa)
  library(plyr)
  library(mvtnorm)
  library(Matrix)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }}

## defining the parameters for the given target multivariate truncated normal
mu<-c(0,1,0)
sig<-matrix(c(1,0.8,0.3,0.8,2,0.4,0.3,0.4,3),3,3)
a<-c(-1,2,3)
b<-c(1,4,5)

## defining the target density function
target<-function(x,mu,sig,a,b){
  m<-length(mu)
  t1<-pmvnorm(a,b,mean=mu,sigma=sig)
  t2<-dmvnorm(x = x,mean = mu,sigma = sig)
  value<-t2/t1
  return(value[1])
}

## defining the function for inverse transformation to generate samples from univariate
## truncated normal distribution
inv.tnorm<-function(x,mu,sig2,a,b){
  t1<-pnorm(a,mean=mu,sd=sqrt(sig2))
  t2<-pnorm(b,mean=mu,sd=sqrt(sig2))
  value<-qnorm(x*(t2-t1)+t1,mean=mu,sd=sqrt(sig2))
  return(value)
}

## defining the proposal density function
prop<-function(x,mu,sig,a,b){
  m<-length(mu)
  q<-rep(NA_real_,m)
  for(i in 1:m){
    q[i]<-dnorm(x[i],mu[i],sqrt(sig[i,i]))/(pnorm(b[i],mean=mu[i],sd=sqrt(sig[i,i]))-pnorm(a[i],mean=mu[i],sd=sqrt(sig[i,i]))) 
  }
  value<-prod(q)
  return(value)
}

## finding K using optim
fq<-function(x,mu,sig,a,b){
  value<-target(x,mu,sig,a,b)/prop(x,mu,sig,a,b)
  return(-value)
}

K<--optim(c(0,3,4),fq,mu=mu,sig=sig,a=a,b=b,method="L-BFGS-B",lower=a,upper=b)$value
K

## defining the function sampler for sampling via accept-reject algorithm 
sampler<-function(mu,sig,a,b){
  m<-length(mu)
  u<-runif(m)
  x<-rep(NA_real_,m)
  for (i in 1:m){
    x[i]<-inv.tnorm(u[i],mu[i],sig[i,i],a[i],b[i])
  }
  q<-prop(x,mu,sig,a,b)
  sample<-rep(NA_real_,m)
  ## sampling u from uniform(0,1) to check for acceptance
  u<-runif(1, min=0, max=1)
  ## initializing the indicator flag=0 to check if the sampled x.star
  ## will be accepted
  flag<-0
  ## if-else to define the next state of the chain based on acceptance
  ## probability
  if(u<=target(x,mu,sig,a,b)/(K*q)){
    flag<-1
    sample<-x
  }
  return(list(sample,flag))
}

## defining the function Expectation for calculating Monte Carlo estimates and standard
## errors
Expectation<-function(n){
  samples<-matrix(NA_real_,3,n)
  k<-0
  j<-1
  sys.t<-system.time(while(j<=n){
    t<-sampler(mu,sig,a,b)
    if(t[[2]]==1){
      samples[,j]<-t[[1]]
      j<-j+1
    }
    else{
      k<-k+1
    }
  })
  sampling.rate<-n/as.numeric(sys.t[3])
  acceptance.rate<-n/(n+k)
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MC.means for storing Monte Carlo estimates
  MC.means<-matrix(NA_real_,3,ml)
  ## Initializing a matrix MC.se for storing Monte Carlo standard errors
  MC.se<-matrix(NA_real_,3,ml)
  ## loop for storing Monte Carlo estimates
  for (j in 1:ml){
    MC.means[,j]<-(rowSums(samples[,1:m[j]]))/m[j]
    MC.se[,j]<-sqrt(((rowSums((samples[,1:m[j]]-MC.means[,j])^2))/(m[j]-1))/m[[j]])
  }
  return(list(samples,MC.means,MC.se,sampling.rate,acceptance.rate))
}

n<-10000
M<-Expectation(n)

## estimate of the expected value
M[[2]][,nrow(M[[2]])]

## Monte Carlo standard error
M[[3]][,nrow(M[[3]])]

## defining the data frame for plotting
f<-data.frame("sample.size"=seq(100,n,by=100),t(M[[2]]),t(M[[3]]))
names(f)<-c("sample.size","X1.m","X2.m","X3.m","X1.se","X2.se","X3.se")
p<-ggplot(data=f)

## plotting the estimates vs. sample size
p1<-p+geom_line(mapping = aes(x=sample.size,y=X1.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X1")
p2<-p+geom_line(mapping = aes(x=sample.size,y=X2.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X2")
p3<-p+geom_line(mapping = aes(x=sample.size,y=X3.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X3")
p4<-p+geom_line(mapping = aes(x=sample.size,y=X1.se),col="tomato")+
  labs(x="Sample size",y="MCse for X1")
p5<-p+geom_line(mapping = aes(x=sample.size,y=X2.se),col="tomato")+
  labs(x="Sample size",y="MCse for X2")
p6<-p+geom_line(mapping = aes(x=sample.size,y=X3.se),col="tomato")+
  labs(x="Sample size",y="MCse for X3")

multiplot(p1, p2, p3, p4, p5, p6,cols=2)

## acceptance rate
M[[5]]

##############################################
## Q3(c)
## samples generated per second
M[[4]]

#########################################################################################
## Q4 (b)

## libraries
{library(ggplot2)
  library(tmvtnorm)
  library(astsa)
  library(plyr)
  library(mvtnorm)
  library(Matrix)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }}

## defining the parameters for the given target multivariate truncated normal
mu<-c(0,1,0)
sig<-matrix(c(1,0.8,0.3,0.8,2,0.4,0.3,0.4,3),3,3)
a<-c(-1,2,3)
b<-c(1,4,5)

## defining the target density function
target<-function(x,mu,sig,a,b){
  m<-length(mu)
  t1<-pmvnorm(a,b,mean=mu,sigma=sig)
  t2<-dmvnorm(x = x,mean = mu,sigma = sig)
  value<-t2/t1
  return(value[1])
}

## defining the proposal density function
prop<-function(x,mu,sig,a,b){
  m<-length(mu)
  q<-rep(NA_real_,m)
  for(i in 1:m){
    q[i]<-dnorm(x[i],mu[i],sqrt(sig[i,i]))/(pnorm(b[i],mean=mu[i],sd=sqrt(sig[i,i]))-pnorm(a[i],mean=mu[i],sd=sqrt(sig[i,i]))) 
  }
  value<-prod(q)
  return(value)
}

## defining the function for inverse transformation to generate samples from univariate
## truncated normal distribution
inv.tnorm<-function(x,mu,sig2,a,b){
  t1<-pnorm(a,mean=mu,sd=sqrt(sig2))
  t2<-pnorm(b,mean=mu,sd=sqrt(sig2))
  value<-qnorm(x*(t2-t1)+t1,mean=mu,sd=sqrt(sig2))
  return(value)
}

## defining the function imp.sampler to generate samples from the proposal and 
## approximate the expected value using basic importance sampling.
imp.sampler<-function(n,mu,sig,a,b){
  m<-length(mu)
  samples<-matrix(NA_real_,3,n)
  w<-rep(NA_real_,n)
  j<-1
  t<-system.time(
    while(j<=n){
      u<-runif(m)
      x<-rep(NA_real_,m)
      for (i in 1:m){
        x[i]<-inv.tnorm(u[i],mu[i],sig[i,i],a[i],b[i])
      }
      w[j]<-(target(x,mu,sig,a,b)/prop(x,mu,sig,a,b))
      samples[,j]<-x*w[j]
      j<-j+1
    })
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix means for storing Monte Carlo estimates
  MC.means<-matrix(NA_real_,3,ml)
  ## Initializing a matrix MC.se for storing Monte Carlo standard errors
  MC.se<-matrix(NA_real_,3,ml)
  ## loop for storing Monte Carlo estimates
  for (j in 1:ml){
    MC.means[,j]<-(rowSums(samples[,1:m[j]]))/m[j]
    MC.se[,j]<-sqrt(((rowSums((samples[,1:m[j]]-MC.means[,j])^2))/(m[j]-1))/m[[j]])
  }
  ## calculating ESS
  ESS.10K<-10000/(mean(w[1:10000]^2)/((mean(w[1:10000]))^2))
  ESS<-n/(mean(w^2)/((mean(w))^2))
  ESS.rate<-ESS/as.numeric(t[3])
  return(list(samples,MC.means,MC.se,ESS.10K,ESS.rate))
}
n<-50000
M<-imp.sampler(n,mu,sig,a,b)

## estimate of the expected value
M[[2]][,nrow(M[[2]])]

## Monte Carlo standard error
M[[3]][,nrow(M[[3]])]

## defining the data frame for plotting
f<-data.frame("sample.size"=seq(100,n,by=100),t(M[[2]]),t(M[[3]]))
names(f)<-c("sample.size","X1.m","X2.m","X3.m","X1.se","X2.se","X3.se")
p<-ggplot(data=f)

## plotting the estimates vs. sample size
p1<-p+geom_line(mapping = aes(x=sample.size,y=X1.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X1")
p2<-p+geom_line(mapping = aes(x=sample.size,y=X2.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X2")
p3<-p+geom_line(mapping = aes(x=sample.size,y=X3.m),col="steelblue")+
  labs(x="Sample size",y="Expectation of X3")
p4<-p+geom_line(mapping = aes(x=sample.size,y=X1.se),col="tomato")+
  labs(x="Sample size",y="MCse for X1")
p5<-p+geom_line(mapping = aes(x=sample.size,y=X2.se),col="tomato")+
  labs(x="Sample size",y="MCse for X2")
p6<-p+geom_line(mapping = aes(x=sample.size,y=X3.se),col="tomato")+
  labs(x="Sample size",y="MCse for X3")

multiplot(p1, p2, p3, p4, p5, p6,cols=2)
#############################################
## Q4 (c)
## Effective samples per 10,000 samples
M[[4]]

## Effective samples generated per second
M[[5]]

#########################################################################################
## Q5 (b)

## libraries
{library(ggplot2)
  library(tmvtnorm)
  library(astsa)
  library(plyr)
  library(mvtnorm)
  library(Matrix)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }}

## sourcing batchmeans function used to calculate MCMC standard
## errors later
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

## defining the parameters for the given target multivariate truncated normal
mu<-c(0,1,0)
sig<-matrix(c(1,0.8,0.3,0.8,2,0.4,0.3,0.4,3),3,3)
a<-c(-1,2,3)
b<-c(1,4,5)

## defining the target density function
target<-function(x,mu,sig,a,b){
  t1<-pmvnorm(a,b,mean=mu,sigma=sig)
  t2<- dmvnorm(x = x,mean = mu,sigma = sig)
  value<-t2/t1
  return(value)
}

## defining the proposal density function
prop<-function(x,mu,sig,a,b){
  m<-length(mu)
  q<-rep(NA_real_,m)
  for(i in 1:m){
    q[i]<-dnorm(x[i],mu[i],sqrt(sig[i]))/(pnorm(b[i],mean=mu[i],sd=sqrt(sig[i]))-pnorm(a[i],mean=mu[i],sd=sqrt(sig[i]))) 
  }
  value<-prod(q)
  return(value)
}

## defining the function for inverse transformation to generate samples from univariate
## truncated normal distribution
inv.tnorm<-function(x,mu,sig2,a,b){
  t1<-pnorm(a,mean=mu,sd=sqrt(sig2))
  t2<-pnorm(b,mean=mu,sd=sqrt(sig2))
  value<-qnorm(x*(t2-t1)+t1,mean=mu,sd=sqrt(sig2))
  return(value)
}

## Metropolis Hastings function with inputs as the variance of proposal
## (tuning parameter) and current state of the MC
MH<-function(variance,current.state){
  ## sampling x.star from truncated normal proposal with mean as current state
  ## of x and specified tuning parameters
  m<-length(mu)
  u<-runif(m)
  x.star<-rep(NA_real_,m)
  for (i in 1:m){
    x.star[i]<-inv.tnorm(u[i],current.state[i],variance[i],a[i],b[i])
  }
  ## defining the acceptance probability for sampled x.star on log
  ## scale
  accept.probab<-(target(x.star,mu,sig,a,b)/target(current.state,mu,sig,a,b))*(prop(x.star,current.state,variance,a,b)/prop(current.state,x.star,variance,a,b))
  ## sampling u from uniform(0,1) to check for acceptance
  u<-runif(1, min=0, max=1)
  ## initializing the indicator flag=0 to check if the sampled x.star
  ## will be accepted
  flag<-0
  ## if-else to define the next state of the chain based on acceptance
  ## probability
  if(u<=accept.probab){
    flag<-1
    next.state<-x.star
  }
  else {next.state<-current.state}
  ## returning the next state and indicator if the sampled value was
  ## accepted
  return(list(next.state,flag))
}

## defining the function Expectation for running Random Walk Metropolis
## algortithm and to calculate Monte Carlo estimates and standard errors
## for x
Expectation<-function(n,start,variance){
  ## Initializing the Markov chain for beta_1
  x<-matrix(NA_real_,3,n)
  ## Defining the initial value for the chain
  x[,1]<-start
  ## Initializing the accept count used to calculate acceptance rate
  ## of x
  accept<-0
  ## loop for RWM updates
  sys.time<-system.time(for(i in 1:(n-1)){
    temp<-MH(variance,x[,i])
    x[,i+1]<-temp[[1]]
    accept<-accept+temp[[2]]
  })
  ## samples obtained from the running the chain for given n
  samples<-data.frame("iterations"=1:n,"X1"=t(x)[,1],"X2"=t(x)[,2],"X3"=t(x)[,3])
  ## calculating the acceptance rate
  acceptance.rate<-accept/n
  ## defining the sequence for number of iterations
  m<-seq(100,n,by=100)
  ## length of the sequence m
  ml<-length(m)
  ## Initializing a matrix MCMC.means for storing MCMC estimates
  MCMC.means<-matrix(NA_real_,3,ml)
  ## Initializing a matrix MCMC.se for storing MCMC standard errors for the 3 components
  ## separately
  MCMC.se<-matrix(NA_real_,3,ml)
  ## loop for storing Monte Carlo estimates and MCMC standard errors in R
  for (j in 1:ml){
    MCMC.means[,j]<-rowSums(x[,1:m[j]])/m[j]
    ## Using the batchmeans function to calculate MCMC standard errors
    MCMC.se[1,j]<-bm(x[1,1:m[j]])$se
    MCMC.se[2,j]<-bm(x[2,1:m[j]])$se
    MCMC.se[3,j]<-bm(x[3,1:m[j]])$se
  }
  ESS.10K<-c(ess(samples[(1:10000),2]),ess(samples[(1:10000),3]),ess(samples[(1:10000),4]))
  ESS<-c(ess(samples[,2]),ess(samples[,3]),ess(samples[,4]))
  ESS.rate<-ESS/as.numeric(sys.time[3])
  return(list("samples"=samples,"iterations"=m, "means"=MCMC.means, "MCMC.se"=MCMC.se, "acceptance.rate"=acceptance.rate,"ESS.10K"=ESS.10K,"ESS.rate"=ESS.rate))
}

n<-100000
start<-matrix(c(0.30,2.66,3.66,0.35,2.71,3.71,0.25,2.61,3.61),3,3)
##start<-c(0,3,4)
start
var.tune<-c(1,2,3)
M<-vector("list",ncol(start))

for (j in 1:ncol(start)){
  M[[j]]<-Expectation(n,start[,j],var.tune)
}

## checking sample variances for tuning
var(M[[1]][[1]][,2])
var(M[[1]][[1]][,3])
var(M[[1]][[1]][,4])

## MCMC estimates for 3 chains
M[[1]][[3]][,nrow(M[[1]][[3]])]
M[[2]][[3]][,nrow(M[[2]][[3]])]
M[[3]][[3]][,nrow(M[[3]][[3]])]

## MCMC standard errors for 3 chains
M[[1]][[4]][,nrow(M[[1]][[4]])]
M[[2]][[4]][,nrow(M[[2]][[4]])]
M[[3]][[4]][,nrow(M[[3]][[4]])]

## checking acf plots for 3 chains
acf(M[[1]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[2]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")
acf(M[[3]][[1]][,-1],lag.max = 50,main="Plot of ACF of samples")

## ESS per 10,000 samples for 3 chains
M[[1]][[6]]
M[[2]][[6]]
M[[3]][[6]]

## ESS rate for 3 chains
M[[1]][[7]]
M[[2]][[7]]
M[[3]][[7]]

## new starting values
M[[1]][[1]][n,-1]
M[[2]][[1]][n,-1]
M[[3]][[1]][n,-1]

## plotting
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

## defining the data frames for plotting
f<-list()
g<-list()
for(j in 1:ncol(start)){
  f[[j]]<-data.frame(M[[j]][[2]],t(M[[j]][[3]]),t(M[[j]][[4]]),j)
  g[[j]]<-data.frame((M[[j]][[1]]),j)
}

f<-do.call(rbind,f)
names(f)<-c("iterations","mean.X1","mean.X2","mean.X3","MCMCse.X1","MCMCse.X2","MCMCse.X3","start.label")

g<-do.call(rbind,g)
names(g)<-c("iterations","samples.X1","samples.X2","samples.X3","start.label")

p<-ggplot(data=f)
q<-ggplot(data=g)

## plotting mean vs. sample size for different starting values
p1<-p+ geom_line(mapping = aes(x=iterations,y=mean.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y="Estimate of the expectation of X1 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X1,ymin=mean.X1-MCMCse.X1,ymax=mean.X1+MCMCse.X1,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
p2<-p+ geom_line(mapping = aes(x=iterations,y=mean.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y="Estimate of the expectation of X2 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X2,ymin=mean.X2-MCMCse.X2,ymax=mean.X2+MCMCse.X2,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
p3<-p+ geom_line(mapping = aes(x=iterations,y=mean.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))+labs(x="Number of samples",y="Estimate of the expectation of X3 with MCMC standard errors", colour="Label for starting values")+thema+
  geom_errorbar(mapping = aes(x=iterations,y=mean.X3,ymin=mean.X3-MCMCse.X3,ymax=mean.X3+MCMCse.X3,group=factor(start.label),colour=factor(start.label)),subset=.((iterations%%1000)==0))
multiplot(p1, p2, p3,cols=1)

## plotting (MCMCse) vs. sample size for different starting values
p4<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X1),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X1",colour="Label for starting values")+thema
p5<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X2),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X2",colour="Label for starting values")+thema
p6<-p+ geom_line(mapping = aes(x=iterations,y=(MCMCse.X3),group=factor(start.label),colour=factor(start.label)))+labs(x="Number of samples",y="MCMCse for X3",colour="Label for starting values")+thema
multiplot(p4, p5, p6,cols=1)

## estimated marginal densities after N/2 and after N
q1<-q+ geom_density(mapping = aes(x=samples.X1),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X1",title="After N/2")+thema
q2<-q+geom_density(mapping = aes(x=samples.X1),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X1",title="After N")+thema

q3<-q+ geom_density(mapping = aes(x=samples.X2),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X2",title="After N/2")+thema
q4<-q+geom_density(mapping = aes(x=samples.X2),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X2",title="After N")+thema

q5<-q+ geom_density(mapping = aes(x=samples.X3),fill="steelblue",subset=.(start.label==1&iterations[1:(length(iterations)/2)]))+
  labs(x="",y="For X3",title="After N/2")+thema
q6<-q+geom_density(mapping = aes(x=samples.X3),fill="tomato",subset=.(start.label==1))+
  labs(x="",y="For X3",title="After N")+thema

multiplot(q1,q3,q5,q2,q4,q6, cols=2)

## estimated density for different starting values
q7<-q+ geom_density(mapping = aes(x=samples.X1,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X1",colour="Label of starting values")+thema
q8<-q+ geom_density(mapping = aes(x=samples.X2,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X2",colour="Label of starting values")+thema
q9<-q+ geom_density(mapping = aes(x=samples.X3,group=factor(start.label),colour=factor(start.label)))+
  labs(x="",y="Estimated marginal density for X3",colour="Label of starting values")+thema
multiplot(q7,q8,q9, cols=1)

##########################################
## Q5 (c)
## Effective samples for every 10,000 samples generated
min(M[[1]][[6]])

## Number of effective samples generated per second
min(M[[1]][[7]])

sim<-function(n,n.cov,n.class,p.vec,mu){
  y <- sample(1:n.class, size = n, prob = p.vec, replace = TRUE)
  y.list <- split(x = y,f = y)
  x.list <- list()
  for(i in 1:n.class){
    x.list[[i]]<- mvrnorm(n = length(y.list[[i]]), mu = mu[,i], Sigma = diag(.1,nrow = n.cov) )
  }
  
  temp2<-lapply(X = 1:n.class,FUN = function(x){
    temp<-data.frame(y.list[[x]],x.list[[x]])
    names(temp)<-c("Y","X1","X2")
    return(temp)
  })
  
  datas<-do.call(what = rbind,args = temp2)
  data <- data.frame(datas)
  invisible(data)
}

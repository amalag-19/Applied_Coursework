#Calculating LDA
my.lda<-function(response,formula,data){
  formula<-paste(formula,-1)
  formulas<-as.formula(formula)
  X.matrix<-model.matrix(object = formulas,data = data)
  Y.matrix<-matrix(response)
  dados<-data.frame("Y"=Y.matrix,X.matrix)
  dados.l<-split(x = dados,f = dados$Y)
  mean.l<-lapply(X = dados.l,FUN = function(x){
    temp<-matrix(rep(colMeans(x),nrow(x)),nrow = nrow(x),byrow = T)
    temp<-data.frame(temp)
    names(temp)<-names(dados)
    return(temp)
  })
  
  dif<- lapply(X = 1:length(dados.l),FUN = function(x){
    dados.l[[x]][,-1]-mean.l[[x]][,-1]
  })
  dif2<- lapply(X = dif,FUN = function(x){
    t(as.matrix(x))%*%as.matrix(x)
  })  
  covM<-matrix(data = 0,nrow = ncol(X.matrix),ncol=ncol(X.matrix))
  for(i in 1:length(dif2)){
    covM<-covM+dif2[[i]]
  }
  CovM<-covM/(nrow(X.matrix)-ncol(X.matrix))
  MeanL<-lapply(X = mean.l,FUN = function(x){
    data.frame(x[1,-1])
  })
  names(MeanL)<-names(dados.l)
  retornar<-list("Covariance"=CovM,"Means"=MeanL)
  return(retornar)
}


predict.my.lda<-function(x,model,prior){
  x<-as.matrix(x)
  Cov<-solve(model$Covariance)
  if(!is.matrix(x)&!is.data.frame(x))x<-matrix(x,nrow=1)
  value<-lapply(1:length(model$Means), function(y){
    Means<-t(as.matrix((model$Means[[y]])))
    priorP<-prior[y]
    x%*%Cov%*%as.matrix(Means)-c(0.5*t(Means)%*%Cov%*%Means+log(priorP))
  })
  
  value<-data.frame(do.call(cbind,value))
  names(value)<-names(model$Means)
  value2 <- apply(X = as.matrix(value),MARGIN = 1,FUN = function(x){
    which.max(x)-1
  })
  value2 <- as.data.frame(value2)
  names(value2)<-"Class"
  return(value2)
}

create.cross.val.data <- function(k,data,...){
  id <- sample(x = nrow(data),replace = F)
  fold.size <- floor(nrow(data)/k)
  final.id <- list()
  for( i in 1:(k-1)){
    final.id[[i]] <- id[seq((fold.size*(i-1)+1),fold.size*i,by=1)]
  }
  final.id[[k]] <- id[seq((fold.size*(k-1)+1),length(id),by=1)]
  names(final.id)<-paste0("Fold",1:k)
  final.data <- list()
  for( i in 1:(k)){
    final.data[[i]] <- data[final.id[[i]],]
  }
  names(final.data)<-paste0("Fold",1:k)
  return(list(Data=final.data,Id=final.id))
}


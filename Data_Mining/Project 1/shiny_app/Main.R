source("untitled folder/untitled folder/QDA.R")
source("untitled folder/untitled folder/logit.R")
source("untitled folder/untitled folder/lda.function.R")
source("untitled folder/untitled folder/transformation_function.R")

mscls.rt <- function(true.cl, pred.cl){
  sum((true.cl!=pred.cl))/length(true.cl)
}

load("untitled folder/cv_data.RData")




cv.fun <- function(cv.data){
  #Parameters
  fold = 5
  alg = 3
  miss.class <- matrix(data = NA_real_,nrow = fold,ncol = alg)
  colnames(miss.class) <- c("LDA","QDA","Logit")
  
  for(i in 1:fold){
    temp <- cv.data$Data
    temp[[i]]<-NULL
    test <- cv.data$Data[[i]]
    train <- do.call(rbind ,temp)
    
    tr <- trsfrm.gen(train.data = train[,-ncol(train)])
    
    tr.train <- cbind(train[,ncol(train)],tr(train[-ncol(train)]))
    tr.test <- cbind(test[,ncol(test)],tr(test[-ncol(train)]))
    colnames(tr.test) <- c(names(test)[ncol(test)],names(test)[1:(ncol(test)-1)])
    colnames(tr.train) <- c(names(train)[ncol(train)],names(train)[1:(ncol(train)-1)])
    ########LDA
    
    lda.model <- my.lda(response = tr.train[,1],formula = "~.",data = data.frame(tr.train[,-1]))
    lda.pred <- predict.my.lda(x = tr.test[,-1],model = lda.model,prior = c(.5,.5))
    miss.class[i,1] <- mscls.rt(tr.test[,1],lda.pred)
    
    #######QDA
    qda.pred<-rep(NA_integer_,nrow(tr.test))
    for(k in 1:nrow(tr.test)){
      qda.pred[k]<-QDA.1(tr.train,as.matrix(tr.test[k,-1]))
    }
    qda.pred
    miss.class[i,2]<-mscls.rt(tr.test[,1],qda.pred)
    
    ########Logistic Regression 
    logit.bs <- logit.m(X =cbind(1,tr.train[,-1]), y = tr.train[,1])
    logit.pred <- predict.logit.cl(logit.bs, X = cbind(1, tr.test[,-1]))                   
    miss.class[i,3] <- mscls.rt(tr.test[,1], logit.pred)
  }
  return(miss.class)
}

miss.class <- cv.fun(cv.data)
mscls.mns <- colMeans(miss.class)




##########PCA
entire.data <- list()
for(i in 1:fold){
  entire.data[[i]] <- data.frame(cv.data$Data[[i]],Fold=i)
}
entire.data<-do.call(rbind,entire.data)

cor(entire.data[,-c(ncol(entire.data)-1,ncol(entire.data))])
pdf(file= "untitled folder/untitled folder/untitled.pdf",width = 100,height = 100)
pairs(entire.data[,-c(ncol(entire.data)-1,ncol(entire.data))])
dev.off()
prince <- princomp(entire.data[,-c(ncol(entire.data)-1,ncol(entire.data))])
prince$loadings
P <- prince$loadings[,1:20]
PX <- as.matrix(entire.data[,-c(ncol(entire.data)-1,ncol(entire.data))])%*%P
pca.data <- cbind(PX, entire.data[,c(ncol(entire.data)-1, ncol(entire.data))])
cv.data.pca <- split(pca.data,f = pca.data[,ncol(pca.data)])
for(i in 1:5){
  cv.data.pca[[i]][ncol(cv.data.pca[[i]])] <- NULL
}                   


cv.data.pca2 <- list()
cv.data.pca2[[2]] <- cv.data.pca
names(cv.data.pca2)<- c("NULL","Data")
colMeans(cv.fun(cv.data.pca2))

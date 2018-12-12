
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(MASS)
library(ggplot2)

shinyServer(function(input, output) {
  source("lda.function.R")
  source("QDA.R")
  source("logit.R")
  source("simulation_function.R")
  output$distPlot <- renderPlot({
    input$go
    sim<-isolate(sim.logistic(n=input$n,n.cov=2,beta=c(input$beta.0,input$beta.1,input$beta.2)))
    ## summary(factor(sim$y))
    ########LDA
    
    lda.model <- my.lda(response = sim[,1],formula = "~.",data = data.frame(sim[,-1]))
    lda.pred <- predict.my.lda(x = sim[,-1],model = lda.model,prior = c(0.5,0.5))
    ##miss.class[i,1] <- mscls.rt(tr.test[,1],lda.pred)
    
    #######QDA
    qda.pred<-rep(NA_integer_,nrow(sim))
    for(k in 1:nrow(sim)){
      qda.pred[k]<-QDA.1(sim,as.matrix(sim[k,-1]))
    }
    ##miss.class[i,2]<-mscls.rt(tr.test[,1],qda.pred)
    
    ########Logistic Regression 
    logit.bs <- logit.m(X =cbind(1,sim[,-1]), y = sim[,1])
    logit.pred <- as.numeric(predict.logit.cl(matrix(logit.bs), X = as.matrix(cbind(1, sim[,-1])))) 
    print(predict.logit(matrix(logit.bs), X = as.matrix(cbind(1, sim[,-1]))))
    ##miss.class[i,3] <- mscls.rt(sim[,1], logit.pred)
    print(nrow(lda.pred))
    df<-data.frame(sim,lda.pred,qda.pred,logit.pred)
    p<-ggplot(data=df)
    
    print(p+geom_point(mapping=aes(x=X2,y=X1,colour=factor(logit.pred),shape=factor(Class),size=factor(qda.pred)))+scale_size_discrete(range = c(3,8))+facet_wrap(~y,nrow = 2))
    
    })

})

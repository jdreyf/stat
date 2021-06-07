source("EM_tools.R",local=TRUE)
source("General_tools.R",local=TRUE)
source("PCA_tools.R",local=TRUE)
source("FISTA_tools.R",local=TRUE)
source("CodeMimi.R",local=TRUE)

library(randomForest)
library(missForest)
library(doParallel)
library(doSNOW)
cl <- makeCluster(20, outfile="") #here, change the number of clusters
registerDoSNOW(cl)

######
# Name: ComparMNAR_Univariate
# Date: 27/12/2018
# Description: For the univariate case (1 missing variable: the first one),this function allow to compare different algorithms and methods to impute and estimate matrices which contain MNAR or MAR missing values.
# The function output is a list containing, for each simulation, the mean squared errors (the prediction error and the total error) for the different algorithms and methods.
# Arguments: 
#Xtrue: the parameter matrix.
#a, b: the logistic regression parameters.
#r: rank of the parameter matrix.
#noise: sigma^2, the added noise to the parameter matrix.
#Ns: number of Monte Carlo simulations in the EM algorithm.
#modmecha: if "logit", the missing-data distribution is a logistic regression, otherwise it is a probit distribution.
#mecha: "MNAR" or "MAR".
#nbsim: number of simulations.
#####

# seems to loop IterEM for FISTA / pred twice: once with ccompt and again with Tt

ComparMNAR_Univariate <- function(Xtrue,a=NULL,b=NULL,r,noise,Ns,modmecha,mecha,nbsim){
  
  nbcol=1
  p<-ncol(Xtrue)
  n<-nrow(Xtrue)
  
  results.list = foreach (ksim = 1:nbsim, .combine = "rbind")  %dopar% {
    source("EM_tools.R",local=TRUE)
    source("General_tools.R",local=TRUE)
    source("PCA_tools.R",local=TRUE)
    source("FISTA_tools.R",local=TRUE)
    source("CodeMimi.R",local=TRUE)
    
    print(paste("it globale",ksim))
    
    conv=FALSE
    while(!conv){
      
      set.seed(ksim)
      X=Xtrue+matrix(data=rnorm(n*p,0,sqrt(noise)),ncol=p)
      
      #Logit or probit distribution
      # uses a, b from environment :-(
      select_prob <- function(x,modmecha){ #probability of selecting coordinate Xij
        if(modmecha=="logit"){
          res=1/(1+exp(-a*(x-b)))}else{ res=pnorm(x)}
        return(res)
      }
      
      # simulate missing values based on mechanism
      if(mecha=="MNAR"){prob <- sapply(X[,1],select_prob,modmecha)} else {prob <- sapply(X[,2],select_prob,modmecha)}
      compt=0
      missing=c()
      for (i in 1:n){
        u<-runif(1)
        compt=compt+(prob[i]>u)
        if(prob[i]>u){missing=c(missing,i)}
      }
      print(compt)
      XNA=X
      XNA[missing,1]=NA
      
      M = 1-is.na(XNA)
      
      # Concatenating the data matrix and the mask
      Y = cbind.data.frame(XNA,M[,1])
      
      X.mean <- as.matrix(ImputeMean(XNA))
      
      # EM with modell
      # theta_hat initialized with mean imputed matrix 
      # "bis" variables are for FISTA; others are for soft thresholding
      ThetaNew=Initialize_theta(ImputeMean0(XNA),r) #on initialise Theta avec une matrice en rang inférieur
      ThetabisNew=ThetaNew
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      
      diff=100
      ccompt<-0
      while(ccompt<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNew,aNew,bNew,M,Ns,noise,algo="soft",lam="Pred",nbcol=1)
        diff=ParamNew$diff
        print(diff)
        ThetaNew=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv1=ParamNew$conv
        ccompt=ccompt+1
      }
      
      
      diff=100
      ccompt2<-0
      while(ccompt2 < 20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNew,abisNew,bbisNew,M,Ns,noise,algo="FISTA",lam="Pred",nbcol=1)
        diff=ParamNew$diff
        ThetabisNew=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv2=ParamNew$conv
        ccompt2=ccompt2+1
      }
      
      Theta=list()
      Theta[[1]]=ThetaNew
      Thetabis=list()
      Thetabis[[1]]=ThetabisNew
      aListNew=list()
      bListNew=list()
      aListNew[[1]]=aNew
      bListNew[[1]]=bNew
      abisListNew=list()
      bbisListNew=list()
      abisListNew[[1]]=abisNew
      bbisListNew[[1]]=bbisNew
      Tt=10
      conv3=c()
      conv4=c()
      pb <- txtProgressBar(min = 0, max = 50, title = 'Blabrf', label = 'fgrtgr')
      for (t in 1:Tt){
        ParamNew <- IterEM(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],M,Ns,noise,algo="soft",lam="Pred",nbcol=1)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        conv3=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],M,Ns,noise,algo="FISTA",lam="Pred",nbcol=1)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        conv4=c(conv4,ParamNew$conv)
        
        setTxtProgressBar(pb = pb, value = t)
      }
      if(sum(conv3)==0){conv3=FALSE}else{conv3=TRUE}
      if(sum(conv4)==0){conv4=FALSE}else{conv4=TRUE}
      
      
      ThetaNew=Theta[[Tt]]
      ThetabisNew=Thetabis[[Tt]]
      
      ThetaNewTot=Initialize_theta(ImputeMean0(XNA),r)
      ThetabisNewTot=ThetaNewTot
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      
      diff=100
      ccompt3=0
      while(ccompt3<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNewTot,aNew,bNew,M,Ns,noise,algo="soft",lam="Tot",nbcol=1)
        diff=ParamNew$diff
        ThetaNewTot=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv5=ParamNew$conv
        ccompt3=ccompt3+1
      }
      
      
      diff=100
      ccompt4=0
      while(ccompt4<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNewTot,abisNew,bbisNew,M,Ns,noise,algo="FISTA",lam="Tot",nbcol=1)
        diff=ParamNew$diff
        ThetabisNewTot=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv6=ParamNew$conv
        ccompt4=ccompt4+1
      }
      
      
      Theta=list()
      Theta[[1]]=ThetaNewTot
      Thetabis=list()
      Thetabis[[1]]=ThetabisNewTot
      aListNew=list()
      bListNew=list()
      aListNew[[1]]=aNew
      bListNew[[1]]=bNew
      abisListNew=list()
      bbisListNew=list()
      abisListNew[[1]]=abisNew
      bbisListNew[[1]]=bbisNew
      Tt=10
      conv7=c()
      conv8=c()
      pb <- txtProgressBar(min = 0, max = 50, title = 'Blabrf', label = 'fgrtgr')
      for (t in 1:Tt){
        ParamNew <- IterEM(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],M,Ns,noise,algo="soft",lam="Tot",nbcol=1)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        conv7=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],M,Ns,noise,algo="FISTA",lam="Tot",nbcol=1)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        conv8=c(conv4,ParamNew$conv)
        
        setTxtProgressBar(pb = pb, value = t)
      }
      if(sum(conv7)==0){conv7=FALSE}else{conv7=TRUE}
      if(sum(conv8)==0){conv8=FALSE}else{conv8=TRUE}
      
      ThetaNewTot=Theta[[Tt]]
      ThetabisNewTot=Thetabis[[Tt]]
      
      paste(c(conv1,conv2,conv3,conv4,conv5,conv6,conv7,conv8))
      if(conv1 & conv2 & conv3 & conv4 & conv5 & conv6 & conv7 & conv8){conv=TRUE}else(print("Pas de convergence"))
      
    }
    
    ###############
    ####### MSE
    ###############
    
    msePred=sapply(list(ThetabisNew[,1:p]*(1-M), ThetabisNewTot[,1:p]*(1-M)), FUN=MSE, X2=X*(1-M))
    mseTrue=sapply(list(ThetabisNew[,1:p], ThetabisNewTot[,1:p]), FUN=MSE, X2=Xtrue)
    
    # names=c("Imputemean","modelsoftPred","modelFISTAPred","modelsoftTot","modelFISTATot","softTot","softPred","softmaskTot","softmaskPred","PCATot","PCAPred","PCAmaskTot","PCAmaskPred","mimiTot","mimiPred","FISTA","FISTAPred","FISTAmask","FISTAmaskPred","randomforest","randomforestmask")
    names <- c("modelFISTAPred", "modelFISTATot")
    
    cbind(msePred, mseTrue, names)
  }
  
  return(results.list)
}

# stopCluster(cl)

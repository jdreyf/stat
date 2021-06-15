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

for (t in 1:Tt){
  ParamNew <- IterEM(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],M,Ns,noise,algo="FISTA",lam="Pred",nbcol=1)
  Thetabis[[t+1]]=ParamNew$ThetaNew
  abisListNew[[t+1]]=ParamNew$a_initNew
  bbisListNew[[t+1]]=ParamNew$b_initNew
  conv4=c(conv4,ParamNew$conv)
  
  setTxtProgressBar(pb = pb, value = t)
}
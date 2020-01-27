firstDgamma<-function(uhatN,conFdata,gamma){
  N=nrow(conFdata)
  fdgamma=rep(0,length(gamma))
  for(i in 1:N){
    gammaw=sum(gamma * c(1,unlist(conFdata[i,])))
    fdgamma=fdgamma+(-uhatN[i]+exp(gammaw)/(1+exp(gammaw)))*c(1,unlist(conFdata[i,]))
  }
  return(fdgamma)
}

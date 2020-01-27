secondDgamma<-function(uhatN,conFdata,gamma){
  N=nrow(conFdata)
  sdgamma=matrix(0,length(gamma),length(gamma))
  for(i in 1:N){
    wimt= c(1,unlist(conFdata[i,]))
    gammaw=sum(gamma * c(1,unlist(conFdata[i,])))
    sdgamma=sdgamma+(exp(gammaw)/(1+exp(gammaw))^2)*outer(wimt,wimt)
  }
  return(sdgamma)
}

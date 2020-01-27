secondDci=function(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group){
  N2=length(unique(group));N=length(uhatN);tmax=ceiling(N/N2)
  diagsDbiN=-(1-uhatN)*exp(theta_bi+as.vector(cihatN2[group]))
  diagsDaiN=exp(logitp_ai+as.vector(cihatN2[(group+N2)]))/(1+exp(logitp_ai+as.vector(cihatN2[(group+N2)])))^2
  if(sum(is.na(diagsDaiN))>0) diagsDaiN=exp(-logitp_ai-as.vector(cihatN2[(group+N2)]))/(1+exp(-logitp_ai-as.vector(cihatN2[(group+N2)])))^2

  c(tapply(diagsDbiN, group, sumn0)+1/sigma[1]^2,
    tapply(diagsDaiN, group, sumn0)+1/sigma[2]^2)
}

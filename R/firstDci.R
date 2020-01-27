firstDci=function(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group){
  N2=length(unique(group));N=length(uhatN)#;tmax=ceiling(N/N2)
  fDbiN=(1-uhatN)*(yFdata-exp(theta_bi+as.vector(cihatN2[group])))
  fDaiN=1/(1+exp(-logitp_ai-as.vector(cihatN2[(group+N2)])))-uhatN
  c(tapply(fDbiN, group, sumn0)+cihatN2[1:N2]/sigma[1]^2,
    tapply(fDaiN, group, sumn0)+cihatN2[((N2+1):(2*N2))]/sigma[2]^2)
}

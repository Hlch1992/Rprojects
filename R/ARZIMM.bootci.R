ARZIMM.bootci=function(obj,nboot=500){
  datalist=obj$datalist
  parabootlist=obj$paralist
  parasetup=obj$parasetup
  tunbootlist=list(lambdabeta=parasetup$lambda,lambdagamma=parasetup$tunpara$lambdaseq2,
                   weight1all=parasetup$initpara$weight1all,weight2all=parasetup$initpara$weight2all,
                   selgamma=parasetup$selectpara$selgamma)

  # create progress bar
  pb <- txtProgressBar(min = 0, max = nboot, style = 3)

  bootres=NULL
  for (boottime in 1:nboot){
    set.seed(boottime)
    yFdataNMboot=bootsim(datalist=datalist,paralist=parabootlist)
    datalisttmp=datalist;datalisttmp$yFdataNM=yFdataNMboot

    boottmp=bootest(datalist=datalisttmp,paralist=parabootlist,tunlist=tunbootlist)

    bootres=rbind(bootres,boottmp)
    setTxtProgressBar(pb, boottime)

  }
  ################## progress bar
  close(pb)

  bootsd=apply(bootres, 2, function(x) sd(x))
  pointest=c(as.vector(beta),as.vector(gamma),sigma)
  chistat=unlist(sapply(1:length(pointest), function(x) ifelse(bootsd[x]==0,0,(pointest[x]/bootsd[x])^2)))
  bootwtpval=pchisq(chistat,1,lower.tail = F)
  obj$bootparapval=list(betapval=matrix(bootwtpval[1:length(as.vector(beta))],ncol = M),
                        gammapval=matrix(bootwtpval[(1+length(as.vector(beta))):(length(as.vector(beta))+length(as.vector(gamma)))],ncol = M),
                        sigmapval=bootwtpval[(1+length(as.vector(beta))+length(as.vector(gamma))):length(bootwtpval)])
  return(obj)

}


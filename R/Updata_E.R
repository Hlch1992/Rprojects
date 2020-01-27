Updata_E=function(yFdata,theta_bi,logitp_ai,sigma,cihatN2,group) {
  ########## parameters
  N=length(yFdata);N2=length(unique(group))#;tmax=ceiling(N/N2)

  ########## set boundary
  BIG=50; SMALL=1.0e-3; converge=0;maxits=100;itse=0; fail=0;

  ########## start value
  while(itse < maxits) {
    itse=itse+1
    cihatN2old=cihatN2

    ###### estimated u hat
    uhatN=as.numeric(yFdata==0)/(1+exp(-logitp_ai-as.vector(cihatN2[group+N2])-
                                         exp(theta_bi+as.vector(cihatN2[group]))))

    ###### updata bi by minimize Si(bi)
    resE=Updateci(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group)
    cihatN2=resE$cihatN2;diagsDci=resE$diagsDci

    eps=sum(abs(cihatN2-cihatN2old))
    ##### print(paste('estbi eps',eps))
    if(eps>BIG) {break
    }else if (eps<SMALL) {converge=1;break}
  }
  ### print(paste('estbi con',converge))
  ###### estimated u hat
  uhatN=as.numeric(yFdata==0)/(1+exp(-logitp_ai-as.vector(cihatN2[group+N2])-
                                       exp(theta_bi+as.vector(cihatN2[group]))))

  if(!exists('diagsDci')) {cihatN2=resE$cihatN2;diagsDci=resE$diagsDci}

  list(cihatN2=cihatN2,civar=1/diagsDci,uhatN=uhatN,conv=converge)
}

Updateci=function(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group){
  ########## parameters
  N=length(yFdata);N2=length(unique(group))#;tmax=ceiling(N/N2)

  ########## set boundary
  BIG=50; SMALL=1.0e-3; converge=0;maxits=200;itsci=0; fail=0;

  while(itsci < maxits) {
    itsci=itsci+1
    cihatN2old=cihatN2

    fDci=firstDci(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group)
    diagsDci=secondDci(yFdata,uhatN,theta_bi,logitp_ai,sigma,cihatN2,group)
    grdd=fDci/diagsDci
    grdd[is.na(grdd)]=0
    cihatN2=cihatN2old-grdd

    eps=sum(abs(cihatN2-cihatN2old))
    ##### print(paste('estbi eps',eps))
    if(eps>BIG) {break
    }else if (eps<SMALL) {converge=1;break}
  }
  if(converge!=1) cihatN2=cihatN2old
  diagsDci[is.na(diagsDci)]=10^7
  list(cihatN2=cihatN2,diagsDci=diagsDci,conv=converge)
}

Updatebi=function(yFdata,linkmuN_bi,sigma,uhatN,bihatN2,group){
  ########## parameters
  N=length(yFdata);N2=length(unique(group));tmax=ceiling(N/N2)

  ########## set boundary
  BIG=1000; SMALL=1.0e-3; converge=0;maxits=200;itsbi=0; fail=0;

  while(itsbi < maxits) {
    itsbi=itsbi+1
    bihatN2old=bihatN2

    fDbi=firstDbi(yFdata,linkmuN_bi,uhatN,bihatN2,sigma,group)
    diagsDbi=secondDbi(linkmuN_bi,uhatN,bihatN2,sigma,group)
    grdd=fDbi/diagsDbi
    grdd[is.na(grdd)]=0
    bihatN2=bihatN2old-grdd

    eps=sum(abs(bihatN2-bihatN2old))
    ##### print(paste('estbi eps',eps))
    if(eps>BIG) {break
    }else if (eps<SMALL) {converge=1;break}
  }
  diagsDbi[is.na(diagsDbi)]=10^7
  list(bihatN2,diagsDbi,converge)
}

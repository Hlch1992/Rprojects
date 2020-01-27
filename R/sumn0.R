sumn0=function(x) {
  if(sum(x)==0) {return(sum(x))
  }else{return(length(x)*mean(x[which(x!=0)]))}
}

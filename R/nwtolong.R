nwtolong=function(nwtable){
  nwtablelong=NULL
  for (j in 1:ncol(nwtable)) {
    for(i in 1: nrow(nwtable)){
      if(nwtable[i,j]!=0) nwtablelong=rbind(nwtablelong,c(colnames(nwtable)[i],rownames(nwtable)[j],nwtable[i,j]))
    }
  }
  nwtablelong=as.data.frame(nwtablelong)
  colnames(nwtablelong)=c('from','to','weights')

  nwtablelong[,3]=as.numeric(as.character(nwtablelong[,3]))

  nwtablelong
}

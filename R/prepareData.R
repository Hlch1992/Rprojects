prepareData=function(phy,Varname,Conname,Covname,IDname,Tname,family){
  fdata=cbind(otu_table(phy),sample_data(phy))

  ##### order fdata by time points and subjects
  fdata[,Tname]=as.numeric(as.character(fdata[,Tname]));fdata=fdata[order(fdata[,IDname],fdata[,Tname]),]

  ##### define subjects' ID
  ID=unique(fdata[,IDname])

  ##### define basic parameters
  N2=length(unique(fdata[,IDname]));M=length(Varname);N=nrow(fdata)
  if(!is.na(Conname[1])) {L=length(Conname)}else{L=0}
  if(!is.na(Covname[1])) {q=length(Covname)}else{q=0}
  #;tmax=max(fdata[,Tname])

  ###### get ind for x and y
  indy=which(fdata[,Tname]!=0)
  indx=match(paste0(fdata[indy,IDname],str_pad(fdata[indy,Tname]-1, 2, pad = "0")),
             paste0(fdata[,IDname],str_pad(fdata[,Tname], 2, pad = "0")))
  ind0=match(fdata[indy,IDname],fdata[which(fdata[,Tname]==0),IDname])
  xFdataNM0=fdata[ind0,Varname]

  ####### define data as y, x, covariates and concomitant
  yFdataNM=fdata[indy,Varname];  xFdataNM=fdata[indx,Varname]
  if(is.na(Conname[1])){conFdataNL=NULL}else{conFdataNL=t(t(fdata[indx,Conname]))}
  if(is.na(Covname[1])){covFdataNq=NULL}else{covFdataNq=t(t(fdata[indx,Covname]))}
  conFdataNL=apply(conFdataNL,2,function(x) as.numeric(as.character(x)))
  if(!is.na(Covname)) covFdataNq=apply(covFdataNq,2,function(x) as.numeric(as.character(x)))

  ####### log transformation x data
  if (family=='Poisson') xFdataNM=log(xFdataNM+1)
  if(!is.na(Covname[1])) {xFdataNMq=cbind(xFdataNM,covFdataNq)}else{xFdataNMq=xFdataNM}

  ####### define group variable for random effect
  group=match(fdata[indx,IDname],unique(fdata[,IDname]))

  datalist=list(yFdataNM=yFdataNM,xFdataNM=xFdataNMq,conFdataNL=conFdataNL,group=group,xFdataNM0=xFdataNM0)

  datalist
}

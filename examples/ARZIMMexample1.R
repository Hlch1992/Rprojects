\dontrun{
data(sampledata)
Varname=colnames(sampledata)[1:20]
Conname=colnames(sampledata)[21:26]
Tname=colnames(sampledata)[27]
IDname=colnames(sampledata)[28]
ARZIMMresult=ARZIMM::ARZIMM(phy=NULL,Varname = Varname,Conname = Conname,fdata = sampledata,
               IDname = IDname,Tname = Tname,bootpara=list(bootpval=TRUE,nboot=5))
}

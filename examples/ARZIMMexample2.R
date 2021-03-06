\dontrun{
  data(phyExample)
  Varname=colnames(otu_table(phyExample))
  Conname=colnames(sample_data(phyExample))[1:6]
  Tname=colnames(sample_data(phyExample))[7]
  IDname=colnames(sample_data(phyExample))[8]

  ARZIMMresult=ARZIMM::ARZIMM(phyExample,Varname = Varname,Conname = Conname,
                 IDname = IDname,Tname = Tname,bootpara=list(bootpval=TRUE,nboot=5))
}

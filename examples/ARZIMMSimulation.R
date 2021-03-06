\dontrun{
  require(ARZIMM)
  data(sampleparameters)
  set.seed(1234)
  ARZIMMSimulateData=simMixTime(baseFdataN2M=parameters$baseFdataN2M,conFdataN2L=parameters$conFdataN2L,
             timeN=parameters$timeN,interceptM=parameters$interceptM,betaMM=parameters$betaMM,
             gammaLM=parameters$gammaLM,sigmaM=parameters$sigmaM,biN2M=parameters$biN2M)
}

filfun=function(phy,relative.abundance.threshold,filnum){
  phy.norm = transformSampleCounts(phy, function(x) x/sum(x))
  ind = rowMeans(otu_table(phy.norm))>relative.abundance.threshold###
  phy=prune_taxa(ind,phy)
  phy = filter_taxa(phy, function(x) sum(x > 3) > filnum, TRUE)
  print(phy)
  return(phy)
}


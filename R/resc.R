resc=function(xvect) {
  ranuni=range(xvect)[2]-range(xvect)[1]
  (xvect-min(xvect))/ranuni+.1
}

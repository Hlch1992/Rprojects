sqrtsingular=function(sqrtmatrix) {
  Xm = eigen(sqrtmatrix)
  Tm = Xm$vectors
  Jm = Diagonal( x=Xm$values )
  Tinv = solve(Tm)
  Tm %*% Jm %*% Tinv
  Jsqrt = Diagonal( x=sqrt( Xm$values ) )
  Jsqrt[is.na(Jsqrt)]=0
  SDsqrt = Tm %*% Jsqrt %*% Tinv
  return(SDsqrt)
}

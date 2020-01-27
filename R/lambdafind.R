lambdafind=function(yFdataNM,xFdataNM,epsilon = c(100,.0001),K = 5) {
  mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
  x=xFdataNM
  lambdapath=matrix(0,nrow = K,ncol = ncol(yFdataNM))
  for(i in 1:ncol(yFdataNM)){
    y=yFdataNM[,i]
    ## Standardize variables: (need to use n instead of (n-1) as denominator)
    sx <- scale(x, scale = apply(x, 2, mysd))
    sy <- as.vector(scale(y, scale = mysd(y)))

    ## Calculate lambda path (first get lambda_max):
    lambda_max = max(abs(colSums(sx*sy)),na.rm = T)/length(y)
    lambdapath[,i] <- round(exp(seq(log(max(lambda_max)*epsilon[1]), log(min(lambda_max)*epsilon[2]),
                                    length.out = K)), digits = 10)
  }
  lambdapath
}

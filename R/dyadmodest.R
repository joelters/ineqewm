dyadmodest <- function(D,X,Y,f,
                       ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL")){
  #create dyads
  n <- length(Y)
  n1 <- n - 1
  XX <- matrix(0,n*(n-1)*0.5,ncol(X)*2 + 2)
  YY <- rep(0,n*(n-1)*0.5)
  cnt <- 0
  for (i in 1:n1){
    j1 <- i + 1
    for (j in j1:n){
      cnt <- cnt + 1
      XX[cnt,] <- c(D[i],as.numeric(X[i,]),D[j],as.numeric(X[j,]))
      YY[cnt] <- f(Y[i],Y[j])
    }
  }
  ML::modest(as.data.frame(XX),YY,ML = ML)
}

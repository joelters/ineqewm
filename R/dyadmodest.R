dyadmodest <- function(X,Y,f,X1 = NULL,welfare = c("ineq","igm"),
                       ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL")){
  #create dyads
  n <- length(Y)
  n1 <- n - 1
  XX <- matrix(0,n*(n-1)*0.5,ncol(X)*2)
  YY <- rep(0,n*(n-1)*0.5)
  cnt <- 0
  if (welfare == "ineq"){
    system.time(
    for (i in 1:n1){
      j1 <- i + 1
      for (j in j1:n){
        cnt <- cnt + 1
        XX[cnt,] <- c(as.numeric(X[i,]),as.numeric(X[j,]))
        YY[cnt] <- f(Y[i],Y[j])
      }
    }
    )
  }
  else if (welfare == "igm"){
    system.time(
      for (i in 1:n1){
        j1 <- i + 1
        for (j in j1:n){
          cnt <- cnt + 1
          XX[cnt,] <- c(as.numeric(X[i,]),as.numeric(X[j,]))
          YY[cnt] <- f(Y[i],Y[j],X1[i],X1[j])
        }
      }
    )
  }
  ML::modest(as.data.frame(XX),YY,ML = ML)
}

dyadFVest <- function(model,
                      Xi,
                      Yi,
                      Xnewi,
                      Ynewi,
                      Xj = NULL,
                      Yj = NULL,
                      Xnewj = NULL,
                      Ynewj = NULL,
                      f,
                      ML = c("Lasso", "Ridge", "RF", "CIF",
                             "XGB", "CB","Logit_lasso", "SL"),
                      shape = c("triangle","square")){
  if (shape == "triangle"){
    n <- length(Ynewi)
    n1 <- n - 1
    XX <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2)
    YY <- rep(0,n*(n-1)*0.5)
    XXnew <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2)
    YYnew <- rep(0,n*(n-1)*0.5)
    cnt <- 0
    for (i in 1:n1){
      j1 <- i + 1
      for (j in j1:n){
        cnt <- cnt + 1
        XXnew[cnt,] <- c(as.numeric(Xnewi[i,]),as.numeric(Xnewi[j,]))
        YYnew[cnt] <- f(Ynewi[i],Ynewi[j])
      }
    }
    ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew),YYnew,ML = ML)
  }
  else if (shape == "square"){
    ni <- length(Ynewi)
    nj <- length(Ynewj)
    XX <- matrix(0,ni*nj,ncol(Xnewi)*2)
    YY <- rep(0,ni*nj)
    XXnew <- matrix(0,ni*nj,ncol(Xnewi)*2)
    YYnew <- rep(0,ni*nj)
    cnt <- 0
    for (i in 1:ni){
      for (j in 1:nj){
        cnt <- cnt + 1
        XXnew[cnt,] <- c(as.numeric(Xnewi[i,]),as.numeric(Xnewj[j,]))
        YYnew[cnt] <- f(Ynewi[i],Ynewj[j])
      }
    }
    ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew),YYnew,ML = ML)
  }
}

dyadFVestab <- function(model,
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
    XX <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2 + 2)
    YY <- rep(0,n*(n-1)*0.5)
    XXnew11 <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2 + 2)
    XXnew10 <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2 + 2)
    XXnew01 <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2 + 2)
    XXnew00 <- matrix(0,n*(n-1)*0.5,ncol(Xnewi)*2 + 2)
    YYnew <- rep(0,n*(n-1)*0.5)
    cnt <- 0
    for (i in 1:n1){
      j1 <- i + 1
      for (j in j1:n){
        cnt <- cnt + 1
        XXnew11[cnt,] <- c(1,as.numeric(Xnewi[i,]),1,as.numeric(Xnewi[j,]))
        XXnew10[cnt,] <- c(1,as.numeric(Xnewi[i,]),0,as.numeric(Xnewi[j,]))
        XXnew01[cnt,] <- c(0,as.numeric(Xnewi[i,]),1,as.numeric(Xnewi[j,]))
        XXnew00[cnt,] <- c(0,as.numeric(Xnewi[i,]),0,as.numeric(Xnewi[j,]))
        YYnew[cnt] <- f(Ynewi[i],Ynewi[j])
      }
    }
    fv11 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew11),YYnew,ML = ML)
    fv10 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew10),YYnew,ML = ML)
    fv01 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew01),YYnew,ML = ML)
    fv00 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew00),YYnew,ML = ML)
    return(list("fv11" = fv11,"fv10" = fv10,"fv01" = fv01,"fv00" = fv00))
  }
  else if (shape == "square"){
    ni <- length(Ynewi)
    nj <- length(Ynewj)
    XX <- matrix(0,ni*nj,ncol(Xnewi)*2 + 2)
    YY <- rep(0,ni*nj)
    XXnew11 <- matrix(0,ni*nj,ncol(Xnewi)*2 + 2)
    XXnew10 <- matrix(0,ni*nj,ncol(Xnewi)*2 + 2)
    XXnew01 <- matrix(0,ni*nj,ncol(Xnewi)*2 + 2)
    XXnew00 <- matrix(0,ni*nj,ncol(Xnewi)*2 + 2)
    YYnew <- rep(0,ni*nj)
    cnt <- 0
    for (i in 1:ni){
      for (j in 1:nj){
        cnt <- cnt + 1
        XXnew11[cnt,] <- c(1,as.numeric(Xnewi[i,]),1,as.numeric(Xnewj[j,]))
        XXnew10[cnt,] <- c(1,as.numeric(Xnewi[i,]),0,as.numeric(Xnewj[j,]))
        XXnew01[cnt,] <- c(0,as.numeric(Xnewi[i,]),1,as.numeric(Xnewj[j,]))
        XXnew00[cnt,] <- c(0,as.numeric(Xnewi[i,]),0,as.numeric(Xnewj[j,]))
        YYnew[cnt] <- f(Ynewi[i],Ynewj[j])
      }
    }
    fv11 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew11),YYnew,ML = ML)
    fv10 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew10),YYnew,ML = ML)
    fv01 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew01),YYnew,ML = ML)
    fv00 <- ML::FVest(model,as.data.frame(XX),YY,as.data.frame(XXnew00),YYnew,ML = ML)
    return(list("fv11" = fv11,"fv10" = fv10,"fv01" = fv01,"fv00" = fv00))
  }
}

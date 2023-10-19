#' Kendall-tau Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param X1 outcome with which to compute Kendall-tau with
#' @param D Treatment assignment.
#' @param t Target for Kendall-tau
#' @param design Whether data comes from an RCT or from observational data
#' @param est_method Whether traditional plug in estimators are to be used
#' or locally robust estimators
#' @param MLps Choice of machine learner for propensity score
#' @param MLigm Choice of machine learner for dyadic regression
#' @param CF Whether to use Cross-fitting
#' @param K Number of folds in Cross-fitting
#'
#' @return A list with the output and a figure.
#' @export
wigm_scores_dyad <- function(Y,X1,D,X,pscore,t,
                        design = c("rct","observational"),
                        est_method = c("PI","LR"),
                        MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
                        MLigm = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
                        CF = TRUE,
                        K = 5,
                        leaveoneout = FALSE,
                        verbose = 1){
  warning("It is assummed that X1 belongs in X")
  D <- as.numeric(D)
  Y <- as.numeric(Y)
  n <- length(Y)
  scores11 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores10 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores01 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores00 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  ############### RCT ####################
  if (design == "rct"){
    n1 <- n-1
    w <- 0
    for(i in 1:n1){
      j1 <- i + 1
      for (j in j1:n){
        a11 <- ((D[i]*D[j])/(p^2))
        a10 <- ((D[i]*(1-D[j]))/(p*(1-p)))
        a01 <- (((1-D[i])*D[j])/((1-p)*p))
        a00 <- (((1-D[i])*(1-D[j]))/((1-p)^2))

        scores11[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a11,
                            i,j)
        scores10[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a10,
                            i,j)
        scores01[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a01,
                            i,j)
        scores00[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a00,
                            i,j)
      }
    }
    return(list(scores11,scores10,scores01,scores00))
  }
  ################ Observational ###########
  if (design == "observational"){
    ############### Plug in ###############
    if (est_method == "PI"){
      ps <- ML::MLest(X,as.numeric(D),ML = MLps,FVs = TRUE)
      ps <- ps$FVs
      if (sum(ps<=0 | ps>=1)!= 0){
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
      }
      n1 <- n-1
      cnt <- 0
      for(i in 1:n1){
        # print(i)
        j1 <- i + 1
        for (j in j1:n){
          cnt <- cnt + 1
          a11 <- ((D[i]*D[j])/(ps[i]*ps[j]))
          a10 <- ((D[i]*(1-D[j]))/(ps[i]*(1-ps[j])))
          a01 <- (((1-D[i])*D[j])/((1-ps[i])*ps[j]))
          a00 <- (((1-D[i])*(1-D[j]))/((1-ps[i])*(1-ps[j])))

          scores11[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a11,
                              i,j)
          scores10[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a10,
                              i,j)
          scores01[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a01,
                              i,j)
          scores00[cnt,] <- c(sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*a00,
                              i,j)
        }
      }
      return(list(scores11,scores10,scores01,scores00))
    }
    ################ Locally robust ###########
    else if (est_method == "LR"){
      id <- 1:n
      if (CF == TRUE){
        #Split data
        ind <- split(seq(n), seq(n) %% K)
        #Loop through blocks of pairs
        cnt <- 0
        cnt2 <- 0
        K1 <- ifelse(leaveoneout == TRUE, n-1, K)
        for (ii in 1:K1){
          ii1 <- ifelse(leaveoneout == TRUE, ii + 1,ii)
          K2 <- ifelse(leaveoneout == TRUE, n, K)
          for (jj in ii1:K2){
            cnt2 <- cnt2 + 1
            ########## Obs not in Ci or Cj ###################
            Xnotij <- X[c(-ind[[ii]],-ind[[jj]]),]
            X1notij <- X1[c(-ind[[ii]],-ind[[jj]])]
            Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
            Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
            ########## Train ps model with obs not in Ci or Cj ############
            mps <- ML::modest(Xnotij,as.numeric(Dnotij),ML = MLps)
            g <- function(yi,yj) sign(yi - yj)
            mreg <- dyadmodest(Xnotij,
                               as.numeric(Ynotij),
                               f = g,
                               X1 = as.numeric(X1notij),
                               ML = MLigm)
            ############## Compute scores evaluating in observations in Ci and Cj ##############
            #If we are in a triangle (Ci = Cj)
            if (ii == jj){
              Yii <- Y[ind[[ii]]]
              Xii <- X[ind[[ii]],]
              X1ii <- X1[ind[[ii]]]
              Dii <- D[ind[[ii]]]
              idii <- id[ind[[ii]]]
              psii <- ML::FVest(mps,X = Xnotij, Y = Dnotij,
                                Xnew = Xii, Ynew = Dii,ML = MLps)
              if (sum(psii <= 0 | psii >= 1)!= 0){
                psii <- (psii > 0 & psii < 1)*psii +
                  0.001*(psii <= 0) + 0.999*(psii >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
              }
              regii <- dyadFVest(mreg,
                                 Xi = Xnotij,
                                 Yi = Ynotij,
                                 Xnewi = Xii,
                                 Ynewi = Yii,
                                 f = g,
                                 ML = MLigm,
                                 shape = "triangle")
              cnt3 <- 0
              nn1 <- length(Yii) - 1
              for (kk in 1:nn1){
                kk1 <- kk + 1
                for (mm in kk1:length(Yii)){
                  cnt <- cnt + 1
                  cnt3 <- cnt3 + 1
                  a11 <- ((Dii[kk]*Dii[mm])/(psii[kk]*psii[mm]))
                  a10 <- ((Dii[kk]*(1-Dii[mm]))/(psii[kk]*(1-psii[mm])))
                  a01 <- (((1-Dii[kk])*Dii[mm])/((1-psii[kk])*psii[mm]))
                  a00 <- (((1-Dii[kk])*(1-Dii[mm]))/((1-psii[kk])*(1-psii[mm])))


                  scores11[cnt,] <- c(regii[cnt3]*sign(X1ii[kk]-X1ii[mm]) +
                                      a11*sign(X1ii[kk]-X1ii[mm])*(g(Yii[kk],Yii[mm]) - regii[cnt3]),
                                      idii[kk],idii[mm])
                  scores10[cnt,] <- c(regii[cnt3]*sign(X1ii[kk]-X1ii[mm]) +
                                      a10*sign(X1ii[kk]-X1ii[mm])*(g(Yii[kk],Yii[mm]) - regii[cnt3]),
                                      idii[kk],idii[mm])
                  scores01[cnt,] <- c(regii[cnt3]*sign(X1ii[kk]-X1ii[mm]) +
                                      a01*sign(X1ii[kk]-X1ii[mm])*(g(Yii[kk],Yii[mm]) - regii[cnt3]),
                                      idii[kk],idii[mm])
                  scores00[cnt,] <- c(regii[cnt3]*sign(X1ii[kk]-X1ii[mm]) +
                                      a00*sign(X1ii[kk]-X1ii[mm])*(g(Yii[kk],Yii[mm]) - regii[cnt3]),
                                      idii[kk],idii[mm])
                }
              }
            }
            #If we are in a square (Ci =/= Cj)
            else if (ii < jj){
              #Obs in Ci and Cj stacked
              Yij <- Y[c(ind[[ii]],ind[[jj]])]
              X1ij <- X1[c(ind[[ii]],ind[[jj]])]
              Xij <- X[c(ind[[ii]],ind[[jj]]),]
              YCi <- Y[ind[[ii]]]
              X1Ci <- X1[ind[[ii]]]
              XCi <- X[ind[[ii]],]
              YCj <- Y[ind[[jj]]]
              X1Cj <- X1[ind[[jj]]]
              XCj <- X[ind[[jj]],]
              Dij <- D[c(ind[[ii]],ind[[jj]])]
              idij <- id[c(ind[[ii]],ind[[jj]])]
              psij <- ML::FVest(mps,Xnotij,Dnotij,Xij,Dij,ML = MLps)
              if (sum(psij <= 0 | psij >= 1)!= 0){
                psij <- (psij > 0 & psij < 1)*psij +
                  0.001*(psij <= 0) + 0.999*(psij >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greater or equal than 1 have been set to 0.999")
              }
              regij <- dyadFVest(mreg,
                                 Xi = Xnotij,
                                 Yi = Ynotij,
                                 Xnewi = XCi,
                                 Ynewi = YCi,
                                 Xj = Xnotij,
                                 Yj = Ynotij,
                                 Xnewj = XCj,
                                 Ynewj = YCj,
                                 f = g,
                                 ML = MLigm,
                                 shape = "square")
              cnt3 <- 0
              for (kk in 1:length(ind[[ii]])){
                mm1 <- length(ind[[ii]]) + 1
                for (mm in mm1:length(Yij)){
                  cnt <- cnt + 1
                  cnt3 <- cnt3 + 1
                  a11 <- ((Dij[kk]*Dij[mm])/(psij[kk]*psij[mm]))
                  a10 <- ((Dij[kk]*(1-Dij[mm]))/(psij[kk]*(1-psij[mm])))
                  a01 <- (((1-Dij[kk])*Dij[mm])/((1-psij[kk])*psij[mm]))
                  a00 <- (((1-Dij[kk])*(1-Dij[mm]))/((1-psij[kk])*(1-psij[mm])))

                  scores11[cnt,] <- c(regij[cnt3]*sign(X1ij[kk]-X1ij[mm]) +
                                      a11*sign(X1ij[kk]-X1ij[mm])*(g(Yij[kk],Yij[mm]) - regij[cnt3]),
                                      idij[kk],idij[mm])
                  scores10[cnt,] <- c(regij[cnt3]*sign(X1ij[kk]-X1ij[mm]) +
                                      a10*sign(X1ij[kk]-X1ij[mm])*(g(Yij[kk],Yij[mm]) - regij[cnt3]),
                                      idij[kk],idij[mm])
                  scores01[cnt,] <- c(regij[cnt3]*sign(X1ij[kk]-X1ij[mm]) +
                                      a01*sign(X1ij[kk]-X1ij[mm])*(g(Yij[kk],Yij[mm]) - regij[cnt3]),
                                      idij[kk],idij[mm])
                  scores00[cnt,] <- c(regij[cnt3]*sign(X1ij[kk]-X1ij[mm]) +
                                      a00*sign(X1ij[kk]-X1ij[mm])*(g(Yij[kk],Yij[mm]) - regij[cnt3]),
                                      idij[kk],idij[mm])

                }
              }
            }
            if (verbose == 2){
              print(paste(round((cnt2/((K1+1)*K1/2))*100,1),"% completed.",sep = ""))
            }
          }
          if (verbose == 1){
            print(paste(round((cnt2/((K1+1)*K1/2))*100,1),"% completed.",sep = ""))
          }
        }
        return(list(scores11,scores10,scores01,scores00))
        ###########################################
      }
      else if (CF == FALSE){
        ############# Estimate nuisance parameters ###################3
        ps <- ML::MLest(X,as.numeric(D),ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        if (sum(ps<= 0 | ps >= 1)!=0){
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
          warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
        }
        n <- length(Y)
        delta <- rep(0,n)
        #Estimation of alpha
        for (k in 1:n){
          a11i <- -((D[k]*D[-k])/((ps[k]^2)*ps[-k]))
          a10i <- -((D[k]*(1-D[-k]))/((ps[k]^2)*(1-ps[-k])))
          a01i <- (((1-D[k])*D[-k])/(((1-ps[k])^2)*ps[-k]))
          a00i <- (((1-D[k])*(1-D[-k]))/(((1-ps[k])^2)*(1-ps[-k])))

          delta11i[k] <- mean((sign(Y[k] - Y[-k])*sign(X1[rr]-X1[-rr]))*(a11i))
          delta10i[k] <- mean((sign(Y[k] - Y[-k])*sign(X1[rr]-X1[-rr]))*(a10i))
          delta01i[k] <- mean((sign(Y[k] - Y[-k])*sign(X1[rr]-X1[-rr]))*(a01i))
          delta00i[k] <- mean((sign(Y[k] - Y[-k])*sign(X1[rr]-X1[-rr]))*(a00i))
          #we are not forgetting j, see next comments
        }
        #alpha11i(x) = alpha11j(x)
        alpha11i <- ML::MLest(X,delta11i,ML = MLalpha,FVs = TRUE)
        alpha11i <- alpha11i$FVs
        #alpha10i(x) = alpha01j(x) NOTE THEY CHANGE
        alpha10i <- ML::MLest(X,delta10i,ML = MLalpha,FVs = TRUE)
        alpha10i <- alpha10i$FVs
        #alpha01i(x) = alpha10j(x) NOTE THEY CHANGE
        alpha01i <- ML::MLest(X,delta01i,ML = MLalpha,FVs = TRUE)
        alpha01i <- alpha01i$FVs
        #alpha11i(x) = alpha11j(x)
        alpha00i <- ML::MLest(X,delta00i,ML = MLalpha,FVs = TRUE)
        alpha00i <- alpha00i$FVs
        #so no need to do everything for j
        ################## Compute Scores #################3
        n1 <- n-1
        cnt <- 0
        for(i in 1:n1){
          # print(i)
          j1 <- i + 1
          for (j in j1:n){
            cnt <- cnt + 1
            a11 <- ((D[i]*D[j])/(ps[i]*ps[j]))
            a10 <- ((D[i]*(1-D[j]))/(ps[i]*(1-ps[j])))
            a01 <- (((1-D[i])*D[j])/((1-ps[i])*ps[j]))
            a00 <- (((1-D[i])*(1-D[j]))/((1-ps[i])*(1-ps[j])))
            #we uses alphai for aplhaj since they are equal (see previous comments),
            #so no typos here
            scores11[cnt,] <- c((sign(Y[i] - Y[j])*sign(X1[i]-X1[j]))*(a11) +
                                  alpha11i[i]*(D[i]-ps[i]) + alpha11i[j]*(D[j]-ps[j]),i,j)
            scores10[cnt,] <- c((sign(Y[i] - Y[j])*sign(X1[i]-X1[j]))*(a10) +
                                  alpha10i[i]*(D[i]-ps[i]) + alpha01i[j]*(D[j]- ps[j]),i,j)
            scores01[cnt,] <- c((sign(Y[i] - Y[j])*sign(X1[i]-X1[j]))*(a01) +
                                  alpha01i[i]*(D[i]- ps[i]) + alpha10i[j]*(D[j]-ps[j]),i,j)
            scores00[cnt,] <- c((sign(Y[i] - Y[j])*sign(X1[i]-X1[j]))*(a00) +
                                  alpha00i[i]*(D[i] - ps[i]) + alpha00i[j]*(D[j] - ps[j]),i,j)
          }
        }
        return(list(scores11,scores10,scores01,scores00))
        ##########################
      }
    }
  }
  ###################
}

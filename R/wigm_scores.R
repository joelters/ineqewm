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
#' @param MLalpha Choice of machine learner for additional nuisance parameter
#' @param CF Whether to use Cross-fitting
#' @param K Number of folds in Cross-fitting
#'
#' @return A list with the output and a figure.
#' @export
wigm_scores <- function(Y,X1,D,X,pscore,t,
                 design = c("rct","observational"),
                 est_method = c("PI","LR"),
                 MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
                 MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
                 CF = TRUE,
                 K = 5){
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
        for (ii in 1:K){
          ii1 <- ii
          for (jj in ii1:K){
            ########## Obs not in Ci or Cj ###################
            Xnotij <- X[c(-ind[[ii]],-ind[[jj]]),]
            X1notij <- X1[c(-ind[[ii]],-ind[[jj]])]
            Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
            Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
            ########## Train ps model with obs not in Ci or Cj ############
            mps <- ML::MLest(Xnotij,as.numeric(Dnotij),ML = MLps)
            psnotij <- mps$FVs
            if (sum(psnotij <= 0 | psnotij >= 1)!= 0){
              psnotij <- (psnotij > 0 & psnotij < 1)*psnotij +
                0.001*(psnotij <= 0) + 0.999*(psnotij >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
            }
            mpsnotij <- mps$model
            ######## Train alpha model with obs not in Ci or Cj ##########
            delta11i <- rep(0,length(Ynotij))
            delta10i <- rep(0,length(Ynotij))
            delta01i <- rep(0,length(Ynotij))
            delta00i <- rep(0,length(Ynotij))
            for (rr in 1:length(Ynotij)){
              a11i <- -((Dnotij[rr]*Dnotij[-rr])/
                          ((psnotij[rr]^2)*psnotij[-rr]))
              a10i <- -((Dnotij[rr]*(1-Dnotij[-rr]))/
                          ((psnotij[rr]^2)*(1-psnotij[-rr])))
              a01i <- (((1-Dnotij[rr])*Dnotij[-rr])/
                         (((1-psnotij[rr])^2)*psnotij[-rr]))
              a00i <- (((1-Dnotij[rr])*(1-Dnotij[-rr]))/
                         (((1-psnotij[rr])^2)*(1-psnotij[-rr])))
              delta11i[rr] <- mean((sign(Ynotij[rr] - Ynotij[-rr])*sign(X1notij[rr]-X1notij[-rr]))*(a11i))
              delta10i[rr] <- mean((sign(Ynotij[rr] - Ynotij[-rr])*sign(X1notij[rr]-X1notij[-rr]))*(a10i))
              delta01i[rr] <- mean((sign(Ynotij[rr] - Ynotij[-rr])*sign(X1notij[rr]-X1notij[-rr]))*(a01i))
              delta00i[rr] <- mean((sign(Ynotij[rr] - Ynotij[-rr])*sign(X1notij[rr]-X1notij[-rr]))*(a00i))
            }
            malphanotij11 <- ML::modest(Xnotij,delta11i,ML = MLalpha)
            malphanotij10 <- ML::modest(Xnotij,delta10i,ML = MLalpha)
            malphanotij01 <- ML::modest(Xnotij,delta01i,ML = MLalpha)
            malphanotij00 <- ML::modest(Xnotij,delta00i,ML = MLalpha)
            ############## Compute scores evaluating in observations in Ci and Cj ##############
            #If we are in a triangle (Ci = Cj)
            if (ii == jj){
              Yii <- Y[ind[[ii]]]
              Xii <- X[ind[[ii]],]
              X1ii <- X1[ind[[ii]]]
              Dii <- D[ind[[ii]]]
              idii <- id[ind[[ii]]]
              psii <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xii,Dii,ML = MLps)
              if (sum(psii <= 0 | psii >= 1)!= 0){
                psii <- (psii > 0 & psii < 1)*psii +
                  0.001*(psii <= 0) + 0.999*(psii >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
              }
              alphaii11i <- ML::FVest(malphanotij11,Xnotij,delta11,
                                      Xii,delta11,ML = MLalpha)
              alphaii10i <- ML::FVest(malphanotij10,Xnotij,delta10,
                                      Xii,delta10,ML = MLalpha)
              alphaii01i <- ML::FVest(malphanotij01,Xnotij,delta01,
                                      Xii,delta01,ML = MLalpha)
              alphaii00i <- ML::FVest(malphanotij00,Xnotij,delta00,
                                      Xii,delta00,ML = MLalpha)
              nn1 <- length(Yii) - 1
              for (kk in 1:nn1){
                kk1 <- kk + 1
                for (mm in kk1:length(Yii)){
                  cnt <- cnt + 1
                  a11 <- ((Dii[kk]*Dii[mm])/(psii[kk]*psii[mm]))
                  a10 <- ((Dii[kk]*(1-Dii[mm]))/(psii[kk]*(1-psii[mm])))
                  a01 <- (((1-Dii[kk])*Dii[mm])/((1-psii[kk])*psii[mm]))
                  a00 <- (((1-Dii[kk])*(1-Dii[mm]))/((1-psii[kk])*(1-psii[mm])))

                  scores11[cnt,] <- c(sign(Yii[kk] - Yii[mm])*sign(X1ii[kk]-X1ii[mm])*(a11) +
                                        alphaii11i[kk]*(Dii[kk]-psii[kk]) + alphaii11i[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                  scores10[cnt,] <- c(sign(Yii[kk] - Yii[mm])*sign(X1ii[kk]-X1ii[mm])*(a10) +
                                        alphaii10i[kk]*(Dii[kk]-psii[kk]) + alphaii01i[mm]*(Dii[mm] - psii[mm]),
                                      idii[kk],idii[mm])
                  scores01[cnt,] <- c(sign(Yii[kk] - Yii[mm])*sign(X1ii[kk]-X1ii[mm])*(a01) +
                                        alphaii01i[kk]*(1-Dii[kk] - 1 +psii[kk]) + alphaii10i[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                  scores00[cnt,] <- c(sign(Yii[kk] - Yii[mm])*sign(X1ii[kk]-X1ii[mm])*(a00) +
                                        alphaii00i[kk]*(Dii[kk] - psii[kk]) + alphaii00i[mm]*(Dii[mm] - psii[mm]),
                                      idii[kk],idii[mm])
                }
              }
            }
            #If we are in a square (Ci =/= Cj)
            else if (ii < jj){
              #Obs in Ci and Cj stacked
              Yij <- Y[c(ind[[ii]],ind[[jj]])]
              Xij <- X[c(ind[[ii]],ind[[jj]]),]
              X1ij <- X1[c(ind[[ii]],ind[[jj]])]
              Dij <- D[c(ind[[ii]],ind[[jj]])]
              idij <- id[c(ind[[ii]],ind[[jj]])]
              psij <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xij,Dij,ML = MLps)
              if (sum(psij<= 0 | psij >= 1)!= 0){
                psij <- (psij > 0 & psij < 1)*psij +
                  0.001*(psij <= 0) + 0.999*(psij >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
              }
              alphaij11i <- ML::FVest(malphanotij11,Xnotij,delta11,
                                      Xij,delta11,ML = MLalpha)
              alphaij10i <- ML::FVest(malphanotij10,Xnotij,delta10,
                                      Xij,delta10,ML = MLalpha)
              alphaij01i <- ML::FVest(malphanotij01,Xnotij,delta01,
                                      Xij,delta01,ML = MLalpha)
              alphaij00i <- ML::FVest(malphanotij00,Xnotij,delta00,
                                      Xij,delta00,ML = MLalpha)
              for (kk in 1:length(ind[[ii]])){
                mm1 <- length(ind[[ii]]) + 1
                for (mm in mm1:length(Yij)){
                  cnt <- cnt + 1
                  a11 <- ((Dij[kk]*Dij[mm])/(psij[kk]*psij[mm]))
                  a10 <- ((Dij[kk]*(1-Dij[mm]))/(psij[kk]*(1-psij[mm])))
                  a01 <- (((1-Dij[kk])*Dij[mm])/((1-psij[kk])*psij[mm]))
                  a00 <- (((1-Dij[kk])*(1-Dij[mm]))/((1-psij[kk])*(1-psij[mm])))

                  scores11[cnt,] <- c(sign(Yij[kk] - Yij[mm])*sign(X1ij[kk]-X1ij[mm])*(a11) +
                                        alphaij11i[kk]*(Dij[kk]-psij[kk]) + alphaij11i[mm]*(Dij[mm]-psij[mm]),
                                      idij[kk],idij[mm])
                  scores10[cnt,] <- c(sign(Yij[kk] - Yij[mm])*sign(X1ij[kk]-X1ij[mm])*(a10) +
                                        alphaij10i[kk]*(Dij[kk]-psij[kk]) + alphaij01i[mm]*(Dij[mm] - psij[mm]),
                                      idij[kk],idij[mm])
                  scores01[cnt,] <- c(sign(Yij[kk] - Yij[mm])*sign(X1ij[kk]-X1ij[mm])*(a01) +
                                        alphaij01i[kk]*(Dij[kk] - psij[kk]) + alphaij10i[mm]*(Dij[mm]-psij[mm]),
                                      idij[kk],idij[mm])
                  scores00[cnt,] <- c(sign(Yij[kk] - Yij[mm])*sign(X1ij[kk]-X1ij[mm])*(a00) +
                                        alphaij00i[kk]*(Dij[kk] - psij[kk]) + alphaij00i[mm]*(Dij[mm] - psij[mm]),
                                      idij[kk],idij[mm])

                }
              }
            }
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

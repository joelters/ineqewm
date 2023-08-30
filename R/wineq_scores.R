#' Standard Gini Social Welfare score computation
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param design Whether data comes from an RCT or from observational data
#' @param est_method Whether traditional plug in estimators are to be used
#' or locally robust estimators
#' @param MLps Choice of machine learner for propensity score
#' @param MLalpha Choice of machine learner for additional nuisance parameter
#' @param CF Whether to use Cross-fitting
#' @param K Number of folds in Cross-fitting
#'
#' @return A matrix with scores and id of the pair
#' @export
wineq_scores <- function(Y,D,X,
                         design = c("rct","observational"),
                         est_method = c("PI","LR"),
                         MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                         MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
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
    p <- mean(D)
    n1 <- n-1
    w1 <- 0
    cnt <- 0
    for(i in 1:n1){
      # print(i)
      j1 <- i + 1
      for (j in j1:n){
        cnt <- cnt + 1
        a11 <- ((D[i]*D[j])/(p^2))
        a10 <- ((D[i]*(1-D[j]))/(p*(1-p)))
        a01 <- (((1-D[i])*D[j])/((1-p)*p))
        a00 <- (((1-D[i])*(1-D[j]))/((1-p)^2))
        scores11[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a11,
                            i,j)
        scores10[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a10,
                            i,j)
        scores01[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a01,
                            i,j)
        scores00[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a00,
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
      if (sum(ps <= 0 | ps >= 1) != 0){
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
      }
      n <- length(Y)
      n1 <- n-1
      w1 <- 0
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
          scores11[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a11,
                              i,j)
          scores10[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a10,
                              i,j)
          scores01[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a01,
                              i,j)
          scores00[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*a00,
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
        n <- length(Y)
        ind <- split(seq(n), seq(n) %% K)
        #Loop through blocks of pairs
        w <- 0
        cnt <- 0
        for (ii in 1:K){
          ii1 <- ii
          for (jj in ii1:K){
            ########## Obs not in Ci or Cj ###################
            Xnotij <- X[c(-ind[[ii]],-ind[[jj]]),]
            Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
            Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
            ########## Train ps model with obs not in Ci or Cj ############
            mps <- ML::MLest(Xnotij,as.numeric(Dnotij),ML = MLps)
            psnotij <- mps$FVs
            if (sum(psnotij <= 0 | psnotij >= 1) != 0){
              psnotij <- (psnotij > 0 & psnotij < 1)*psnotij +
                0.001*(psnotij <= 0) + 0.999*(psnotij >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
            }

            mpsnotij <- mps$model
            ######## Train alpha model with obs not in Ci or Cj ##########
            delta11 <- rep(0,length(Ynotij))
            delta10 <- rep(0,length(Ynotij))
            delta01 <- rep(0,length(Ynotij))
            delta00 <- rep(0,length(Ynotij))
            for (rr in 1:length(Ynotij)){
              a11 <- ((Dnotij[rr]*Dnotij[-rr])/
                        ((psnotij[rr]^2)*psnotij[-rr]))
              a10 <- ((Dnotij[rr]*(1-Dnotij[-rr]))/
                        ((psnotij[rr]^2)*(1-psnotij[-rr])))
              a01 <- (((1-Dnotij[rr])*Dnotij[-rr])/
                        (((1-psnotij[rr])^2)*psnotij[-rr]))
              a00 <- (((1-Dnotij[rr])*(1-Dnotij[-rr]))/
                        (((1-psnotij[rr])^2)*(1-psnotij[-rr])))
              delta11[rr] <- mean(0.5*(Ynotij[rr] + Ynotij[-rr] -
                                         abs(Ynotij[rr] - Ynotij[-rr]))*(-a11))
              delta10[rr] <- mean(0.5*(Ynotij[rr] + Ynotij[-rr] -
                                         abs(Ynotij[rr] - Ynotij[-rr]))*(-a10))
              delta01[rr] <- mean(0.5*(Ynotij[rr] + Ynotij[-rr] -
                                         abs(Ynotij[rr] - Ynotij[-rr]))*(a01))
              delta00[rr] <- mean(0.5*(Ynotij[rr] + Ynotij[-rr] -
                                         abs(Ynotij[rr] - Ynotij[-rr]))*(a00))
            }
            malphanotij11 <- ML::modest(Xnotij,delta11,ML = MLalpha)
            malphanotij10 <- ML::modest(Xnotij,delta10,ML = MLalpha)
            malphanotij01 <- ML::modest(Xnotij,delta01,ML = MLalpha)
            malphanotij00 <- ML::modest(Xnotij,delta00,ML = MLalpha)
            ############## Compute welfare evaluating in observations in Ci and Cj ##############
            #If we are in a triangle (Ci = Cj)
            if (ii == jj){
              Yii <- Y[ind[[ii]]]
              Xii <- X[ind[[ii]],]
              Dii <- D[ind[[ii]]]
              idii <- id[ind[[ii]]]
              psii <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xii,Dii,ML = MLps)
              if (sum(psii <= 0 | psii >= 1) != 0){
                psii <- (psii > 0 & psii < 1)*psii +
                  0.001*(psii <= 0) + 0.999*(psii >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
              }
              alphaii11 <- ML::FVest(malphanotij11,Xnotij,delta11,
                                     Xii,delta11,ML = MLalpha)
              alphaii10 <- ML::FVest(malphanotij10,Xnotij,delta10,
                                     Xii,delta10,ML = MLalpha)
              alphaii01 <- ML::FVest(malphanotij01,Xnotij,delta01,
                                     Xii,delta01,ML = MLalpha)
              alphaii00 <- ML::FVest(malphanotij00,Xnotij,delta00,
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
                  scores11[cnt,] <- c(0.5*(Yii[kk] + Yii[mm] - abs(Yii[kk] - Yii[mm]))*(a11) +
                                        alphaii11[kk]*(Dii[kk]-psii[kk]) + alphaii11[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                  scores10[cnt,] <- c(0.5*(Yii[kk] + Yii[mm] - abs(Yii[kk] - Yii[mm]))*(a10) +
                                        alphaii10[kk]*(Dii[kk]-psii[kk]) + alphaii10[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                  scores01[cnt,] <- c(0.5*(Yii[kk] + Yii[mm] - abs(Yii[kk] - Yii[mm]))*(a01) +
                                        alphaii01[kk]*(Dii[kk]-psii[kk]) + alphaii01[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                  scores00[cnt,] <- c(0.5*(Yii[kk] + Yii[mm] - abs(Yii[kk] - Yii[mm]))*(a00) +
                                        alphaii00[kk]*(Dii[kk]-psii[kk]) + alphaii00[mm]*(Dii[mm]-psii[mm]),
                                      idii[kk],idii[mm])
                }
              }
            }
            #If we are in a square (Ci =/= Cj)
            else if (ii < jj){
              #Obs in Ci and Cj stacked
              Yij <- Y[c(ind[[ii]],ind[[jj]])]
              Xij <- X[c(ind[[ii]],ind[[jj]]),]
              Dij <- D[c(ind[[ii]],ind[[jj]])]
              idij <- id[c(ind[[ii]],ind[[jj]])]
              psij <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xij,Dij,ML = MLps)
              if (sum(psij <= 0 | psij >= 1)!= 0){
                psij <- (psij > 0 & psij < 1)*psij +
                  0.001*(psij <= 0) + 0.999*(psij >= 1)
                warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
              }
              alphaij11 <- ML::FVest(malphanotij11,Xnotij,delta11,
                                     Xij,delta11,ML = MLalpha)
              alphaij10 <- ML::FVest(malphanotij10,Xnotij,delta10,
                                     Xij,delta10,ML = MLalpha)
              alphaij01 <- ML::FVest(malphanotij01,Xnotij,delta01,
                                     Xij,delta01,ML = MLalpha)
              alphaij00 <- ML::FVest(malphanotij00,Xnotij,delta00,
                                     Xij,delta00,ML = MLalpha)
              for (kk in 1:length(ind[[ii]])){
                mm1 <- length(ind[[ii]]) + 1
                for (mm in mm1:length(Yij)){
                  cnt <- cnt + 1
                  a11 <- ((Dij[kk]*Dij[mm])/(psij[kk]*psij[mm]))
                  a10 <- ((Dij[kk]*(1-Dij[mm]))/(psij[kk]*(1-psij[mm])))
                  a01 <- (((1-Dij[kk])*Dij[mm])/((1-psij[kk])*psij[mm]))
                  a00 <- (((1-Dij[kk])*(1-Dij[mm]))/((1-psij[kk])*(1-psij[mm])))
                  scores11[cnt,] <- c(0.5*(Yij[kk] + Yij[mm] - abs(Yij[kk] - Yij[mm]))*(a11) +
                                        alphaij11[kk]*(Dij[kk]-psij[kk]) + alphaij11[mm]*(Dij[mm]-psij[mm]),
                                      idij[kk],idij[mm])
                  scores10[cnt,] <- c(0.5*(Yij[kk] + Yij[mm] - abs(Yij[kk] - Yij[mm]))*(a10) +
                                        alphaij10[kk]*(Dij[kk]-psij[kk]) + alphaij10[mm]*(Dij[mm]-psij[mm]),
                                      idij[kk],idij[mm])
                  scores01[cnt,] <- c(0.5*(Yij[kk] + Yij[mm] - abs(Yij[kk] - Yij[mm]))*(a01) +
                                        alphaij01[kk]*(Dij[kk]-psij[kk]) + alphaij01[mm]*(Dij[mm]-psij[mm]),
                                      idij[kk],idij[mm])
                  scores00[cnt,] <- c(0.5*(Yij[kk] + Yij[mm] - abs(Yij[kk] - Yij[mm]))*(a00) +
                                        alphaij00[kk]*(Dij[kk]-psij[kk]) + alphaij00[mm]*(Dij[mm]-psij[mm]),
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
        if (sum(ps <= 0 | ps >= 1)!=0){
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
          warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
        }
        n <- length(Y)
        delta11 <- rep(0,n)
        delta10 <- rep(0,n)
        delta01 <- rep(0,n)
        delta00 <- rep(0,n)
        #Estimation of alpha
        for (i in 1:n){
          a11 <- ((D[i]*D[-i])/((ps[i]^2)*ps[-i]))
          a10 <- ((D[i]*(1-D[-i]))/((ps[i]^2)*(1-ps[-i])))
          a01 <- (((1-D[i])*D[-i])/(((1-ps[i])^2)*ps[-i]))
          a00 <- (((1-D[i])*(1-D[-i]))/(((1-ps[i])^2)*(1-ps[-i])))
          delta11[i] <- mean(0.5*(Y[i] + Y[-i] - abs(Y[i] - Y[-i]))*(-a11))
          delta10[i] <- mean(0.5*(Y[i] + Y[-i] - abs(Y[i] - Y[-i]))*(-a10))
          delta01[i] <- mean(0.5*(Y[i] + Y[-i] - abs(Y[i] - Y[-i]))*(a01))
          delta00[i] <- mean(0.5*(Y[i] + Y[-i] - abs(Y[i] - Y[-i]))*(a00))
        }
        alpha11 <- ML::MLest(X,delta11,ML = MLalpha,FVs = TRUE)
        alpha11 <- alpha$FVs
        alpha10 <- ML::MLest(X,delta10,ML = MLalpha,FVs = TRUE)
        alpha10 <- alpha$FVs
        alpha01 <- ML::MLest(X,delta01,ML = MLalpha,FVs = TRUE)
        alpha01 <- alpha$FVs
        alpha00 <- ML::MLest(X,delta00,ML = MLalpha,FVs = TRUE)
        alpha00 <- alpha$FVs
        ################## Compute Welfare #################
        n1 <- n-1
        w <- 0
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
            scores11[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a11) +
                                  alpha11[i]*(D[i]-ps[i]) + alpha11[j]*(D[j]-ps[j]),i,j)
            scores10[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a10) +
                                  alpha10[i]*(D[i]-ps[i]) + alpha10[j]*(D[j]-ps[j]),i,j)
            scores01[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a01) +
                                  alpha01[i]*(D[i]-ps[i]) + alpha01[j]*(D[j]-ps[j]),i,j)
            scores00[cnt,] <- c(0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a00) +
                                  alpha00[i]*(D[i]-ps[i]) + alpha00[j]*(D[j]-ps[j]),i,j)
          }
        }
        return(list(scores11,scores10,scores01,scores00))
        ##########################
      }
    }
  }
  ###################
}

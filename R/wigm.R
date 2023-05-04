#' Kendall-tau Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param rule Treatment rule.
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
wigm <- function(Y,X1,D,X,rule,t,
                  design = c("rct","observational"),
                  est_method = c("PI","LR"),
                  MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  CF = TRUE,
                  K = 5){
  D <- as.numeric(D)
  Y <- as.numeric(Y)
  ############### RCT ####################
  if (design == "rct"){
    p <- mean(D)
    n <- length(Y)
    n1 <- n-1
    w <- 0
    for(i in 1:n1){
      # print(i)
      j1 <- i + 1
      for (j in j1:n){
        a11 <- ((D[i]*D[j])/(p^2))*(rule[i]*rule[j])
        a10 <- ((D[i]*(1-D[j]))/(p*(1-p)))*(rule[i]*(1-rule[j]))
        a01 <- (((1-D[i])*D[j])/((1-p)*p))*((1-rule[i])*rule[j])
        a00 <- (((1-D[i])*(1-D[j]))/((1-p)^2))*((1-rule[i])*(1-rule[j]))
        w <- w + sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*
          (a11+a10+a01+a00)
      }
    }
    WT <- -abs((2/(n*(n-1)))*w-t)
    mu <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
    return(list("Welfare" = WT, "Mean" = mu))
  }
  ################ Observational ###########
  if (design == "observational"){
    ############### Plug in ###############
    if (est_method == "PI"){
      ps <- ML::MLest(X,as.numeric(D),ML = MLps,FVs = TRUE)
      ps <- ps$FVs
      ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps == 0) + 0.999*(ps >= 1)
      n <- length(Y)
      n1 <- n-1
      w <- 0
      for(i in 1:n1){
        # print(i)
        j1 <- i + 1
        for (j in j1:n){
          a11 <- ((D[i]*D[j])/(ps[i]*ps[j]))*(rule[i]*rule[j])
          a10 <- ((D[i]*(1-D[j]))/(ps[i]*(1-ps[j])))*(rule[i]*(1-rule[j]))
          a01 <- (((1-D[i])*D[j])/((1-ps[i])*ps[j]))*((1-rule[i])*rule[j])
          a00 <- (((1-D[i])*(1-D[j]))/((1-ps[i])*(1-ps[j])))*((1-rule[i])*(1-rule[j]))
          w <- w + sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*
            (a11+a10+a01+a00)
        }
      }
      WT <- -abs((2/(n*(n-1)))*w-t)
      mu <- mean(((Y*D)/ps - Y*(1-D)/(1-ps))*rule + Y*(1-D)/(1-ps))
      return(list("Welfare" = WT, "Mean" = mu))
    }
    ################ Locally robust ###########
    else if (est_method == "LR"){
      if (CF == TRUE){
        #Split data
        n <- length(Y)
        ind <- split(seq(n), seq(n) %% K)
        #Loop through blocks of pairs
        w <- 0
        for (ii in 1:K){
          ii1 <- ii
          for (jj in ii1:K){
            ########## Obs not in Ci or Cj ###################
            Xnotij <- X[c(-ind[[ii]],-ind[[jj]]),]
            X1notij <- X1[c(-ind[[ii]],-ind[[jj]])]
            Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
            Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
            rulenotij <- rule[c(-ind[[ii]],-ind[[jj]])]
            ########## Train ps model with obs not in Ci or Cj ############
            mps <- ML::MLest(Xnotij,as.numeric(Dnotij),ML = MLps)
            psnotij <- mps$FVs
            psnotij <- (psnotij > 0 & psnotij < 1)*psnotij +
              0.001*(psnotij == 0) + 0.999*(psnotij >= 1)
            mpsnotij <- mps$model
            ######## Train alpha model with obs not in Ci or Cj ##########
            delta <- rep(0,length(Ynotij))
            for (rr in 1:length(Ynotij)){
              a11 <- ((Dnotij[rr]*Dnotij[-rr])/
                        ((psnotij[rr]^2)*psnotij[-rr]))*(rulenotij[rr]*rulenotij[-rr])
              a10 <- ((Dnotij[rr]*(1-Dnotij[-rr]))/
                        ((psnotij[rr]^2)*(1-psnotij[-rr])))*(rulenotij[rr]*(1-rulenotij[-rr]))
              a01 <- (((1-Dnotij[rr])*Dnotij[-rr])/
                        (((1-psnotij[rr])^2)*psnotij[-rr]))*((1-rulenotij[rr])*rulenotij[-rr])
              a00 <- (((1-Dnotij[rr])*(1-Dnotij[-rr]))/
                        (((1-psnotij[rr])^2)*(1-psnotij[-rr])))*((1-rulenotij[rr])*(1-rulenotij[-rr]))
              delta[rr] <- mean(sign(Ynotij[rr] - Ynotij[-rr])*sign(X1notij[rr]-X1notij[-rr])*
                                  (a11+a10+a01+a00))

            }
            malphanotij <- ML::modest(Xnotij,delta,ML = MLalpha)
            ############## Compute welfare evaluating in observations in Ci and Cj ##############
            #If we are in a triangle (Ci = Cj)
            if (ii == jj){
              Yii <- Y[ind[[ii]]]
              Xii <- X[ind[[ii]],]
              X1ii <- X1[ind[[ii]]]
              Dii <- D[ind[[ii]]]
              ruleii <- rule[ind[[ii]]]
              psii <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xii,Dii,ML = MLps)
              psii <- (psii > 0 & psii < 1)*psii +
                0.001*(psii == 0) + 0.999*(psii >= 1)
              alphaii <- ML::FVest(malphanotij,Xnotij,delta,
                                   Xii,delta,ML = MLalpha)
              nn1 <- length(Yii) - 1
              for (kk in 1:nn1){
                kk1 <- kk + 1
                for (mm in kk1:length(Yii)){
                  a11 <- ((Dii[kk]*Dii[mm])/(psii[kk]*psii[mm]))*(ruleii[kk]*ruleii[mm])
                  a10 <- ((Dii[kk]*(1-Dii[mm]))/(psii[kk]*(1-psii[mm])))*(ruleii[kk]*(1-ruleii[mm]))
                  a01 <- (((1-Dii[kk])*Dii[mm])/((1-psii[kk])*psii[mm]))*((1-ruleii[kk])*ruleii[mm])
                  a00 <- (((1-Dii[kk])*(1-Dii[mm]))/((1-psii[kk])*(1-psii[mm])))*((1-ruleii[kk])*(1-ruleii[mm]))
                  w <- w +  sign(Yii[kk] - Yii[mm])*sign(X1ii[kk]-X1ii[mm])*(a11+a10+a01+a00) +
                    alphaii[kk]*(Dii[kk]-psii[kk]) + alphaii[mm]*(Dii[mm]-psii[mm])
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
              ruleij <- rule[c(ind[[ii]],ind[[jj]])]
              psij <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xij,Dij,ML = MLps)
              psij <- (psij > 0 & psij < 1)*psij +
                0.001*(psij == 0) + 0.999*(psij >= 1)
              alphaij <- ML::FVest(malphanotij,Xnotij,delta,
                                   Xij,delta,ML = MLalpha)
              for (kk in 1:length(ind[[ii]])){
                mm1 <- length(ind[[ii]]) + 1
                for (mm in mm1:length(Yij)){
                  a11 <- ((Dij[kk]*Dij[mm])/(psij[kk]*psij[mm]))*(ruleij[kk]*ruleij[mm])
                  a10 <- ((Dij[kk]*(1-Dij[mm]))/(psij[kk]*(1-psij[mm])))*(ruleij[kk]*(1-ruleij[mm]))
                  a01 <- (((1-Dij[kk])*Dij[mm])/((1-psij[kk])*psij[mm]))*((1-ruleij[kk])*ruleij[mm])
                  a00 <- (((1-Dij[kk])*(1-Dij[mm]))/((1-psij[kk])*(1-psij[mm])))*((1-ruleij[kk])*(1-ruleij[mm]))
                  w <- w + sign(Yij[kk] - Yij[mm])*sign(X1ij[kk]-X1ij[mm])*(a11+a10+a01+a00) +
                    alphaij[kk]*(Dij[kk]-psij[kk]) + alphaij[mm]*(Dij[mm]-psij[mm])

                }
              }
            }
            WT <- -abs((2/(n*(n-1)))*w-t)
            mu <- wutil(Y = Y,D = D,X = X,rule = rule,
                        design = design,
                        est_method = est_method,
                        MLps = MLps,
                        MLreg = MLalpha,
                        CF = CF, K = K)
            mu <- mu$Welfare
            return(list("Welfare" = WT, "Mean" = mu))
            ###########################################
          }
        }
      }
      else if (CF == FALSE){
        ############# Estimate nuisance parameters ###################3
        ps <- ML::MLest(X,as.numeric(D),ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps == 0) + 0.999*(ps >= 1)
        n <- length(Y)
        delta <- rep(0,n)
        #Estimation of alpha
        for (i in 1:n){
          a11 <- ((D[i]*D[-i])/((ps[i]^2)*ps[-i]))*(rule[i]*rule[-i])
          a10 <- ((D[i]*(1-D[-i]))/((ps[i]^2)*(1-ps[-i])))*(rule[i]*(1-rule[-i]))
          a01 <- (((1-D[i])*D[-i])/(((1-ps[i])^2)*ps[-i]))*((1-rule[i])*rule[-i])
          a00 <- (((1-D[i])*(1-D[-i]))/(((1-ps[i])^2)*(1-ps[-i])))*((1-rule[i])*(1-rule[-i]))
          delta[i] <- mean(sign(Y[i] - Y[-i])*sign(X1[i]-X1[-i])*(a11+a10+a01+a00))
        }
        alpha <- ML::MLest(X,delta,ML = MLalpha,FVs = TRUE)
        alpha <- alpha$FVs
        ################## Compute Welfare #################3
        n1 <- n-1
        w <- 0
        for(i in 1:n1){
          # print(i)
          j1 <- i + 1
          for (j in j1:n){
            a11 <- ((D[i]*D[j])/(ps[i]*ps[j]))*(rule[i]*rule[j])
            a10 <- ((D[i]*(1-D[j]))/(ps[i]*(1-ps[j])))*(rule[i]*(1-rule[j]))
            a01 <- (((1-D[i])*D[j])/((1-ps[i])*ps[j]))*((1-rule[i])*rule[j])
            a00 <- (((1-D[i])*(1-D[j]))/((1-ps[i])*(1-ps[j])))*((1-rule[i])*(1-rule[j]))
            w <- w + sign(Y[i] - Y[j])*sign(X1[i]-X1[j])*(a11+a10+a01+a00) +
              alpha[i]*(D[i]-ps[i]) + alpha[j]*(D[j]-ps[j])
          }
        }
        WT <--abs((2/(n*(n-1)))*w-t)
        mu <- wutil(Y = Y,D = D,X = X,rule = rule,
                    design = design,
                    est_method = est_method,
                    MLps = MLps,
                    MLreg = MLalpha,
                    CF = CF, K = K)
        mu <- mu$Welfare
        return(list("Welfare" = WT, "Mean" = mu))
        ##########################
      }
    }
  }
  ###################
}

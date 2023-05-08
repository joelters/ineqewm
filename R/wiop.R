#' Inequality of Opportunity Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X Tibble with circumstances.
#' @param rule Treatment rule.
#' @param design Whether data comes from an RCT or from observational data
#' @param est_method Whether traditional plug in estimators are to be used
#' or locally robust estimators
#' @param MLiop Choice of Machine Learner for outcome regression
#' @param MLps Choice of Machine Learner for propensity score
#'
#' @return A list with the output and a figure.
#' @export
wiop <- function(Y,D,X,rule,
                 design = c("rct","observational"),
                 est_method = c("PI","LR"),
                 MLiop = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                 MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                 MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                 CF = TRUE,
                 K = 5){
  D <- as.numeric(D)
  Y <- as.numeric(Y)
  if (est_method == "PI"){
    n <- length(Y)
    DX <- as_tibble(cbind(D = D,X))
    mdx <- ML::MLest(DX,Y,ML = MLiop,FVs = TRUE)
    D1X <- as_tibble(cbind(D = rep(1,nrow(X)),X))
    D0X <- as_tibble(cbind(D = rep(0,nrow(X)),X))
    fv1 <- ML::FVest(mdx$model,DX,Y,D1X,Y,ML = MLiop)
    fv0 <- ML::FVest(mdx$model,DX,Y,D0X,Y,ML = MLiop)
    fvt <- fv1*rule + fv0*(1-rule)
    n1 <- n-1
    w <- 0
    mu <- 0
    for(i in 1:n1){
      # print(i)
      j1 <- i + 1
      for (j in j1:n){
        mu <- mu + 0.5*(fvt[i] + fvt[j])
        w <- w + 0.5*(- abs(fvt[i]-fvt[j]))
      }
    }
    WT <- (2/(n*(n-1)))*(w+mu)
    mu <- (2/(n*(n-1)))*(mu)
    iop <- 1 - WT/mu
  }
  else if (est_method == "LR"){
    if (CF == TRUE){
      #Split data
      n <- length(Y)
      ind <- split(seq(n), seq(n) %% K)
      #Loop through blocks of pairs
      jt <- 0
      for (ii in 1:K){
        ii1 <- ii
        for (jj in ii1:K){
          ########## Obs not in Ci or Cj ###################
          Xnotij <- X[c(-ind[[ii]],-ind[[jj]]),]
          Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
          DXnotij <- as_tibble(cbind(D = Dnotij,Xnotij))
          Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
          rulenotij <- rule[c(-ind[[ii]],-ind[[jj]])]
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
          ########## Train fvs iop model with obs not in Ci or Cj ########
          mdx <- ML::MLest(DXnotij,Ynotij,ML = MLiop,FVs = TRUE)
          mdx <- mdx$model
          ############## Compute welfare evaluating in observations in Ci and Cj ##############
          #If we are in a triangle (Ci = Cj)
          if (ii == jj){
            Yii <- Y[ind[[ii]]]
            Xii <- X[ind[[ii]],]
            Dii <- D[ind[[ii]]]
            DXii <- as_tibble(cbind(D = Dii,Xii))
            D1Xii <- as_tibble(cbind(D = rep(1,nrow(Xii)),Xii))
            D0Xii <- as_tibble(cbind(D = rep(0,nrow(Xii)),Xii))
            ruleii <- rule[ind[[ii]]]
            psii <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xii,Dii,ML = MLps)
            if (sum(psii <= 0 | psii >= 1)!= 0){
              psii <- (psii > 0 & psii < 1)*psii +
              0.001*(psii <= 0) + 0.999*(psii >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
            }
            fv1ii <- ML::FVest(mdx,DXnotij,Ynotij,D1Xii,Yii,ML = MLiop)
            fv0ii <- ML::FVest(mdx,DXnotij,Ynotij,D0Xii,Yii,ML = MLiop)
            fvdxii <- ML::FVest(mdx,DXnotij,Ynotij,DXii,Yii,ML = MLiop)
            fvtii <- fv1ii*ruleii + fv0ii*(1-ruleii)
            wii <- ruleii*Dii/psii + (1-ruleii)*(1-Dii)/(1-psii)
            a1 <- 1/2
            c11 <- 1
            c12 <- 1
            c21 <- 1
            c22 <- -1
            n1 <- n-1
            nn1 <- length(Yii) - 1
            for (kk in 1:nn1){
              kk1 <- kk + 1
              for (mm in kk1:length(Yii)){
                m <- 0.5*(fvtii[kk] + fvtii[mm] - abs(fvtii[kk]-fvtii[mm]))
                a2 <- ((fvtii[kk]>fvtii[mm]) - (fvtii[kk] < fvtii[mm]))/2
                phi1 <- a1*(c11*wii[kk]*(Yii[kk] - fvdxii[kk]) + c12*wii[mm]*(Yii[mm] - fvdxii[mm]))
                phi2 <- a2*(c21*wii[kk]*(Yii[kk] - fvdxii[kk]) + c22*wii[mm]*(Yii[mm] - fvdxii[mm]))
                jt <- jt + (m + phi1 + phi2)
              }
            }
          }
          #If we are in a square (Ci =/= Cj)
          else if (ii < jj){
            #Obs in Ci and Cj stacked
            Yij <- Y[c(ind[[ii]],ind[[jj]])]
            Xij <- X[c(ind[[ii]],ind[[jj]]),]
            Dij <- D[c(ind[[ii]],ind[[jj]])]
            DXij <- as_tibble(cbind(D = Dij,Xij))
            D1Xij <- as_tibble(cbind(D = rep(1,nrow(Xij)),Xij))
            D0Xij <- as_tibble(cbind(D = rep(0,nrow(Xij)),Xij))
            ruleij <- rule[c(ind[[ii]],ind[[jj]])]
            psij <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xij,Dij,ML = MLps)
            if (sum(psij <= 0 | psij >= 1)!= 0){
              psij <- (psij > 0 & psij < 1)*psij +
              0.001*(psij <= 0) + 0.999*(psij >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
            }
            fv1ij <- ML::FVest(mdx,DXnotij,Ynotij,D1Xij,Yij,ML = MLiop)
            fv0ij <- ML::FVest(mdx,DXnotij,Ynotij,D0Xij,Yij,ML = MLiop)
            fvdxij <- ML::FVest(mdx,DXnotij,Ynotij,DXij,Yij,ML = MLiop)
            fvtij <- fv1ij*ruleij + fv0ij*(1-ruleij)
            wij <- ruleij*Dij/psij + (1-ruleij)*(1-Dij)/(1-psij)
            a1 <- 1/2
            c11 <- 1
            c12 <- 1
            c21 <- 1
            c22 <- -1
            for (kk in 1:length(ind[[ii]])){
              mm1 <- length(ind[[ii]]) + 1
              for (mm in mm1:length(Yij)){
                m <- 0.5*(fvtij[kk] + fvtij[mm] - abs(fvtij[kk]-fvtij[mm]))
                a2 <- ((fvtij[kk]>fvtij[mm]) - (fvtij[kk] < fvtij[mm]))/2
                phi1 <- a1*(c11*wij[kk]*(Yij[kk] - fvdxij[kk]) + c12*wij[mm]*(Yij[mm] - fvdxij[mm]))
                phi2 <- a2*(c21*wij[kk]*(Yij[kk] - fvdxij[kk]) + c22*wij[mm]*(Yij[mm] - fvdxij[mm]))
                jt <- jt + (m + phi1 + phi2)
              }
            }
          }
        }
      }
      WT <- (2/(n*(n-1)))*jt
      mu <- wutil(Y = Y,D = D,X = X,rule = rule,
                  design = design,
                  est_method = est_method,
                  MLps = MLps,
                  MLreg = MLalpha,
                  CF = CF, K = K)
      mu <- mu$Welfare
      iop <- 1 - WT/mu
      return(list("Welfare" = WT, "Mean" = mu, "IOp" = iop))
          ###########################################
    }
    else if (CF == FALSE){
      n <- length(Y)
      DX <- as_tibble(cbind(D = D,X))
      mdx <- ML::MLest(DX,Y,ML = MLiop,FVs = TRUE)
      D1X <- as_tibble(cbind(D = rep(1,nrow(X)),X))
      D0X <- as_tibble(cbind(D = rep(0,nrow(X)),X))
      fv1 <- ML::FVest(mdx$model,DX,Y,D1X,Y,ML = MLiop)
      fv0 <- ML::FVest(mdx$model,DX,Y,D0X,Y,ML = MLiop)
      fvdx <- mdx$FVs
      fvt <- fv1*rule + fv0*(1-rule)
      ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
      ps <- ps$FVs
      if (sum(ps<=0 | ps >= 1)!=0){
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
      }
      w <- rule*D/ps + (1-rule)*(1-D)/(1-ps)
      a1 <- 1/2
      c11 <- 1
      c12 <- 1
      c21 <- 1
      c22 <- -1
      n1 <- n-1
      jt <- 0

      for(i in 1:n1){
        # print(i)
        j1 <- i + 1
        for (j in j1:n){
          m <- 0.5*(fvt[i] + fvt[j] - abs(fvt[i]-fvt[j]))
          a2 <- ((fvt[i]>fvt[j]) - (fvt[i] < fvt[j]))/2
          phi1 <- a1*(c11*w[i]*(Y[i] - fvdx[i]) + c12*w[j]*(Y[j] - fvdx[j]))
          phi2 <- a2*(c21*w[i]*(Y[i] - fvdx[i]) + c22*w[j]*(Y[j] - fvdx[j]))
          jt <- jt + (m + phi1 + phi2)
        }
      }
      WT <- (2/(n*(n-1)))*jt
      mu <- wutil(Y = Y,D = D,X = X,rule = rule,
                  design = design,
                  est_method = est_method,
                  MLps = MLps,
                  MLreg = MLalpha,
                  CF = CF, K = K)
      mu <- mu$Welfare
      iop <- 1 - WT/mu
    }
  }
  return(list("Welfare" = WT, "Mean" = mu, "IOp" = iop))
}

#' Inequality of Opportunity Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X Tibble with circumstances.
#' @param design Whether data comes from an RCT or from observational data
#' @param est_method Whether traditional plug in estimators are to be used
#' or locally robust estimators
#' @param MLiop Choice of Machine Learner for outcome regression
#' @param MLps Choice of Machine Learner for propensity score
#'
#' @return A list with the output and a figure.
#' @export
wiop_scores <- function(Y,D,X,
                 design = c("rct","observational"),
                 est_method = c("PI","LR"),
                 MLiop = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
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
  if (est_method == "PI"){
    DX <- dplyr::as_tibble(cbind(D = D,X))
    mdx <- ML::MLest(DX,Y,ML = MLiop,FVs = TRUE)
    D1X <- dplyr::as_tibble(cbind(D = rep(1,nrow(X)),X))
    D0X <- dplyr::as_tibble(cbind(D = rep(0,nrow(X)),X))
    fv1 <- ML::FVest(mdx$model,DX,Y,D1X,Y,ML = MLiop)
    fv0 <- ML::FVest(mdx$model,DX,Y,D0X,Y,ML = MLiop)
    n1 <- n-1
    cnt <- 0
    for(i in 1:n1){
      j1 <- i + 1
      for (j in j1:n){
        cnt <- cnt + 1
        scores11[cnt,] <- c(0.5*(fv1[i] + fv1[j] - abs(fv1[i] - fv1[j])),
                            i,j)
        scores10[cnt,] <- c(0.5*(fv1[i] + fv0[j] - abs(fv1[i] - fv0[j])),
                            i,j)
        scores01[cnt,] <- c(0.5*(fv0[i] + fv1[j] - abs(fv0[i] - fv1[j])),
                            i,j)
        scores00[cnt,] <- c(0.5*(fv0[i] + fv0[j] - abs(fv0[i] - fv0[j])),
                            i,j)
      }
    }
    return(list(scores11,scores10,scores01,scores00))
  }
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
          Dnotij <- D[c(-ind[[ii]],-ind[[jj]])]
          DXnotij <- dplyr::as_tibble(cbind(D = Dnotij,Xnotij))
          Ynotij <- Y[c(-ind[[ii]],-ind[[jj]])]
          ########## Train ps model with obs not in Ci or Cj ############
          mps <- ML::MLest(Xnotij,as.numeric(Dnotij),ML = MLps)
          psnotij <- mps$FVs
          if (sum(psnotij <= 0 | psnotij >= 1)!= 0){
            psnotij <- (psnotij > 0 & psnotij < 1)*psnotij +
              0.001*(psnotij <= 0) + 0.999*(psnotij >= 1)
            warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greater or equal than 1 have been set to 0.999")
          }
          mpsnotij <- mps$model
          ########## Train fvs iop model with obs not in Ci or Cj ########
          mdx <- ML::MLest(DXnotij,Ynotij,ML = MLiop,FVs = TRUE)
          mdx <- mdx$model
          ############## Compute scores evaluating in observations in Ci and Cj ##############
          #If we are in a triangle (Ci = Cj)
          if (ii == jj){
            Yii <- Y[ind[[ii]]]
            Xii <- X[ind[[ii]],]
            Dii <- D[ind[[ii]]]
            DXii <- dplyr::as_tibble(cbind(D = Dii,Xii))
            D1Xii <- dplyr::as_tibble(cbind(D = rep(1,nrow(Xii)),Xii))
            D0Xii <- dplyr::as_tibble(cbind(D = rep(0,nrow(Xii)),Xii))
            idii <- id[ind[[ii]]]
            psii <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xii,Dii,ML = MLps)
            if (sum(psii <= 0 | psii >= 1)!= 0){
              psii <- (psii > 0 & psii < 1)*psii +
                0.001*(psii <= 0) + 0.999*(psii >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greater or equal than 1 have been set to 0.999")
            }
            fv1ii <- ML::FVest(mdx,DXnotij,Ynotij,D1Xii,Yii,ML = MLiop)
            fv0ii <- ML::FVest(mdx,DXnotij,Ynotij,D0Xii,Yii,ML = MLiop)
            fvdxii <- ML::FVest(mdx,DXnotij,Ynotij,DXii,Yii,ML = MLiop)
            nn1 <- length(Yii) - 1
            for (kk in 1:nn1){
              kk1 <- kk + 1
              for (mm in kk1:length(Yii)){
                cnt <- cnt + 1
                sgn11 <- (fv1ii[kk] > fv1ii[mm]) - (fv1ii[kk] < fv1ii[mm])
                sgn10 <- (fv1ii[kk] > fv0ii[mm]) - (fv1ii[kk] < fv0ii[mm])
                sgn01 <- (fv0ii[kk] > fv1ii[mm]) - (fv0ii[kk] < fv1ii[mm])
                sgn00 <- (fv0ii[kk] > fv0ii[mm]) - (fv0ii[kk] < fv0ii[mm])

                scores11[cnt,] <- c(0.5*(fv1ii[kk] + fv1ii[mm] - abs(fv1ii[kk] - fv1ii[mm]) +
                                         Dii[kk]*(1-sgn11)*(Yii[kk] - fvdxii[kk])/psii[kk] +
                                         Dii[mm]*(1-sgn11)*(Yii[mm] - fvdxii[mm])/psii[mm]),
                                    idii[kk],idii[mm])
                scores10[cnt,] <- c(0.5*(fv1ii[kk] + fv0ii[mm] - abs(fv1ii[kk] - fv0ii[mm]) +
                                           Dii[kk]*(1-sgn10)*(Yii[kk] - fvdxii[kk])/psii[kk] +
                                           (1-Dii[mm])*(1-sgn10)*(Yii[mm] - fvdxii[mm])/(1-psii[mm])),
                                    idii[kk],idii[mm])
                scores01[cnt,] <- c(0.5*(fv0ii[kk] + fv1ii[mm] - abs(fv0ii[kk] - fv1ii[mm]) +
                                           (1-Dii[kk])*(1-sgn01)*(Yii[kk] - fvdxii[kk])/(1-psii[kk]) +
                                           Dii[mm]*(1-sgn01)*(Yii[mm] - fvdxii[mm])/psii[mm]),
                                    idii[kk],idii[mm])
                scores00[cnt,] <- c(0.5*(fv0ii[kk] + fv0ii[mm] - abs(fv0ii[kk] - fv0ii[mm]) +
                                           (1-Dii[kk])*(1-sgn00)*(Yii[kk] - fvdxii[kk])/(1-psii[kk]) +
                                           (1-Dii[mm])*(1-sgn00)*(Yii[mm] - fvdxii[mm])/(1-psii[mm])),
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
            DXij <- dplyr::as_tibble(cbind(D = Dij,Xij))
            D1Xij <- dplyr::as_tibble(cbind(D = rep(1,nrow(Xij)),Xij))
            D0Xij <- dplyr::as_tibble(cbind(D = rep(0,nrow(Xij)),Xij))
            idij <- id[c(ind[[ii]],ind[[jj]])]
            psij <- ML::FVest(mpsnotij,Xnotij,Dnotij,Xij,Dij,ML = MLps)
            if (sum(psij <= 0 | psij >= 1)!= 0){
              psij <- (psij > 0 & psij < 1)*psij +
                0.001*(psij <= 0) + 0.999*(psij >= 1)
              warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greater or equal than 1 have been set to 0.999")
            }
            fv1ij <- ML::FVest(mdx,DXnotij,Ynotij,D1Xij,Yij,ML = MLiop)
            fv0ij <- ML::FVest(mdx,DXnotij,Ynotij,D0Xij,Yij,ML = MLiop)
            fvdxij <- ML::FVest(mdx,DXnotij,Ynotij,DXij,Yij,ML = MLiop)
            for (kk in 1:length(ind[[ii]])){
              mm1 <- length(ind[[ii]]) + 1
              for (mm in mm1:length(Yij)){
                cnt <- cnt + 1
                sgn11 <- (fv1ij[kk] > fv1ij[mm]) - (fv1ij[kk] < fv1ij[mm])
                sgn10 <- (fv1ij[kk] > fv0ij[mm]) - (fv1ij[kk] < fv0ij[mm])
                sgn01 <- (fv0ij[kk] > fv1ij[mm]) - (fv0ij[kk] < fv1ij[mm])
                sgn00 <- (fv0ij[kk] > fv0ij[mm]) - (fv0ij[kk] < fv0ij[mm])

                scores11[cnt,] <- c(0.5*(fv1ij[kk] + fv1ij[mm] - abs(fv1ij[kk] - fv1ij[mm]) +
                                           Dij[kk]*(1-sgn11)*(Yij[kk] - fvdxij[kk])/psij[kk] +
                                           Dij[mm]*(1-sgn11)*(Yij[mm] - fvdxij[mm])/psij[mm]),
                                    idij[kk],idij[mm])
                scores10[cnt,] <- c(0.5*(fv1ij[kk] + fv0ij[mm] - abs(fv1ij[kk] - fv0ij[mm]) +
                                           Dij[kk]*(1-sgn10)*(Yij[kk] - fvdxij[kk])/psij[kk] +
                                           (1-Dij[mm])*(1-sgn10)*(Yij[mm] - fvdxij[mm])/(1-psij[mm])),
                                    idij[kk],idij[mm])
                scores01[cnt,] <- c(0.5*(fv0ij[kk] + fv1ij[mm] - abs(fv0ij[kk] - fv1ij[mm]) +
                                           (1-Dij[kk])*(1-sgn01)*(Yij[kk] - fvdxij[kk])/(1-psij[kk]) +
                                           Dij[mm]*(1-sgn01)*(Yij[mm] - fvdxij[mm])/psij[mm]),
                                    idij[kk],idij[mm])
                scores00[cnt,] <- c(0.5*(fv0ij[kk] + fv0ij[mm] - abs(fv0ij[kk] - fv0ij[mm]) +
                                           (1-Dij[kk])*(1-sgn00)*(Yij[kk] - fvdxij[kk])/(1-psij[kk]) +
                                           (1-Dij[mm])*(1-sgn00)*(Yij[mm] - fvdxij[mm])/(1-psij[mm])),
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
      n <- length(Y)
      DX <- dplyr::as_tibble(cbind(D = D,X))
      mdx <- ML::MLest(DX,Y,ML = MLiop,FVs = TRUE)
      D1X <- dplyr::as_tibble(cbind(D = rep(1,nrow(X)),X))
      D0X <- dplyr::as_tibble(cbind(D = rep(0,nrow(X)),X))
      fv1 <- ML::FVest(mdx$model,DX,Y,D1X,Y,ML = MLiop)
      fv0 <- ML::FVest(mdx$model,DX,Y,D0X,Y,ML = MLiop)
      fvdx <- mdx$FVs
      ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
      ps <- ps$FVs
      if (sum(ps<=0 | ps >= 1)!=0){
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greater or equal than 1 have been set to 0.999")
      }
      n1 <- n-1
      for(i in 1:n1){
        j1 <- i + 1
        for (j in j1:n){
          sgn11 <- (fv1[i] > fv1[j]) - (fv1[i] < fv1[j])
          sgn10 <- (fv1[i] > fv0[j]) - (fv1[i] < fv0[j])
          sgn01 <- (fv0[i] > fv1[j]) - (fv0[i] < fv1[j])
          sgn00 <- (fv0[i] > fv0[j]) - (fv0[i] < fv0[j])

          scores11[cnt,] <- c(0.5*(fv1[i] + fv1[j] - abs(fv1[i] - fv1[j]) +
                                     D[kk]*(1-sgn11)*(Y[i] - fvdx[i])/ps[i] +
                                     D[mm]*(1-sgn11)*(Y[j] - fvdx[j])/ps[j]),
                              i,j)
          scores10[cnt,] <- c(0.5*(fv1[i] + fv0[j] - abs(fv1[i] - fv0[j]) +
                                     D[kk]*(1-sgn10)*(Y[i] - fvdx[i])/ps[i] +
                                     (1-D[mm])*(1-sgn10)*(Y[j] - fvdx[j])/(1-ps[j])),
                              i,j)
          scores01[cnt,] <- c(0.5*(fv0[i] + fv1[j] - abs(fv0[i] - fv1[j]) +
                                     (1-D[kk])*(1-sgn01)*(Y[i] - fvdx[i])/(1-ps[i]) +
                                     D[mm]*(1-sgn01)*(Y[j] - fvdx[j])/ps[j]),
                              i,j)
          scores00[cnt,] <- c(0.5*(fv0[i] + fv0[j] - abs(fv0[i] - fv0[j]) +
                                     (1-D[kk])*(1-sgn00)*(Y[i] - fvdx[i])/(1-ps[i]) +
                                     (1-D[mm])*(1-sgn00)*(Y[j] - fvdx[j])/(1-ps[j])),
                              i,j)
        }
      }
      return(list(scores11,scores10,scores01,scores00))
    }
  }
}

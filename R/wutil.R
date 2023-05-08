#' Utilitarian Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X tibble with covariates
#' @param rule Treatment rule.
#' @param design Whether data comes from an RCT or from observational data
#' @param est_method Whether traditional plug in estimators are to be used
#' or locally robust estimators
#' @param MLps Choice of machine learner for propensity score
#' @param MLreg Choice of machine learner for outcome regression
#' @param CF Whether to use Cross-fitting
#' @param K number of folds in cross-fitting
#'
#' @return A list with the output and a figure.
#' @export
wutil <- function(Y,D,X,rule,
                  design = c("rct","observational"),
                  est_method = c("PI","LR"),
                  MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  MLreg = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  CF = TRUE,
                  K = 5){
D <- as.numeric(D)
Y <- as.numeric(Y)
n <- length(Y)
############### RCT ####################
  if (design == "rct"){
    p <- mean(D)
    WT <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  }
################ Observational ###########
  else if (design == "observational"){
    ############### Plug in ###############
    if (est_method == "PI"){
      ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
      ps <- ps$FVs
      if (sum(ps <= 0 | ps >= 1)!= 0){
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
      }
      WT <- mean(Y*(1-D)/(1-ps) + ((Y*D)/ps - Y*(1-D)/(1-ps))*rule)
    }
    ################ Locally robust ###########
    else if (est_method == "LR"){
      if (CF == TRUE){
        fv1 <- rep(0,n)
        fv0 <- rep(0,n)
        ps <- rep(0,n)
        ind <- split(seq(n), seq(n) %% K)
        for (i in 1:K){
          mps <- ML::modest(X[-ind[[i]],],D[-ind[[i]]],ML = MLps)
          ps[ind[[i]]] <- ML::FVest(mps,X[-ind[[i]],],D[-ind[[i]]],
                                     X[ind[[i]],],D[ind[[i]]],ML = MLps)
          XD <- dplyr::filter(as_tibble(data.frame(Y = Y[-ind[[i]]],
                                                   D = D[-ind[[i]]],
                                                   X[-ind[[i]],])))
          X1 <- XD %>% filter(D==1) %>% dplyr::select(-c("D","Y"))
          X0 <- XD %>% filter(D==0) %>% dplyr::select(-c("D","Y"))
          Y1 <- XD %>% filter(D==1) %>% dplyr::select(Y)
          Y0 <- XD %>% filter(D==0) %>% dplyr::select(Y)

          mfv1 <- ML::modest(X1,Y1$Y,ML = MLreg)
          mfv0 <- ML::modest(X0,Y0$Y,ML = MLreg)

          fv1[ind[[i]]] <- ML::FVest(mfv1, X1,Y1$Y,
                                      X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
          fv0[ind[[i]]] <- ML::FVest(mfv0, X0,Y0$Y,
                                      X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
        }
        if (sum(ps <= 0 | ps >= 1)!= 0){
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
          warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
        }
        WT <- Y*(1-D)/(1-ps) +
          (fv1 - fv0 + (D/ps)*(Y-fv1) - ((1-D)/(1-ps))*(Y-fv0))*rule
        WT <- mean(WT)
        }
      else if (CF == FALSE){
          ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
          ps <- ps$FVs
          if (sum(ps <= 0 | ps >= 1)!= 0){
            ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
            warning("There are estimated propensity scores outside (0,1).
                  Values lower or equal than zero have been set to 0.001
                  and values greter or equal than 1 have been set to 0.999")
          }
          mfv1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
          fv1 <- ML::FVest(mfv1,X[D==1,],Y[D==1],X,Y,ML = MLreg)
          mfv0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
          fv0 <- ML::FVest(mfv0,X[D==0,],Y[D==0],X,Y,ML = MLreg)
          WT <- Y*(1-D)/(1-ps) +
            (fv1 - fv0 + (D/ps)*(Y-fv1) - ((1-D)/(1-ps))*(Y-fv0))*rule
          WT <- mean(WT)
        }
      }
  }
###################
return(list("Welfare" = WT))
}

#' Standard Gini Social Welfare Function maximization
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
#'
#' @return A list with the output and a figure.
#' @export
wineq <- function(Y,D,X,rule,
                  design = c("rct","observational"),
                  est_method = c("PI","LR"),
                  MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                  CF = FALSE){
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
        w <- w + 0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a11+a10+a01+a00)
      }
    }
    WT <- (2/(n*(n-1)))*w
    mu <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
    G <- 1 - WT/mu
    return(list("Welfare" = WT, "Mean" = mu, "Gini" = G))
  }
  if (design == "observational"){
    if (est_method == "PI"){
      ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
      ps <- ps$FVs
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
          w <- w + 0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a11+a10+a01+a00)
        }
      }
      WT <- (2/(n*(n-1)))*w
      mu <- mean(((Y*D)/ps - Y*(1-D)/(1-ps))*rule + Y*(1-D)/(1-ps))
      G <- 1 - WT/mu
      return(list("Welfare" = WT, "Mean" = mu, "Gini" = G))
    }
    else if (est_method == "LR"){
      if (CF == TRUE){
        stop("LR with CF not coded yet")
      }
      else if (CF == FALSE){
        ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        n <- length(Y)
        delta <- rep(0,n)
        #Estimation of alpha
        for (i in 1:n){
          a11 <- ((D[i]*D[-i])/((ps[i]^2)*ps[-i]))*(rule[i]*rule[-i])
          a10 <- ((D[i]*(1-D[-i]))/((ps[i]^2)*(1-ps[-i])))*(rule[i]*(1-rule[-i]))
          a01 <- (((1-D[i])*D[-i])/(((1-ps[i])^2)*ps[-i]))*((1-rule[i])*rule[-i])
          a00 <- (((1-D[i])*(1-D[-i]))/(((1-ps[i])^2)*(1-ps[-i])))*((1-rule[i])*(1-rule[-i]))
          delta[i] <- mean(0.5*(Y[i] + Y[-i] - abs(Y[i] - Y[-i]))*(a00+a01-a10-a11))
        }
        alpha <- ML::MLest(X,delta,ML = MLalpha,FVs = TRUE)
        alpha <- alpha$FVs
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
            w <- w + 0.5*(Y[i] + Y[j] - abs(Y[i] - Y[j]))*(a11+a10+a01+a00) +
              alpha[i]*(D[i]-ps[i]) + alpha[j]*(D[j]-ps[j])
          }
        }
        WT <- (2/(n*(n-1)))*w
        mu <- mean(((Y*D)/ps - Y*(1-D)/(1-ps))*rule + Y*(1-D)/(1-ps))
        G <- 1 - WT/mu
        return(list("Welfare" = WT, "Mean" = mu, "Gini" = G))
      }
    }
  }
}

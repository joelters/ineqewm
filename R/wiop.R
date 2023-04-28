#' Inequality of Opportunity Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X Tibble with circumstances.
#' @param rule Treatment rule.
#' @param ML Choice of Machine Learner
#'
#' @return A list with the output and a figure.
#' @export
wiop <- function(Y,D,X,rule,
         ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL")){
  p <- mean(D)
  n <- length(Y)
  m1 <- ML::MLest(X,Y*D,ML = ML, FVs = TRUE)
  m0 <- ML::MLest(X,Y*(1-D),ML = ML, FVs = TRUE)
  fvs1 <- m1$FVs
  fvs0 <- m0$FVs
  fvst <- rule*fvs1/p + (1-rule)*fvs0/(1-p)
  n1 <- n-1
  w <- 0
  for(i in 1:n1){
    # print(i)
    j1 <- i + 1
    for (j in j1:n){
      w <- w + 0.5*(fvst[i] + fvst[j] - abs(fvst[i]-fvst[j]))
    }
  }
  WT <- (2/(n*(n-1)))*w
  mu <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  iop <- DescTools::Gini(fvst)
  return(list("Welfare" = WT, "Mean" = mu, "IOp" = iop))
}

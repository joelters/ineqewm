wutil <- function(Y,D,rule){
  p <- mean(D)
  n <- length(Y)
  WT <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  G <- Gini(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  return(list("Welfare" = WT, "Gini" = G))
}

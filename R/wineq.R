wineq <- function(Y,D,rule){
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
  G <- Gini(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  return(list("Welfare" = WT, "Mean" = mu, "Gini" = G))
}

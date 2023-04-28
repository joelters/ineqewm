#' Kendall tau Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param X1 variable to compute Kendall-tau with
#' @param D Treatment assignment.
#' @param t target Kendall-tau
#' @param rule Treatment rule.
#'
#' @return A list with the output and a figure.
#' @export
wigm <- function(Y,X1,D,t,rule){
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
}

potential_scores_igm <- function(Y1,Y0,X){
  n <- length(Y1)

  scores11 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores10 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores01 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))
  scores00 <- cbind(rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)),rep(0,0.5*n*(n-1)))

  cnt <- 0
  n1 <- n - 1
  for (i in 1:n1){
    j1 <- i + 1
    for (j in j1:n){
      cnt <- cnt + 1
      scores11[cnt,] <- c(sign(Y1[i] - Y1[j])*sign(X[i] - X[j]),i,j)
      scores10[cnt,] <- c(sign(Y1[i] - Y0[j])*sign(X[i] - X[j]),i,j)
      scores01[cnt,] <- c(sign(Y0[i] - Y1[j])*sign(X[i] - X[j]),i,j)
      scores00[cnt,] <- c(sign(Y0[i] - Y0[j])*sign(X[i] - X[j]),i,j)
    }
  }
  return(list(scores11,scores10,scores01,scores00))
}

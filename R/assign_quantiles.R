assign_quantiles <- function(X,quants = 4){
  if (quants == 0){
    X
  }
  else {
    a <- 0
    for (i in 1:quants){
      a <- a + quantile(X,i/quants)*(X >= quantile(X,i/quants - 1/quants) &
                                   X <= quantile(X,i/quants))
    }
    return(a)
  }
}

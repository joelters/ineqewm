#' Utilitarian Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param rule Treatment rule.
#'
#' @return A list with the output and a figure.
#' @export
wutil <- function(Y,D,rule){
  p <- mean(D)
  n <- length(Y)
  WT <- mean(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  G <- DescTools::Gini(((Y*D)/p - Y*(1-D)/(1-p))*rule + Y*(1-D)/(1-p))
  return(list("Welfare" = WT))
}

#' True ewm simulation
#'
#' @param Y1 Potential outcome under treatment
#' @param Y0 Potential outcome under no treatment
#' @param fv1 Potential fitted values under treatment
#' @param fv0 Potential fitted values under no treatment
#' @param rule Treatment rule.
#' @param X1 outcome with which to compute Kendall-tau with
#' @param t Target for Kendall-tau
#' @param WF Welfare function
#'
#' @return A list with the output and a figure.
#' @export
wsims <- function(Y1,Y0,fv1,fv0,rule,
                  X1,t,
                  WF = c("ineq","IOp","IGM","util")){
  YT <- Y1*rule + Y0*(1-rule)
  fvT <- fv1*rule + fv0*(1-rule)
  if (WF == "util"){
    WT <- mean(YT)
  }
  else if (WF == "ineq"){
    WT <- mean(YT)*(1-Gini(YT))
  }
  else if (WF == "IOp"){
    WT <- mean(fvT)*(1-Gini(fvT))
  }
  else if (WF == "IGM"){
    WT <- abs(cor.test(YT,X1, method="kendall")$estimate-t)
  }
}

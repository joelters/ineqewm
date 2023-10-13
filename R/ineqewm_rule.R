#' Social Welfare Function maximization
#'
#' @param scores11 Scores when (i,j) both treated
#' @param scores10 Scores when i treated and j not treated
#' @param scores01 Scores when i not treated and j treated
#' @param scores11 Scores when (i,j) both not treated
#' @param rule rule for which we want to compute welfare.
#' @param welfare Welfare function
#' @param t Target for Kendall-tau
#' @return Estimated welfare and Kendall-tau for IGM case.
#' @export
ineqewm_rule <-function(scores11 = NULL,
                    scores10 = NULL,
                    scores01 = NULL,
                    scores00 = NULL,
                    rule = NULL,
                    welfare = c("ineq","iop","igm","util"),
                    t = 0){
  if (welfare == "util"){
    s11_1 <- scores11[,1]
    s00_1 <- scores00[,1]

    s11_2 <- scores11[,2]
    s00_2 <- scores00[,2]
    W <- (1/n)*collapse::fsum(s11_1*rule[s11_2] + s00_1*(1-rule[s00_2]))
    return("Welfare" = W)
  }
  else if (welfare == "igm"){
    s11_1 <- scores11[,1]
    s10_1 <- scores10[,1]
    s01_1 <- scores01[,1]
    s00_1 <- scores00[,1]

    s11_2 <- scores11[,2]
    s10_2 <- scores10[,2]
    s01_2 <- scores01[,2]
    s00_2 <- scores00[,2]

    s11_3 <- scores11[,3]
    s10_3 <- scores10[,3]
    s01_3 <- scores01[,3]
    s00_3 <- scores00[,3]

    tau <- (2/(n*(n-1)))*collapse::fsum(s11_1*rule[s11_2]*rule[s11_3] +
                                         s10_1*rule[s10_2]*(1-rule[s10_3]) +
                                         s01_1*(1-rule[s01_2])*rule[s01_3] +
                                         s00_1*(1-rule[s00_2])*(1-rule[s00_3]))
    W <- -abs(tau - t)
    return("Welfare" = W, "tau" = tau)
  }
  else {
    s11_1 <- scores11[,1]
    s10_1 <- scores10[,1]
    s01_1 <- scores01[,1]
    s00_1 <- scores00[,1]

    s11_2 <- scores11[,2]
    s10_2 <- scores10[,2]
    s01_2 <- scores01[,2]
    s00_2 <- scores00[,2]

    s11_3 <- scores11[,3]
    s10_3 <- scores10[,3]
    s01_3 <- scores01[,3]
    s00_3 <- scores00[,3]

    W <- (2/(n*(n-1)))*collapse::fsum(s11_1*rule[s11_2]*rule[s11_3] +
                                         s10_1*rule[s10_2]*(1-rule[s10_3]) +
                                         s01_1*(1-rule[s01_2])*rule[s01_3] +
                                         s00_1*(1-rule[s00_2])*(1-rule[s00_3]))
    return("Welfare" = W)
  }


























}

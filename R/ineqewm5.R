#' Social Welfare Function maximization
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X Tibble with circumstances.
#' @param targetX Variables to use for treatment allocation
#' @param rule Family of treatment rules to search from
#' @param WF Welfare function
#' @param tigm Target for Kendall-tau
#' @param X1 variable to compute Kendall-tau with
#' @param quants how many quantiles to split continuous vars in targetX
#' @param ML Choice of Machine Learner (for IOp)
#'
#' @return A list with the output and a figure.
#' @export
ineqewm5 <-function(scores,
                    targetX,
                    rule = c("tree","linscore"),
                    approx = FALSE,
                    approxmatrix = NULL){
  ############# Welfare when no one is treated ###################
  n <- nrow(targetX)
  rl0 <- rep(0,n)
  WT0 <- (2/(n*(n-1)))*sum(scores[[1]][,1]*rl0[scores[[1]][,2]]*rl0[scores[[1]][,3]] +
                             scores[[2]][,1]*rl0[scores[[1]][,2]]*(1-rl0[scores[[1]][,3]]) +
                             scores[[3]][,1]*(1-rl0[scores[[1]][,2]])*rl0[scores[[1]][,3]] +
                             scores[[4]][,1]*(1-rl0[scores[[1]][,2]])*(1-rl0[scores[[1]][,3]]))

  ##################### TReeees ###########
  s <- NULL
  if (approx == TRUE){
    ns <- approxmatrix
    for (i in 1:ncol(targetX)){
      # s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
      #                   seq(1/(ns[i]+1), ns[i]/(ns[i] + 1), 1/(ns[i]+1)),
      #                   names = FALSE)),1)
      s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
                                      seq(0, 1, 1/(ns[i]+1)),
                                      names = FALSE)),1)
    }
  }
  else{
    for (i in 1:ncol(targetX)){
      s[[i]] <- sort(unique(targetX[,i,drop = TRUE]))
    }
  }

  res <- list(Welfare = -9999)

  for (ii in 1:ncol(targetX)){
    print(ii)
    for (jj in 1:length(s[[ii]])){
      leaves0 <- targetX[,ii,drop = TRUE] <= s[[ii]][jj]
      for (kk in 1:ncol(targetX)){
        for (ll in 1:length(s[[kk]])){
          leaves1 <- (targetX[,ii,drop = TRUE] <= s[[ii]][jj])*
            (targetX[,kk,drop = TRUE] <= s[[kk]][ll]) +
            2*(targetX[,ii,drop = TRUE] <= s[[ii]][jj])*
            (targetX[,kk,drop = TRUE] > s[[kk]][ll])
          for (pp in 1:ncol(targetX)){
            for (tt in 1:length(s[[pp]])){
              leaves2 <- 3*(targetX[,ii,drop = TRUE] > s[[ii]][jj])*
                (targetX[,pp,drop = TRUE] <= s[[pp]][tt]) +
                4*(targetX[,ii,drop = TRUE] > s[[ii]][jj])*
                (targetX[,pp,drop = TRUE] > s[[pp]][tt])
              rules <- data.frame(rule_12_34 = (leaves1 == 1 | leaves2 == 3),
                                  rule_12_43 = (leaves1 == 1 | leaves2 == 4),
                                  rule_21_34 = (leaves1 == 2 | leaves2 == 3),
                                  rule_21_43 = (leaves1 == 2 | leaves2 == 4))


              restr <- apply(rules, 2, function(rl){
              (2/(n*(n-1)))*sum(scores[[1]][,1]*rl[scores[[1]][,2]]*rl[scores[[1]][,3]] +
              scores[[2]][,1]*rl[scores[[1]][,2]]*(1-rl[scores[[1]][,3]]) +
              scores[[3]][,1]*(1-rl[scores[[1]][,2]])*rl[scores[[1]][,3]] +
              scores[[4]][,1]*(1-rl[scores[[1]][,2]])*(1-rl[scores[[1]][,3]]))
              })

              res2 <- list(Welfare = max(restr),
                           node0 = colnames(targetX[,ii]),
                           thres0 = s[[ii]][jj],
                           node1 = colnames(targetX[,kk]),
                           thres1 = s[[kk]][ll],
                           node2 = colnames(targetX[,pp]),
                           thres2 = s[[pp]][tt],
                           rule = names(restr)[which.max(restr)])
              if (res$Welfare <= res2$Welfare){
                res <- res2
              }
            }
          }
        }
      }
    }
  }

  tn1 <- data.table::fifelse(substr(res$rule,6,7) == "12","Treat","Do not treat")
  tn2 <- data.table::fifelse(substr(res$rule,6,7) == "21","Do not Treat","Treat")
  tn3 <- data.table::fifelse(substr(res$rule,9,10) == "34","Treat","Do not treat")
  tn4 <- data.table::fifelse(substr(res$rule,9,10) == "43","Do not Treat","Treat")

  tree <- data.tree::Node$new(paste(res$node0, "<=", res$thres0))
  node1 <- tree$AddChild(paste(res$node1, "<=", res$thres1))
  trnode1_1 <- node1$AddChild(tn1)
  trnode1_2 <- node1$AddChild(tn2)
  node2 <- tree$AddChild(paste(res$node2, "<=", res$thres2))
  trnode2_3 <- node2$AddChild(tn3)
  trnode2_4 <- node2$AddChild(tn4)

  # plot(tree)

  WTg <- (res$Welfare - WT0)/WT0
  return(list(c(Welfare = res$Welfare,
                Welfare0 = WT0,
                Welfare_gain = WTg),tree = tree))
}

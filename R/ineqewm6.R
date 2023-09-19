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
ineqewm6 <-function(scores11,scores10,scores01,scores00,
                    welfare = c("ineq","iop","igm","util"),
                    t,
                    targetX,
                    rule = c("tree","linear_rule"),
                    approx = FALSE,
                    approxmatrix = NULL){
  n <- nrow(targetX)

  ##################### Trees ###########
  if (rule == "tree"){
    s <- NULL
    ns <- rep(0,ncol(targetX))
    for (i in 1:ncol(targetX)){
      s[[i]] <- sort(unique(targetX[,i,drop = TRUE]))
      ns[i] <- length(s[[i]])
    }
    # if (approx == TRUE){
    #   ns <- approxmatrix
    #   for (i in 1:ncol(targetX)){
    #     # s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
    #     #                   seq(1/(ns[i]+1), ns[i]/(ns[i] + 1), 1/(ns[i]+1)),
    #     #                   names = FALSE)),1)
    #     s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
    #                                     seq(1/(ns[i]+1), ns[i]/(ns[i]+1), 1/(ns[i]+1)),
    #                                     names = FALSE)),1)
    #   }
    # }
    # else{
    #   for (i in 1:ncol(targetX)){
    #     s[[i]] <- sort(unique(targetX[,i,drop = TRUE]))
    #     ns[i] <- length(s[[i]])
    #   }
    # }

    s <- lapply(s,function(u){
      if(length(u) < max(sapply(s,function(u) length(u)))){
        c(u,rep(0,max(sapply(s,function(u) length(u)))- length(u)))
      } else{u}
    })
    s <- do.call(cbind,s)

    res <- Rlrpltree(scores11, scores10, scores01, scores00,
        welfare, t, targetX, depth = 2, s, ns)

    rls <- c("rl_ttnt","rl_tttn",
             "rl_tntn","rl_tnnt","rl_tntt","rl_tnnn",
             "rl_ntnt","rl_nttn","rl_nttt","rl_ntnn",
             "rl_nnnt","rl_nntn")

    res2 <- list(Welfare = res[1,7],
             node0 = colnames(targetX[,res[1,1]]),
             thres0 = s[res[1,2],res[1,1]],
             node1 = colnames(targetX[,res[1,3]]),
             thres1 = s[res[1,4],res[1,3]],
             node2 = colnames(targetX[,res[1,5]]),
             thres2 = s[res[1,6],res[1,5]],
             rule = rls[res[1,8]])

    tn1 <- data.table::fifelse(substr(res2$rule,4,4) == "t","Treat","Do not treat")
    tn2 <- data.table::fifelse(substr(res2$rule,5,5) == "t","Treat","Do not treat")
    tn3 <- data.table::fifelse(substr(res2$rule,6,6) == "t","Treat","Do not treat")
    tn4 <- data.table::fifelse(substr(res2$rule,7,7) == "t","Treat","Do not treat")

    tree <- data.tree::Node$new(paste(res2$node0, "<=", res2$thres0))
    node1 <- tree$AddChild(paste(res2$node1, "<=", res2$thres1))
    trnode1_1 <- node1$AddChild(tn1)
    trnode1_2 <- node1$AddChild(tn2)
    node2 <- tree$AddChild(paste(res2$node2, "<=", res2$thres2))
    trnode2_3 <- node2$AddChild(tn3)
    trnode2_4 <- node2$AddChild(tn4)

    # plot(tree)
    WT0 <- res[1,9]
    WTg <- (res2$Welfare - WT0)/WT0
    return(list(Welfare = c(Welfare = res2$Welfare,
                  Welfare0 = WT0,
                  Welfare_gain = WTg),tree = tree,
                tnode = res[,10]))
  }
  if (rule == "linear_rule"){
    X1 <- targetX[,1]
    X2 <- targetX[,2]
    if (ncol(targetX) > 2){
      warning("Linear rule only uses the first two columns of targetX")
    }
    res <- linear_rules(X1 = X1,
                        X2 = X2,
                        scores11 = scores11,
                        scores10 = scores10,
                        scores01 = scores01,
                        scores00 = scores00,
                        welfare = welfare,
                        t = t)
    W0 <- res$W0
    Wg <- (res$W - W0)/W0
    return(list(Welfare = c(Welfare = res$W,
                            Welfare0 = W0,
                            Welfare_gain = Wg),
                            coefficients = res$L,
                            plot = res$plot))
  }
}

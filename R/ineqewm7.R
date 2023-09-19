#' Social Welfare Function maximization
#'
#' @param scores11 Scores when (i,j) both treated
#' @param scores10 Scores when i treated and j not treated
#' @param scores01 Scores when i not treated and j treated
#' @param scores11 Scores when (i,j) both not treated
#' @param welfare Welfare function
#' @param t Target for Kendall-tau
#' @param targetX Variables to use for treatment allocation
#' @param rule Family of treatment rules to search from
#' @param quants vector indicating how many quantiles to split each var
#' in targetX. If 0 the var is not split, i.e. quants = c(0,4) would
#' leave the first var (column) in targetX unchanged and split the second
#' one into quartiles.
#' @param verbose level of progress messaging. For trees: 0 none, 1 information
#' about the root node, 2 information about root node and first subnode,
#' 3 information about root node and both subnodes. For linear rules: 0 none,
#' 1 information about first coordinate, 2 information about both coordinates.
#' @return A list with the output and a figure.
#' @export
ineqewm7 <-function(scores11 = NULL,
                    scores10 = NULL,
                    scores01 = NULL,
                    scores00 = NULL,
                    welfare = c("ineq","iop","igm","util"),
                    t = 0,
                    targetX,
                    rule = c("tree","linear_rule"),
                    depth = 2,
                    quants = NULL,
                    verbose = 1){
  n <- nrow(targetX)

  ##################### Trees ###########
  if (rule == "tree"){
    if (!is.null(quants)){
      targetXq <- sapply(1:ncol(targetX), function(u){
        assign_quantiles(unlist(targetX[,u]),quants = quants[u])
      })
      s <- NULL
      ns <- rep(0,ncol(targetXq))
      for (i in 1:ncol(targetXq)){
        s[[i]] <- sort(unique(targetXq[,i,drop = TRUE]))
        #we take out quantile 1 (quantile 0 is already not computed by assign_quantiles)
        ns[i] <- length(s[[i]]) - 1
      }
      s <- lapply(s,function(u) u[1:(length(u)-1)])
      s <- lapply(s,function(u){
        if(length(u) < max(sapply(s,function(u) length(u)))){
          c(u,rep(0,max(sapply(s,function(u) length(u)))- length(u)))
        } else{u}
      })
      s <- do.call(cbind,s)
    }
    else{
      s <- NULL
      ns <- rep(0,ncol(targetX))
      for (i in 1:ncol(targetX)){
        s[[i]] <- sort(unique(targetX[,i,drop = TRUE]))
        # we take out lowest and higest value
        ns[i] <- length(s[[i]]) - 2
      }
      # we take out lowest and higest value
      s <- lapply(s,function(u) u[2:(length(u)-1)])
      s <- lapply(s,function(u){
        if(length(u) < max(sapply(s,function(u) length(u)))){
          c(u,rep(0,max(sapply(s,function(u) length(u)))- length(u)))
        } else{u}
      })
      s <- do.call(cbind,s)
    }

    res <- Rlrpltree(scores11 = scores11,
                     scores10 = scores10,
                     scores01 = scores01,
                     scores00 = scores00,
                     welfare = welfare,
                     t = t,
                     targetX = targetX,
                     depth = depth,
                     s = s,
                     ns = ns,
                     verbose = verbose)
    # Optimal not to treat anyone
    if (res[1,8] == 0){
      return(list(Welfare = c(Welfare = res[1,7],
                              Welfare0 = res[1,9],
                              Welfare_gain = 0),tree = "Treat no one"))
    }
    #Optimal to treat everyone
    else if (res[1,8] == 99){
      return(list(Welfare = c(Welfare = res[1,7],
                              Welfare0 = res[1,9],
                              Welfare_gain = (res[1,7] - res[1,9])/res[1,9]),tree = "Treat everyone"))
    }
    # Optimal tree has depth 1
    # treat if lower or equal to
    else if (res[1,8] == 9912){
      res2 <- list(Welfare = res[1,7],
                   node0 = colnames(targetX[,res[1,1]]),
                   thres0 = s[res[1,2],res[1,1]],
                   rule = "rl_tn")
      tn1 <- data.table::fifelse(substr(res2$rule,4,4) == "t","Treat","Do not treat")
      tn2 <- data.table::fifelse(substr(res2$rule,5,5) == "t","Treat","Do not treat")

      tree <- data.tree::Node$new(paste(res2$node0, "<=", res2$thres0))
      node1 <- tree$AddChild(tn1)
      node2 <- tree$AddChild(tn2)

      WT0 <- res[1,9]
      WTg <- (res2$Welfare - WT0)/WT0
      return(list(Welfare = c(Welfare = res2$Welfare,
                              Welfare0 = WT0,
                              Welfare_gain = WTg),tree = tree,
                  tnode = res[,10]))
    }
    #treat if strictly greater than
    else if (res[1,8] == 9921){
      res2 <- list(Welfare = res[1,7],
                   node0 = colnames(targetX[,res[1,1]]),
                   thres0 = s[res[1,2],res[1,1]],
                   rule = "rl_nt")
      tn1 <- data.table::fifelse(substr(res2$rule,4,4) == "t","Treat","Do not treat")
      tn2 <- data.table::fifelse(substr(res2$rule,5,5) == "t","Treat","Do not treat")

      tree <- data.tree::Node$new(paste(res2$node0, "<=", res2$thres0))
      node1 <- tree$AddChild(tn1)
      node2 <- tree$AddChild(tn2)

      WT0 <- res[1,9]
      WTg <- (res2$Welfare - WT0)/WT0
      return(list(Welfare = c(Welfare = res2$Welfare,
                              Welfare0 = WT0,
                              Welfare_gain = WTg),tree = tree,
                  tnode = res[,10]))
    }
    # optimal tree has depth 2
    else{
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

      WT0 <- res[1,9]
      WTg <- (res2$Welfare - WT0)/WT0
      return(list(Welfare = c(Welfare = res2$Welfare,
                              Welfare0 = WT0,
                              Welfare_gain = WTg),tree = tree,
                  tnode = res[,10]))
    }
  }
  if (rule == "linear_rule"){
    if (!is.null(quants)){
      names_targetX <- names(targetX)
      targetX <- sapply(1:ncol(targetX), function(u){
        assign_quantiles(unlist(targetX[,u]),quants = quants[u])})
      targetX <- dplyr::as_tibble(targetX)
      names(targetX) <- names_targetX

    }
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
                        t = t,
                        verbose = verbose)
    W0 <- res$W0
    Wg <- (res$W - W0)/W0
    return(list(Welfare = c(Welfare = res$W,
                            Welfare0 = W0,
                            Welfare_gain = Wg),
                coefficients = res$L,
                plot = res$plot))
  }
}

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
ineqewm4 <-function(Y,
                    D,
                    X,
                    targetX,
                    rule = c("tree","linscore"),
                    approx = FALSE,
                    approxmatrix = NULL,
                    WF = c("ineq","IOp","IGM","util"),
                    tigm = 0,
                    X1 = NULL,
                    quants,
                    design = c("rct","observational"),
                    est_method = c("PI","LR"),
                    MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                    MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                    CF = FALSE,
                    K = 5,
                    MLiop = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                    parallel = FALSE){
  ############# Welfare when no one is treated ###################
  rl0 <- rep(0,length(Y))
  if (WF == "ineq"){
    WT0all <- wineq(Y = Y,
                    D = D,
                    X = X,
                    rule = rl0,
                    design = design,
                    est_method = est_method,
                    MLps = MLps,
                    MLalpha = MLalpha,
                    CF = CF,
                    K = K)
    WT0 <- WT0all$Welfare
    mu0 <- WT0all$Mean
    G0 <- WT0all$Gini
  }
  else if (WF == "IOp"){
    WT0all <- wiop(Y = Y,
                   D = D,
                   X = X,
                   rule = rl0,
                   design = design,
                   est_method = est_method,
                   MLiop = MLiop,
                   MLps = MLps,
                   CF = CF,
                   K = K)
    WT0 <- WT0all$Welfare
    mu0 <- WT0all$Mean
    IOp0 <- WT0all$IOp
  }
  else if (WF == "IGM"){
    WT0all <- wigm(Y = Y,
                   X1 = X1,
                   D = D,
                   X = X,
                   rule = rl0,
                   t = tigm,
                   design = design,
                   est_method = est_method,
                   MLps = MLps,
                   MLalpha = MLalpha,
                   CF = CF,
                   K = K)
    WT0 <- WT0all$Welfare
    mu0 <- WT0all$Mean
    G0 <- NA
  }
  else if (WF == "util"){
    WT0all <- wutil(Y = Y,
                    D = D,
                    X = X,
                    rule = rl0,
                    design = design,
                    est_method = est_method,
                    MLps = MLps,
                    MLreg = MLalpha,
                    CF = CF,
                    K = K)
    WT0 <- WT0all$Welfare
  }

  ##################### TReeees ###########
  s <- NULL
  if (approx == TRUE){
    ns <- approxmatrix
    for (i in 1:p){
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
    for (jj in 1:length(s[[ii]])){
      leaves0 <- targetX[,ii,drop = TRUE] <= s[[ii]][jj]
      for (kk in 1:ncol(targetX)){
        for (ll in 1:length(s[[kk]])){
          leaves1 <- (targetX[,ii,drop = TRUE] <= s[[ii]][jj])*
            (targetX[,kk,drop = TRUE] <= s[[kk]][ll]) +
            2*(targetX[,ii,drop = TRUE] <= s[[ii]][jj])*
            (targetX[,-ii,drop = TRUE] > s[[kk]][ll])
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
                if (WF == "ineq"){
                  WTall <- wineq(Y = Y,
                                 D = D,
                                 X = X,
                                 rule = rl,
                                 design = design,
                                 est_method = est_method,
                                 MLps = MLps,
                                 MLalpha = MLalpha,
                                 CF = CF,
                                 K = K)
                  WT <- WTall$Welfare
                  muT <- WTall$Mean
                  GT <- WTall$Gini
                }
                else if (WF == "IOp"){
                  WTall <- wiop(Y = Y,
                                D = D,
                                X = X,
                                rule = rl,
                                design = design,
                                est_method = est_method,
                                MLiop = MLiop,
                                MLps = MLps,
                                CF = CF,
                                K = K)
                  WT <- WTall$Welfare
                  muT <- WTall$Mean
                  GT <- WTall$IOp
                }
                else if (WF == "IGM"){
                  WTall <- wigm(Y = Y,
                                X1 = X1,
                                D = D,
                                X = X,
                                rule = rl,
                                t = tigm,
                                design = design,
                                est_method = est_method,
                                MLps = MLps,
                                MLalpha = MLalpha,
                                CF = CF,
                                K = K)
                  WT <- WTall$Welfare
                  muT <- WTall$Mean
                  GT <- NA
                }
                else if (WF == "util"){
                  WTall <- wutil(Y = Y,
                                 D = D,
                                 X = X,
                                 rule = rl,
                                 design = design,
                                 est_method = est_method,
                                 MLps = MLps,
                                 MLreg = MLalpha,
                                 CF = CF,
                                 K = K)
                  WT <- WTall$Welfare
                  muT <- NA
                  GT <- NA
                }
                data.frame(WT = WT, muT = muT, GT= GT)
              }) #here should be welfare computation


              res2 <- list(Welfare = max(restr$WT),
                           Mean = restr$muT[which.maX(restr$WT)],
                           Gini = restr$GT[which.maX(restr$WT)],
                           node0 = colnames(targetX[,ii]),
                           thres0 = s[[ii]][jj],
                           node1 = colnames(targetX[,kk]),
                           thres1 = s[[kk]][ll],
                           node2 = colnames(targetX[,pp]),
                           thres2 = s[[pp]][tt],
                           rule = names(which.max(restr)))
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
  tn2 <- data.table::fifelse(substr(res$rule,6,7) == "12","Do not Treat","Treat")
  tn3 <- data.table::fifelse(substr(res$rule,9,10) == "34","Treat","Do not treat")
  tn4 <- data.table::fifelse(substr(res$rule,9,10) == "34","Do not Treat","Treat")

  tree <- data.tree::Node$new(paste(res$node0, "<=", res$thres0))
  node1 <- tree$AddChild(paste(res$node1, "<=", res$thres1))
  trnode1_1 <- node1$AddChild(tn1)
  trnode1_2 <- node1$AddChild(tn2)
  node2 <- tree$AddChild(paste(res$node2, "<=", res$thres2))
  trnode2_3 <- node2$AddChild(tn3)
  trnode2_4 <- node2$AddChild(tn4)

  # plot(tree)





  WTall <- lapply(aux,function(u) u[[1]])
  WTall <- do.call(rbind,WTall)
  res <- WTall[which.max(WTall$WT),]
  WTmax <- res$WT
  WTg <- (res$Welfare - WT0)/WT0
  muTmax <- res$muT
  GTmax <- res$GT
  rlstar <- lapply(aux,function(u) u[[2]])
  rlstar <- rlstar[[which.max(WTall$WT)]]
  if (WF == "ineq"){
    res <- c("W0" = WT0,"Mean0" = mu0, "Gini0" = G0, "W*" = res$Welfare, "Wg*" = WTg,
             "MeanT" = res$Mean, "GiniT" = res$Gini)
  }
  else if (WF == "util" | WF == "IGM"){
    res <- c("W0" = WT0, "W*" = res$Welfare, "Wg*" = WTg)
  }
  else if (WF == "IOp"){
    res <- c("W0" = WT0, "Mean0" = mu0, "IOp0" = IOp0,
             "W*" = res$Welfare, "Wg*" = WTg,
             "MeanT" = res$Mean, "GiniT" = res$Gini)
  }

  return(list(res = res, tree = tree))
}

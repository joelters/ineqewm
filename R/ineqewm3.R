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
ineqewm3 <-function(Y,
                    D,
                    X,
                    targetX,
                    rule = c("lexi","monot"),
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
    WT0 <- wigm(Y = Y,
                X1 = X1,
                D = D,
                X = X,
                rule = rl0,
                design = design,
                est_method = est_method,
                MLps = MLps,
                MLalpha = MLaplha,
                CF = CF,
                K = K)
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
  ##################### Prepare targets to create rules ###########
  cns <- unlist(lapply(targetX,is.double))
  cns_names <- colnames(targetX[,cns])
  discr_names <- colnames(targetX[,!cns])
  if (ncol(targetX[,cns]) > 0){
    aux <- dplyr::as_tibble(lapply(targetX[,cns], function(u){
      as.numeric(ggplot2::cut_number(u, n = quants, labels = as.character(1:quants)))
    }))
    colnames(aux) <- paste("Q",cns_names,sep="")
    targetnames <- c(paste("Q",cns_names,sep=""),discr_names)
    aa <- dplyr::as_tibble(cbind(targetX,aux)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(targetnames))) %>%
      dplyr::mutate(group_id = dplyr::cur_group_id()) #group id for each unique combination
  }
  else if (ncol(targetX[,cns]) == 0){
    targetnames <- colnames(targetX)
    aa <- targetX %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(targetnames))) %>%
      dplyr::mutate(group_id = dplyr::cur_group_id()) #group id for each unique combination
  }
  aasum <- dplyr::summarise(aa,count = dplyr::n())
  if (parallel == FALSE){
    aux <- lapply(1:nrow(aasum),function(u){
      if(rule == "lexi"){
        rl <- aa$group_id <= sort(unique(aa$group_id))[u]
      }
      else if (rule == "monot"){
        rl <- (data.frame(aa[,3])[,] <= as.numeric(aasum[u,1]) &
                 data.frame(aa[,2])[,] <= as.numeric(aasum[u,2]))
      }
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
        WT <- wigm(Y = Y,
                       X1 = X1,
                       D = D,
                       X = X,
                       rule = rl,
                       t = tigm,
                       design = design,
                       est_method = est_method,
                       MLps = MLps,
                       MLalpha = MLaplha,
                       CF = CF,
                       K = K)
        muT <- NA
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
      list(data.frame(WT = WT, muT = muT, GT= GT),rl)
    })
  }
  else if (parallel == TRUE){
    n.cores <- parallel::detectCores()
    clust <- parallel::makeCluster(n.cores)
    parallel::clusterExport(clust, c("aa","aasum","rule",
                  "Y","X","D", "design",
                  "est_method","MLps","MLalpha",
                  "CF","K","X1","tigm",
                  "wineq","wiop","wigm","wutil"),
                  envir=environment())
    aux <- parallel::parLapply(clust, 1:nrow(aasum),function(u){
      if(rule == "lexi"){
        rl <- aa$group_id <= sort(unique(aa$group_id))[u]
      }
      else if (rule == "monot"){
        rl <- (data.frame(aa[,3])[,] <= as.numeric(aasum[u,1]) &
                 data.frame(aa[,2])[,] <= as.numeric(aasum[u,2]))
      }
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
        WT <- wigm(Y = Y,
                   X1 = X1,
                   D = D,
                   X = X,
                   rule = rl,
                   t = tigm,
                   design = design,
                   est_method = est_method,
                   MLps = MLps,
                   MLalpha = MLaplha,
                   CF = CF,
                   K = K)
        muT <- NA
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
      list(data.frame(WT = WT, muT = muT, GT= GT),rl)
    })
    parallel::stopCluster(clust)
  }
  WTall <- lapply(aux,function(u) u[[1]])
  WTall <- do.call(rbind,WTall)
  res <- WTall[which.max(WTall$WT),]
  WTmax <- res$WT
  WTg <- (WTmax - WT0)/WT0
  muTmax <- res$muT
  GTmax <- res$GT
  rlstar <- lapply(aux,function(u) u[[2]])
  rlstar <- rlstar[[which.max(WTall$WT)]]
  if (WF == "ineq"){
    res <- c("W0" = WT0,"Mean0" = mu0, "Gini0" = G0, "W*" = WTmax, "Wg*" = WTg,
             "MeanT" = muTmax, "GiniT" = GTmax)
  }
  else if (WF == "util" | WF == "IGM"){
    res <- c("W0" = WT0, "W*" = WTmax, "Wg*" = WTg)
  }
  else if (WF == "IOp"){
    res <- c("W0" = WT0, "Mean0" = mu0, "IOp0" = IOp0,
             "W*" = WTmax, "Wg*" = WTg,
             "MeanT" = muTmax, "IOpT" = GTmax)
  }
  ############ PLOT ############################
  rlmax <- sort(unique(aa$group_id))[which.max(WTall$WT)]
  aux <- dplyr::select(aa[aa$group_id == rlmax,],dplyr::all_of(targetnames))
  group_id <- NULL
  count <- NULL
  dfpl <- aa %>% dplyr::select(dplyr::all_of(c(targetnames,"group_id"))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(targetnames))) %>%
    dplyr::summarise(count = dplyr::n(), group_id = mean(group_id))
  dfpl <- data.frame(dfpl)

  p <- ggplot2::ggplot(dfpl, ggplot2::aes(x = dfpl[,1],
                                          y = dfpl[,2],
                                          size = count),
                       environment = environment()) +
    ggplot2::geom_point() +
    ggplot2::xlab(colnames(dfpl)[1]) +
    ggplot2::ylab(colnames(dfpl)[2]) +
    ggplot2::ggtitle(WF) +
    {if(rule == "lexi")ggplot2::annotate('rect', xmin = -0.5,
                                         xmax= dfpl[dfpl$group_id == rlmax,1] - 1, ymin=-0.5, ymax = max(dfpl[,2]),
                                         alpha=.3, fill='red')} +
    {if(rule == "lexi")ggplot2::annotate('rect', xmin=dfpl[dfpl$group_id == rlmax,1] - 1,
                                         xmax=dfpl[dfpl$group_id == rlmax,1],
                                         ymin=-0.5, ymax=dfpl[dfpl$group_id == rlmax,2], alpha=.3, fill='red')} +
    {if(rule == "monot")ggplot2::annotate('rect', xmin = 0.5,
                                          xmax= dfpl[rlmax,1],
                                          ymin=-0.5,
                                          ymax = dfpl[rlmax,2],
                                          alpha=.3, fill='red')}
  #############################
  return(list(res = res, rulemax = rlstar, plot = p))
}

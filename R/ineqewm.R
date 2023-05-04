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
ineqewm <-function(Y,
                   D,
                   X,
                   targetX,
                   rule = c("lexi","monot"),
                   WF = c("ineq","IOp","IGM","util"),tigm, X1,
                   quants,
                   design = c("rct","observational"),
                   est_method = c("PI","LR"),
                   MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                   MLalpha = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL"),
                   CF = FALSE, K = 5,
                   MLiop = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL")){
  if(rule == "lexi"){
    cns <- unlist(lapply(targetX,is.double))
    cns_names <- colnames(targetX[,cns])
    discr_names <- colnames(targetX[,!cns])
    aux <- dplyr::as_tibble(lapply(targetX[,cns], function(u){
      as.numeric(ggplot2::cut_number(u, n = quants, labels = as.character(1:quants)))
    }))
    colnames(aux) <- paste("Q",cns_names,sep="")
    targetnames <- c(paste("Q",cns_names,sep=""),discr_names)
    aa <- dplyr::as_tibble(cbind(targetX,aux)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(targetnames))) %>%
      dplyr::mutate(group_id = dplyr::cur_group_id()) #group id for each unique combination
    rl0 <- rep(0,length(aa$group_id))
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
    WT <- rep(0,length(aa$group_id))
    muT <- rep(0,length(aa$group_id))
    GT <- rep(0,length(aa$group_id))
    for (rr in 1:length(unique(aa$group_id))){
      # print(rr)
      # print(rr/length(unique(aa$group_id)))
      rl <- aa$group_id <= sort(unique(aa$group_id))[rr]
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
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$Gini
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
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$IOp
      }
      else if (WF == "IGM"){
        WT[rr] <- wigm(Y = Y,
                       X1 = X1,
                       D = D,
                       X = X,
                       rule = rl,
                       design = design,
                       est_method = est_method,
                       MLps = MLps,
                       MLalpha = MLaplha,
                       CF = CF,
                       K = K)
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
        WT[rr] <- WTall$Welfare
      }
      if (WT[rr] >= max(WT)){rlstar <- rl}
    }
    WTmax <- max(WT)
    WTg <- (WTmax - WT0)/WT0
    muTmax <- muT[which.max(WT)]
    GTmax <- GT[which.max(WT)]
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
    rlmax <- sort(unique(aa$group_id))[which.max(WT)]
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
      ggplot2::annotate('rect', xmin = -0.5,
               xmax= dfpl[dfpl$group_id == rlmax,1] - 1, ymin=-0.5, ymax = max(dfpl[,2]),
               alpha=.3, fill='red') +
      ggplot2::annotate('rect', xmin=dfpl[dfpl$group_id == rlmax,1] - 1,
               xmax=dfpl[dfpl$group_id == rlmax,1],
               ymin=-0.5, ymax=dfpl[dfpl$group_id == rlmax,2], alpha=.3, fill='red')
  }
  if(rule == "monot"){
    cns <- unlist(lapply(targetX,is.double))
    cns_names <- colnames(targetX[,cns])
    discr_names <- colnames(targetX[,!cns])
    aux <- dplyr::as_tibble(lapply(targetX[,cns], function(u){
      as.numeric(ggplot2::cut_number(u, n = quants, labels = as.character(1:quants)))
    }))
    colnames(aux) <- paste("Q",cns_names,sep="")
    targetnames <- c(paste("Q",cns_names,sep=""),discr_names)
    aa <- dplyr::as_tibble(cbind(targetX,aux)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(targetnames)))
    aasum <- dplyr::summarise(aa,count = dplyr::n())
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
    WT <- rep(0,length(nrow(aasum)))
    muT <- rep(0,length(nrow(aasum)))
    GT <- rep(0,length(nrow(aasum)))
    for (rr in 1:nrow(aasum)){
      rl <- (data.frame(aa[,3])[,] <= as.numeric(aasum[rr,1]) &
               data.frame(aa[,2])[,] <= as.numeric(aasum[rr,2]))
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
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$Gini
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
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$IOp
      }
      else if (WF == "IGM"){
        WT[rr] <- wigm(Y = Y,
                       X1 = X1,
                       D = D,
                       X = X,
                       rule = rl,
                       design = design,
                       est_method = est_method,
                       MLps = MLps,
                       MLalpha = MLaplha,
                       CF = CF,
                       K = K)
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
        WT[rr] <- WTall$Welfare
      }
      if (WT[rr] >= max(WT)){rlstar <- rl}
    }
    WTmax <- max(WT)
    WTg <- (WTmax - WT0)/WT0
    muTmax <- muT[which.max(WT)]
    GTmax <- GT[which.max(WT)]
    if (WF == "ineq"){
      res <- c("W0" = WT0,"Mean0" = mu0, "Gini0" = G0, "W*" = WTmax, "Wg*" = WTg,
               "MeanT" = muTmax, "GiniT" = GTmax)
    } else if (WF == "util" | WF == "IGM"){
      res <- c("W0" = WT0, "W*" = WTmax, "Wg*" = WTg)
    } else if (WF == "IOp"){
      res <- c("W0" = WT0, "Mean0" = mu0, "IOp0" = IOp0,
               "W*" = WTmax, "Wg*" = WTg,
               "MeanT" = muTmax, "IOpT" = GTmax)
    }
    rlmax <- which.max(WT)
    aux <- dplyr::select(aasum[rlmax,],dplyr::all_of(targetnames))

    dfpl <- data.frame(aasum)

    p <- ggplot2::ggplot(dfpl, ggplot2::aes(x = dfpl[,1],
                          y = dfpl[,2],
                          size = count),
                environment = environment()) +
      ggplot2::geom_point() +
      ggplot2::xlab(colnames(dfpl)[1]) +
      ggplot2::ylab(colnames(dfpl)[2]) +
      ggplot2::ggtitle(WF) +
      ggplot2::annotate('rect', xmin = 0.5,
               xmax= dfpl[rlmax,1],
               ymin=-0.5,
               ymax = dfpl[rlmax,2],
               alpha=.3, fill='red')
  }
  return(list(res = res, rulemax = rlstar, plot = p))
}

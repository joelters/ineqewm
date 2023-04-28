ineqewm <-function(Y,D,X,targetX,rule = c("lexi","monot"),
         WF = c("ineq","IOp","IGM","util"),tigm, X1,
         quants,
         ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "SL")){
  if(rule == "lexi"){
    cns <- unlist(lapply(targetX,is.double))
    cns_names <- colnames(targetX[,cns])
    discr_names <- colnames(targetX[,!cns])
    aux <- as_tibble(lapply(targetX[,cns], function(u){
      as.numeric(cut_number(u, n = quants, labels = as.character(1:quants)))
    }))
    colnames(aux) <- paste("Q",cns_names,sep="")
    targetnames <- c(paste("Q",cns_names,sep=""),discr_names)
    aa <- as_tibble(cbind(targetX,aux)) %>%
      group_by(across(all_of(targetnames))) %>%
      mutate(group_id = cur_group_id())
    rl0 <- rep(0,length(aa$group_id))
    if (WF == "ineq"){
      WT0all <- Wineqipw(Y,D,rl0)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Mean
      G0 <- WT0all$Gini
    } else if (WF == "IOp"){
      WT0all <- WIOpPI(Y,D,X,rl0, ML = ML)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Mean
      IOp0 <- WT0all$IOp
    } else if (WF == "IGM"){
      WT0 <- Wintergmipw(Y = Y,X1 = X1,D = D,t = tigm,rule = rl0)
    } else if (WF == "util"){
      WT0all <- Wutil(Y,D,rl0)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Welfare
      G0 <- WT0all$Gini
    }
    WT <- rep(0,length(aa$group_id))
    muT <- rep(0,length(aa$group_id))
    GT <- rep(0,length(aa$group_id))
    for (rr in 1:length(unique(aa$group_id))){
      # print(rr)
      # print(rr/length(unique(aa$group_id)))
      rl <- aa$group_id <= sort(unique(aa$group_id))[rr]
      if (WF == "ineq"){
        WTall <- Wineqipw(Y,D,rl)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$Gini
      } else if (WF == "IOp"){
        WTall <- WIOpPI(Y,D,X,rl, ML = ML)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$IOp
      } else if (WF == "IGM"){
        WT[rr] <- Wintergmipw(Y = Y,X1 = X1,D = D,X = X, t = tigm,rl)
      } else if (WF == "util"){
        WTall <- Wutil(Y,D,rl)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Welfare
        GT[rr] <- WTall$Gini
      }
    }
    WTmax <- max(WT)
    WTg <- (WTmax - WT0)/WT0
    muTmax <- muT[which.max(WT)]
    GTmax <- GT[which.max(WT)]
    if (WF != "IOp"){
      res <- c("W0" = WT0,"Mean0" = mu0, "Gini0" = G0, "W*" = WTmax, "Wg*" = WTg,
               "MeanT" = muTmax, "GiniT" = GTmax)
    } else{res <- c("W0" = WT0, "Mean0" = mu0, "IOp0" = IOp0,
                    "W*" = WTmax, "Wg*" = WTg,
                    "MeanT" = muTmax, "IOpT" = GTmax)}
    rlmax <- sort(unique(aa$group_id))[which.max(WT)]
    aux <- dplyr::select(aa[aa$group_id == rlmax,],all_of(targetnames))

    dfpl <- aa %>% dplyr::select(all_of(c(targetnames,"group_id"))) %>%
      dplyr::group_by(across(all_of(targetnames))) %>%
      dplyr::summarise(count = n(), group_id = mean(group_id))
    dfpl <- data.frame(dfpl)

    p <- ggplot(dfpl, aes(x = dfpl[,1],
                          y = dfpl[,2],
                          size = count),
                environment = environment()) +
      geom_point() +
      xlab(colnames(dfpl)[1]) +
      ylab(colnames(dfpl)[2]) +
      ggtitle(WF) +
      annotate('rect', xmin = -0.5,
               xmax= dfpl[dfpl$group_id == rlmax,1] - 1, ymin=-0.5, ymax = max(dfpl[,2]),
               alpha=.3, fill='red') +
      annotate('rect', xmin=dfpl[dfpl$group_id == rlmax,1] - 1,
               xmax=dfpl[dfpl$group_id == rlmax,1],
               ymin=-0.5, ymax=dfpl[dfpl$group_id == rlmax,2], alpha=.3, fill='red')
  }
  if(rule == "monot"){
    cns <- unlist(lapply(targetX,is.double))
    cns_names <- colnames(targetX[,cns])
    discr_names <- colnames(targetX[,!cns])
    aux <- as_tibble(lapply(targetX[,cns], function(u){
      as.numeric(cut_number(u, n = quants, labels = as.character(1:quants)))
    }))
    colnames(aux) <- paste("Q",cns_names,sep="")
    targetnames <- c(paste("Q",cns_names,sep=""),discr_names)
    aa <- as_tibble(cbind(targetX,aux)) %>%
      group_by(across(all_of(targetnames)))
    aasum <- dplyr::summarise(aa,count = n())
    rl0 <- rep(0,length(Y))
    if (WF == "ineq"){
      WT0all <- Wineqipw(Y,D,rl0)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Mean
      G0 <- WT0all$Gini
    } else if (WF == "IOp"){
      WT0all <- WIOpPI(Y,D,X,rl0, ML = ML)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Mean
      IOp0 <- WT0all$IOp
    } else if (WF == "IGM"){
      WT0 <- Wintergmipw(Y = Y,X1 = X1,D = D,t = tigm,rule = rl0)
    } else if (WF == "util"){
      WT0all <- Wutil(Y,D,rl0)
      WT0 <- WT0all$Welfare
      mu0 <- WT0all$Welfare
      G0 <- WT0all$Gini
    }
    WT <- rep(0,length(nrow(aasum)))
    muT <- rep(0,length(nrow(aasum)))
    GT <- rep(0,length(nrow(aasum)))
    for (rr in 1:nrow(aasum)){
      # print(rr)
      # print(rr/nrow(aasum))
      rl <- (data.frame(aa[,3])[,] <= as.numeric(aasum[rr,1]) &
               data.frame(aa[,2])[,] <= as.numeric(aasum[rr,2]))
      if (WF == "ineq"){
        WTall <- Wineqipw(Y,D,rl)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$Gini
      } else if (WF == "IOp"){
        WTall <- WIOpPI(Y,D,X,rl, ML = ML)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Mean
        GT[rr] <- WTall$IOp
      } else if (WF == "IGM"){
        WT[rr] <- Wintergmipw(Y = Y,X1 = X1,D = D,X = X, t = tigm,rl)
      } else if (WF == "util"){
        WTall <- Wutil(Y,D,rl)
        WT[rr] <- WTall$Welfare
        muT[rr] <- WTall$Welfare
        GT[rr] <- WTall$Gini
      }
    }
    WTmax <- max(WT)
    WTg <- (WTmax - WT0)/WT0
    muTmax <- muT[which.max(WT)]
    GTmax <- GT[which.max(WT)]
    if (WF != "IOp"){
      res <- c("W0" = WT0,"Mean0" = mu0, "Gini0" = G0, "W*" = WTmax, "Wg*" = WTg,
               "MeanT" = muTmax, "GiniT" = GTmax)
    } else{res <- c("W0" = WT0, "Mean0" = mu0, "IOp0" = IOp0,
                    "W*" = WTmax, "Wg*" = WTg,
                    "MeanT" = muTmax, "IOpT" = GTmax)}
    rlmax <- which.max(WT)
    aux <- dplyr::select(aasum[rlmax,],all_of(targetnames))

    dfpl <- data.frame(aasum)

    p <- ggplot(dfpl, aes(x = dfpl[,1],
                          y = dfpl[,2],
                          size = count),
                environment = environment()) +
      geom_point() +
      xlab(colnames(dfpl)[1]) +
      ylab(colnames(dfpl)[2]) +
      ggtitle(WF) +
      annotate('rect', xmin = 0.5,
               xmax= dfpl[rlmax,1],
               ymin=-0.5,
               ymax = dfpl[rlmax,2],
               alpha=.3, fill='red')
  }
  return(list(res = res, plot = p))
}

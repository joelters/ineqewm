#' Causal estimators
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment.
#' @param X Tibble with circumstances.
#' @param design Whether data comes from a Randomized control trial or
#' from observational data
#' @param method IPW or AIPW
#' @param estimand Population over which we want the treatment effect (only
#' relevant for observational studies)
#' @param MLreg machine learner choice for outcome regressions
#' @param MLps machine learner choice for propensity scores
#' @param CF whether to use Cross-Fitting
#' @param K Number of folds in cross-fitting
#' @param csplot Plot histograms of propensity score by treatment status
#' @param cvps Whether machine learning employed for propensity score should
#' be chosen by cross-validation
#' @param cvreg Whether machine learning employed for outcome regression should
#' be chosen by cross-validation
#' @param Kcv number of folds for cross-validation
#' @param trimm two-dimensional vector with lower and upper bound for trimming
#' low and high propensity scores
#'
#' @return A list with the output and a figure.
#' @export
causaljt <- function(Y,
                     D,
                     X,
                     design = c("rct","observational"),
                     method = c("meandiff","ols","directmethod", "ipw","lr"),
                     estimand = c("ate","att"),
                     MLreg = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Logit_lasso","SL"),
                     MLps = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB","Logit_lasso", "SL"),
                     CF = TRUE,
                     K = 5,
                     csplot = FALSE,
                     cvps = FALSE,
                     cvreg = FALSE,
                     Kcv = 5,
                     trimm = NULL){
D <- as.numeric(D)
Y <- as.numeric(Y)
if(design == "rct"){
  if (method == "meandiff"){
    ate <- mean(Y[D == 1]) - mean(Y[D == 0])
    se <- sqrt(var(Y[D == 1])/sum(D == 1) +
                    var(Y[D == 0])/sum(D == 0))
    pval <- 2*(1-pnorm(abs(ate/se)))
    return(list(results = data.frame(estimate = ate,
                                     se = se,
                                     pval = pval)))
  }
  else if (method == "ols"){
    stop("Figure out exactly demeaning and interactions and so on")
  }
  else if (method == "lr"){
    p <- mean(D)
    if (CF == TRUE){
      n <- length(Y)
      if (cvreg == TRUE){
        cvregm <- ML::MLcv(X,Y,ML = MLreg,
                           Kcv = Kcv)
        MLreg <- cvregm$mlbest
        MLregrmse <- cvregm$rmse
      }
      fv1 <- rep(0,n)
      fv0 <- rep(0,n)
      ind <- split(seq(n), seq(n) %% K)
      for (i in 1:K){
        XD <- dplyr::filter(dplyr::as_tibble(data.frame(Y = Y[-ind[[i]]],
                                                 D = D[-ind[[i]]],
                                                 X[-ind[[i]],])))
        X1 <- XD %>% dplyr::filter(D==1) %>% dplyr::select(-c("D","Y"))
        X0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(-c("D","Y"))
        Y1 <- XD %>% dplyr::filter(D==1) %>% dplyr::select(Y)
        Y0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(Y)

        mfv1 <- ML::modest(X1,Y1$Y,ML = MLreg)
        mfv0 <- ML::modest(X0,Y0$Y,ML = MLreg)

        fv1[ind[[i]]] <- ML::FVest(mfv1, X1,Y1$Y,
                                   X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
        fv0[ind[[i]]] <- ML::FVest(mfv0, X0,Y0$Y,
                                   X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
      }
      ate_sc <- fv1 - fv0 + (D/p)*(Y-fv1) - ((1-D)/(1-p))*(Y-fv0)
      ate <- mean(ate_sc)
      se <- sd(ate_sc)/sqrt(n)
      pval <- 2*(1-pnorm(abs(ate/se)))
      return(list(results = data.frame(estimate = ate, se = se, pval = pval,
                                       MLreg = MLreg,rmsereg = MLregrmse)))
    } else if (CF == FALSE){
      if (cvreg == FALSE){
        m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
        fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                         X,Y,ML = MLreg)
        m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
        fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                         X,Y,ML = MLreg)
        MLregrmse <- NULL
      }
      else if (cvreg == TRUE){
        cvregm <- ML::MLcv(X,Y,ML = MLreg,
                           Kcv = Kcv)
        MLreg <- cvregm$mlbest
        MLregrmse <- cvregm$rmse
        m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
        fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                         X,Y,ML = MLreg)
        m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
        fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                         X,Y,ML = MLreg)
      }
      ate_sc <- fv1 - fv0 + (D/ps)*(Y-fv1) - ((1-D)/(1-ps))*(Y-fv0)
      ate <- mean(ate_sc)
      return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                       MLreg = MLreg,rmsereg = MLregrmse)))
    }
  }
  else {stop("For RCT designs method has to be either meandiff,
             ols or lr")}

}
if(design == "observational"){
  if (estimand == "ate"){
    if (method == "meandiff"){
      ate <- mean(Y[D == 1]) - mean(Y[D == 0])
      se <- sqrt(var(Y[D == 1])/sum(D == 1) +
                   var(Y[D == 0])/sum(D == 0))
      pval <- 2*(1-pnorm(abs(ate/se)))
      return(list(results = data.frame(estimate = ate,
                                       se = se,
                                       pval = pval)))
    }
    else if (method == "ols"){
      stop("make sure you know the correct ols with observational data")
    }
    else if (method == "directmethod"){
      if (cvreg == FALSE){
        m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
        fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                         X,Y,ML = MLreg)
        m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
        fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                         X,Y,ML = MLreg)
        MLregrmse <- NA
      }
      else if (cvreg == TRUE){
        cvregm <- ML::MLcv(X,Y,ML = MLreg,
                           Kcv = Kcv)
        MLreg <- cvregm$mlbest
        MLregrmse <- cvregm$rmse
        m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
        fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                         X,Y,ML = MLreg)
        m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
        fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                         X,Y,ML = MLreg)
      }
      ate_sc <- fv1 - fv0
      ate <- mean(ate_sc)
      return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                       MLreg = MLreg,rmsereg = MLregrmse)))
    }
    else if (method == "ipw"){
      if (cvps == FALSE){
        ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        MLpsrmse <- NULL
      } else if (cvps == TRUE){
        cvpsm <- ML::MLcv(X,D,ML = MLps,
                      Kcv = Kcv)
        MLps <- cvpsm$mlbest
        MLpsrmse <- cvpsm$rmse
        ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
      }
      if (!is.null(trimm)){
        trps <- ps >= trimm[1] & ps <= trimm[2]
        D <- D[trps]
        Y <- Y[trps]
        ps <- ps[trps]
      }
      ate = mean(D*Y/ps - (1-D)*Y/(1-ps))
      if (csplot == TRUE){
        c1 <- rgb(173,216,230,max = 255, alpha = 180, names = "lt.blue")
        c2 <- rgb(255,192,203, max = 255, alpha = 180, names = "lt.pink")
        par(mar = c(5, 5, 5, 5) + 0.3)
        hist(ps[D==1], freq = TRUE, col = c1, axes = FALSE, xlab = "", ylab = "", main = "")
        axis(side = 1, xlim = c(0,1))
        axis(side = 4, ylab = "")
        mtext(side = 4, text = "D=1", line = 2.5, col = "blue")
        par(new=TRUE)
        hist(ps[D==0], freq = TRUE, axes = FALSE, col = c2, xlab = "", ylab = "", main = "Common Support")
        axis(side = 2)
        mtext(side = 2, text = "D=0", line = 2.5, col = "pink")
        p <- recordPlot()
        return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                MLps = MLps,rmseps = MLpsrmse), csplot = p))
      }
      else{return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                   MLps = MLps,rmseps = MLpsrmse)))}
  }
    else if (method == "lr"){
      if (CF == TRUE){
        n <- length(Y)
        if(cvps == TRUE){
          cvpsm <- ML::MLcv(X,D,ML = MLps,
                        Kcv = Kcv)
          MLps <- cvpsm$mlbest
          MLpsrmse <- cvpsm$rmse
        }
        if (cvreg == TRUE){
          cvregm <- ML::MLcv(X,Y,ML = MLreg,
                         Kcv = Kcv)
          MLreg <- cvregm$mlbest
          MLregrmse <- cvregm$rmse
        }
        fv1 <- rep(0,n)
        fv0 <- rep(0,n)
        ps <- rep(0,n)
        ind <- split(seq(n), seq(n) %% K)
        for (i in 1:K){
          mps <- ML::modest(X[-ind[[i]],],D[-ind[[i]]],ML = MLps)
          ps[ind[[i]]] <- ML::FVest(mps,X[-ind[[i]],],D[-ind[[i]]],
                                 X[ind[[i]],],D[ind[[i]]],ML = MLps)
          XD <- dplyr::filter(dplyr::as_tibble(data.frame(Y = Y[-ind[[i]]],
                                                   D = D[-ind[[i]]],
                                                   X[-ind[[i]],])))
          X1 <- XD %>% dplyr::filter(D==1) %>% dplyr::select(-c("D","Y"))
          X0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(-c("D","Y"))
          Y1 <- XD %>% dplyr::filter(D==1) %>% dplyr::select(Y)
          Y0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(Y)

          mfv1 <- ML::modest(X1,Y1$Y,ML = MLreg)
          mfv0 <- ML::modest(X0,Y0$Y,ML = MLreg)

          fv1[ind[[i]]] <- ML::FVest(mfv1, X1,Y1$Y,
                                  X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
          fv0[ind[[i]]] <- ML::FVest(mfv0, X0,Y0$Y,
                                  X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
        }
        if (!is.null(trimm)){
          trps <- ps >= trimm[1] & ps <= trimm[2]
          D <- D[trps]
          Y <- Y[trps]
          ps <- ps[trps]
          fv1 <- fv1[trps]
          fv0 <- fv0[trps]
        }
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        ate_sc <- fv1 - fv0 + (D/ps)*(Y-fv1) - ((1-D)/(1-ps))*(Y-fv0)
        ate <- mean(ate_sc)
        se <- sd(ate_sc)/sqrt(n)
        pval <- 2*(1-pnorm(abs(ate/se)))
        if (csplot == TRUE){
          c1 <- rgb(173,216,230,max = 255, alpha = 180, names = "lt.blue")
          c2 <- rgb(255,192,203, max = 255, alpha = 180, names = "lt.pink")
          par(mar = c(5, 5, 5, 5) + 0.3)
          hist(ps[D==1], freq = TRUE, col = c1, axes = FALSE, xlab = "", ylab = "", main = "")
          axis(side = 1, xlim = c(0,1))
          axis(side = 4, ylab = "")
          mtext(side = 4, text = "D=1", line = 2.5, col = "blue")
          par(new=TRUE)
          hist(ps[D==0], freq = TRUE, axes = FALSE, col = c2, xlab = "", ylab = "", main = "Common Support")
          axis(side = 2)
          mtext(side = 2, text = "D=0", line = 2.5, col = "pink")
          p <- recordPlot()
          return(list(results = data.frame(estimate = ate, se = se, pval = pval,
                                  MLps = MLps,
                                  MLreg = MLreg), csplot = p))
        }
        else{return(list(results = data.frame(estimate = ate, se = se, pval = pval,
                                     MLps = MLps,
                                     MLreg = MLreg)))}
      }
      else if (CF == FALSE){
        if (cvps == FALSE){
          ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
          ps <- ps$FVs
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
          MLpsrmse <- NULL
        } else if (cvps == TRUE){
          cvpsm <- ML::MLcv(X,D,ML = MLps,
                        Kcv = Kcv)
          MLps <- cvpsm$mlbest
          MLpsrmse <- cvpsm$rmse
          ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
          ps <- ps$FVs
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        }
        if (cvreg == FALSE){
          m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
          fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                           X,Y,ML = MLreg)
          m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
          fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                           X,Y,ML = MLreg)
          MLregrmse <- NULL
        }
        else if (cvreg == TRUE){
          cvregm <- ML::MLcv(X,Y,ML = MLreg,
                         Kcv = Kcv)
          MLreg <- cvregm$mlbest
          MLregrmse <- cvregm$rmse
          m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
          fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                           X,Y,ML = MLreg)
          m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
          fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                           X,Y,ML = MLreg)
        }
        if (!is.null(trimm)){
          trps <- ps >= trimm[1] & ps <= trimm[2]
          D <- D[trps]
          Y <- Y[trps]
          ps <- ps[trps]
          fv1 <- fv1[trps]
          fv0 <- fv0[trps]
        }
        ate_sc <- fv1 - fv0 + (D/ps)*(Y-fv1) - ((1-D)/(1-ps))*(Y-fv0)
        ate <- mean(ate_sc)
        if (csplot == TRUE){
          c1 <- rgb(173,216,230,max = 255, alpha = 180, names = "lt.blue")
          c2 <- rgb(255,192,203, max = 255, alpha = 180, names = "lt.pink")
          par(mar = c(5, 5, 5, 5) + 0.3)
          hist(ps[D==1], freq = TRUE, col = c1, axes = FALSE, xlab = "", ylab = "", main = "")
          axis(side = 1, xlim = c(0,1))
          axis(side = 4, ylab = "")
          mtext(side = 4, text = "D=1", line = 2.5, col = "blue")
          par(new=TRUE)
          hist(ps[D==0], freq = TRUE, axes = FALSE, col = c2, xlab = "", ylab = "", main = "Common Support")
          axis(side = 2)
          mtext(side = 2, text = "D=0", line = 2.5, col = "pink")
          p <- recordPlot()
          return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                  MLps = MLps,rmseps = MLpsrmse,
                                  MLreg = MLreg,rmsereg = MLregrmse), csplot = p))
        }
        else{return(list(results = data.frame(estimate = ate, se = NA, pval = NA,
                                     MLps = MLps,rmseps = MLpsrmse,
                                     MLreg = MLreg,rmsereg = MLregrmse)))}
        }
    }
  }
  if (estimand == "att"){
    if (method == "meandiff"){
      stop("is there meandiff with att?")
    }
    else if (method == "ols"){
      stop("figure out att with ols")
    }
    else if (method == "directmethod"){
      stop("figure out att with directmethod")
    }
    else if (method == "ipw"){
      if (cvps == FALSE){
        ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        MLpsrmse <- NULL
      } else if (cvps == TRUE){
        cvpsm <- ML::MLcv(X,D,ML = MLps,
                      Kcv = Kcv)
        MLps <- cvpsm$mlbest
        MLpsrmse <- cvpsm$rmse
        ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
        ps <- ps$FVs
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
      }
      if (!is.null(trimm)){
        trps <- ps >= trimm[1] & ps <= trimm[2]
        D <- D[trps]
        Y <- Y[trps]
        ps <- ps[trps]
      }
      att = mean(((D-ps)*Y)/((1-ps)*mean(D)))
      return(list(results = data.frame(estimate = att, se = NA, pval = NA,
                              MLps = MLps,rmseps = MLpsrmse)))
    }
    else if (method == "lr"){
      if (CF == TRUE){
        n <- length(Y)
        if(cvps == TRUE){
          cvpsm <- ML::MLcv(X,D,ML = MLps,
                        Kcv = Kcv)
          MLps <- cvpsm$mlbest
          MLpsrmse <- cvpsm$rmse
        }
        if (cvreg == TRUE){
          cvregm <- ML::MLcv(X,Y,ML = MLreg,
                         Kcv = Kcv)
          MLreg <- cvregm$mlbest
          MLregrmse <- cvregm$rmse
        }
        fv1 <- rep(0,n)
        fv0 <- rep(0,n)
        ps <- rep(0,n)
        ind <- split(seq(n), seq(n) %% K)
        for (i in 1:K){
          mps <- ML::modest(X[-ind[[i]],],D[-ind[[i]]],ML = MLps)
          ps[ind[[i]]] <- ML::FVest(mps,X[-ind[[i]],],D[-ind[[i]]],
                                     X[ind[[i]],],D[ind[[i]]],ML = MLps)
          XD <- dplyr::filter(dplyr::as_tibble(data.frame(Y = Y[-ind[[i]]],
                                                   D = D[-ind[[i]]],
                                                   X[-ind[[i]],])))
          X0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(-c("D","Y"))
          Y0 <- XD %>% dplyr::filter(D==0) %>% dplyr::select(Y)

          mfv0 <- ML::modest(X0,Y0$Y,ML = MLreg)
          fv0[ind[[i]]] <- ML::FVest(mfv0, X0,Y0$Y,
                                      X[ind[[i]],],Y[ind[[i]]],ML = MLreg)
        }
        if (!is.null(trimm)){
          trps <- ps >= trimm[1] & ps <= trimm[2]
          D <- D[trps]
          Y <- Y[trps]
          ps <- ps[trps]
          fv1 <- fv1[trps]
          fv0 <- fv0[trps]
        }
        ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        att_sc <- ((D-ps)*(Y-fv0))/((1-ps)*mean(D))
        att <- mean(att_sc)
        se <- sd(att_sc)/sqrt(n)
        pval <- 2*(1-pnorm(abs(att/se)))
        return(list(results = data.frame(estimate = att, se = se, pval = pval,
                                MLps = MLps,rmseps = MLpsrmse,
                                MLreg = MLreg,rmsereg = MLregrmse)))
      }
      else if (CF == FALSE){
        if (cvps == FALSE){
          ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
          ps <- ps$FVs
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
          MLpsrmse <- NULL
        } else if (cvps == TRUE){
          cvpsm <- ML::MLcv(X,D,ML = MLps,
                        Kcv = Kcv)
          MLps <- cvpsm$mlbest
          MLpsrmse <- cvpsm$rmse
          ps <- ML::MLest(X,D,ML = MLps,FVs = TRUE)
          ps <- ps$FVs
          ps <- (ps > 0 & ps < 1)*ps + 0.001*(ps <= 0) + 0.999*(ps >= 1)
        }
        if (cvreg == FALSE){
          m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
          fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                           X,Y,ML = MLreg)
          m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
          fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                           X,Y,ML = MLreg)
          MLregrmse <- NULL
        }
        else if (cvreg == TRUE){
          cvregm <- ML::MLcv(X,Y,ML = MLreg,
                         Kcv = Kcv)
          MLreg <- cvregm$mlbest
          MLregrmse <- cvregm$rmse
          m1 <- ML::modest(X[D==1,],Y[D==1],ML = MLreg)
          fv1 <- ML::FVest(m1,X[D==1,],Y[D==1],
                           X,Y,ML = MLreg)
          m0 <- ML::modest(X[D==0,],Y[D==0],ML = MLreg)
          fv0 <- ML::FVest(m1,X[D==0,],Y[D==0],
                           X,Y,ML = MLreg)
        }
        if (!is.null(trimm)){
          trps <- ps >= trimm[1] & ps <= trimm[2]
          D <- D[trps]
          Y <- Y[trps]
          ps <- ps[trps]
          fv1 <- fv1[trps]
          fv0 <- fv0[trps]
        }
        att_sc <- ((D-ps)*(Y-fv0))/((1-ps)*mean(D))
        att <- mean(att_sc)
        return(list(results = data.frame(estimate = att, se = NA, pval = pval,
                                MLps = MLps,rmseps = MLpsrmse,
                                MLreg = MLreg,rmsereg = MLregrmse)))
      }
    }
  }
}
}


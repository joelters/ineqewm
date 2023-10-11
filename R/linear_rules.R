#I follow the implementation in Kitagawa, T., & Tetenov, A. (2021).
#Equality-minded treatment choice. Journal of Business & Economic Statistics, 39(2), 561-574.

linear_rules <- function(X1, X2,
                         scores11 = NULL,
                         scores10 = NULL,
                         scores01 = NULL,
                         scores00 = NULL,
                         welfare = c("ineq","iop","igm","util"),
                         t,
                         verbose = c(0,1,2)){
  targetnames <- c(names(X1),names(X2))
  X1 = as.numeric(as.matrix(X1))
  X2 = as.numeric(as.matrix(X2))
  n <- length(X1)
  XX <- data.frame(X1 = X1, X2 = X2) %>%
    dplyr::group_by(X1,X2) %>%
    dplyr::summarize(cnt = dplyr::n()) %>%
    dplyr::ungroup()
  nupairs <- nrow(XX)

  CXX <- as.matrix(cbind(rep(1,n),X1,X2))

  # 1) hyperplanes based only on x1
  X1u <- sort(unique(X1))
  nX1u <- length(X1u)
  X1range <- X1u[nX1u] - X1u[1]
  dX1min <- min(X1u[2:nX1u] - X1u[1:(nX1u-1)])
  eps1 <- dX1min/2

  W <- -9999
  if (verbose == 0){
    if (welfare == "util") {
      s11_1 <- scores11[,1]
      s00_1 <- scores00[,1]

      s11_2 <- scores11[,2]
      s00_2 <- scores00[,2]

      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (welfare == "igm"){
              Waux <- -abs(Waux - t)
            }
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (1/n)*collapse::fsum(s11_1*rl0[s11_2] + s00_1*(1-rl0[s00_2]))
      if (welfare == "igm"){
        W0 <- -abs(W0 - t)
      }
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
    else if (welfare == "igm") {
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          Waux <- -abs(Waux - t)

          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      W0 <- -abs(W0 - t)

      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
  }
  else if (verbose == 1){
    if (welfare == "util") {
      s11_1 <- scores11[,1]
      s00_1 <- scores00[,1]

      s11_2 <- scores11[,2]
      s00_2 <- scores00[,2]

      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (welfare == "igm"){
              Waux <- -abs(Waux - t)
            }
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
        print(paste("Computing Linear rules, approximately ",
                    round(100*ia/nupairs1),
                    "% complete. Current maximum welfare is ",
                    round(W,3),sep =""))
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (1/n)*collapse::fsum(s11_1*rl0[s11_2] + s00_1*(1-rl0[s00_2]))
      if (welfare == "igm"){
        W0 <- -abs(W0 - t)
      }
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
    else if (welfare == "igm") {
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          Waux <- -abs(Waux - t)

          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
        print(paste("Computing Linear rules, approximately ",
                    round(100*ia/nupairs1),
                    "% complete. Current maximum welfare is ",
                    round(W,3),sep =""))
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      W0 <- -abs(W0 - t)

      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
        }
        print(paste("Computing Linear rules, approximately ",
                    round(100*ia/nupairs1),
                    "% complete. Current maximum welfare is ",
                    round(W,3),sep =""))
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
  }
  else if (verbose == 2){
    if (welfare == "util") {
      s11_1 <- scores11[,1]
      s00_1 <- scores00[,1]

      s11_2 <- scores11[,2]
      s00_2 <- scores00[,2]

      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (1/n)*collapse::fsum(s11_1*rl[s11_2] + s00_1*(1-rl[s00_2]))
            if (welfare == "igm"){
              Waux <- -abs(Waux - t)
            }
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
          print(paste("Computing Linear rules: first coordinate approximately ",
                      round(100*ia/nupairs1),
                      "% complete. Second coordinate approximately ",
                      round(100*ib/(nupairs1-ib1)),
                      "% complete. Current maximum welfare is ",
                      round(W,3),sep =""))
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (1/n)*collapse::fsum(s11_1*rl0[s11_2] + s00_1*(1-rl0[s00_2]))
      if (welfare == "igm"){
        W0 <- -abs(W0 - t)
      }
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
    else if (welfare == "igm") {
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          Waux <- -abs(Waux - t)

          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            Waux <- -abs(Waux - t)

            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
          print(paste("Computing Linear rules: first coordinate approximately ",
                      round(100*ia/nupairs1),
                      "% complete. Second coordinate approximately ",
                      round(100*ib/(nupairs1-ib1)),
                      "% complete. Current maximum welfare is ",
                      round(W,3),sep =""))
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      W0 <- -abs(W0 - t)

      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
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
      for (i1 in 1:nX1u){
        x1 <- X1u[i1] - eps1
        for (uu in c(-1,1)){
          sgn1 <- uu
          L <- c(sgn1*x1, -sgn1, 0)
          rl <- CXX %*% L <= 0
          Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                 s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                 s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                 s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
          if (Waux >= W){
            W <- Waux
            Lstar <- L
          }
        }
      }

      # 2 all hyperplanes that separate points with same x1, see notes on matlab code in
      # equality minded treatment choice
      X2u <- sort(unique(X2))
      nX2u <- length(X2u)
      dX2min <- min(X2u[2:nX2u] - X2u[1:(nX2u-1)])
      eps2 <- ((dX2min/dX1min)/X1range)/8

      # Cycle through all x_a in X
      nupairs1 <- nupairs - 1
      for (ia in 1:nupairs1){
        x1a <- XX$X1[ia]
        x2a <- XX$X2[ia]
        x2ap <- x2a + eps2
        x2am <- x2a - eps2
        # Cycle through all x_b in X with x1_b > x1_a
        ib1 <- ia + 1
        for (ib in ib1:nupairs){
          x1b <- XX$X1[ib]
          if (x1b > x1a){
            x2b <- XX$X2[ib]
            x2bp <- x2b + eps2
            x2bm <- x2b - eps2
            # consider 6 treatment rules

            # above x_a, above x_b, (resting on x_a and  x_b) treatment assignment decreasing in x2
            #(i.e. below hyperplane)
            L <- c(x1a*x2bp - x1b*x2ap, x2ap-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl = (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bm - x1b*x2ap, x2ap-x2bm, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # above x_a , below x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment decreasing in x2
            L <- c(x1a*x2bp - x1b*x2am, x2am-x2bp, x1b-x1a)
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }

            # below x_a, above x_b, treatment assignment increasing in x2
            L <- -L
            rl <- (CXX %*% L<=0)
            Waux <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                                   s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                                   s01_1*(1-rl[s01_2])*rl[s01_3] +
                                                   s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
            if (Waux >= W){
              W <- Waux
              Lstar <- L
            }
          }
          print(paste("Computing Linear rules: first coordinate approximately ",
                      round(100*ia/nupairs1),
                      "% complete. Second coordinate approximately ",
                      round(100*ib/(nupairs1-ib1)),
                      "% complete. Current maximum welfare is ",
                      round(W,3),sep =""))
        }
      }
      # In the first loop we already compute welfare when no one is treated,
      # probably more efficient way to get it from there but i just compute
      # it here again
      L0 <- c(-(X1u[1] - eps1), 1, 0)
      rl0 <- CXX %*% L0 <= 0
      W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl[s11_2]*rl[s11_3] +
                                           s10_1*rl[s10_2]*(1-rl[s10_3]) +
                                           s01_1*(1-rl[s01_2])*rl[s01_3] +
                                           s00_1*(1-rl[s00_2])*(1-rl[s00_3]))
      #Plot
      hyperplane <- function(x1,L) {
        x2 <- (-Lstar[1] - Lstar[2]*x1)/Lstar[3]
        return(x2)
      }
      pl <- ggplot2::ggplot(XX, ggplot2::aes(x = X1, y = X2, size = cnt)) +
        ggplot2::geom_point() +
        ggplot2::stat_function(fun = function(u) hyperplane(u,L = L),
                               geom = "area", fill = "red", alpha = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(targetnames[1]) + ggplot2::ylab(targetnames[2])

      return(list(W = W, W0 = W0, L = Lstar, plot = pl))
    }
  }
}


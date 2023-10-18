#IMPORTANT: order of observations in target X has to be the same as in scores
# i.e. the individuals ij in score Gamma_ij have to be the same as those
# in targetX[i] and targetX[j], Yes is left
Rlrpltree <- function(scores11 = NULL,
                      scores10 = NULL,
                      scores01 = NULL,
                      scores00 = NULL,
                      welfare = c("ineq","iop","igm","util"),
                      t,
                      targetX, depth = 2, s, ns,
                      verbose = c(0,1,2,3)){
  if (welfare == "util"){
    s11_1 <- scores11[,1]
    s00_1 <- scores00[,1]

    s11_2 <- scores11[,2]
    s00_2 <- scores00[,2]

    # 0) treat everyone vs treat no one
    n <- nrow(targetX)
    rl0 <- rep(0,n)
    WW0 <- c(-9999,0)
    WW1 <- c(1,1,-9999,12)
    WW2 <- c(1,1,1,1,1,1,-9999,1)
    W0 <- (1/n)*collapse::fsum(s11_1*rl0[s11_2] + s00_1*(1-rl0[s00_2]))

    rl1 <- rep(1,n)
    W1 <- (1/n)*collapse::fsum(s11_1*rl1[s11_2] + s00_1*(1-rl1[s00_2]))
    WW0 <- rbind(c(W0,0),c(W1,1),WW0)[which.max(c(c(W0,0)[1],c(W1,1)[1],WW0[1])),]
    if (W0 > WW2[7]){
      WW2[7] <- W0
      WW2[8] <- 0
      tnodemax <- rep(0,n)
    }
    if (W1 > WW2[7]){
      WW2[7] <- W1
      # we are going to indicate treat all by 99
      WW2[8] <- 99
      tnodemax <- rep(1,n)
    }

    if (depth == 0){
      return(WW0)
    }

    # 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
    for (ii in 1:ncol(targetX)){
      for (jj in 1:ns[ii]){
        rl12 <- (targetX[,ii] <= s[jj,ii])
        W12 <- (1/n)*collapse::fsum(s11_1*rl12[s11_2] + s00_1*(1-rl12[s00_2]))

        rl21 <- (targetX[,ii] > s[jj,ii])
        W21 <- (1/n)*collapse::fsum(s11_1*rl21[s11_2] + s00_1*(1-rl21[s00_2]))
        WW1 <- rbind(c(ii,jj,W12,12),
                     c(ii,jj,W21,21),
                     WW1)[which.max(c(c(ii,jj,W12,12)[3],c(ii,jj,W21,21)[3],WW1[3])),]
        if (W12 > WW2[7]){
          WW2[7] <- W12
          # we are going to indicate rule just one split treat if lower or equal by 9912
          WW2[8] <- 9912
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (treated 1) or right node (not treated 2)
          tnodemax <- rl12 + (1-rl12)*2
          # save rule
          rulemax <- rl12
        }
        if (W21 > WW2[7]){
          WW2[7] <- W21
          # we are going to indicate rule just one split treat if strictly greater by 9921
          WW2[8] <- 9921
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (not treated 1) or right node (treated 2)
          tnodemax <- rl21 + (1-rl21)*2
          # save rule
          rulemax <- rl21
        }
        # 2 Depth 2 tree
        if (depth == 2){
          #sub-trees
          if (verbose == 3){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    print(paste("Root node is ",
                                names(targetX)[ii],
                                " (",round(100*jj/ns[ii]),"% splits completed). ",
                                "Subnode 1 is ",
                                names(targetX)[kk],
                                " (",round(100*ll/ns[kk]),"% splits completed). ",
                                "Subnode 2 is ",
                                names(targetX)[pp],
                                " (",round(100*tt/ns[pp]),"% splits completed). ",
                                "Current Welfare is ",round(WW2[7],2),sep = ""))

                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (1/n)*collapse::fsum(s11_1*r1234[s11_2] + s00_1*(1-r1234[s00_2]))
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 2){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                print(paste("Root node is ",
                            names(targetX)[ii],
                            " (",round(100*jj/ns[ii]),"% splits completed). ",
                            "Subnode 1 is ",
                            names(targetX)[kk],
                            " (",round(100*ll/ns[kk]),"% splits completed). ",
                            "Current Welfare is ",round(WW2[7],2),sep = ""))
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){

                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (1/n)*collapse::fsum(s11_1*r1234[s11_2] + s00_1*(1-r1234[s00_2]))
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 1){
            print(paste("Root node is ",
                        names(targetX)[ii],
                        " (",round(100*jj/ns[ii]),"% splits completed). ",
                        "Current Welfare is ",round(WW2[7],2),sep = ""))
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (1/n)*collapse::fsum(s11_1*r1234[s11_2] + s00_1*(1-r1234[s00_2]))
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 0){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (1/n)*collapse::fsum(s11_1*r1234[s11_2] + s00_1*(1-r1234[s00_2]))
                      if (welfare == "igm"){
                        W1234 <- -abs(W1234 - t)
                      }
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (depth == 2) {return(cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax))}
    else {return(c(WW1,W0))}
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
    # 0) treat everyone vs treat no one
    n <- nrow(targetX)
    rl0 <- rep(0,n)
    WW0 <- c(-9999,0)
    WW1 <- c(1,1,-9999,12)
    WW2 <- c(1,1,1,1,1,1,-9999,1)
    tau0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl0[s11_2]*rl0[s11_3] +
                                         s10_1*rl0[s10_2]*(1-rl0[s10_3]) +
                                         s01_1*(1-rl0[s01_2])*rl0[s01_3] +
                                         s00_1*(1-rl0[s00_2])*(1-rl0[s00_3]))
    W0 <- -abs(tau0 - t)

    rl1 <- rep(1,n)
    tau1 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl1[s11_2]*rl1[s11_3] +
                                         s10_1*rl1[s10_2]*(1-rl1[s10_3]) +
                                         s01_1*(1-rl1[s01_2])*rl1[s01_3] +
                                         s00_1*(1-rl1[s00_2])*(1-rl1[s00_3]))
    W1 <- -abs(tau1 - t)

    WW0 <- rbind(c(W0,0),c(W1,1),WW0)[which.max(c(c(W0,0)[1],c(W1,1)[1],WW0[1])),]
    if (W0 > WW2[7]){
      WW2[7] <- W0
      WW2[8] <- 0
      tnodemax <- rep(0,n)
      taumax <- tau0
    }
    if (W1 > WW2[7]){
      WW2[7] <- W1
      # we are going to indicate treat all by 99
      WW2[8] <- 99
      tnodemax <- rep(1,n)
      taumax <- tau1
    }

    if (depth == 0){
      return(c(WW0,"tau" = taumax))
    }

    # 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
    for (ii in 1:ncol(targetX)){
      for (jj in 1:ns[ii]){
        rl12 <- (targetX[,ii] <= s[jj,ii])
        tau <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl12[s11_2]*rl12[s11_3] +
                                              s10_1*rl12[s10_2]*(1-rl12[s10_3]) +
                                              s01_1*(1-rl12[s01_2])*rl12[s01_3] +
                                              s00_1*(1-rl12[s00_2])*(1-rl12[s00_3]))
        W12 <- -abs(tau - t)


        rl21 <- (targetX[,ii] > s[jj,ii])
        tau <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl21[s11_2]*rl21[s11_3] +
                                              s10_1*rl21[s10_2]*(1-rl21[s10_3]) +
                                              s01_1*(1-rl21[s01_2])*rl21[s01_3] +
                                              s00_1*(1-rl21[s00_2])*(1-rl21[s00_3]))
        W21 <- -abs(tau - t)
        WW1 <- rbind(c(ii,jj,W12,12),
                     c(ii,jj,W21,21),
                     WW1)[which.max(c(c(ii,jj,W12,12)[3],c(ii,jj,W21,21)[3],WW1[3])),]
        if (W12 > WW2[7]){
          WW2[7] <- W12
          # we are going to indicate rule just one split treat if lower or equal by 9912
          WW2[8] <- 9912
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (treated 1) or right node (not treated 2)
          tnodemax <- rl12 + (1-rl12)*2
          # save rule
          rulemax <- rl12
          # save Kendall-tau
          taumax <- tau
        }
        if (W21 > WW2[7]){
          WW2[7] <- W21
          # we are going to indicate rule just one split treat if strictly greater by 9921
          WW2[8] <- 9921
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (not treated 1) or right node (treated 2)
          tnodemax <- rl21 + (1-rl21)*2
          # save rule
          rulemax <- rl21
          # save Kendall-tau
          taumax <- tau
        }
        # 2 Depth 2 tree
        if (depth == 2){
          #sub-trees
          if (verbose == 3){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    print(paste("Root node is ",
                                names(targetX)[ii],
                                " (",round(100*jj/ns[ii]),"% splits completed). ",
                                "Subnode 1 is ",
                                names(targetX)[kk],
                                " (",round(100*ll/ns[kk]),"% splits completed). ",
                                "Subnode 2 is ",
                                names(targetX)[pp],
                                " (",round(100*tt/ns[pp]),"% splits completed). ",
                                "Current Welfare is ",round(WW2[7],4),sep = ""))

                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      tau <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      W1234 <- -abs(tau - t)

                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                        # save Kendall-tau
                        taumax <- tau
                        if (W1234 == 0 & depth == 2){
                          # save Kendall-tau
                          taumax <- tau
                          jt <- cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax)
                          jt[2,1:3] <- c(taumax,tau0,tau1)
                          return(jt)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 2){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                print(paste("Root node is ",
                            names(targetX)[ii],
                            " (",round(100*jj/ns[ii]),"% splits completed). ",
                            "Subnode 1 is ",
                            names(targetX)[kk],
                            " (",round(100*ll/ns[kk]),"% splits completed). ",
                            "Current Welfare is ",round(WW2[7],4),sep = ""))
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){

                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      tau <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      W1234 <- -abs(tau - t)

                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                        # save Kendall-tau
                        taumax <- tau
                        if (W1234 == 0 & depth == 2){
                          # save Kendall-tau
                          taumax <- tau
                          jt <- cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax)
                          jt[2,1:3] <- c(taumax,tau0,tau1)
                          return(jt)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 1){
            print(paste("Root node is ",
                        names(targetX)[ii],
                        " (",round(100*jj/ns[ii]),"% splits completed). ",
                        "Current Welfare is ",round(WW2[7],4),sep = ""))
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      tau <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      W1234 <- -abs(tau - t)

                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                        # save Kendall-tau
                        taumax <- tau
                        if (W1234 == 0 & depth == 2){
                          # save Kendall-tau
                          taumax <- tau
                          jt <- cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax)
                          jt[2,1:3] <- c(taumax,tau0,tau1)
                          return(jt)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 0){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      tau <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      W1234 <- -abs(tau - t)

                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                        # save Kendall-tau
                        taumax <- tau
                        if (W1234 == 0 & depth == 2){
                          # save Kendall-tau
                          taumax <- tau
                          jt <- cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax)
                          jt[2,1:3] <- c(taumax,tau0,tau1)
                          return(jt)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (depth == 2) {
      jt <- cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax)
      jt[2,1:3] <- c(taumax,tau0,tau1)
      return(jt)
    }
    else {return(c(WW1,W0,W1,taumax))}
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
    # 0) treat everyone vs treat no one
    n <- nrow(targetX)
    rl0 <- rep(0,n)
    WW0 <- c(-9999,0)
    WW1 <- c(1,1,-9999,12)
    WW2 <- c(1,1,1,1,1,1,-9999,1)
    W0 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl0[s11_2]*rl0[s11_3] +
                               s10_1*rl0[s10_2]*(1-rl0[s10_3]) +
                               s01_1*(1-rl0[s01_2])*rl0[s01_3] +
                               s00_1*(1-rl0[s00_2])*(1-rl0[s00_3]))

    rl1 <- rep(1,n)
    W1 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl1[s11_2]*rl1[s11_3] +
                                         s10_1*rl1[s10_2]*(1-rl1[s10_3]) +
                                         s01_1*(1-rl1[s01_2])*rl1[s01_3] +
                                         s00_1*(1-rl1[s00_2])*(1-rl1[s00_3]))

    WW0 <- rbind(c(W0,0),c(W1,1),WW0)[which.max(c(c(W0,0)[1],c(W1,1)[1],WW0[1])),]
    if (W0 > WW2[7]){
      WW2[7] <- W0
      WW2[8] <- 0
      tnodemax <- rep(0,n)
      rulemax <- rl0
    }
    if (W1 > WW2[7]){
      WW2[7] <- W1
      # we are going to indicate treat all by 99
      WW2[8] <- 99
      tnodemax <- rep(1,n)
      rulemax <- rl1
    }

    if (depth == 0){
      return(WW0)
    }

    # 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
    for (ii in 1:ncol(targetX)){
      for (jj in 1:ns[ii]){
        rl12 <- (targetX[,ii] <= s[jj,ii])
        W12 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl12[s11_2]*rl12[s11_3] +
                                              s10_1*rl12[s10_2]*(1-rl12[s10_3]) +
                                              s01_1*(1-rl12[s01_2])*rl12[s01_3] +
                                              s00_1*(1-rl12[s00_2])*(1-rl12[s00_3]))

        rl21 <- (targetX[,ii] > s[jj,ii])
        W21 <- (2/(n*(n-1)))*collapse::fsum(s11_1*rl21[s11_2]*rl21[s11_3] +
                                              s10_1*rl21[s10_2]*(1-rl21[s10_3]) +
                                              s01_1*(1-rl21[s01_2])*rl21[s01_3] +
                                              s00_1*(1-rl21[s00_2])*(1-rl21[s00_3]))
        WW1 <- rbind(c(ii,jj,W12,12),
                     c(ii,jj,W21,21),
                     WW1)[which.max(c(c(ii,jj,W12,12)[3],c(ii,jj,W21,21)[3],WW1[3])),]
        if (W12 > WW2[7]){
          WW2[7] <- W12
          # we are going to indicate rule just one split treat if lower or equal by 9912
          WW2[8] <- 9912
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (treated 1) or right node (not treated 2)
          tnodemax <- rl12 + (1-rl12)*2
          # save rule
          rulemax <- rl12
        }
        if (W21 > WW2[7]){
          WW2[7] <- W21
          # we are going to indicate rule just one split treat if strictly greater by 9921
          WW2[8] <- 9921
          #split variable
          WW2[1] <- ii
          # split point
          WW2[2] <- jj
          # whether in left node (not treated 1) or right node (treated 2)
          tnodemax <- rl21 + (1-rl21)*2
          # save rule
          rulemax <- rl21
        }
        # 2 Depth 2 tree
        if (depth == 2){
          #sub-trees
          if (verbose == 3){
            for (kk in 1:ncol(targetX)){
            for (ll in 1:ns[kk]){
              for (pp in 1:ncol(targetX)){
                for (tt in 1:ns[pp]){
                  print(paste("Root node is ",
                              names(targetX)[ii],
                              " (",round(100*jj/ns[ii]),"% splits completed). ",
                              "Subnode 1 is ",
                              names(targetX)[kk],
                              " (",round(100*ll/ns[kk]),"% splits completed). ",
                              "Subnode 2 is ",
                              names(targetX)[pp],
                              " (",round(100*tt/ns[pp]),"% splits completed). ",
                              "Current Welfare is ",round(WW2[7],2),sep = ""))

                  tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                    2*rl12*(targetX[,kk] > s[ll,kk]) +
                    3*rl21*(targetX[,pp] <= s[tt,pp]) +
                    4*rl21*(targetX[,pp] > s[tt,pp])

                  rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                  rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                  rl_tntn <- (tnode == 1 | tnode == 3)
                  rl_tnnt <- (tnode == 1 | tnode == 4)
                  rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                  rl_tnnn <- (tnode == 1)

                  rl_ntnt <- (tnode == 2 | tnode == 4)
                  rl_nttn <- (tnode == 2 | tnode == 3)
                  rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                  rl_ntnn <- (tnode == 2)

                  rl_nnnt <- (tnode == 4)
                  rl_nntn <- (tnode == 3)

                  rls <- cbind(rl_ttnt,rl_tttn,
                               rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                               rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                               rl_nnnt,rl_nntn)

                  for (rr in 1:ncol(rls)){
                    r1234 <- rls[,rr]
                    W1234 <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                            s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                            s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                            s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))

                    if (W1234 > WW2[7]){
                      # ii is splitting variable at node 0 and jj the split point
                      # kk is splitting variable at first sub-tree and ll the split point
                      # pp is splitting variable at second sub-tree and tt the split point
                      # rr is rule (ttnt, etc.)
                      WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                      tnodemax <- tnode
                      # save rule
                      rulemax <- r1234
                    }
                  }
                }
              }
            }
            }
          }
          else if (verbose == 2){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                print(paste("Root node is ",
                            names(targetX)[ii],
                            " (",round(100*jj/ns[ii]),"% splits completed). ",
                            "Subnode 1 is ",
                            names(targetX)[kk],
                            " (",round(100*ll/ns[kk]),"% splits completed). ",
                            "Current Welfare is ",round(WW2[7],2),sep = ""))
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){

                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))

                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        #save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 1){
            print(paste("Root node is ",
                        names(targetX)[ii],
                        " (",round(100*jj/ns[ii]),"% splits completed). ",
                        "Current Welfare is ",round(WW2[7],2),sep = ""))
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        #save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
          else if (verbose == 0){
            for (kk in 1:ncol(targetX)){
              for (ll in 1:ns[kk]){
                for (pp in 1:ncol(targetX)){
                  for (tt in 1:ns[pp]){
                    tnode <- rl12*(targetX[,kk] <= s[ll,kk]) +
                      2*rl12*(targetX[,kk] > s[ll,kk]) +
                      3*rl21*(targetX[,pp] <= s[tt,pp]) +
                      4*rl21*(targetX[,pp] > s[tt,pp])

                    rl_ttnt <- (tnode == 1 | tnode == 2 | tnode == 4)
                    rl_tttn <- (tnode == 1 | tnode == 2 | tnode == 3)

                    rl_tntn <- (tnode == 1 | tnode == 3)
                    rl_tnnt <- (tnode == 1 | tnode == 4)
                    rl_tntt <- (tnode == 1 | tnode == 3 | tnode == 4)
                    rl_tnnn <- (tnode == 1)

                    rl_ntnt <- (tnode == 2 | tnode == 4)
                    rl_nttn <- (tnode == 2 | tnode == 3)
                    rl_nttt <- (tnode == 2 | tnode == 3 | tnode == 4)
                    rl_ntnn <- (tnode == 2)

                    rl_nnnt <- (tnode == 4)
                    rl_nntn <- (tnode == 3)

                    rls <- cbind(rl_ttnt,rl_tttn,
                                 rl_tntn,rl_tnnt,rl_tntt,rl_tnnn,
                                 rl_ntnt,rl_nttn,rl_nttt,rl_ntnn,
                                 rl_nnnt,rl_nntn)

                    for (rr in 1:ncol(rls)){
                      r1234 <- rls[,rr]
                      W1234 <- (2/(n*(n-1)))*collapse::fsum(s11_1* r1234[s11_2]* r1234[s11_3] +
                                                              s10_1* r1234[s10_2]*(1- r1234[s10_3]) +
                                                              s01_1*(1- r1234[s01_2])* r1234[s01_3] +
                                                              s00_1*(1- r1234[s00_2])*(1- r1234[s00_3]))
                      if (W1234 > WW2[7]){
                        # ii is splitting variable at node 0 and jj the split point
                        # kk is splitting variable at first sub-tree and ll the split point
                        # pp is splitting variable at second sub-tree and tt the split point
                        # rr is rule (ttnt, etc.)
                        WW2 <- c(ii,jj,kk,ll,pp,tt,W1234,rr)
                        tnodemax <- tnode
                        # save rule
                        rulemax <- r1234
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (depth == 2) {return(cbind(rbind(c(WW2,W0,W1),matrix(rep(0,(n-1)*10),n-1,10)),tnodemax,rulemax))}
    else {return(c(WW1,W0))}
  }
}


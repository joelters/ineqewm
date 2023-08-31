#IMPORTANT: order of observations in target X has to be the same as in scores
# i.e. the individuals ij in score Gamma_ij have to be the same as those
# in targetX[i] and targetX[j]
lrpltree <- function(scores11, scores10, scores01, scores00,
                     targetX, depth = 2, s, ns){
  # 0) treat everyone vs treat no one
  n <- nrow(targetX)
  rl0 <- rep(0,n)
  WW0 <- c(-9999,0)
  W0 <- (2/(n*(n-1)))*sum(scores11[,1]*rl0[scores11[,2]]*rl0[scores11[,3]] +
                             scores10[,1]*rl0[scores10[,2]]*(1-rl0[scores10[,3]]) +
                             scores01[,1]*(1-rl0[scores01[,2]])*rl0[scores01[,3]] +
                             scores00[,1]*(1-rl0[scores00[,2]])*(1-rl0[scores00[,3]]))

  rl1 <- rep(1,n)
  W1 <- (2/(n*(n-1)))*sum(scores11[,1]*rl1[scores11[,2]]*rl1[scores11[,3]] +
                             scores10[,1]*rl1[scores10[,2]]*(1-rl1[scores10[,3]]) +
                             scores01[,1]*(1-rl1[scores01[,2]])*rl1[scores01[,3]] +
                             scores00[,1]*(1-rl1[scores00[,2]])*(1-rl1[scores00[,3]]))
  WW0 <- rbind(c(W0,0),c(W1,1),WW0)[which.max(c(c(W0,0)[1],c(W1,1)[1],WW0[1])),]

  if (depth == 0){
    return(WW0)
  }

  # 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
  WW1 <- c(1,1,-9999,12)
  for (ii in 1:ncol(targetX)){
    print(ii)
    for (jj in 1:ns[ii]){
      rl12 <- (targetX[,ii] <= s[jj,ii])
      W12 <- (2/(n*(n-1)))*sum(scores11[,1]*rl12[scores11[,2]]*rl12[scores11[,3]] +
                                scores10[,1]*rl12[scores10[,2]]*(1-rl12[scores10[,3]]) +
                                scores01[,1]*(1-rl12[scores01[,2]])*rl12[scores01[,3]] +
                                scores00[,1]*(1-rl12[scores00[,2]])*(1-rl12[scores00[,3]]))

      rl21 <- (targetX[,ii] > s[jj,ii])
      W21 <- (2/(n*(n-1)))*sum(scores11[,1]*rl21[scores11[,2]]*rl21[scores11[,3]] +
                                 scores10[,1]*rl21[scores10[,2]]*(1-rl21[scores10[,3]]) +
                                 scores01[,1]*(1-rl21[scores01[,2]])*rl21[scores01[,3]] +
                                 scores00[,1]*(1-rl21[scores00[,2]])*(1-rl21[scores00[,3]]))
      WW1 <- rbind(c(ii,jj,W12,12),
                   c(ii,jj,W21,21),
                   WW1)[which.max(c(c(ii,jj,W12,12)[3],c(ii,jj,W21,21)[3],WW1[3])),]
      # 2 Depth 2 tree
      if (depth == 2){
        #sub-trees
        WW2 <- c(1,1,1,1,1,1,-9999,1)
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
                  W1234 <- (2/(n*(n-1)))*sum(scores11[,1]*rls[,rr][scores11[,2]]*rls[,rr][scores11[,3]] +
                                               scores10[,1]*rls[,rr][scores10[,2]]*(1-rls[,rr][scores10[,3]]) +
                                               scores01[,1]*(1-rls[,rr][scores01[,2]])*rls[,rr][scores01[,3]] +
                                               scores00[,1]*(1-rls[,rr][scores00[,2]])*(1-rls[,rr][scores00[,3]]))
                  WW2 <- rbind(c(ii,jj,kk,ll,pp,tt,W1234,rr),
                               WW2)[which.max(c(c(ii,jj,kk,ll,pp,tt,W1234,rr)[7],WW2[7])),]
                  if (WW2[7] == W1234){
                    tnodemax <- tnode
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (depth == 2) {return(cbind(rbind(c(WW2,W0),matrix(rep(0,(n-1)*9),n-1,9)),tnodemax))}
  else {return(c(WW1,W0))}
}


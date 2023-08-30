

dta <- dplyr::as_tibble(iris)
dta$Species <- as.integer(dta$Species)
Y <- dta$Petal.Length
X <- dplyr::select(dta,-Petal.Length)
targetX <- dplyr::select(X,c("Sepal.Length","Species","Sepal.Width"))
p <- ncol(targetX)
D <- runif(nrow(dta)) <= 0.5
# ns <- c(nrow(unique(targetX[,1])),nrow(unique(targetX[,2])),nrow(unique(targetX[,3])))
ns <- c(5,3,5)

# if approx
s <- NULL
for (i in 1:p){
  # s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
  #                   seq(1/(ns[i]+1), ns[i]/(ns[i] + 1), 1/(ns[i]+1)),
  #                   names = FALSE)),1)
  s[[i]] <- round(unique(quantile(targetX[,i,drop = TRUE],
                                  seq(0, 1, 1/(ns[i]+1)),
                                  names = FALSE)),1)
}

# ow
# s <- NULL
# for (i in 1:p){
#   s[[i]] <- sort(unique(targetX[,i,drop = TRUE]))
# }


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
            restr <- apply(rules, 2, mean) #here should be welfare computation
            res2 <- list(Welfare = max(restr),
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

tn1 <- fifelse(substr(res$rule,6,7) == "12","Treat","Do not treat")
tn2 <- fifelse(substr(res$rule,6,7) == "12","Do not Treat","Treat")
tn3 <- fifelse(substr(res$rule,9,10) == "34","Treat","Do not treat")
tn4 <- fifelse(substr(res$rule,9,10) == "34","Do not Treat","Treat")

tree <- Node$new(paste(res$node0, "<=", res$thres0))
  node1 <- tree$AddChild(paste(res$node1, "<=", res$thres1))
    trnode1_1 <- node1$AddChild(tn1)
    trnode1_2 <- node1$AddChild(tn2)
  node2 <- tree$AddChild(paste(res$node2, "<=", res$thres2))
    trnode2_3 <- node2$AddChild(tn3)
    trnode2_4 <- node2$AddChild(tn4)

plot(tree)







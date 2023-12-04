rm(list=ls())
gc()
setwd("C:/Users/david/Desktop/AdaBDP-main")
#installing required packages and libraries
library(ggplot2)
library(cowplot)
library(showtext)
library(gridExtra)
library(ggrepel)
library(rrcov)
library(MASS)
library(Matrix)
library(robustbase)
library(rlist)
library(mvoutlier)
library(mvtnorm)
library(Rcpp)
sourceCpp("block.cpp")
##########################
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

Bdp_D <- function(x, Sigma=NULL){
  n_1 <- nrow(x)
  p_1 <- ncol(x)
  list_1 <- NULL
  if (is.null(Sigma)){
    if (p_1 %% 2 == 0){
      p_2 <- p_1 / 2
      even_cols <- seq(2, p_2 * 2, 2)
      odd_cols <- seq(1, (p_2 - 1) * 2 + 1, 2)
      list_1 <- lapply(1:p_2, function(i){
        cov(x[, c(odd_cols[i], even_cols[i])])
      })
    }
    else{
      p_2 <- (p_1-1)/2
      p_3 <- p_2 + 1
      even_cols <- seq(2, p_2 * 2, 2)
      odd_cols <- seq(1, (p_2 - 1) * 2 + 1, 2)
      list_1 <- lapply(1:p_2, function(i){
        cov(x[, c(odd_cols[i], even_cols[i])])
      })
      list_1[[p_3]] <- matrix(c(var(x[, p_1]),0,0,0),2,2)
    }
  }
  else{
    if (p_1 %% 2 == 0) {
      p_2 <- p_1 / 2
      even_cols <- seq(2, p_2 * 2, 2)
      odd_cols <- seq(1, (p_2 - 1) * 2 + 1, 2)
      list_1 <- lapply(1:p_2, function(i){
        Sigma[c(odd_cols[i], even_cols[i]), c(odd_cols[i], even_cols[i])]
      })
    }
    else{
      p_2 <- (p_1-1)/2
      p_3 <- p_2 + 1
      even_cols <- seq(2, p_2 * 2, 2)
      odd_cols <- seq(1, (p_2 - 1) * 2 + 1, 2)
      list_1 <- lapply(1:p_2, function(i){
        Sigma[c(odd_cols[i], even_cols[i]), c(odd_cols[i], even_cols[i])]
      })
      list_1[[p_3]] <- matrix(c(Sigma[p_1, p_1],0,0,0),2,2)
    }
  }
  D_block <- as.matrix(bdiag_m(list_1))[1:p_1, 1:p_1]
}

BDP <- function(X, m1=NULL){
  if(is.null(m1)){
    m1 <- 100
  }
  p <- ncol(X) #number of varaibles in data
  N <- nrow(X) #number of observations in data
  h <- floor(N/2) + 1 #subset size
  #------------------------------------------------------------------------------------------
  # Reorder the components of X by correlation structure
  #------------------------------------------------------------------------------------------
  ro_x <- rep(0, p)
  ro_x[1] <- 1
  corr_x <- cor(X)
  for (i in 1:(p-1))
    try({
      od <- which(abs(corr_x[ro_x[i], ])==max(abs(corr_x[ro_x[i], -ro_x[1:i]])))[1]
      ro_x[i+1] <- ifelse(od %in% ro_x, i+1, od)
    }, silent = TRUE)
  X <- X[ ,ro_x]
  ###########################
  H0 <- matrix(0, h, m1)
  H0 <- apply(H0, 2, function(X){sort(sample(N, h))}) #initial subset for BDP
  H0 <- t(H0) #dimension m1*h
  
  H0_LTS <- matrix(0, m1, h)  #save the best BDP subset 
  detD_Bdp <- matrix(0, m1, 1) #save the best BDP detD
  dis_Bdp <- matrix(0, m1, N)
  
  for(i in 1:m1){
    Y <- X[H0[i, ], ]
    Ybar <- apply(Y, 2, mean)
    S1 <- Bdp_D(Y)
    temp <- t(t(X) - Ybar)
    dis <- sapply(1:N, function(j) blockdis(temp[j,], S1)[1])
    nn <- sort(order(dis)[1:h]) 
    crit <- 100 
    k <- 1
    while(crit !=0  & k < 10){
      Y <- X[nn, ]
      Ybar <- apply(Y, 2, mean)
      S2 <- Bdp_D(Y)
      temp <- t(t(X) - Ybar)
      detD <- blockdis(temp[1,], S2)[2]
      dis <- sapply(1:N, function(j) blockdis(temp[j,], S2)[1])
      nn2 <- sort(order(dis)[1:h])
      crit <- sum(abs(nn2-nn))
      nn <- nn2
      k <- k + 1
    }
    H0_LTS[i, ] <- nn
    detD_Bdp[i] <- detD
    dis_Bdp[i,] <- dis
  }
  loc_Bdp <- which.min(detD_Bdp)
  list(HBdp=H0_LTS[loc_Bdp, ], dis=dis_Bdp[loc_Bdp, ], X=X)
}####end function BDP

##MDP
MDP<-function(X, m1=NULL){
  if(is.null(m1)){
    m1 <- 100
  }
  p=ncol(X) #number of varaibles in data
  N=nrow(X) #number of observations in data
  h=floor(N/2)+1
  H0=matrix(0,2,m1)
  H0=apply(H0,2,function(X){sort(sample(N,2))})    #initial subset for MDP
  H0=t(H0) #dimension m1*k
  
  H0_LTS=matrix(0,m1,h )  #save the best MDP subset 
  detD_mdp=matrix(0,m1,1) #save the best MDP detD
  dis_mdp=matrix(0,m1,N)
  
  for(i in 1:m1){
    Y=X[H0[i,],]
    Ybar=apply(Y,2,mean)
    D=apply(Y, 2, var)
    detD=prod(D)
    dis=matrix(0,N,1)
    for (j in 1:N){
      temp2=as.matrix(X[j,]-Ybar)
      dis[j]=t(temp2/D)%*%temp2
    }
    nn=sort(order(dis)[1:h])
    crit=100
    
    k=1
    while(crit!=0 & k<10){
      Y=X[nn,]
      Ybar=apply(Y,2,mean)
      D1=apply(Y, 2, var)
      detD=prod(D1)
      dis=matrix(0,N,1)
      for (j in 1:N){
        temp2=as.matrix(X[j,]-Ybar)
        dis[j]=t(temp2/D1)%*%temp2
      }
      nn2=sort(order(dis)[1:h])
      crit=sum(abs(nn2-nn))
      nn=nn2
      k=k+1
    }
    
    H0_LTS[i,]=nn
    detD_mdp[i]=detD
    dis_mdp[i,]=dis
  }
  
  loc_mdp=which.min(detD_mdp)
  list(Hmdp=H0_LTS[loc_mdp,],dis=dis_mdp[loc_mdp,])
  
}#####end function MDP

# FDR control based mirror
outlier_det <- function(X, alpha, options){
  X_seq <- sort(abs(X))
  if(options=='+'){
    alpha_hat <- sapply(X_seq, function(x){(1 + length(X[X <= (-x)]))/max(1, length(X[X >= x]))})
  }else{
    alpha_hat <- sapply(X_seq, function(x){(length(X[X <= (-x)]))/max(1, length(X[X >= x]))})
  }
  threshold <- min(X_seq[which(alpha_hat <= alpha)])
  delta <- which(X >= threshold)
  a_seq <- rep(0, length(X))
  if(length(delta) > 0) a_seq[delta] <- 1
  return(list(a_seq, threshold))
}

set.seed(210)
stockdata <- as.data.frame(read.table('moneyflow.csv', head=T, sep=',', fileEncoding = 'UTF-8'))

X <- as.matrix(stockdata)
n <- nrow(X)
p <- ncol(X)
X_null <- X

alpha <- c(0.1, 0.2) #FDR level
c.rate <- c(0.15, 0.2, 0.25)
FDP <- array(0, dim=c(2, length(alpha), length(c.rate))) 
TDP <- FDP

X <- X_null
cc <- 2
ou <- rep(0, n)
n_ou <- sum(ou)
while(n_ou==0){
  ou <- rbinom(n, size=1, prob=c.rate[cc]) #randomly select contaminated indices
  n_ou <- sum(ou)
}
ou_loc <- which(ou==1)

for (ol in 1:n_ou) {
  X[ou_loc[ol], ] <- (sample(c(-1, 1), p, replace = TRUE))*X[ou_loc[ol], ]
}

MDP_output <- MDP(X)
X_H <- X[MDP_output$Hmdp, ]
n_H <- dim(X_H)[1]
T_H <- apply(X_H, 2, mean)

V_H <- diag(diag(cov(X_H)))
mah_distH <- mahalanobis(X, T_H, V_H, inverted = FALSE)
mah_distH <- mah_distH*p/median(mah_distH)

BDP_output <- BDP(X)
X <- BDP_output$X
X_H <- X[BDP_output$HBdp, ]
n_H <- dim(X_H)[1]
T_H <- apply(X_H, 2, mean)

V_H <- Bdp_D(X_H)
mah_distHB <- mahalanobis(X, T_H, V_H, inverted = FALSE)
mah_distHB <- mah_distHB*p/median(mah_distHB)

for (aa in 1:1){ #1:2 for three FDR levels
  
  result <- outlier_det(mah_distH-(n-1)/(n-3)*p, alpha[aa], '+')
  
  resultB <- outlier_det(mah_distHB-(n-1)/(n-3)*p, alpha[aa], '+')
  
  FDP[2, aa, cc] <- length(setdiff(which(result[[1]]==1), which(ou==1)))/max(1, sum(result[[1]]))
  TDP[2, aa, cc] <- 1-length(setdiff(which(ou==1), which(result[[1]]==1)))/n_ou
  
  FDP[1, aa, cc] <- length(setdiff(which(resultB[[1]]==1), which(ou==1)))/max(1, sum(resultB[[1]]))
  TDP[1, aa, cc] <- 1-length(setdiff(which(ou==1), which(resultB[[1]]==1)))/n_ou
}
# FDP
# TDP
ou_AdaBDP <- which(resultB[[1]]==1)
ou_AdaMDP <- which(result[[1]]==1)
cutoff_AdaBDP <- resultB[[2]]+(n-1)/(n-3)*p
#[1] 665.5876
cutoff_AdaMDP <- result[[2]]+(n-1)/(n-3)*p
#[1] 683.1158

df1 <- as.data.frame(t(rbind(c(1:n), mah_distHB)))
df1.highlight <- df1[ou_loc, ]
df2 <- as.data.frame(t(rbind(c(1:n), mah_distH)))
df2.highlight <- df2[ou_loc, ]

pic1 <- ggplot(df1, aes(V1, mah_distHB)) +
  geom_point(size = 3, shape = 1, fill = "transparent") +
  geom_point(data = df1.highlight, aes(V1, mah_distHB), size = 4, col = "black") +
  geom_hline(yintercept = cutoff_AdaBDP, linetype = "dotted", col = "darkcyan", size = 1) +
  annotate("text", x = 1, y = cutoff_AdaBDP, label = "cutoff: 665.59", vjust = -0.5, size = 5) +
  theme(
    text = element_text(family = "serif"),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20)
  ) +
  labs(x = "Index", y = "Modified distance", title = "method AdaBDP")

pic2 <- ggplot(df2, aes(V1, mah_distH)) +
  geom_point(size = 3, shape = 1, fill = "transparent") +
  geom_point(data = df2.highlight, aes(V1, mah_distH), size = 4, col = "black") +
  geom_hline(yintercept = cutoff_AdaMDP, linetype = "dotted", col = "darkcyan", size = 1) +
  annotate("text", x = 1, y = cutoff_AdaMDP, label = "cutoff: 683.12", vjust = -0.5, size = 5) +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20)
  ) +
  labs(x = "Index", y = "Modified distance", title = "method AdaMDP")

plot_rd <- plot_grid(pic1, pic2, ncol = 1, align = "v")
ggsave("plot_rd.pdf", plot_rd, device = cairo_pdf, width = 16.18, height = 10)

set.seed(210)
ou_Null <- array(0, dim=c(2, length(alpha))) 
X <- X_null
MDP_output <- MDP(X)
X_H <- X[MDP_output$Hmdp, ]
n_H <- dim(X_H)[1]
T_H <- apply(X_H, 2, mean)

V_H <- diag(diag(cov(X_H)))
mah_distH <- mahalanobis(X, T_H, V_H, inverted = FALSE)
mah_distH <- mah_distH*p/median(mah_distH)

BDP_output <- BDP(X)
X <- BDP_output$X
X_H <- X[BDP_output$HBdp, ]
n_H <- dim(X_H)[1]
T_H <- apply(X_H, 2, mean)

V_H <- Bdp_D(X_H)
mah_distHB <- mahalanobis(X, T_H, V_H, inverted = FALSE)
mah_distHB <- mah_distHB*p/median(mah_distHB)

for (aa in 1:2){ #1:2 for two FDR levels
  
  result <- outlier_det(mah_distH-(n-1)/(n-3)*p, alpha[aa], '+')
  
  resultB <- outlier_det(mah_distHB-(n-1)/(n-3)*p, alpha[aa], '+')
  
  ou_Null[2, aa] <- length((which(result[[1]]==1)))
  
  ou_Null[1, aa] <- length((which(resultB[[1]]==1)))
  
}
ou_Null
# Warning messages:
#   1: In min(X_seq[which(alpha_hat <= alpha)]) :
#   no non-missing arguments to min; returning Inf
# 2: In min(X_seq[which(alpha_hat <= alpha)]) :
#   no non-missing arguments to min; returning Inf
# > ou_Null
# [,1] [,2]
# [1,]    0    8
# [2,]    0    6
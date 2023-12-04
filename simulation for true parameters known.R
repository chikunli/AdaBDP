#Compute the FDR and TDR by methods: AdaBDP&AdaMDP
rm(list=ls())
gc()
setwd("C:/Users/david/Desktop/AdaBDP")

Simfunction <- function(ns){
  #---------------------------------------------------------------------------------------------------------
  # Function BDP (min the block-diagonal matrix determinant & min the modified mah.distance)
  #---------------------------------------------------------------------------------------------------------
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
      # S1 <- cov(Y)##########
      S1 <- Bdp_D(Y)
      temp <- t(t(X) - Ybar)
      dis <- sapply(1:N, function(j) blockdis(temp[j,], S1)[1])
      nn <- sort(order(dis)[1:h]) #
      crit <- 100 #
      k <- 1
      while(crit !=0  & k < 10){
        Y <- X[nn, ]
        Ybar <- apply(Y, 2, mean)
        # S2 <- cov(Y)#########
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
    # trR2_mdp=matrix(0,m1,1) #save the best MDP trR2
    dis_mdp=matrix(0,m1,N)
    
    for(i in 1:m1){
      Y=X[H0[i,],]
      Ybar=apply(Y,2,mean)
      # S1=cov(Y)# maybe there is a faster way here
      # D=diag(S1)
      D=apply(Y, 2, var)
      detD=prod(D)
      dis=matrix(0,N,1)
      for (j in 1:N){
        temp2=as.matrix(X[j,]-Ybar)
        dis[j]=t(temp2/D)%*%temp2
      }
      nn=sort(order(dis)[1:h]) #
      crit=100 #
      
      k=1
      while(crit!=0 & k<10){
        Y=X[nn,]
        Ybar=apply(Y,2,mean)
        # S2=cov(Y)# maybe there is a faster way here
        # D1=diag(S2)
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
      # ER=cor(X[nn,])
      # trR2=sum(diag(ER%*%ER))-p^2/h
      
      H0_LTS[i,]=nn
      detD_mdp[i]=detD
      # trR2_mdp[i]=trR2
      dis_mdp[i,]=dis
      
    }
    
    loc_mdp=which.min(detD_mdp)
    # list(Hmdp=H0_LTS[loc_mdp,],dis=dis_mdp[loc_mdp,],trR2=trR2_mdp[loc_mdp])
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
  
  #installing required packages and libraries
  library(rrcov)
  library(rlist)
  library(MASS)
  library(Matrix)
  library(robustbase)
  library(mvoutlier)
  library(mvtnorm)
  library(Rcpp)
  sourceCpp("block.cpp")
  
  #Model a#(AR multivariate normal)
  model_a <- function(n, p, rho, Sigma){
    Sigma <- diag(p)
    # for (j in 1:p){
    #   for (k in 1:p){
    #     Sigma[j, k] <- rho^(abs(j-k))
    #   }
    # }
    
    # for (j in 1:p){
    #   for (k in 1:p){
    #     if ((j+k)%%4 == 3 & abs(j-k)==1) {
    #       Sigma[j, k] <- rho
    #     }
    #   }
    # }
    
    for (j in 1:p){
      for (k in 1:p){
        if (abs(j-k)==1) {
          random_value <- runif(1, -0.5, -0.45)
          Sigma[j, k] <- random_value
          Sigma[k, j] <- random_value 
        }
      }
    }
    Sigma <- solve(Sigma)
    Sigma <- diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    
    X <- mvrnorm(n = n, rep(0, p), Sigma)
    list(X=X, Sigma=Sigma)
  }
  
  #Model b#(t)
  model_b <- function(n, p, rho, v, Sigma){
    Sigma <- diag(p)
    v <- 8
    # for (j in 1:p){
    #   for (k in 1:p){
    #     Sigma[j, k] <- rho^(abs(j-k))
    #   }
    # }
    
    # for (j in 1:p){
    #   for (k in 1:p){
    #     if ((j+k)%%4 == 3 & abs(j-k)==1) {
    #       Sigma[j, k] <- rho
    #     }
    #   }
    # }
    
    for (j in 1:p){
      for (k in 1:p){
        if (abs(j-k)==1) {
          random_value <- runif(1, -0.5, -0.45)
          Sigma[j, k] <- random_value
          Sigma[k, j] <- random_value 
        }
      }
    }
    Sigma <- solve(Sigma)
    Sigma <- diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    
    LL <- t(chol(Sigma))
    X <- matrix(rt(n * p, df = v), nrow = n, byrow = TRUE)
    LL <- chol(diag(p)*((v-2)/v)) %*% LL
    X <- t(LL %*% t(X))
    list(X=X, Sigma=Sigma)
  }
  
  #Model c#(Uniform&normal mixture)
  model_c <- function(n, p, rho, Sigma){
    Sigma <- diag(p)
    # 
    # for (j in 1:p){
    #   for (k in 1:p){
    #     Sigma[j, k] <- rho^(abs(j-k))
    #   }
    # }
    
    # for (j in 1:p){
    #   for (k in 1:p){
    #     if ((j+k)%%4 == 3 & abs(j-k)==1) {
    #       Sigma[j, k] <- rho
    #     }
    #   }
    # }
    
    for (j in 1:p){
      for (k in 1:p){
        if (abs(j-k)==1) {
          random_value <- runif(1, -0.5, -0.45)
          Sigma[j, k] <- random_value
          Sigma[k, j] <- random_value 
        }
      }
    }
    Sigma <- solve(Sigma)
    Sigma <- diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    
    LL <- t(chol(Sigma))
    w1 <- sqrt(sqrt(120)/12)
    w2 <- sqrt(1-sqrt(120)/12)
    X1 <- matrix(runif(n*p, min=-sqrt(3), max=sqrt(3)), nrow = n, byrow = TRUE)
    X2 <- matrix(rnorm(n*p, 0, 1), nrow = n, byrow = TRUE)
    X <- w1*X1 + w2*X2
    X <- t(LL %*% t(X))
    list(X=X, Sigma=Sigma)
  }
  
  #Model d#(Laplace)
  model_d <- function(n, p, rho, Sigma){
    Sigma <- diag(p)
    # 
    # for (j in 1:p){
    #   for (k in 1:p){
    #     Sigma[j, k] <- rho^(abs(j-k))
    #   }
    # }
    
    # for (j in 1:p){
    #   for (k in 1:p){
    #     if ((j+k)%%4 == 3 & abs(j-k)==1) {
    #       Sigma[j, k] <- rho
    #     }
    #   }
    # }
    
    for (j in 1:p){
      for (k in 1:p){
        if (abs(j-k)==1) {
          random_value <- runif(1, -0.5, -0.45)
          Sigma[j, k] <- random_value
          Sigma[k, j] <- random_value 
        }
      }
    }
    Sigma <- solve(Sigma)
    Sigma <- diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    
    LL <- t(chol(Sigma))
    u <- runif(n*p, min=0, max=1)
    X <- matrix(ifelse(u>=0.5, -log(2*(1-u))/sqrt(2), log(2*u)/sqrt(2)), nrow = n, byrow = TRUE)
    X <- t(LL %*% t(X))
    list(X=X, Sigma=Sigma)
  }
  
  n <- 300 #sample size
  # pd <- seq(300, 600, by=50) #dimension
  p <- 600
  # rhos <- seq(-0.9, 0.9, by=0.1)
  rhos <- 0.8
  alpha <- c(0.1, 0.2) #FDR level
  c.rate <- c(0.15, 0.2, 0.25) #contamination ratio
  kd <- seq(8, 10, by=1/3) #outlier strength
  m_spa <- c(1, 0.1) #sparse parameter
  
  #Method_result(with true \mu and Sigma known): 1-AdaBDP, 2-AdaMDP
  # FDP <- array(0, dim=c(2, length(alpha), length(c.rate), length(m_spa), length(rhos)))
  FDP <- array(0, dim=c(2, length(alpha), length(c.rate), length(m_spa), length(kd)))
  #2-FDP level; 3-contamination ratio; 4-sparse parameter; 5-rho
  TDP <- FDP
  
  for (pp in 1:7){ #1:19 for different rho
    # model_output <- model_a(n, p, rhos[pp], Sigma)
    # k <- kd[1]
    # model_output <- model_b(n, p, rhos[pp], Sigma)
    # k <- kd[1]
    # model_output <- model_c(n, p, rhos[pp], Sigma)
    # k <- kd[1]
    model_output <- model_d(n, p, rho = rhos, Sigma)
    k <- kd[pp]
    
    for (cc in 2:2){ #1:3 for three contamination ratios
      X <- as.matrix(model_output$X)
      n <- nrow(X)
      p <- ncol(X)
      ou <- rep(0, n)
      for (mm in 1:1){ #1:2 for two sparse settings
        n_ou <- sum(ou)
        while(n_ou==0){
          ou <- rbinom(n, size=1, prob=c.rate[cc]) #randomly select contaminated indices
          n_ou <- sum(ou)
        }
        n_ou <- sum(ou)
        m_s <- floor(p^m_spa[mm])
        ou_loc <- which(ou==1)
        bi <- runif(p)*(2*rbinom(p, 1, 0.5)-1)
        bi[sample(1:p, (p-m_s))] <- 0 #sparse setting
        bi <- bi/sqrt(sum(bi^2))
        
        X[ou_loc, ] <- mvrnorm(n = n_ou, k*bi, as.matrix(model_output$Sigma))
        
        mah_dist0 <- mahalanobis(X, rep(0, p), diag(p), inverted = FALSE)
        
        D_Sigma <- Bdp_D(X, Sigma=as.matrix(model_output$Sigma))
        mah_distB0 <- mahalanobis(X, rep(0, p), D_Sigma, inverted = FALSE)
        
        for (aa in 1:2){ #1:2 for two FDR levels
          result_0 <- outlier_det(mah_dist0-(n-1)/(n-3)*p, alpha[aa], '+')
          result_B0 <- outlier_det(mah_distB0-(n-1)/(n-3)*p, alpha[aa], '+')
          
          FDP[2, aa, cc, mm, pp] <- length(setdiff(which(result_0[[1]]==1), which(ou==1)))/max(1, sum(result_0[[1]]))
          TDP[2, aa, cc, mm, pp] <- 1-length(setdiff(which(ou==1), which(result_0[[1]]==1)))/n_ou
          
          FDP[1, aa, cc, mm, pp] <- length(setdiff(which(result_B0[[1]]==1), which(ou==1)))/max(1, sum(result_B0[[1]]))
          TDP[1, aa, cc, mm, pp] <- 1-length(setdiff(which(ou==1), which(result_B0[[1]]==1)))/n_ou
        }
      }
    }
  }
  list(FDP = FDP, TDP = TDP)
} ####end function Simfunction

library(parallel)
ns <- 500# 500
#Calculate the number of cores for parallel computing
no_cores <- 5
#Initiate cluster
cl <- makeCluster(no_cores)
clusterSetRNGStream(cl, iseed = 123)
tic <- proc.time()
result <- parLapply(cl, 1:ns, Simfunction)
proc.time()-tic
stopCluster(cl)
#Compute the average FDR & TDR
FDR <- array(0, dim=c(2, 2, 3, 2, 7)) 
rownames(FDR) <- c("AdaBDP", "AdaMDP")
TDR <- FDR

for(i in 1:ns){
  FDR <- result[[i]]$FDP + FDR
  TDR <- result[[i]]$TDP + TDR
}
FDR <- FDR/ns
TDR <- TDR/ns
#---------------------------------------------------------------------------------------------------------
# Save Model a-d results
#---------------------------------------------------------------------------------------------------------
# save(FDR, TDR, result, file="Model a_true_newcontami_sparsek10rho0.5.RData")
# result_a <- load("C:/Users/david/Desktop/AdaBDP/Model a_true_newcontami_sparsek10rho0.5.RData")
# result_aFDR <- eval(parse(text = "FDR"))
# result_aTDR <- eval(parse(text = "TDR"))
# round((result_aFDR)*100, 2)
# round((result_aTDR)*100, 2)
# 
# save(FDR, TDR, result, file="Model b_true_newcontami_sparsek10rho0.5.RData")
# result_b <- load("C:/Users/david/Desktop/AdaBDP/Model b_true_newcontami_sparsek10rho0.5.RData")
# result_bFDR <- eval(parse(text = "FDR"))
# result_bTDR <- eval(parse(text = "TDR"))
# round((result_bFDR)*100, 2)
# round((result_bTDR)*100, 2)

# save(FDR, TDR, result, file="Model c_true_newcontami_sparsek10rho0.5.RData")
# result_c <- load("C:/Users/david/Desktop/AdaBDP/Model c_true_newcontami_sparsek10rho0.5.RData")
# result_cFDR <- eval(parse(text = "FDR"))
# result_cTDR <- eval(parse(text = "TDR"))
# round((result_cFDR)*100, 2)
# round((result_cTDR)*100, 2)
# # 
save(FDR, TDR, result, file="Model d_true_newcontami_sparsediffk.RData")
result_d <- load("C:/Users/david/Desktop/AdaBDP/Model d_true_newcontami_sparsediffk.RData")
result_dFDR <- eval(parse(text = "FDR"))
result_dTDR <- eval(parse(text = "TDR"))
round((result_dFDR)*100, 2)
round((result_dTDR)*100, 2)
# 




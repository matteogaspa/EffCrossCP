library(MASS)
library(glmnet)
library(randomForest)

# Utils -----
f_e_avg <- function(p){
  ## given a vector of p-values output the min (over K) of the averages obtained using the first K p-values
  # p: vector of p-values
  K   <- length(p)  # number of p-values
  mus <- cumsum(p)/(1:K)
  return(min(mus))
}

f_u_avg <- function(p, u){
  ## given a vector of p-value output 1/(2-u)A(p)
  # p: vector of p-values
  # u: realization from a U(0,1)
  return(mean(p)/(2-u))
}

f_eu_avg <- function(p, u){
  ## given a vector of p-value output 1/(2-u)p_1, min_k (1/k) \sum_k p_k
  # p: vector of p-values
  # u: realization from a U(0,1)
  a   <- p[1]/(2-u)
  K   <- length(p)  
  mus <- cumsum(p)/(1:K)
  b   <- min(mus)
  return(min(a, b))
}

set_cc <- function(p, y, alpha){
  ## given a vector of "p-values" (one for each point y) output the cross-conf set at level alpha
  # p: vector of p-values
  # y: grid of values
  # alpha: conf. level
  n      <- length(p)
  lower <- upper <- NULL
  
  i <- 1
  while(i < n) {
    cond <- I(p[i] > alpha)
    if(cond){
      lower <- c(lower, y[i])
      j <- i
      while(j < n & cond){
        j <- j+1
        cond <- I(p[j] > alpha)
      }
      i <- j
      upper <- c(upper, y[i])
    }
    i <- i+1
  }
  if(is.null(lower)){
    return(NA)
  }else{
    return(list(set = cbind(lower, upper)))
  }
}

cov_int <- function(set, y_test){
  ## given a prediction set and a test point says if y_test in set
  # set: cross conformal pred. set (a matrix containing upper and lowers because can be an union of disjoint sets)
  if(sum(is.na(set))>0){
    return(0)
  }else{
    inds <- apply(set, 1, function(x) I(x[1] <= y_test && y_test <= x[2]))
    return(sum(inds)) 
  }
}

len_int <- function(set){
  # given a prediction set compute the length of the set
  if(sum(is.na(set))>0){
    return(0)
  }else{
    return(sum(set[,2] - set[,1]))
  }
}



cc_rf <- function(y, X, x_test, K, alpha = 0.1, ntree = 250, n_grid = 1000, grid_factor = 1){
  ## given training dataset it outputs the modified/e-modified/u-modified/standard cross conformal interval using rf
  # y: vector of responses
  # X: matrix of covariates
  # x_test: covariates for the test point
  # alpha: miscoverage rate
  # ntree: number of trees
  # mtry: number of variables randomly sampled
  # n_grid: number of points in the grid
  # grid_factor: number that premultiply (max and min)  
  
  # set hyper-parameters
  n      <- length(y)
  n_test <- NROW(x_test)
  grid   <- seq(grid_factor*min(y), grid_factor*max(y), length = n_grid)
  idx    <- sample(1:n)
  S_k    <- split(idx, cut(1:n, breaks = K, labels=FALSE))  # list of indices
  m      <- n/K
  
  
  # define the n_grid x K matrix containing the p-values for points y (one for each fold)
  p_vals <- array(NA, dim = c(n_grid, K, n_test))
  
  for(k in 1:K){
    idx_test  <- S_k[[k]]
    idx_train <- setdiff(1:n, idx_test)
    mu_hat    <- randomForest(X[idx_train,], y[idx_train], ntree = ntree)
    pred_hat  <- predict(mu_hat, as.matrix(X[idx_test,]))
    pred_new  <- predict(mu_hat, as.matrix(x_test))
    for(i in 1:n_grid){
      for(j in 1:n_test){
        inds          <- I(abs(y[idx_test] - pred_hat) >= c(abs(grid[i] - pred_new[j])))
        p_vals[i,k,j] <- (1 + sum(inds))/(m + 1)
      }
    }
  }
  
  # modified cross-conformal, e-modified cross-conformal, u-modified cross-conformal and the standard cross-conformal by Vovk
  pv_cc   <- apply(p_vals, c(1,3), mean)
  pv_ecc  <- apply(p_vals, c(1,3), f_e_avg)
  pv_ucc  <- pv_cc
  for(l in 1:n_test){
    u           <- runif(1)
    pv_ucc[,l]  <- apply(p_vals[,,l], 1, function(x) f_u_avg(x, u))
  }
  pv_eucc <- pv_cc
  for(l in 1:n_test){
    u            <- runif(1)
    pv_eucc[,l]  <- apply(p_vals[,,l], 1, function(x) f_eu_avg(x, u))
  }
  pv_ccs  <- apply(p_vals, c(1,3), function(x) (1 + sum(x*(m+1) - 1))/(n+1))
  
  int_cc   <- apply(pv_cc, 2, function(x) set_cc(x, grid, alpha + (1-alpha)*(K-1)/(K+n)))
  int_cce  <- apply(pv_ecc, 2, function(x) set_cc(x, grid, alpha + (1-alpha)*(K-1)/(K+n)))
  int_ccu  <- apply(pv_ucc, 2, function(x) set_cc(x, grid, alpha + (1-alpha)*(K-1)/(K+n)))
  int_cceu <- apply(pv_eucc, 2, function(x) set_cc(x, grid, alpha + (1-alpha)*(K-1)/(K+n)))
  int_ccs  <- apply(pv_ccs, 2, function(x) set_cc(x, grid, alpha))
  
  return(list(p_vals = p_vals, ys = grid, int_cc = int_cc, int_cce = int_cce, int_ccu = int_ccu, int_cceu = int_cceu, int_ccs = int_ccs)) 
}




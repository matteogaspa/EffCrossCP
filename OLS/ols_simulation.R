rm(list = ls())
library(conformalInference)
library(showtext)
source("utils.R")



# Replication of the simulation study in Barber et al. (2021)
dgp_data <- function(n, p, scale_fct = sqrt(10)){
  beta      <- rnorm(p)
  beta_norm <- scale_fct * beta/sqrt(sum(beta^2))
  X         <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y         <- X %*% beta_norm + rnorm(n)
  return(list(X = X, y = y))
}

# set hyper-parameters
n     <- 100                   
ps    <- seq(5, 200, by = 5)
B     <- 1000                # number of replications
alpha <- 0.1

n_met   <- 9
res_cov <- res_len <- array(NA, c(length(ps), B, n_met))

set.seed(5)
for(i in 1:length(ps)){
  p <- ps[i]
  for(j in 1:B){
    dati   <- dgp_data(n+1, p)
    y      <- dati$y[1:n]
    X      <- dati$X[1:n,]
    y_test <- dati$y[n+1]
    x_test <- dati$X[n+1,]
    
    # cross-conformal
    grid_factor <- 1
    n_grid      <- 1000
    if(75 <= p && p <= 85){
      # algorithm is unstable
      grid_factor <- 1000
      n_grid      <- 300000
    }
    cross <- cc_ols(y, X, x_test, K = 5, alpha, grid_factor = grid_factor, n_grid = n_grid)
    
    res_cov[i,j,1] <- cov_int(cross$int_cc, y_test)
    res_cov[i,j,2] <- cov_int(cross$int_cce, y_test)
    res_cov[i,j,3] <- cov_int(cross$int_ccu, y_test)
    res_cov[i,j,4] <- cov_int(cross$int_cceu, y_test)
    res_cov[i,j,5] <- cov_int(cross$int_ccs, y_test)
    
    res_len[i,j,1] <- len_int(cross$int_cc)
    res_len[i,j,2] <- len_int(cross$int_cce)
    res_len[i,j,3] <- len_int(cross$int_ccu)
    res_len[i,j,4] <- len_int(cross$int_cceu)
    res_len[i,j,5] <- len_int(cross$int_ccs)
    
    # split conformal 
    splt           <- conformal.pred.split(X, y, x_test, ols_pseudo, ols_pseudo_predict, alpha, 0.5)
    splt_int       <- matrix(c(splt$lo, splt$up), ncol = 2)
    res_cov[i,j,6] <- cov_int(splt_int, y_test)
    res_len[i,j,6] <- len_int(splt_int)
    
    # full conformal 
    full           <- conformal.pred(X, y, x_test, ols_pseudo, ols_pseudo_predict, alpha, num.grid.pts = 999, grid.factor = 1)
    full_int       <- matrix(c(full$lo, full$up), ncol = 2)
    res_cov[i,j,7] <- cov_int(full_int, y_test)
    res_len[i,j,7] <- len_int(full_int)
    
    # jackknife +
    jack           <- jackp_ols(y, X, x_test, alpha)
    jack_int       <- matrix(jack, ncol = 2)
    res_cov[i,j,8] <- cov_int(jack_int, y_test)
    res_len[i,j,8] <- len_int(jack_int)
    
    # split conformal 2*alpha
    splt2          <- conformal.pred.split(X, y, x_test, ols_pseudo, ols_pseudo_predict, 2*alpha, 0.5)
    splt_int2      <- matrix(c(splt2$lo, splt2$up), ncol = 2)
    res_cov[i,j,9] <- cov_int(splt_int2, y_test)
    res_len[i,j,9] <- len_int(splt_int2)
    
    if(j %% 10 == 0) cat("p = ", i, "; B = ", j, "\n")
  }
}


cov_p <- apply(res_cov, c(1,3), mean)
len_p <- apply(res_len, c(1,3), mean)

save(res_cov, res_len, file = "result_ols_n.RData")


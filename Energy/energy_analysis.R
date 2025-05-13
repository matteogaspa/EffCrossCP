rm(list = ls())
library(tidyr)
library(dplyr)
library(conformalInference)
source("utils2.R")
load("energy.RData")

# set hyper parameters
n <- NROW(y)
p <- NCOL(y)
n_train <- 30000
n_test  <- n-n_train
K       <- 10
alpha   <- 0.1
ntree   <- 100
funs1 <- rf.funs(ntree=ntree)

# Simulation study 
B <- 100

cov_rf_fn  <- len_rf_fn <- NULL

cov_res_rf <- len_res_rf <- matrix(NA, ncol = 10, nrow = B)

for(b in 1:B){
  set.seed(b+123)
  
  train <- sample(1:NROW(y), n_train, replace = F)
  test  <- sample(setdiff(1:NROW(y), train), n_test)
  
  ytrain <- y[train]
  Xtrain <- X[train,]
  
  ytest <- y[test]
  Xtest <- X[test,]
  
  # Cross-conf
  cr_rf <- cc_rf(y = ytrain, X = Xtrain, x_test = Xtest, K = K, alpha = alpha, ntree = ntree)
  
  # split
  splt_rf <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = Xtest, train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = alpha)
  
  # split2
  splt2_rf <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = Xtest, train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = 2*alpha)
  
  # bonf 
  splt3_rf <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = Xtest, train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = alpha, rho = (K-1)/K)
  
  cov_rf <- len_rf <- matrix(NA, ncol = 10, nrow = n_test)
  
  for(i in 1:n_test){
    # mod-cc
    cov_rf[i,1]    <- cov_int(cr_rf$int_cc[[i]]$set, ytest[i])
    len_rf[i,1]    <- len_int(cr_rf$int_cc[[i]]$set)
    # e-cc
    cov_rf[i,2]    <- cov_int(cr_rf$int_cce[[i]]$set, ytest[i])
    len_rf[i,2]    <- len_int(cr_rf$int_cce[[i]]$set)
    # u-cc
    cov_rf[i,3]    <- cov_int(cr_rf$int_ccu[[i]]$set, ytest[i])
    len_rf[i,3]    <- len_int(cr_rf$int_ccu[[i]]$set)
    # eu-cc
    cov_rf[i,4]    <- cov_int(cr_rf$int_cceu[[i]]$set, ytest[i])
    len_rf[i,4]    <- len_int(cr_rf$int_cceu[[i]]$set)
    # cc
    cov_rf[i,5]    <- cov_int(cr_rf$int_ccs[[i]]$set, ytest[i])
    len_rf[i,5]    <- len_int(cr_rf$int_ccs[[i]]$set)
    # bon
    cov_rf[i,6]    <- cov_int(cr_rf$int_bon[[i]]$set, ytest[i])
    len_rf[i,6]    <- len_int(cr_rf$int_bon[[i]]$set)
    # bonc
    cov_rf[i,7]    <- cov_int(cr_rf$int_bonc[[i]]$set, ytest[i])
    len_rf[i,7]    <- len_int(cr_rf$int_bonc[[i]]$set)
    # splt
    cov_rf[i,8]    <- cov_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2), ytest[i])
    len_rf[i,8]    <- len_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2))
    #splt2
    cov_rf[i,9]    <- cov_int(matrix(c(splt2_rf$lo[i], splt2_rf$up[i]), ncol = 2), ytest[i])
    len_rf[i,9]    <- len_int(matrix(c(splt2_rf$lo[i], splt2_rf$up[i]), ncol = 2))
    #splt3
    cov_rf[i,10]    <- cov_int(matrix(c(splt3_rf$lo[i], splt3_rf$up[i]), ncol = 2), ytest[i])
    len_rf[i,10]    <- len_int(matrix(c(splt3_rf$lo[i], splt3_rf$up[i]), ncol = 2))
  }
  cov_rf_fn    <- rbind(cov_rf_fn, cov_rf)
  len_rf_fn    <- rbind(len_rf_fn, len_rf)
  
  cov_res_rf[b,]    <- colMeans(cov_rf)
  len_res_rf[b,]    <- colMeans(len_rf)
  
  cat(b,"\n")
}

save(cov_res_rf, cov_rf_fn, len_res_rf, len_rf_fn, file = "energy_r.RData")
load("energy_r.RData")


tab <- rbind(
  #c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "full"),
  apply(len_res_rf, 2, mean),
  apply(len_res_rf, 2, sd),
  apply(len_res_rf, 2, median),
  apply(len_res_rf, 2, min),
  apply(len_res_rf, 2, max),
  apply(cov_res_rf, 2, mean)
)
tab <- tab[,c(1:5,8,10,9,6:7)]
xtable::xtable(tab, digits = 2)


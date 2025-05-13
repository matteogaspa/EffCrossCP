rm(list = ls())
library(conformalInference)
library(tidyr)
library(ggplot2)
library(MASS)
source("utils2.R")

# Load dataset
data("Boston")
df <- Boston

y <- df$medv
X <- df[,1:13]


# set hyper parameters -----
n_train <- 200
n       <- NROW(y)
n_test  <- n-n_train
K       <- 5
alpha   <- 0.1
lambda  <- 0.01
ntree   <- 20

# Simulation study 
B <- 100

cov_ols_fn <- len_ols_fn <- NULL 
cov_res_ols <- len_res_ols <- matrix(NA, ncol = 7, nrow = B)

for(b in 1:B){
  set.seed(b+45)
  
  train <- sample(1:NROW(y), n_train, replace = F)
  test  <- setdiff(1:NROW(y), train)
  
  ytrain <- y[train]
  Xtrain <- X[train,]
  
  ytest <- y[test]
  Xtest <- X[test,]
  
  # Cross-conf
  cr_ols <- cc_ols(y = ytrain, X = as.matrix(Xtrain), x_test = as.matrix(Xtest), K = K, alpha = alpha)

  # split
  splt_ols <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = ols_pseudo, predict.fun = ols_pseudo_predict, alpha = alpha)

  # full
  full_ols <- conformal.pred(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = ols_pseudo, predict.fun = ols_pseudo_predict, alpha = alpha, num.grid.pts = 200)

  cov_ols <- len_ols <- matrix(NA, nrow = n_test, ncol = 7)
  
  for(i in 1:n_test){
    # mod-cc
    cov_ols[i,1] <- cov_int(cr_ols$int_cc[[i]]$set, ytest[i])
    len_ols[i,1] <- len_int(cr_ols$int_cc[[i]]$set)
    # e-cc
    cov_ols[i,2] <- cov_int(cr_ols$int_cce[[i]]$set, ytest[i])
    len_ols[i,2] <- len_int(cr_ols$int_cce[[i]]$set)
    # u-cc
    cov_ols[i,3] <- cov_int(cr_ols$int_ccu[[i]]$set, ytest[i])
    len_ols[i,3] <- len_int(cr_ols$int_ccu[[i]]$set)
    # eu-cc
    cov_ols[i,4] <- cov_int(cr_ols$int_cceu[[i]]$set, ytest[i])
    len_ols[i,4] <- len_int(cr_ols$int_cceu[[i]]$set)
    # cc
    cov_ols[i,5] <- cov_int(cr_ols$int_ccs[[i]]$set, ytest[i])
    len_ols[i,5] <- len_int(cr_ols$int_ccs[[i]]$set)
    # splt
    cov_ols[i,6] <- cov_int(matrix(c(splt_ols$lo[i], splt_ols$up[i]), ncol = 2), ytest[i])
    len_ols[i,6] <- len_int(matrix(c(splt_ols$lo[i], splt_ols$up[i]), ncol = 2))
    # full
    cov_ols[i,7] <- cov_int(matrix(c(full_ols$lo[i], full_ols$up[i]), ncol = 2), ytest[i])
    len_ols[i,7] <- len_int(matrix(c(full_ols$lo[i], full_ols$up[i]), ncol = 2))
  }
  cov_ols_fn <- rbind(cov_ols_fn, cov_ols)
  len_ols_fn <- rbind(len_ols_fn, len_ols)

  cov_res_ols[b,] <- colMeans(cov_ols)
  len_res_ols[b,] <- colMeans(len_ols)

  cat(b,"\n")
}


save(cov_res_ols, len_res_ols, cov_ols_fn, len_ols_fn, file = "boston.RData")
load("boston.RData")

tab <- rbind(
  #c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "full"),
  apply(len_res_ols, 2, mean),
  apply(len_res_ols, 2, sd),
  apply(len_res_ols, 2, median),
  apply(len_res_ols, 2, min),
  apply(len_res_ols, 2, max),
  apply(cov_res_ols, 2, mean)
)

rownames(tab) <- c("Mean", "Sd", "Median", "Min", "Max", "Coverage")
xtable::xtable(tab, digits = 3)


rm(list = ls())
library(conformalInference)
library(tidyr)
library(ggplot2)
source("utils2.R")

# Load dataset
df <- read.csv("abalone.csv", header = T, stringsAsFactors = T)

y <- df$Rings
X1 <- model.matrix(~df$Sex)   # cat. vars
X2 <- df[,2:8]                # num. vars
X <- cbind(X1, X2)            # covariates


# set hyper parameters -----
n_train <- 4000
n       <- NROW(y)
n_test  <- n-n_train
K       <- 10
alpha   <- 0.1
ntree   <- 25
funs1   <- rf.funs(ntree = ntree)

# Simulation study 
B <- 100

cov_rf_fn  <- len_rf_fn <- NULL 
cov_res_rf <- len_res_rf <- matrix(NA, ncol = 7, nrow = B)

for(b in 1:B){
  set.seed(b+100)
  
  train <- sample(1:NROW(y), n_train, replace = F)
  test  <- setdiff(1:NROW(y), train)
  
  ytrain <- y[train]
  Xtrain <- X[train,]
  
  ytest <- y[test]
  Xtest <- X[test,]
  
  # Cross-conf
  cr_rf <- cc_rf(y = ytrain, X = as.matrix(Xtrain), x_test = as.matrix(Xtest), K = K, alpha = alpha, ntree = ntree)
  
  # split
  splt_rf <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = alpha)
  
  # split 2alpha
  splt_rf2 <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = 2*alpha)
  
  cov_rf <- len_rf <- matrix(NA, nrow = n_test, ncol = 7)
  
  for(i in 1:n_test){
    # mod-cc
    cov_rf[i,1] <- cov_int(cr_rf$int_cc[[i]]$set, ytest[i])
    len_rf[i,1] <- len_int(cr_rf$int_cc[[i]]$set)
    # e-cc
    cov_rf[i,2] <- cov_int(cr_rf$int_cce[[i]]$set, ytest[i])
    len_rf[i,2] <- len_int(cr_rf$int_cce[[i]]$set)
    # u-cc
    cov_rf[i,3] <- cov_int(cr_rf$int_ccu[[i]]$set, ytest[i])
    len_rf[i,3] <- len_int(cr_rf$int_ccu[[i]]$set)
    # eu-cc
    cov_rf[i,4] <- cov_int(cr_rf$int_cceu[[i]]$set, ytest[i])
    len_rf[i,4] <- len_int(cr_rf$int_cceu[[i]]$set)
    # cc
    cov_rf[i,5] <- cov_int(cr_rf$int_ccs[[i]]$set, ytest[i])
    len_rf[i,5] <- len_int(cr_rf$int_ccs[[i]]$set)
    # splt
    cov_rf[i,6] <- cov_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2), ytest[i])
    len_rf[i,6] <- len_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2))
    # splt2
    cov_rf[i,7] <- cov_int(matrix(c(splt_rf2$lo[i], splt_rf2$up[i]), ncol = 2), ytest[i])
    len_rf[i,7] <- len_int(matrix(c(splt_rf2$lo[i], splt_rf2$up[i]), ncol = 2))
  }
  cov_rf_fn <- rbind(cov_rf_fn, cov_rf)
  len_rf_fn <- rbind(len_rf_fn, len_rf)
  
  cov_res_rf[b,] <- colMeans(cov_rf)
  len_res_rf[b,] <- colMeans(len_rf)
  
  cat(b,"\n")
}

save(cov_res_rf, len_res_rf, file = "abalone100.Rdata")

tab <- rbind(
  #c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "split(2alpha)"),
  apply(len_res_rf, 2, mean),
  apply(len_res_rf, 2, sd),
  apply(len_res_rf, 2, median),
  apply(len_res_rf, 2, min),
  apply(len_res_rf, 2, max),
  apply(cov_res_rf, 2, mean)
)

rownames(tab) <- c("Mean", "Sd", "Median", "Min", "Max", "Coverage")
xtable::xtable(tab, digits = 3)

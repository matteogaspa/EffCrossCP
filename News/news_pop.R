rm(list = ls())
library(conformalInference)
library(tidyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(gridExtra)
source("utils2.R")

# Load dataset
df <- read.csv("news_pop.csv", sep = ",", header = T, stringsAsFactors = FALSE, na.strings = "?")
head(df)

n <- nrow(df)
p <- ncol(df)

y <- log(df$shares)
X <- df[,1:(p-1)]
X <- X[,-c(4,5,19)]

# set hyper parameters
n <- NROW(y)
n_train <- 10000
n_test  <- 2500
K       <- 10
alpha   <- 0.1
lambda  <- 0.2
ntree   <- 200
funs1 <- lasso.funs(lambda=lambda,standardize=F)
funs2 <- rf.funs(ntree=ntree)

# Simulation study 
B <- 100

cov_lm_fn    <- len_lm_fn    <- NULL
cov_lasso_fn <- len_lasso_fn <- NULL 
cov_rf_fn    <- len_rf_fn    <- NULL

cov_res_lm <- cov_res_lasso <- cov_res_rf <- len_res_lm <- len_res_lasso <- len_res_rf <- matrix(NA, ncol = 7, nrow = B)

for(b in 1:B){
  set.seed(b)
  
  train <- sample(1:NROW(y), n_train, replace = F)
  test  <- sample(setdiff(1:NROW(y), train), n_test)
  
  ytrain <- y[train]
  Xtrain <- X[train,]
  
  ytest <- y[test]
  Xtest <- X[test,]
  
  # Cross-conf
  cr_lm    <- cc_ols(y = ytrain, X = as.matrix(Xtrain), x_test = as.matrix(Xtest), K = K, alpha = alpha)
  cr_lasso <- cc_lasso(y = ytrain, X = Xtrain, x_test = Xtest, K = K, alpha = alpha, lambda = lambda)
  cr_rf    <- cc_rf(y = ytrain, X = Xtrain, x_test = Xtest, K = K, alpha = alpha, ntree = ntree)
  
  # split
  splt_lm    <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = ols_pseudo, predict.fun = ols_pseudo_predict, alpha = alpha)
  splt_lasso <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = alpha)
  splt_rf    <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs2$train.fun, predict.fun = funs2$predict.fun, alpha = alpha)
  
  # split2
  splt2_lm    <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = ols_pseudo, predict.fun = ols_pseudo_predict, alpha = 2*alpha)
  splt2_lasso <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs1$train.fun, predict.fun = funs1$predict.fun, alpha = 2*alpha)
  splt2_rf    <- conformal.pred.split(x = Xtrain, y = ytrain, x0 = as.matrix(Xtest), train.fun = funs2$train.fun, predict.fun = funs2$predict.fun, alpha = 2*alpha)
  
  cov_lm <- cov_lasso <- cov_rf <- len_lm <- len_lasso <- len_rf <- matrix(NA, ncol = 7, nrow = n_test)
  
  for(i in 1:n_test){
    # mod-cc
    cov_lm[i,1]    <- cov_int(cr_lm$int_cc[[i]]$set, ytest[i])
    cov_lasso[i,1] <- cov_int(cr_lasso$int_cc[[i]]$set, ytest[i])
    cov_rf[i,1]    <- cov_int(cr_rf$int_cc[[i]]$set, ytest[i])
    len_lm[i,1]    <- len_int(cr_lm$int_cc[[i]]$set)
    len_lasso[i,1] <- len_int(cr_lasso$int_cc[[i]]$set)
    len_rf[i,1]    <- len_int(cr_rf$int_cc[[i]]$set)
    # e-cc
    cov_lm[i,2]    <- cov_int(cr_lm$int_cce[[i]]$set, ytest[i])
    cov_lasso[i,2] <- cov_int(cr_lasso$int_cce[[i]]$set, ytest[i])
    cov_rf[i,2]    <- cov_int(cr_rf$int_cce[[i]]$set, ytest[i])
    len_lm[i,2]    <- len_int(cr_lm$int_cce[[i]]$set)
    len_lasso[i,2] <- len_int(cr_lasso$int_cce[[i]]$set)
    len_rf[i,2]    <- len_int(cr_rf$int_cce[[i]]$set)
    # u-cc
    cov_lm[i,3]    <- cov_int(cr_lm$int_ccu[[i]]$set, ytest[i])
    cov_lasso[i,3] <- cov_int(cr_lasso$int_ccu[[i]]$set, ytest[i])
    cov_rf[i,3]    <- cov_int(cr_rf$int_ccu[[i]]$set, ytest[i])
    len_lm[i,3]    <- len_int(cr_lm$int_ccu[[i]]$set)
    len_lasso[i,3] <- len_int(cr_lasso$int_ccu[[i]]$set)
    len_rf[i,3]    <- len_int(cr_rf$int_ccu[[i]]$set)
    # eu-cc
    cov_lm[i,4]    <- cov_int(cr_lm$int_cceu[[i]]$set, ytest[i])
    cov_lasso[i,4] <- cov_int(cr_lasso$int_cceu[[i]]$set, ytest[i])
    cov_rf[i,4]    <- cov_int(cr_rf$int_cceu[[i]]$set, ytest[i])
    len_lm[i,4]    <- len_int(cr_lm$int_cceu[[i]]$set)
    len_lasso[i,4] <- len_int(cr_lasso$int_cceu[[i]]$set)
    len_rf[i,4]    <- len_int(cr_rf$int_cceu[[i]]$set)
    # cc
    cov_lm[i,5]    <- cov_int(cr_lm$int_ccs[[i]]$set, ytest[i])
    cov_lasso[i,5] <- cov_int(cr_lasso$int_ccs[[i]]$set, ytest[i])
    cov_rf[i,5]    <- cov_int(cr_rf$int_ccs[[i]]$set, ytest[i])
    len_lm[i,5]    <- len_int(cr_lm$int_ccs[[i]]$set)
    len_lasso[i,5] <- len_int(cr_lasso$int_ccs[[i]]$set)
    len_rf[i,5]    <- len_int(cr_rf$int_ccs[[i]]$set)
    # splt
    cov_lm[i,6]    <- cov_int(matrix(c(splt_lm$lo[i], splt_lm$up[i]), ncol = 2), ytest[i])
    cov_lasso[i,6] <- cov_int(matrix(c(splt_lasso$lo[i], splt_lasso$up[i]), ncol = 2), ytest[i])
    cov_rf[i,6]    <- cov_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2), ytest[i])
    len_lm[i,6]    <- len_int(matrix(c(splt_lm$lo[i], splt_lm$up[i]), ncol = 2))
    len_lasso[i,6] <- len_int(matrix(c(splt_lasso$lo[i], splt_lasso$up[i]), ncol = 2))
    len_rf[i,6]    <- len_int(matrix(c(splt_rf$lo[i], splt_rf$up[i]), ncol = 2))
    #splt2
    cov_lm[i,7]    <- cov_int(matrix(c(splt2_lm$lo[i], splt2_lm$up[i]), ncol = 2), ytest[i])
    cov_lasso[i,7] <- cov_int(matrix(c(splt2_lasso$lo[i], splt2_lasso$up[i]), ncol = 2), ytest[i])
    cov_rf[i,7]    <- cov_int(matrix(c(splt2_rf$lo[i], splt2_rf$up[i]), ncol = 2), ytest[i])
    len_lm[i,7]    <- len_int(matrix(c(splt2_lm$lo[i], splt2_lm$up[i]), ncol = 2))
    len_lasso[i,7] <- len_int(matrix(c(splt2_lasso$lo[i], splt2_lasso$up[i]), ncol = 2))
    len_rf[i,7]    <- len_int(matrix(c(splt2_rf$lo[i], splt2_rf$up[i]), ncol = 2))
  }
  cov_lm_fn    <- rbind(cov_lm_fn, cov_lm)
  cov_lasso_fn <- rbind(cov_lasso_fn, cov_lasso)
  cov_rf_fn    <- rbind(cov_rf_fn, cov_rf)
  len_lm_fn    <- rbind(len_lm_fn, len_lm)
  len_lasso_fn <- rbind(len_lasso_fn, len_lasso)
  len_rf_fn    <- rbind(len_rf_fn, len_rf)
  
  cov_res_lm[b,]    <- colMeans(cov_lm)
  cov_res_lasso[b,] <- colMeans(cov_lasso)
  cov_res_rf[b,]    <- colMeans(cov_rf)
  len_res_lm[b,]    <- colMeans(len_lm)
  len_res_lasso[b,] <- colMeans(len_lasso)
  len_res_rf[b,]    <- colMeans(len_rf)
  
  cat(b,"\n")
}

save(cov_res_lm, cov_res_lasso, cov_res_rf, cov_lm_fn, cov_lasso_fn, cov_rf_fn,
     len_res_lm, len_res_lasso, len_res_rf, len_lm_fn, len_lasso_fn, len_rf_fn,
     file = "news.RData")
load("news.RData")

colnames(cov_res_lm) <- colnames(cov_res_lasso) <- colnames(cov_res_rf) <- colnames(len_res_lm) <- colnames(len_res_lasso) <- colnames(len_res_rf) <- c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "split(2α)")
cov_res_lm    <- as.data.frame(cov_res_lm)
cov_res_lasso <- as.data.frame(cov_res_lasso)
cov_res_rf    <- as.data.frame(cov_res_rf)
len_res_lm    <- as.data.frame(len_res_lm)
len_res_lasso <- as.data.frame(len_res_lasso)
len_res_rf    <- as.data.frame(len_res_rf)

len_res_lm_long        <- pivot_longer(len_res_lm, cols = everything(), names_to = "Col", values_to = "Val")
len_res_lm_long$Col    <- factor(len_res_lm_long$Col, levels = colnames(len_res_lm[,]))
len_res_lasso_long     <- pivot_longer(len_res_lasso, cols = everything(), names_to = "Col", values_to = "Val")
len_res_lasso_long$Col <- factor(len_res_lasso_long$Col, levels = colnames(len_res_lasso[,]))
len_res_rf_long        <- pivot_longer(len_res_rf, cols = everything(), names_to = "Col", values_to = "Val")
len_res_rf_long$Col    <- factor(len_res_rf_long$Col, levels = colnames(len_res_rf[,]))

p0 <- ggplot(len_res_lm_long, aes(x = Col, y = Val, fill = Col)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "LM", x = "", y = "Size") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") + ylim(1.8, 3)


p1 <- ggplot(len_res_lasso_long, aes(x = Col, y = Val, fill = Col)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Lasso", x = "", y = "Size") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") + ylim(1.8, 3)

p2 <- ggplot(len_res_rf_long, aes(x = Col, y = Val, fill = Col)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Random Forest", x = "", y = "Size") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") + ylim(1.8, 3)

df_cov <- data.frame(
  method = c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "split(2α)"),
  LM     = colMeans(cov_res_lm),
  LM1    = apply(cov_res_lm, 2, min),
  LM2    = apply(cov_res_lm, 2, max),
  Lasso  = colMeans(cov_res_lasso),
  Lasso1 = apply(cov_res_lasso, 2, min),
  Lasso2 = apply(cov_res_lasso, 2, max),
  RF     = colMeans(cov_res_rf),
  RF1    = apply(cov_res_rf, 2, min),
  RF2    = apply(cov_res_rf, 2, max)
)

#df_cov <- data.frame(
#  method = c("mod-cross", "e-mod-cross", "u-mod-cross", "eu-mod-cross", "cross", "split", "split(2α)"),
#  LM     = colMeans(cov_res_lm),
#  LM1    = apply(cov_res_lm, 2, sd),
#  Lasso  = colMeans(cov_res_lasso),
#  Lasso1 = apply(cov_res_lasso, 2, sd),
#  RF     = colMeans(cov_res_rf),
#  RF1    = apply(cov_res_rf, 2, sd)
#)

df_long <- pivot_longer(df_cov, cols = c(LM, Lasso, RF), names_to = "Model", values_to = "Value")
df_long$method <- factor(df_long$method, levels = df_cov$method)
p3 <- ggplot(df_long, aes(x = method, y = Value, color = Model, shape = Model)) +
  geom_point(size = 3) +  
  theme_minimal() +
  labs(title = "Empirical Coverage", x = "", y = "Coverage", color = NULL, shape = NULL) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.5, 0.1),  legend.justification = "center", legend.background = element_rect(fill = "white")) +
  scale_color_manual(values = c("LM" = "seagreen3", "Lasso" = "darkblue", "RF" = "orange2" )) +  
  scale_shape_manual(values = c("Lasso" = 15, "RF" = 18, "LM" = 19)) + 
  ylim(0.7,1)

showtext_auto()
plot_news <-(p0 + p1 + p2)
ggsave("news_plot.pdf", plot_news, width = 9, height = 4.5)


xtable::xtable(t(df_cov), digits = 3)


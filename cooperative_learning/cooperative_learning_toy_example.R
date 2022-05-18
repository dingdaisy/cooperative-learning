library(caret)
library(glmnet)
source("cooperative_regression_function.R")

simulate_data <- function(n=2000, px=500, pz=500, p_imp=30, sigma=39, 
                          sy=1, sx=3, sz=3, u_std=1, factor_strength=7,
                          train_frac=0.1, nfolds = 10){
  
  # Simulate data based on the factor model
  x = matrix(rnorm(n*px), n, px)
  z = matrix(rnorm(n*pz), n, pz)
  U = matrix(rep(0, n*p_imp), n, p_imp)
  
  for (m in seq(p_imp)){
    u = rnorm(n, sd = u_std)
    x[, m] = x[, m] + sx*u
    z[, m] = z[, m] + sz*u
    U[, m] = U[, m] + sy*u
  }
  x = scale(x, T, F)
  z = scale(z, T, F)
  
  beta_U = c(rep(factor_strength, p_imp))
  mu_all = U %*% beta_U
  y = mu_all + sigma * rnorm(n) 
  
  snr = var(mu_all) / var(y-mu_all)
  cat("", fill=T)
  cat(c("snr =",snr),fill=T)
  cat("",fill=T)
  
  # Split training and test sets
  smp_size_train = floor(train_frac * nrow(x)) 
  train_ind = sort(sample(seq_len(nrow(x)), size = smp_size_train))
  test_ind = setdiff(seq_len(nrow(x)), train_ind)
  
  colnames(x) = seq(ncol(x))
  colnames(z) = seq(ncol(z))
  
  train_X_raw <- x[train_ind, ]
  test_X_raw <- x[test_ind, ]
  train_Z_raw <- z[train_ind, ]
  test_Z_raw <- z[test_ind, ]
  train_y <- y[train_ind, ]
  test_y <- y[test_ind, ]
  
  preprocess_values_train = preProcess(train_X_raw, method = c("center", "scale"))
  train_X = predict(preprocess_values_train, train_X_raw)
  test_X = predict(preprocess_values_train, test_X_raw)
  
  preprocess_values_train_Z = preProcess(train_Z_raw, method = c("center", "scale"))
  train_Z = predict(preprocess_values_train_Z, train_Z_raw)
  test_Z = predict(preprocess_values_train_Z, test_Z_raw)
  
  foldid = sample(rep_len(1:nfolds, dim(train_X)[1]))
  
  return(data=list(x=x, z=z, y=y, train_X=train_X, train_Z=train_Z, train_y=train_y,
                   test_X=test_X, test_Z=test_Z, test_y=test_y, foldid=foldid))
}

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

data = simulate_data()
alpha = 0.5
coop_fit = coop_cv_new(data$train_X,data$train_Z,data$train_y,alpha=alpha,
                       foldid=data$foldid,nfolds=max(data$foldid),
                       pf_values=rep(1,ncol(data$train_X)+ncol(data$train_Z)))
coop_pred_test = cbind(data$test_X, data$test_Z) %*% coop_fit$best_fit$beta + (coop_fit$best_fit$a0*2)
calc_mse(coop_pred_test, data$test_y)

early_fusion = cv.glmnet(cbind(data$train_X,data$train_Z), data$train_y, 
                         standardize = F, foldid = data$foldid)
early_pred_test = predict(early_fusion, cbind(data$test_X, data$test_Z), s="lambda.min")
calc_mse(early_pred_test, data$test_y)



library(caret)
library(glmnet)
library(multiview)

source("./../../cooperative_regression_function.R")
load('Data.rda')

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

y = DOS
prot = Proteomics
met = Metabolomics

G1_index = match(unique(Id), Id)
prot_G1 = prot[G1_index,]
met_G1 = met[G1_index,]

x_fil = as.matrix(prot_G1)
z_fil = as.matrix(met_G1)
y = DOS[G1_index]

fit_mode = "min"
simN = 10
nfolds = 10
train_frac = 0.75
sim_seed = 321
val_frac = 0.4
set.seed(sim_seed)

alphalist = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,3)

err_train_coop = err_test_coop = support_coop = 
  diff_err = matrix(NA, simN, length(alphalist))
new_train_coop = new_cv_coop = new_test_coop = new_support_coop = matrix(NA, simN, length(alphalist))
new_train_coop_no_pf = new_support_coop_no_pf = new_test_coop_no_pf = new_cv_coop_no_pf = matrix(NA, simN, length(alphalist))

new_err_test_coop_validation = new_coop_support_validation = new_alpha_select = rep(NA, simN)

err_train_lasso  = err_test_lasso = lasso_support = 
  X_err_train_lasso = X_err_test_lasso = X_lasso_support = 
  Z_err_train_lasso = Z_err_test_lasso = Z_lasso_support = rep(NA, simN)

coop_selected_by_cv = support_by_cv = alpha_by_cv = rep(NA, simN)
coop_selected_by_cv_no_pf = support_by_cv_no_pf = alpha_by_cv_no_pf = rep(NA, simN)
err_null = rep(NA, simN)
err_fuse = rep(NA, simN)
support_fuse_late = rep(NA, simN)

for (ii in 1:simN){
  cat(ii)
  
  smp_size_train = floor(train_frac * nrow(x_fil)) 
  train_ind = sort(sample(seq_len(nrow(x_fil)), size = smp_size_train))
  test_ind = setdiff(seq_len(nrow(x_fil)), train_ind)
  
  train_X_raw <- x_fil[train_ind, ]
  test_X_raw <- x_fil[test_ind, ]
  train_Z_raw <- z_fil[train_ind, ]
  test_Z_raw <- z_fil[test_ind, ]
  
  # Filter by variance
  X_ind = which(log(apply(train_X_raw,2,var)) > 12)
  Z_ind = which(log(apply(train_Z_raw,2,var)) > 25)
  train_X_raw = train_X_raw[,X_ind]
  train_Z_raw = train_Z_raw[,Z_ind]
  test_X_raw = test_X_raw[,X_ind]
  test_Z_raw = test_Z_raw[,Z_ind]
  
  preprocess_values_train = preProcess(train_X_raw, method = c("center", "scale"))
  train_X = predict(preprocess_values_train, train_X_raw)
  test_X = predict(preprocess_values_train, test_X_raw)
  
  preprocess_values_train_Z = preProcess(train_Z_raw, method = c("center", "scale"))
  train_Z = predict(preprocess_values_train_Z, train_Z_raw)
  test_Z = predict(preprocess_values_train_Z, test_Z_raw)
  
  train_y <- y[train_ind]
  test_y <- y[test_ind]
  
  foldid = sample(rep_len(1:nfolds, dim(train_X)[1]))

  # Null model
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  # Separate models
  print("Only X")
  X_lasso_fit = cv.glmnet(train_X, train_y, standardize = F, foldid=foldid)
  
  if (fit_mode == "min"){
    X_yhat_lasso_train = predict(X_lasso_fit, train_X, s = "lambda.min")
    X_yhat_lasso_test = predict(X_lasso_fit, test_X, s = "lambda.min")
    index_X = which(X_lasso_fit$lambda == X_lasso_fit$lambda.min)
  } else if (fit_mode == "1se"){
    X_yhat_lasso_train = predict(X_lasso_fit, train_X, s = "lambda.1se")
    X_yhat_lasso_test = predict(X_lasso_fit, test_X, s = "lambda.1se")
    index_X = which(X_lasso_fit$lambda == X_lasso_fit$lambda.1se)
  }
  
  X_train_e = calc_mse(X_yhat_lasso_train, train_y)
  X_err_train_lasso[ii] = X_train_e

  X_test_e = calc_mse(X_yhat_lasso_test, test_y)
  X_err_test_lasso[ii] = X_test_e
  
  X_lasso_support[ii] = X_lasso_fit$nzero[index_X] 
  print(X_test_e)
  
  print("Only Z")
  Z_lasso_fit = cv.glmnet(train_Z, train_y, standardize = F, foldid=foldid)
  
  if (fit_mode == "min"){
    Z_yhat_lasso_train = predict(Z_lasso_fit, train_Z, s = "lambda.min")
    Z_yhat_lasso_test = predict(Z_lasso_fit, test_Z, s = "lambda.min")
    index_Z = which(Z_lasso_fit$lambda == Z_lasso_fit$lambda.min)
  } else if (fit_mode == "1se"){
    Z_yhat_lasso_train = predict(Z_lasso_fit, train_Z, s = "lambda.1se")
    Z_yhat_lasso_test = predict(Z_lasso_fit, test_Z, s = "lambda.1se")
    index_Z = which(Z_lasso_fit$lambda == Z_lasso_fit$lambda.1se)
  }
  
  Z_train_e = calc_mse(Z_yhat_lasso_train, train_y)
  Z_err_train_lasso[ii] = Z_train_e
  Z_test_e = calc_mse(Z_yhat_lasso_test, test_y)
  Z_err_test_lasso[ii] = Z_test_e
  Z_lasso_support[ii] = Z_lasso_fit$nzero[index_Z] 
  print(Z_test_e)
  
  # Late fusion
  print("Late Fusion")
  second_stage_smp = floor(val_frac * nrow(train_X))
  val_ind = sort(sample(seq_len(nrow(train_X)), size = second_stage_smp))
  train_late_ind = setdiff(seq_len(nrow(train_X)), val_ind)
  val_X = train_X[val_ind,]
  train_X_late = train_X[train_late_ind,]
  val_Z = train_Z[val_ind,]
  train_Z_late = train_Z[train_late_ind,]
  val_y = train_y[val_ind]
  train_y_late = train_y[train_late_ind]
  X_lasso_fit_late = cv.glmnet(train_X_late, train_y_late, standardize = F)
  Z_lasso_fit_late = cv.glmnet(train_Z_late, train_y_late, standardize = F)
  
  if (fit_mode == "min"){
    X_yhat_lasso_late_val = predict(X_lasso_fit_late, val_X, s = "lambda.min")
    X_yhat_lasso_late_test = predict(X_lasso_fit_late, test_X, s = "lambda.min")
    Z_yhat_lasso_late_val = predict(Z_lasso_fit_late, val_Z, s = "lambda.min")
    Z_yhat_lasso_late_test = predict(Z_lasso_fit_late, test_Z, s = "lambda.min")
    late_index_X = which(X_lasso_fit_late$lambda == X_lasso_fit_late$lambda.min)
    late_index_Z = which(Z_lasso_fit_late$lambda == Z_lasso_fit_late$lambda.min)
  } else if (fit_mode == "1se"){
    X_yhat_lasso_late_val = predict(X_lasso_fit_late, val_X, s = "lambda.1se")
    X_yhat_lasso_late_test = predict(X_lasso_fit_late, test_X, s = "lambda.1se")
    Z_yhat_lasso_late_val = predict(Z_lasso_fit_late, val_Z, s = "lambda.1se")
    Z_yhat_lasso_late_test = predict(Z_lasso_fit_late, test_Z, s = "lambda.1se")
    late_index_X = which(X_lasso_fit_late$lambda == X_lasso_fit_late$lambda.1se)
    late_index_Z = which(Z_lasso_fit_late$lambda == Z_lasso_fit_late$lambda.1se)
  }
  
  fuse_data = data.frame(y=val_y, X_pred=as.vector(X_yhat_lasso_late_val), Z_pred=as.vector(Z_yhat_lasso_late_val))
  fit_fuse = lm(y ~ X_pred + Z_pred, data=fuse_data)
  fuse_pred_test = predict(fit_fuse, data.frame(X_pred=as.vector(X_yhat_lasso_late_test), 
                                                Z_pred=as.vector(Z_yhat_lasso_late_test)))
  err_fuse[ii] = calc_mse(fuse_pred_test, test_y)
  print(err_fuse[ii])
  support_fuse_late[ii] = X_lasso_fit_late$nzero[late_index_X] + 
    Z_lasso_fit_late$nzero[late_index_Z]
  
  # Early fusion
  print("Early Fusion")
  fit_lasso_cv = cv.glmnet(cbind(train_X,train_Z), train_y, standardize = F, foldid = foldid)
  
  if (fit_mode == "min"){
    yhat_lasso_train = predict(fit_lasso_cv, cbind(train_X,train_Z), s="lambda.min")
    yhat_lasso_test = predict(fit_lasso_cv, cbind(test_X, test_Z), s="lambda.min")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.min)
  } else if (fit_mode == "1se"){
    yhat_lasso_train = predict(fit_lasso_cv, cbind(train_X,train_Z), s="lambda.1se")
    yhat_lasso_test = predict(fit_lasso_cv, cbind(test_X, test_Z), s="lambda.1se")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se)
  }
  
  train_e = calc_mse(yhat_lasso_train, train_y)
  err_train_lasso[ii] = train_e
  lasso_support[ii] = fit_lasso_cv$nzero[index_early]
  test_e = calc_mse(yhat_lasso_test, test_y)
  err_test_lasso[ii] = test_e
  print(test_e)
  print(lasso_support[ii])
  
  # Cooperative Regression
  print("Cooperative Regression")
  cvm_min = rep(0, times = length(alphalist))
  test_MSE_min = rep(0, times = length(alphalist))
  support_min = rep(0, times = length(alphalist))
  
  cvm_min_no_pf = rep(0, times = length(alphalist))
  test_MSE_min_no_pf = rep(0, times = length(alphalist))
  support_min_no_pf = rep(0, times = length(alphalist))
  
  for (j in 1:length(alphalist)){
    alpha=alphalist[j]
    print(alpha)
    
    # Cooperative Regression
    pf_null = rep(1, ncol(train_X) + ncol(train_Z))
    full_fit_no_pf = cv.multiview(list(x=train_X,z=train_Z),train_y,
                     rho=alpha,foldid=foldid, penalty.factor=pf_null)
    
    yhat_coop_new_no_pf = predict(full_fit_no_pf, list(x=train_X,z=train_Z), s="lambda.min")
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, train_y)
    
    new_cv_coop_no_pf[ii, j] = min(full_fit_no_pf$cvm)
    new_support_coop_no_pf[ii, j] = full_fit_no_pf$nzero[which.min(full_fit_no_pf$cvm)]
    yhat_coop_new_test_no_pf = predict(full_fit_no_pf, list(x=test_X,z=test_Z), s="lambda.min")
    test_e_coop_new_no_pf = calc_mse(yhat_coop_new_test_no_pf, test_y)
    new_test_coop_no_pf[ii,j] = test_e_coop_new_no_pf
    print("Cooperative Regression")
    print(new_test_coop_no_pf[ii,j])
    
    cvm_min_no_pf[j] = min(full_fit_no_pf$cvm)
    test_MSE_min_no_pf[j] = test_e_coop_new_no_pf
    support_min_no_pf[j] = new_support_coop_no_pf[ii, j]
  }
  
  coop_selected_by_cv_no_pf[ii] = test_MSE_min_no_pf[which.min(cvm_min_no_pf)]
  support_by_cv_no_pf[ii] = support_min_no_pf[which.min(cvm_min_no_pf)]
  alpha_by_cv_no_pf[ii] = alphalist[which.min(cvm_min_no_pf)]
  print("Cooperative regression selected by cv")
  print(coop_selected_by_cv_no_pf[ii])
}

test_err = cbind(X_err_test_lasso, Z_err_test_lasso, err_test_lasso, err_fuse, 
                 coop_selected_by_cv_no_pf, coop_selected_by_cv,
                 new_test_coop_no_pf, new_test_coop)

support_df = cbind(X_lasso_support, Z_lasso_support,
                   support_fuse_late, lasso_support,
                   support_by_cv_no_pf, support_by_cv, 
                   new_support_coop_no_pf, 
                   new_support_coop)

out = rbind(colMeans(test_err), apply(test_err,2,sd)/sqrt(10))
test_err_all_diff = test_err - test_err[,3]
out_diff = rbind(colMeans(test_err_all_diff), apply(test_err_all_diff,2,sd)/sqrt(10))
colMeans(support_df[,1:6])

saveRDS(test_err, paste0("labor_test_err.rds"))
saveRDS(support_df, paste0("labor_support.rds"))
saveRDS(alpha_by_cv_no_pf, paste0("alpha_no_pf.rds"))

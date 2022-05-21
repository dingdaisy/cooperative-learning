source("cooperative_learning_multi_view.R")

library(RColorBrewer)
library(caret)
library(latex2exp)

#Simulation
fit_mode_final = "comb"
n = 10000
px1 = 300 #dimensionality
px2 = 300
px3 = 300
p_imp = 60
sigma = 25 #level of snr, higher more noisy
sy = 2
sx1 = 2 #factor strength
sx2 = 2
sx3 = 2
u_std = 1 #std of u
factor_strength = 2
train_frac = 0.02
val_frac = 0
test_frac = 0.98
nfolds = 10
sim_seed = 33
fit_mode_coop_adap = 'min'
fit_mode_coop = 'min'
simN = 10
alphalist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)

mse_coop = err_coop_mu = err_coop_y = val_coop = support_coop = matrix(NA, simN, length(alphalist))
coop_cv_x = coop_cv_z = matrix(NA, simN, length(alphalist))
new_train_coop = new_test_coop = new_cv_coop = new_support_coop = matrix(NA, simN, length(alphalist))
new_train_coop_no_pf = new_cv_coop_no_pf = new_support_coop_no_pf = new_test_coop_no_pf = matrix(NA, simN, length(alphalist))

mse_lasso = train_lasso_cv = test_lasso_cv = err_lasso = val_lasso = lasso_support = lasso_support_cv = rep(NA, simN)
snr_list = rep(NA, simN)
err_null = rep(NA, simN)

train_lasso_cv_X1 = train_lasso_cv_X2 = train_lasso_cv_X3 = 
  lasso_support_cv_X1 = lasso_support_cv_X2 = lasso_support_cv_X3 =
  test_lasso_cv_X1 = test_lasso_cv_X2 = test_lasso_cv_X3 = rep(NA, simN)
err_fuse_late = support_fuse_late = rep(NA, simN)

coop_selected_by_cv = support_by_cv = alpha_by_cv = rep(NA, simN)
coop_selected_by_cv_no_pf = support_by_cv_no_pf = alpha_by_cv_no_pf = rep(NA, simN)
imin_alpha = imin_alpha_val = imin_alpha_test = rep(NA, length(alphalist))

fit_pf = list()
fit_no_pf = list()

set.seed(sim_seed)
for (ii in 1:simN){
  cat(ii)
  
  x1 = matrix(rnorm(n*px1), n, px1)
  x2 = matrix(rnorm(n*px2), n, px2)
  x3 = matrix(rnorm(n*px3), n, px3)
  U = matrix(rep(0, n*p_imp), n, p_imp)
  
  for (m in seq(p_imp)){
    u = rnorm(n, sd = u_std)
    x1[, m] = x1[, m] + sx1*u
    x2[, m] = x2[, m] + sx2*u
    x3[, m] = x3[, m] + sx3*u
    U[, m] = U[, m] + sy*u
  }
  x1 = scale(x1, T, F)
  x2 = scale(x2, T, F)
  x3 = scale(x3, T, F)
  
  #Only a subset of features are useful
  beta_U = c(rep(factor_strength, p_imp))
  mu_all = U %*% beta_U
  y = mu_all + sigma * rnorm(n)
  
  snr = var(mu_all) / var(y-mu_all)
  cat("", fill=T)
  cat(c("snr=",snr),fill=T)
  cat("",fill=T)
  
  #split training, validation and test sets
  smp_size_train = floor(train_frac * nrow(x1)) 
  smp_size_val = floor(val_frac * nrow(x1))
  train_ind = sort(sample(seq_len(nrow(x1)), size = smp_size_train))
  ind_no_train = setdiff(seq_len(nrow(x1)), train_ind)
  val_ind = sort(sample(ind_no_train, size = smp_size_val))
  test_ind = setdiff(ind_no_train, val_ind)
  
  colnames(x1) = seq(ncol(x1))
  colnames(x2) = seq(ncol(x2))
  colnames(x3) = seq(ncol(x3))
  
  train_x1_raw <- x1[train_ind, ]
  val_x1_raw <- x1[val_ind,]
  test_x1_raw <- x1[test_ind, ]
  train_x2_raw <- x2[train_ind, ]
  val_x2_raw <- x2[val_ind,]
  test_x2_raw <- x2[test_ind, ]
  train_x3_raw <- x3[train_ind, ]
  val_x3_raw <- x3[val_ind,]
  test_x3_raw <- x3[test_ind, ]
  train_y <- y[train_ind, ]
  val_y <- y[val_ind, ]
  test_y <- y[test_ind, ]
  
  preprocess_values_train_X1 = preProcess(train_x1_raw, method = c("center", "scale"))
  train_X1 = predict(preprocess_values_train_X1, train_x1_raw)
  test_X1 = predict(preprocess_values_train_X1, test_x1_raw)
  
  preprocess_values_train_X2 = preProcess(train_x2_raw, method = c("center", "scale"))
  train_X2 = predict(preprocess_values_train_X2, train_x2_raw)
  test_X2 = predict(preprocess_values_train_X2, test_x2_raw)
  
  preprocess_values_train_X3 = preProcess(train_x3_raw, method = c("center", "scale"))
  train_X3 = predict(preprocess_values_train_X3, train_x3_raw)
  test_X3 = predict(preprocess_values_train_X3, test_x3_raw)
  
  mu_train = mu_all[train_ind, ]
  mu_test= mu_all[test_ind, ]
  
  snr_list[ii] = snr
  foldid = sample(rep_len(1:nfolds, dim(train_X1)[1]))
  
  #Null
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  #Lasso
  print("Early Fusion")
  fit_mode = "min"
  fit_lasso_cv = cv.glmnet(cbind(train_X1,train_X2, train_X3), train_y, standardize = F, foldid = foldid)
  
  if (fit_mode == "min"){
    yhat_lasso_cv = predict(fit_lasso_cv, cbind(train_X1,train_X2, train_X3), s="lambda.min")
    yhat_lasso_cv_test = predict(fit_lasso_cv, cbind(test_X1, test_X2, test_X3), s="lambda.min")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.min)
  } else if (fit_mode == "1se"){
    yhat_lasso_cv = predict(fit_lasso_cv, cbind(train_X1,train_X2, train_X3), s="lambda.1se")
    yhat_lasso_cv_test = predict(fit_lasso_cv, cbind(test_X1, test_X2, test_X3), s="lambda.1se")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se)
  }
  
  train_e_cv = calc_mse(yhat_lasso_cv, mu_train)
  train_lasso_cv[ii] = train_e_cv
  lasso_support_cv[ii] = fit_lasso_cv$nzero[index_early]
  test_e_cv = calc_mse(yhat_lasso_cv_test, mu_test)
  test_lasso_cv[ii] = test_e_cv
  print(test_e_cv)
  
  #Lasso X1
  print("Separate X1")
  fit_lasso_cv_X1 = cv.glmnet(train_X1, train_y, standardize = F, foldid = foldid)
  
  if (fit_mode == "min"){
    yhat_lasso_cv_X1 = predict(fit_lasso_cv_X1, train_X1, s="lambda.min")
    yhat_lasso_cv_test_X1 = predict(fit_lasso_cv_X1, test_X1, s="lambda.min")
    index_X1 = which(fit_lasso_cv_X1$lambda == fit_lasso_cv_X1$lambda.min)
  } else if (fit_mode == "1se"){
    yhat_lasso_cv_X1 = predict(fit_lasso_cv_X1, train_X1, s="lambda.1se")
    yhat_lasso_cv_test_X1 = predict(fit_lasso_cv_X1, test_X1, s="lambda.1se")
    index_X1 = which(fit_lasso_cv_X1$lambda == fit_lasso_cv_X1$lambda.1se)
  }
  
  train_e_cv_X1 = calc_mse(yhat_lasso_cv_X1, mu_train)
  train_lasso_cv_X1[ii] = train_e_cv_X1
  lasso_support_cv_X1[ii] = fit_lasso_cv_X1$nzero[index_X1]
  test_e_cv_X1 = calc_mse(yhat_lasso_cv_test_X1, mu_test)
  test_lasso_cv_X1[ii] = test_e_cv_X1
  print(test_e_cv_X1)
  
  #Lasso X2
  print("Separate X2")
  fit_lasso_cv_X2 = cv.glmnet(train_X2, train_y, standardize = F, foldid = foldid)
  
  if (fit_mode == "min"){
    yhat_lasso_cv_X2 = predict(fit_lasso_cv_X2, train_X2, s="lambda.min")
    yhat_lasso_cv_test_X2 = predict(fit_lasso_cv_X2, test_X2, s="lambda.min")
    index_X2 = which(fit_lasso_cv_X2$lambda == fit_lasso_cv_X2$lambda.min)
  } else if (fit_mode == "1se"){
    yhat_lasso_cv_X2 = predict(fit_lasso_cv_X2, train_X2, s="lambda.1se")
    yhat_lasso_cv_test_X2 = predict(fit_lasso_cv_X2, test_X2, s="lambda.1se")
    index_X2 = which(fit_lasso_cv_X2$lambda == fit_lasso_cv_X2$lambda.1se)
  }
  
  train_e_cv_X2 = calc_mse(yhat_lasso_cv_X2, mu_train)
  train_lasso_cv_X2[ii] = train_e_cv_X2
  lasso_support_cv_X2[ii] = fit_lasso_cv_X2$nzero[index_X2]
  test_e_cv_X2 = calc_mse(yhat_lasso_cv_test_X2, mu_test)
  test_lasso_cv_X2[ii] = test_e_cv_X2
  print(test_e_cv_X2)
  
  #Lasso X3
  print("Separate X3")
  fit_lasso_cv_X3 = cv.glmnet(train_X3, train_y, standardize = F, foldid = foldid)
  
  if (fit_mode == "min"){
    yhat_lasso_cv_X3 = predict(fit_lasso_cv_X3, train_X3, s="lambda.min")
    yhat_lasso_cv_test_X3 = predict(fit_lasso_cv_X3, test_X3, s="lambda.min")
    index_X3 = which(fit_lasso_cv_X3$lambda == fit_lasso_cv_X3$lambda.min)
  } else if (fit_mode == "1se"){
    yhat_lasso_cv_X3 = predict(fit_lasso_cv_X3, train_X3, s="lambda.1se")
    yhat_lasso_cv_test_X3 = predict(fit_lasso_cv_X3, test_X3, s="lambda.1se")
    index_X3 = which(fit_lasso_cv_X3$lambda == fit_lasso_cv_X3$lambda.1se)
  }
  
  train_e_cv_X3 = calc_mse(yhat_lasso_cv_X3, mu_train)
  train_lasso_cv_X3[ii] = train_e_cv_X3
  lasso_support_cv_X3[ii] = fit_lasso_cv_X3$nzero[index_X3]
  test_e_cv_X3 = calc_mse(yhat_lasso_cv_test_X3, mu_test)
  test_lasso_cv_X3[ii] = test_e_cv_X3
  print(test_e_cv_X3)
  
  #Late fusion
  print("Late Fusion")
  val_frac = 0.3
  second_stage_smp = floor(val_frac * nrow(train_X1))
  val_ind = sort(sample(seq_len(nrow(train_X1)), size = second_stage_smp))
  train_late_ind = setdiff(seq_len(nrow(train_X1)), val_ind)
  val_X1 = train_X1[val_ind,]
  train_X1_late = train_X1[train_late_ind,]
  val_X2 = train_X2[val_ind,]
  train_X2_late = train_X2[train_late_ind,]
  val_X3 = train_X3[val_ind,]
  train_X3_late = train_X3[train_late_ind,]
  val_y = train_y[val_ind]
  train_y_late = train_y[train_late_ind]
  X1_lasso_fit_late = cv.glmnet(train_X1_late, train_y_late, standardize = F)
  X2_lasso_fit_late = cv.glmnet(train_X2_late, train_y_late, standardize = F)
  X3_lasso_fit_late = cv.glmnet(train_X3_late, train_y_late, standardize = F)
  
  if (fit_mode == "min"){
    X1_yhat_lasso_late_val = predict(X1_lasso_fit_late, val_X1, s = "lambda.min")
    X1_yhat_lasso_late_test = predict(X1_lasso_fit_late, test_X1, s = "lambda.min")
    X2_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.min")
    X2_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.min")
    X3_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.min")
    X3_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.min")
    late_index_X1 = which(X1_lasso_fit_late$lambda == X1_lasso_fit_late$lambda.min)
    late_index_X2 = which(X2_lasso_fit_late$lambda == X2_lasso_fit_late$lambda.min)
    late_index_X3 = which(X3_lasso_fit_late$lambda == X3_lasso_fit_late$lambda.min)
  } else if (fit_mode == "1se"){
    X1_yhat_lasso_late_val = predict(X1_lasso_fit_late, val_X1, s = "lambda.1se")
    X1_yhat_lasso_late_test = predict(X1_lasso_fit_late, test_X1, s = "lambda.1se")
    X2_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.1se")
    X2_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.1se")
    X3_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.1se")
    X3_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.1se")
    late_index_X1 = which(X1_lasso_fit_late$lambda == X1_lasso_fit_late$lambda.1se)
    late_index_X2 = which(X2_lasso_fit_late$lambda == X2_lasso_fit_late$lambda.1se)
    late_index_X3 = which(X3_lasso_fit_late$lambda == X3_lasso_fit_late$lambda.1se)
  }
  
  fuse_data = data.frame(y=val_y, X1_pred=as.vector(X1_yhat_lasso_late_val), 
                         X2_pred=as.vector(X2_yhat_lasso_late_val),
                         X3_pred=as.vector(X3_yhat_lasso_late_val))
  fit_fuse = lm(y ~ X1_pred + X2_pred + X3_pred, data=fuse_data)
  fuse_pred_test = predict(fit_fuse, data.frame(X1_pred=as.vector(X1_yhat_lasso_late_test), 
                                                X2_pred=as.vector(X2_yhat_lasso_late_test),
                                                X3_pred=as.vector(X3_yhat_lasso_late_test)))
  err_fuse_late[ii] = calc_mse(fuse_pred_test, mu_test)
  print(err_fuse_late[ii])
  
  support_fuse_late[ii] = X1_lasso_fit_late$nzero[late_index_X1] + 
    X2_lasso_fit_late$nzero[late_index_X2] + X3_lasso_fit_late$nzero[late_index_X3]
  
  #Cooperative Regression
  
  print("Cooperative Regression")
  cvm_min = rep(0, times = length(alphalist))
  test_MSE_min = rep(0, times = length(alphalist))
  support_min = rep(0, times = length(alphalist))
  
  cvm_min_no_pf = rep(0, times = length(alphalist))
  test_MSE_min_no_pf = rep(0, times = length(alphalist))
  support_min_no_pf = rep(0, times = length(alphalist))
  
  for(j in 1:length(alphalist)){
    alpha = alphalist[j]
    print(alpha)
    
    fit_mode = fit_mode_coop
    pf_null = rep(1, ncol(train_X1) + ncol(train_X2) + ncol(train_X3))
    full_fit_no_pf = coop_cv_multi(train_X1,train_X2,train_X3,train_y,
                                   alpha=alpha,foldid=foldid,
                                   nfolds=max(foldid),
                                   pf_values=pf_null,
                                   fit_mode = fit_mode)
    
    yhat_coop_new_no_pf = cbind(train_X1,train_X2,train_X3) %*% full_fit_no_pf$best_fit_coef + (full_fit_no_pf$best_fit_intercept)
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, mu_train)
    
    if (fit_mode == 'min'){
      new_cv_coop_no_pf[ii, j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_min]
      cvm_min_no_pf[j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_min]
      new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[full_fit_no_pf$ind_min]
    } else if (fit_mode == '1se'){
      new_cv_coop_no_pf[ii, j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_1se]
      cvm_min_no_pf[j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_1se]
      new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[full_fit_no_pf$ind_1se]
    }
    
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, mu_train)
    yhat_coop_new_test_no_pf = cbind(test_X1, test_X2, test_X3) %*% full_fit_no_pf$best_fit_coef + (full_fit_no_pf$best_fit_intercept)
    test_e_coop_new_no_pf = calc_mse(yhat_coop_new_test_no_pf, mu_test)
    new_test_coop_no_pf[ii,j] = test_e_coop_new_no_pf
    print("Full (no penalty factor)")
    print(new_test_coop_no_pf[ii,j])
    
    test_MSE_min_no_pf[j] = test_e_coop_new_no_pf
    support_min_no_pf[j] = new_support_coop_no_pf[ii, j]
    
    #Iterative Coop Learning
    fit_mode = fit_mode_coop_adap
    coop_fit_iter = coop_regression_iter_multi(train_X1,train_X2,train_X3,train_y,
                               alpha=alpha,foldid=foldid,
                               max_iteration=5,
                               thr=0.05)
    yhat_coop_iter = (train_X1%*%coop_fit_iter$thetax1 + coop_fit_iter$intercept_x1) +
      (train_X2%*%coop_fit_iter$thetax2 + coop_fit_iter$intercept_x2) +
      (train_X3%*%coop_fit_iter$thetax3 + coop_fit_iter$intercept_x3)
    support_coop[ii, j] = coop_fit_iter$n_thetax1 + coop_fit_iter$n_thetax2 + coop_fit_iter$n_thetax3
    
    yhat_coop_test = (test_X1%*%coop_fit_iter$thetax1 + coop_fit_iter$intercept_x1) +
      (test_X2%*%coop_fit_iter$thetax2 + coop_fit_iter$intercept_x2) +
      (test_X3%*%coop_fit_iter$thetax3 + coop_fit_iter$intercept_x3)
    test_e_coop = calc_mse(yhat_coop_test, mu_test)
    err_coop_mu[ii,j] = test_e_coop
    err_coop_y[ii,j] = calc_mse(yhat_coop_test, test_y)
    print("Iterative")
    print(err_coop_mu[ii,j])
    
    #Iterative Coop Learning
    lambda_x1 = coop_fit_iter$lam_x1
    lambda_x2 = coop_fit_iter$lam_x2
    lambda_x3 = coop_fit_iter$lam_x3

    nx1 = ncol(train_X1)
    nx2 = ncol(train_X2)
    nx3 = ncol(train_X3)
    
    adjust_factor = (nx1 * lambda_x1 + nx2 * lambda_x2 + + nx3 * lambda_x3) / (nx1 + nx2 + nx3)
    pfs = c(rep(lambda_x1, ncol(train_X1)),rep(lambda_x2, ncol(train_X2)),rep(lambda_x3, ncol(train_X3)))
    full_fit = coop_cv_multi(train_X1,train_X2,train_X3,train_y,
                           alpha=alpha,foldid=foldid,
                           nfolds=max(foldid),
                           pf_values=pfs,
                           fit_mode = fit_mode)
    
    yhat_coop_new = cbind(train_X1,train_X2,train_X3) %*% full_fit$best_fit_coef + (full_fit$best_fit_intercept)
    new_train_coop[ii,j] = calc_mse(yhat_coop_new, mu_train)
    
    if (fit_mode == 'min'){
      new_cv_coop[ii, j] = full_fit$cvm[full_fit$ind_min]
      cvm_min[j] = full_fit$cvm[full_fit$ind_min]
      new_support_coop[ii, j] = full_fit$support[full_fit$ind_min]
    } else if (fit_mode == '1se'){
      new_cv_coop[ii, j] = full_fit$cvm[full_fit$ind_1se]
      cvm_min[j] = full_fit$cvm[full_fit$ind_1se]
      new_support_coop[ii, j] = full_fit$support[full_fit$ind_1se]
    }
    
    new_train_coop[ii,j] = calc_mse(yhat_coop_new, mu_train)
    yhat_coop_new_test = cbind(test_X1, test_X2, test_X3) %*% full_fit$best_fit_coef + (full_fit$best_fit_intercept)
    test_e_coop_new = calc_mse(yhat_coop_new_test, mu_test)
    new_test_coop[ii,j] = test_e_coop_new
    print("Full (with penalty factor)")
    print(new_test_coop[ii,j])
    
    test_MSE_min[j] = test_e_coop_new
    support_min[j] = new_support_coop[ii, j]
  }
  
  coop_selected_by_cv_no_pf[ii] = test_MSE_min_no_pf[which.min(cvm_min_no_pf)]
  support_by_cv_no_pf[ii] = support_min_no_pf[which.min(cvm_min_no_pf)]
  alpha_by_cv_no_pf[ii] = alphalist[which.min(cvm_min_no_pf)]
  print("Full selected by cv, without pf")
  print(coop_selected_by_cv_no_pf[ii])
  print(alpha_by_cv_no_pf[ii])
  
  coop_selected_by_cv[ii] = test_MSE_min[which.min(cvm_min)]
  support_by_cv[ii] = support_min[which.min(cvm_min)]
  alpha_by_cv[ii] = alphalist[which.min(cvm_min)]
  print("Full selected by cv, with pf")
  print(coop_selected_by_cv[ii])
  print(alpha_by_cv[ii])
}

snr_avg = mean(snr_list)

sim1_filename = paste(sim_seed, paste0("pimp", p_imp),
                      paste0("px1", px1), paste0("px2", px2), paste0("px3", px3),
                      paste0("sigma", sigma),
                      paste0("factorstr", sy), 
                      paste0("factorx1_", sx1), 
                      paste0("factorx2_", sx2), 
                      paste0("factorx3_", sx3), 
                      paste0("ustd", u_std),  
                      paste0("factor", factor_strength), 
                      paste0("SNR", as.integer(snr_avg)),
                      paste0("mode_", fit_mode_final),
                      sep = "_")
mse_df_file = paste("test_mse", sim1_filename, sep = "_")
mse_diff_df_file = paste("diff_test_mse", sim1_filename, sep = "_")
support_df_file = paste("support", sim1_filename, sep = "_")
alpha_df_file = paste("alpha", sim1_filename, sep = "_") 
plot_name = paste("plot", sim1_filename, sep = "_", ".pdf")

test_mse_df <- data.frame(cbind(test_lasso_cv_X1, test_lasso_cv_X2, test_lasso_cv_X3, test_lasso_cv, 
                                err_fuse_late, coop_selected_by_cv_no_pf, coop_selected_by_cv, 
                                new_test_coop_no_pf, new_test_coop))
#test_mse_df = readRDS("test_mse_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

lasso_diff_df = data.frame(cbind(test_lasso_cv_X1, test_lasso_cv_X2, test_lasso_cv_X3, test_lasso_cv,
                                 err_fuse_late,
                                 coop_selected_by_cv_no_pf, 
                                 coop_selected_by_cv,
                                 new_test_coop_no_pf,
                                 new_test_coop) - test_lasso_cv)
#lasso_diff_df = readRDS("diff_test_mse_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

support_df = data.frame(cbind(lasso_support_cv_X1, lasso_support_cv_X2, lasso_support_cv_X3,
                              lasso_support_cv, support_fuse_late,
                              support_by_cv_no_pf, support_by_cv, 
                              new_support_coop_no_pf,
                              new_support_coop))
#support_df = readRDS("support_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

alpha_df = data.frame(cbind(alpha_by_cv, alpha_by_cv_no_pf))
#alpha_df = readRDS("alpha_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

alpha_axis = #as.character(round(alphalist,1)) 
  sapply(as.character(round(alphalist,1)),
         function(x) paste0("(alpha=",x,")"))

saveRDS(test_mse_df, paste0(mse_df_file, ".rds"))
saveRDS(lasso_diff_df, paste0(mse_diff_df_file, ".rds"))
saveRDS(support_df, paste0(support_df_file, ".rds"))
saveRDS(alpha_df, paste0(alpha_df_file, ".rds"))





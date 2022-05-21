source("../cooperative_regression_function.R", chdir=TRUE)

library(RColorBrewer)
library(caret)
library(latex2exp)

#Simulation
n = 10000
px = 500 #dimensionality
pz = 500
px_imp = 30 #level of sparsity
pz_imp = 30
p_imp = 30
sigma = 33 #level of snr, higher more noisy
sy = 1
sx = 2
sz = 2
u_std = 1 #std of u
factor_strength = 8
train_frac = 0.02
val_frac = 0
test_frac = 0.98
nfolds = 10
sim_seed = 33
simN = 10
alphalist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)

mse_coop = err_coop_mu = err_coop_y = val_coop = support_coop = support_coop_1se = matrix(NA, simN, length(alphalist))
coop_cv_x = coop_cv_z = matrix(NA, simN, length(alphalist))
new_train_coop = new_test_coop = new_cv_coop = new_support_coop = new_support_coop_1se = matrix(NA, simN, length(alphalist))
new_train_coop_no_pf = new_cv_coop_no_pf = new_support_coop_no_pf = new_support_coop_no_pf_1se = new_test_coop_no_pf = matrix(NA, simN, length(alphalist))

mse_lasso  = train_lasso_cv = test_lasso_cv = err_lasso = val_lasso = lasso_support = lasso_support_1se = lasso_support_cv = lasso_support_cv_1se =rep(NA, simN)
snr_list = rep(NA, simN)
err_null = rep(NA, simN)

train_lasso_cv_Z = train_lasso_cv_X = lasso_support_cv_Z = lasso_support_cv_Z_1se = lasso_support_cv_X = lasso_support_cv_X_1se = test_lasso_cv_Z = test_lasso_cv_X = rep(NA, simN)
err_fuse_late = support_fuse_late = support_fuse_late_1se = rep(NA, simN)

coop_selected_by_cv = support_by_cv = support_by_cv_1se = alpha_by_cv = rep(NA, simN)
coop_selected_by_cv_no_pf = support_by_cv_no_pf = alpha_by_cv_no_pf = support_by_cv_no_pf_1se = rep(NA, simN)
imin_alpha = imin_alpha_val = imin_alpha_test = rep(NA, length(alphalist))

fit_pf = list()
fit_no_pf = list()

set.seed(sim_seed)
for (ii in 1:simN){
  cat(ii)

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
  
  #Only a subset of features are useful
  beta_U = c(rep(factor_strength, p_imp))
  mu_all = U %*% beta_U
  y = mu_all + sigma * rnorm(n)
  
  snr = var(mu_all) / var(y-mu_all)
  cat("", fill=T)
  cat(c("snr=",snr),fill=T)
  cat("",fill=T)
  
  #split training, validation and test sets
  smp_size_train <- floor(train_frac * nrow(x)) 
  smp_size_val <- floor(val_frac * nrow(x))
  train_ind <- sort(sample(seq_len(nrow(x)), size = smp_size_train))
  ind_no_train = setdiff(seq_len(nrow(x)), train_ind)
  val_ind = sort(sample(ind_no_train, size = smp_size_val))
  test_ind = setdiff(ind_no_train, val_ind)
  
  colnames(x) = seq(ncol(x))
  colnames(z) = seq(ncol(z))
  
  train_X_raw <- x[train_ind, ]
  val_X_raw <- x[val_ind,]
  test_X_raw <- x[test_ind, ]
  train_Z_raw <- z[train_ind, ]
  val_Z_raw <- z[val_ind,]
  test_Z_raw <- z[test_ind, ]
  train_y <- y[train_ind, ]
  val_y <- y[val_ind, ]
  test_y <- y[test_ind, ]
  
  preprocess_values_train = preProcess(train_X_raw, method = c("center", "scale"))
  train_X = predict(preprocess_values_train, train_X_raw)
  test_X = predict(preprocess_values_train, test_X_raw)
  
  preprocess_values_train_Z = preProcess(train_Z_raw, method = c("center", "scale"))
  train_Z = predict(preprocess_values_train_Z, train_Z_raw)
  test_Z = predict(preprocess_values_train_Z, test_Z_raw)
  
  mu_train = mu_all[train_ind, ]
  mu_test= mu_all[test_ind, ]
  
  snr_list[ii] = snr
  foldid = sample(rep_len(1:nfolds, dim(train_X)[1]))
  
  #Null
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  #Lasso
  print("Early Fusion")
  fit_lasso_cv = cv.glmnet(cbind(train_X,train_Z), train_y, standardize = F, foldid = foldid)
  yhat_lasso_cv = predict(fit_lasso_cv, cbind(train_X,train_Z), s="lambda.min")
  train_e_cv = calc_mse(yhat_lasso_cv, mu_train)
  train_lasso_cv[ii] = train_e_cv
  imin = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.min)
  imin_1se = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se)
  lasso_support_cv[ii] = fit_lasso_cv$nzero[imin_1se]
  lasso_support_cv_1se[ii] = fit_lasso_cv$nzero[imin_1se]
  
  yhat_lasso_cv_test = predict(fit_lasso_cv, cbind(test_X, test_Z), s="lambda.min")
  test_e_cv = calc_mse(yhat_lasso_cv_test, mu_test)
  test_lasso_cv[ii] = test_e_cv
  print(test_e_cv)
  
  #Lasso X
  print("Separate X")
  fit_lasso_cv_X = cv.glmnet(train_X, train_y, standardize = F, foldid = foldid)
  yhat_lasso_cv_X = predict(fit_lasso_cv_X, train_X, s="lambda.min")
  train_e_cv_X = calc_mse(yhat_lasso_cv_X, mu_train)
  train_lasso_cv_X[ii] = train_e_cv_X
  imin_X = which(fit_lasso_cv_X$lambda == fit_lasso_cv_X$lambda.min)
  lasso_support_cv_X[ii] = fit_lasso_cv_X$nzero[imin_X]
  
  yhat_lasso_cv_test_X = predict(fit_lasso_cv_X, test_X, s="lambda.min")
  test_e_cv_X = calc_mse(yhat_lasso_cv_test_X, mu_test)
  test_lasso_cv_X[ii] = test_e_cv_X
  print(test_e_cv_X)
  
  #Lasso Z
  print("Separate Z")
  fit_lasso_cv_Z = cv.glmnet(train_Z, train_y, standardize = F, foldid = foldid)
  yhat_lasso_cv_Z = predict(fit_lasso_cv_Z, train_Z, s="lambda.min")
  train_e_cv_Z = calc_mse(yhat_lasso_cv_Z, mu_train)
  train_lasso_cv_Z[ii] = train_e_cv_Z
  imin_Z = which(fit_lasso_cv_Z$lambda == fit_lasso_cv_Z$lambda.min)
  lasso_support_cv_Z[ii] = fit_lasso_cv_Z$nzero[imin_Z]
  
  yhat_lasso_cv_test_Z = predict(fit_lasso_cv_Z, test_Z, s="lambda.min")
  test_e_cv_Z = calc_mse(yhat_lasso_cv_test_Z, mu_test)
  test_lasso_cv_Z[ii] = test_e_cv_Z
  print(test_e_cv_Z)
  
  #Late fusion
  print("Late Fusion")
  val_frac = 0.3
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
  
  X_yhat_lasso_late_val = predict(X_lasso_fit_late, val_X, s = "lambda.min")
  X_yhat_lasso_late_test = predict(X_lasso_fit_late, test_X, s = "lambda.min")
  Z_yhat_lasso_late_val = predict(Z_lasso_fit_late, val_Z, s = "lambda.min")
  Z_yhat_lasso_late_test = predict(Z_lasso_fit_late, test_Z, s = "lambda.min")
  
  fuse_data = data.frame(y=val_y, X_pred=as.vector(X_yhat_lasso_late_val), Z_pred=as.vector(Z_yhat_lasso_late_val))
  fit_fuse = lm(y ~ X_pred + Z_pred, data=fuse_data)
  fuse_pred_test = predict(fit_fuse, data.frame(X_pred=as.vector(X_yhat_lasso_late_test), 
                                                Z_pred=as.vector(Z_yhat_lasso_late_test)))
  err_fuse_late[ii] = calc_mse(fuse_pred_test, mu_test)
  print(err_fuse_late[ii])
  support_fuse_late[ii] = X_lasso_fit_late$nzero[which(X_lasso_fit_late$lambda == X_lasso_fit_late$lambda.min)] + 
    Z_lasso_fit_late$nzero[which(Z_lasso_fit_late$lambda == Z_lasso_fit_late$lambda.min)]
  
  #Cooperative Regression
  print("Cooperative Regression")
  cvm_min = rep(0, times = length(alphalist))
  test_MSE_min = rep(0, times = length(alphalist))
  support_min = rep(0, times = length(alphalist))
  
  cvm_min_no_pf = rep(0, times = length(alphalist))
  test_MSE_min_no_pf = rep(0, times = length(alphalist))
  support_min_no_pf = rep(0, times = length(alphalist))
  
  if (ii == 1){
    fit_null = coop_cv_new(train_X,train_Z,train_y,alpha=0,
                           foldid=foldid,nfolds=max(foldid),
                           pf_values=rep(1, ncol(train_X)+ncol(train_Z)))
  }
  
  for(j in 1:length(alphalist)){
    alpha=alphalist[j]
    print(alpha)
    
    coop_fit = coop_regression_iter_order(train_X,train_Z,train_y,alpha,
                                          foldid=foldid,
                                          max_iteration=5,
                                          thr=0.05)
    if (coop_fit$order == "xz"){
      yhat_coop = (train_X%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) +
        (train_Z%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      train_e_coop = calc_mse(yhat_coop, mu_train)
      mse_coop[ii,j] = train_e_coop
      support_coop[ii, j] = sum(coop_fit$fit$thetax != 0) + sum(coop_fit$fit$thetaz != 0)
      
      yhat_coop_test = (test_X%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) + 
        (test_Z%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      test_e_coop = calc_mse(yhat_coop_test, mu_test)
      err_coop_mu[ii,j] = test_e_coop
      err_coop_y[ii,j] = calc_mse(yhat_coop_test, test_y)
      print("Iterative")
      print(err_coop_mu[ii,j])
      
      #New Cooperative Regression
      lambda_x = coop_fit$fit$lam_x
      lambda_z = coop_fit$fit$lam_z
      nx = ncol(train_X)
      nz = ncol(train_Z)
      adjust_factor = (nx * lambda_x + nz * lambda_z) / (nx + nz)
      pfs = c(rep(lambda_x, ncol(train_X)),rep(lambda_z, ncol(train_Z)))
      full_fit = coop_cv_new(train_X,train_Z,train_y,
                             alpha=alpha,foldid=foldid,
                             nfolds=max(foldid),
                             pf_values=pfs)
      if (ii == 1){
        fit_pf[[j]] = full_fit}
      
      yhat_coop_new = cbind(train_X,train_Z) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      new_train_coop[ii,j] = calc_mse(yhat_coop_new, mu_train)
      
      new_cv_coop[ii, j] = min(full_fit$cvm)
      new_support_coop[ii, j] = full_fit$support[which.min(full_fit$cvm)]
      yhat_coop_new_test = cbind(test_X, test_Z) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      test_e_coop_new = calc_mse(yhat_coop_new_test, mu_test)
      new_test_coop[ii,j] = test_e_coop_new
      print("Full")
      print(new_test_coop[ii,j])
      
      cvm_min[j] = min(full_fit$cvm)
      test_MSE_min[j] = test_e_coop_new
      support_min[j] = new_support_coop[ii, j]
    }
    
    else if (coop_fit$order == "zx"){
      
      yhat_coop = (train_Z%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) +
        (train_X%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      train_e_coop = calc_mse(yhat_coop, mu_train)
      mse_coop[ii,j] = train_e_coop
      support_coop[ii, j] = sum(coop_fit$fit$thetax != 0) + sum(coop_fit$fit$thetaz != 0)
      
      yhat_coop_test = (test_Z%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) + 
        (test_X%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      test_e_coop = calc_mse(yhat_coop_test, mu_test)
      err_coop_mu[ii,j] = test_e_coop
      err_coop_y[ii,j] = calc_mse(yhat_coop_test, test_y)
      print("Iterative")
      print(err_coop_mu[ii,j])
      
      #New Cooperative Regression, with penalty factor
      lambda_z = coop_fit$fit$lam_x
      lambda_x = coop_fit$fit$lam_z
      
      nx = ncol(train_X)
      nz = ncol(train_Z)
      adjust_factor = (nx * lambda_x + nz * lambda_z) / (nx + nz)
      pfs = c(rep(lambda_z, ncol(train_Z)), rep(lambda_x, ncol(train_X)))
      full_fit = coop_cv_new(train_Z,train_X,train_y,
                             alpha=alpha,foldid=foldid,
                             nfolds=max(foldid),
                             pf_values=pfs)
      
      if (ii == 1){
        fit_pf[[j]] = full_fit}
      
      yhat_coop_new = cbind(train_Z,train_X) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      new_train_coop[ii,j] = calc_mse(yhat_coop_new, mu_train)
      
      new_cv_coop[ii, j] = min(full_fit$cvm)
      new_support_coop[ii, j] = full_fit$support[which.min(full_fit$cvm)]
      yhat_coop_new_test = cbind(test_Z,test_X) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      test_e_coop_new = calc_mse(yhat_coop_new_test, mu_test)
      new_test_coop[ii,j] = calc_mse(yhat_coop_new_test, mu_test)
      
      print("Full")
      print(new_test_coop[ii,j])
      cvm_min[j] = min(full_fit$cvm)
      test_MSE_min[j] = test_e_coop_new
      support_min[j] = new_support_coop[ii, j]
    }
    
    pf_null = rep(1, ncol(train_X) + ncol(train_Z))
    full_fit_no_pf = coop_cv_new(train_X,train_Z,train_y,
                                 alpha=alpha,foldid=foldid,
                                 nfolds=max(foldid),
                                 pf_values=pf_null)
    if (ii == 1){
      fit_no_pf[[j]] = full_fit_no_pf}
    
    yhat_coop_new_no_pf = cbind(train_X,train_Z) %*% full_fit_no_pf$best_fit$beta + (full_fit_no_pf$best_fit$a0*2)
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, mu_train)
    
    new_cv_coop_no_pf[ii, j] = min(full_fit_no_pf$cvm)
    new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[which.min(full_fit_no_pf$cvm)]
    yhat_coop_new_test_no_pf = cbind(test_X, test_Z) %*% full_fit_no_pf$best_fit$beta + (full_fit_no_pf$best_fit$a0*2)
    test_e_coop_new_no_pf = calc_mse(yhat_coop_new_test_no_pf, mu_test)
    new_test_coop_no_pf[ii,j] = test_e_coop_new_no_pf
    print("Full (no penalty factor)")
    print(new_test_coop_no_pf[ii,j])
    
    cvm_min_no_pf[j] = min(full_fit_no_pf$cvm)
    test_MSE_min_no_pf[j] = test_e_coop_new_no_pf
    support_min_no_pf[j] = new_support_coop_no_pf[ii, j]
  }
  
  coop_selected_by_cv[ii] = test_MSE_min[which.min(cvm_min)]
  support_by_cv[ii] = support_min[which.min(cvm_min)]
  alpha_by_cv[ii] = alphalist[which.min(cvm_min)]
  print("Full selected by cv, with pf")
  print(coop_selected_by_cv[ii])
  
  coop_selected_by_cv_no_pf[ii] = test_MSE_min_no_pf[which.min(cvm_min_no_pf)]
  support_by_cv_no_pf[ii] = support_min_no_pf[which.min(cvm_min_no_pf)]
  alpha_by_cv_no_pf[ii] = alphalist[which.min(cvm_min_no_pf)]
  print("Full selected by cv, without pf")
  print(coop_selected_by_cv_no_pf[ii])
}

snr_avg = mean(snr_list)

sim1_filename = paste(sim_seed, paste0("pimp", p_imp),
                      paste0("px", px), paste0("pz", pz),
                      paste0("sigma", sigma),
                      paste0("factorstr", sy), 
                      paste0("ustd", u_std),  
                      paste0("factor", factor_strength), 
                      paste0("SNR", as.integer(snr_avg)),
                      sep = "_")
mse_df_file = paste("test_mse", sim1_filename, sep = "_")
mse_diff_df_file = paste("diff_test_mse", sim1_filename, sep = "_")
support_df_file = paste("support", sim1_filename, sep = "_")
alpha_df_file = paste("alpha", sim1_filename, sep = "_") 
plot_name = paste("plot", sim1_filename, sep = "_", ".pdf")


test_mse_df <- data.frame(cbind(test_lasso_cv_X, test_lasso_cv_Z, test_lasso_cv, 
                                err_fuse_late, coop_selected_by_cv_no_pf,
                                new_test_coop_no_pf, coop_selected_by_cv, new_test_coop))
saveRDS(test_mse_df, paste0(mse_df_file, ".rds"))
#test_mse_df = readRDS("test_mse_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

lasso_diff_df = data.frame(cbind(test_lasso_cv_X, test_lasso_cv_Z, test_lasso_cv,
                                 err_fuse_late,
                                 coop_selected_by_cv_no_pf, 
                                 new_test_coop_no_pf,
                                 coop_selected_by_cv,
                                 new_test_coop) - test_lasso_cv)
saveRDS(lasso_diff_df, paste0(mse_diff_df_file, ".rds"))
#lasso_diff_df = readRDS("diff_test_mse_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

support_df = data.frame(cbind(lasso_support_cv_X, lasso_support_cv_Z,
                              lasso_support_cv, support_fuse_late,
                              support_by_cv_no_pf, new_support_coop_no_pf,
                              support_by_cv, new_support_coop))
saveRDS(support_df, paste0(support_df_file, ".rds"))
#support_df = readRDS("support_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

alpha_df = data.frame(cbind(alpha_by_cv, alpha_by_cv_no_pf))
saveRDS(alpha_df, paste0(alpha_df_file, ".rds"))
#alpha_df = readRDS("alpha_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1.rds")

alpha_axis = #as.character(round(alphalist,1)) 
  sapply(as.character(round(alphalist,1)),
         function(x) paste0("(alpha=",x,")"))

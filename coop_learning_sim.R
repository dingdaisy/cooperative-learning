source("cooperative_regression_function.R")
source("coop_learning_procedure.R")
library(RColorBrewer)
library(caret)
library(latex2exp)

#Simulation
sim_setting = "sim1"
fitMethod = "iter"
#setwd('./../12/')

#Sim 1
if (sim_setting == "sim1"){
  n = 10000
  px = 500 #dimensionality
  pz = 500
  px_imp = 30 #level of sparsity
  pz_imp = 30
  sigma = 35 #level of snr, higher more noisy
  s = 2 #factor strength
  u_std = 1 #std of u
  beta_strength_x = 2 #level of signal in x
  beta_strength_z = 2 #level of signal in z
  train_frac = 0.02
  val_frac = 0
  test_frac = 0.98
  nfolds = 10
}

simN = 10
alphalist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)

mse_coop = err_coop_mu = err_coop_y = val_coop = support_coop = matrix(NA, simN, length(alphalist))
coop_cv_x = coop_cv_z = matrix(NA, simN, length(alphalist))
new_train_coop = new_test_coop = new_cv_coop = new_support_coop = matrix(NA, simN, length(alphalist))
new_train_coop_no_pf = new_cv_coop_no_pf = new_support_coop_no_pf = new_test_coop_no_pf = matrix(NA, simN, length(alphalist))

mse_lasso  = train_lasso_cv = test_lasso_cv = err_lasso = val_lasso = lasso_support = lasso_support_cv = rep(NA, simN)
snr_list = rep(NA, simN)
err_null = rep(NA, simN)

train_lasso_cv_Z = train_lasso_cv_X = lasso_support_cv_Z = lasso_support_cv_X = test_lasso_cv_Z = test_lasso_cv_X = rep(NA, simN)
err_fuse_late = support_fuse_late = rep(NA, simN)

coop_selected_by_cv = support_by_cv = alpha_by_cv = rep(NA, simN)
coop_selected_by_cv_no_pf = support_by_cv_no_pf = alpha_by_cv_no_pf = rep(NA, simN)
imin_alpha = imin_alpha_val = imin_alpha_test = rep(NA, length(alphalist))

fit_pf = list()
fit_no_pf = list()

sim_seed = 3
set.seed(sim_seed)
for (ii in 1:simN){
  cat(ii)
  
  #Simulation 1
  if (sim_setting == "sim1"){
    x = matrix(rnorm(n*px), n, px)
    z = matrix(rnorm(n*pz), n, pz)
    
    for (m in seq(px_imp)){
      u = rnorm(n, sd = u_std)
      x[, m] = x[, m] + s*u
      z[, m] = z[, m] + s*u
    }
    x = scale(x, T, F)
    z = scale(z, T, F)
    
    #Only a subset of features are useful
    beta = c(rep(beta_strength_x,px_imp), rep(0,px-px_imp), 
             rep(beta_strength_z,pz_imp), rep(0,pz-pz_imp))
    mu_all = cbind(x,z) %*% beta
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
    
    #train_X = t(apply(train_X_raw,1,
    #            function(i)(i-colMeans(train_X_raw))/apply(train_X_raw, 2, sd)))
    #train_X0 = scale(train_X_raw, T, T)
    
    preprocess_values_train = preProcess(train_X_raw, method = c("center", "scale"))
    train_X = predict(preprocess_values_train, train_X_raw)
    test_X = predict(preprocess_values_train, test_X_raw)
    
    preprocess_values_train_Z = preProcess(train_Z_raw, method = c("center", "scale"))
    train_Z = predict(preprocess_values_train_Z, train_Z_raw)
    test_Z = predict(preprocess_values_train_Z, test_Z_raw)
    
    mu_train = mu_all[train_ind, ]
    mu_test= mu_all[test_ind, ]
  }
  
  snr_list[ii] = snr
  foldid = sample(rep_len(1:nfolds, dim(train_X)[1]))
  
  #Null
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  #Lasso
  print("Early Fusion")
  #cross-validation
  fit_lasso_cv = cv.glmnet(cbind(train_X,train_Z), train_y, standardize = F, foldid = foldid)
  yhat_lasso_cv = predict(fit_lasso_cv, cbind(train_X,train_Z), s="lambda.min")
  train_e_cv = calc_mse(yhat_lasso_cv, mu_train)
  train_lasso_cv[ii] = train_e_cv
  imin = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.min)
  lasso_support_cv[ii] = fit_lasso_cv$nzero[imin]
  
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
  #print(support_fuse_late[ii])
  
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
      print("xz")
      print(paste("lambda_x", lambda_x))
      print(paste("lambda_z", lambda_z))
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
      
      print("zx")
      print(paste("lambda_z", lambda_z))
      print(paste("lambda_x", lambda_x))
      
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

if (sim_setting == "sim1"){
  sim1_filename = paste(sim_seed, paste0("pximp", px_imp), paste0("pzimp", pz_imp),
                        paste0("px", px), paste0("pz", pz),
                        paste0("sigma", sigma),
                        paste0("factorstr", s), 
                        paste0("ustd", u_std),  
                        paste0("beta_xstr", beta_strength_x), 
                        paste0("betazstr", beta_strength_z), 
                        paste0("SNR", as.integer(snr_avg)),
                        #paste0("fitMethod", fitMethod),
                        sim_setting, sep = "_")
  mse_df_file = paste("test_mse", sim1_filename, sep = "_")
  mse_diff_df_file = paste("diff_test_mse", sim1_filename, sep = "_")
  support_df_file = paste("support", sim1_filename, sep = "_")
  alpha_df_file = paste("alpha", sim1_filename, sep = "_") 
  plot_name = paste("plot", sim1_filename, sep = "_", ".pdf")
}

if (sim_setting == "sim1"){
  if (s > 3) corr_level = "High"
  if (s > 0 && s <= 3) corr_level = "Medium"
  if (s == 0) corr_level = "No"
  
  if (beta_strength_x == beta_strength_z) contri = "Both X and Z contribute |"
  if (beta_strength_z == 0) contri = "Only X contributes |"
  if (beta_strength_z != 0 && beta_strength_x > beta_strength_z) contri = "X contributes more |"
  if (beta_strength_x != 0 && beta_strength_x < beta_strength_z) contri = "Z contributes more |"
}

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

#plot_name = "plot_3_pximp30_pzimp30_px500_pz500_sigma60_factorstr2_ustd1_beta_xstr2_betazstr2_SNR0_sim1_.pdf"
pdf(file=plot_name, width=18, height=5)
par(mfrow=c(1,4), las=2, cex.main=2.33, cex.axis = 2.15, cex.lab=1.95, mar=c(10.3,6.6,6.1,0.9))
p1 = boxplot(unlist(cbind(test_mse_df[,1:5], test_mse_df[,15]))
        ~sort(rep(1:6, simN)), 
        #ylim=c(3000,5000),
        names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "      Coop","Adap Coop"), 
        xlab="",ylab="", xaxt = "n")#, main="Test MSE")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Test MSE", line=1.8, cex.lab=2.75)

boxplot(unlist(cbind(lasso_diff_df[,1:5],lasso_diff_df[,15]))~
          sort(rep(1:6, simN)), 
        #ylim=c(-1300,800),
        names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "      Coop","Adap Coop"), 
        xlab="",ylab="", main="", xaxt = "n")
mtext("Relative to early fusion", side=3, line=0.6, las=1, cex=1.6)
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Test MSE Difference", line=3, cex.lab=2.75)
abline(h=0, lty=2)

boxplot(unlist(cbind(support_df[,1:5], support_df[,15]))~
          sort(rep(1:6, simN)),
        ylim=c(0,300),
        names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "      Coop","Adap Coop"), 
        xlab="", ylab="", main="", xaxt = "n")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Number of Features Selected", line=1.8, cex.lab=2.75)

counts = sapply(alphalist, function(x) sum(alpha_df$alpha_by_cv_no_pf == x))
counts_ada = sapply(alphalist, function(x) sum(alpha_df$alpha_by_cv == x))
barplot(rbind(counts,counts_ada),
        beside=TRUE,
        main="",
        #xlab="alpha", #ylab="Count",
        names=c(as.character(round(alphalist,1))),#as.character(round(alphalist,1))),
        ylim=c(0,10), cex.lab=2.0, 
        density=c(0,50),
        angle=c(0,36),
        col=c("#a4aaa4", "black"))#"#89c589"))
legend("topright", c("Coop","Adap Coop"), cex=2.15, 
       col=c("#a4aaa4", "black"), density=c(0,66),angle=c(0,36))
#title(main=paste0(TeX(sprintf("$\\alpha$")), "selected by CV (counts)"), line=1.8, cex.lab=2.75)
title(main=expression(bold(paste("Selected ", alpha, " by CV (Counts)"))), line=2.2, cex.lab=2.75)
mtext(TeX(sprintf("$\\alpha$")), side=1, las=1, line=5, cex=1.72)
dev.off()

pdf(file=paste0("path_",plot_name), width=10, height=5)
par(mfrow=c(1,2), las=2, cex.main=1.06, cex.axis = 0.8, mar=c(6,3.9,3.9,0.9))
cl = brewer.pal(length(alphalist), "Paired")

plot(0,0, xlim = range(log(fit_null$lambda)),
     ylim=range(fit_null$cvup,fit_null$cvlo),
     type = "n",
     xlab=expression(Log(lambda)),
     ylab="CV Error",
     main="Coop")

for(j in 1:length(alphalist)){
  object = fit_no_pf[[j]]
  lines(x=log(object$lambda), y=object$cvm, ylim=range(object$cvup,object$cvlo),
        pch=20, type="b", 
        cex=0.6,
        col=cl[j])
}

legend(legend=as.character(round(alphalist,1)), 
       x = "topleft",
       col=cl,
       pch=20,
       lty=rep(1,length(alphalist)), 
       ncol=1)

plot(0,0, xlim = range(log(fit_null$lambda)),
     ylim=range(fit_null$cvup,fit_null$cvlo),
     type = "n",
     xlab=expression(Log(lambda)),
     ylab="CV Error",
     main="Adaptive Coop")

for(j in 1:length(alphalist)){
  object = fit_pf[[j]]
  lines(x=log(object$lambda), y=object$cvm, ylim=range(object$cvup,object$cvlo),
        pch=20, type="b", 
        cex=0.6,
        col=cl[j])
}
dev.off()



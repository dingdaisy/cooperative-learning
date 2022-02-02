source("cooperative_regression_function.R")
source("coop_learning_procedure.R")
library(stringr)
library(plyr)
library(caret)

convert_stage = function(stage_str){
  if (stage_str == "Stage I" | stage_str == "Stage IA" | stage_str == "Stage IB")
    stage_int = 1
  else if (stage_str == "Stage II" | stage_str == "Stage IIA" | stage_str == "Stage IIB")
    stage_int  = 2
  else if (stage_str == "Stage III" | stage_str == "Stage IIIA" | stage_str == "Stage IIIB")
    stage_int  = 3
  else if (stage_str == "Stage IV" | stage_str == "Stage IVA" | stage_str == "Stage IVB")
    stage_int  = 4
  else
    stage_int = NA
  return(stage_int)
}

load("./tcga/TCGA_KIRC_MET.RData")
df = t(Dataset_TCGA_MET)
X_rna = df[,751:(751+19533-1)]
X_mic = df[,20284:21015]

stage = Clinical$ajcc_pathologic_tumor_stage
stage_int = sapply(seq(length(stage)), function(i) convert_stage(stage[i]))
fil_ind = which(is.na(stage_int))

x = X_rna[-fil_ind,]
z = X_mic[-fil_ind,]
y = stage_int[-fil_ind]

x_fil = as.matrix(x)
z_fil = as.matrix(z)

simN = 10
nfolds = 5
train_frac = 0.7
sim_seed = 1 
set.seed(sim_seed)

alphalist = c(0,0.2,0.4,0.6,0.8,1,2,3,9)

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

fit_pf = list()
fit_no_pf = list()

for (ii in 1:simN){
  cat(ii)
  
  smp_size_train = floor(train_frac * nrow(x)) 
  train_ind = sort(sample(seq_len(nrow(x)), size = smp_size_train))
  test_ind = setdiff(seq_len(nrow(x)), train_ind)
  
  train_X_raw <- x_fil[train_ind, ]
  test_X_raw <- x_fil[test_ind, ]
  train_Z_raw <- z_fil[train_ind, ]
  test_Z_raw <- z_fil[test_ind, ]
  
  X_ind = which(log(apply(train_X_raw,2,var)) > 0.9)
  Z_ind = which(log(apply(train_Z_raw,2,var)) > -3)
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

  #null
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  #separate model
  print("Only X")
  X_lasso_fit = cv.glmnet(train_X, train_y, standardize = F, foldid=foldid)
  X_yhat_lasso_train = predict(X_lasso_fit, train_X, s = "lambda.min")
  X_train_e = calc_mse(X_yhat_lasso_train, train_y)
  X_err_train_lasso[ii] = X_train_e
  
  X_yhat_lasso_test = predict(X_lasso_fit, test_X, s = "lambda.min")
  X_test_e = calc_mse(X_yhat_lasso_test, test_y)
  X_err_test_lasso[ii] = X_test_e
  
  imin = which(X_lasso_fit$lambda == X_lasso_fit$lambda.min)
  X_lasso_support[ii] = X_lasso_fit$nzero[imin] 
  print(X_test_e)
  
  print("Only Z")
  Z_lasso_fit = cv.glmnet(train_Z, train_y, standardize = F, foldid=foldid)
  Z_yhat_lasso_train = predict(Z_lasso_fit, train_Z, s = "lambda.min")
  Z_train_e = calc_mse(Z_yhat_lasso_train, train_y)
  Z_err_train_lasso[ii] = Z_train_e
  
  Z_yhat_lasso_test = predict(Z_lasso_fit, test_Z, s = "lambda.min")
  Z_test_e = calc_mse(Z_yhat_lasso_test, test_y)
  Z_err_test_lasso[ii] = Z_test_e
  imin = which(Z_lasso_fit$lambda == Z_lasso_fit$lambda.min)
  Z_lasso_support[ii] = Z_lasso_fit$nzero[imin] 
  print(Z_test_e)
  
  #late fusion
  print("Late Fusion")
  val_frac = 0.4
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
  err_fuse[ii] = calc_mse(fuse_pred_test, test_y)
  print(err_fuse[ii])
  support_fuse_late[ii] = X_lasso_fit_late$nzero[which(X_lasso_fit_late$lambda == X_lasso_fit_late$lambda.min)] + 
    Z_lasso_fit_late$nzero[which(Z_lasso_fit_late$lambda == Z_lasso_fit_late$lambda.min)]
  
  #lasso
  print("Lasso")
  lasso_fit = cv.glmnet(cbind(train_X,train_Z), train_y, standardize = F, foldid=foldid)
  yhat_lasso_train = predict(lasso_fit, cbind(train_X,train_Z), 
                             s = "lambda.min")
  train_e = calc_mse(yhat_lasso_train, train_y)
  err_train_lasso[ii] = train_e
  
  yhat_lasso_test = predict(lasso_fit, cbind(test_X, test_Z), s = "lambda.min")
  test_e = calc_mse(yhat_lasso_test, test_y)
  err_test_lasso[ii] = test_e
  imin = which(lasso_fit$lambda == lasso_fit$lambda.min)
  lasso_support[ii] = lasso_fit$nzero[imin] 
  print(test_e)
  
  #cooperative regression
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
  
  for (j in 1:length(alphalist)){
    alpha=alphalist[j]
    print(alpha)
    
    coop_fit = coop_regression_iter_order(train_X,train_Z,train_y,alpha,
                                          foldid=foldid,
                                          max_iteration=5,
                                          thr=0.05)
    if (coop_fit$order == "xz"){
      yhat_coop = (train_X%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) +
        (train_Z%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      train_e_coop = calc_mse(yhat_coop, train_y)
      err_train_coop[ii,j] = train_e_coop
      support_coop[ii, j] = sum(coop_fit$fit$thetax != 0) + sum(coop_fit$fit$thetaz != 0)
      
      yhat_coop_test = (test_X%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) + 
        (test_Z%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      err_test_coop[ii,j] = calc_mse(yhat_coop_test, test_y)
      print("Iterative")
      print(err_test_coop[ii,j])
      
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
      new_train_coop[ii,j] = calc_mse(yhat_coop_new, train_y)
      
      new_cv_coop[ii, j] = min(full_fit$cvm)
      new_support_coop[ii, j] = full_fit$support[which.min(full_fit$cvm)]
      yhat_coop_new_test = cbind(test_X, test_Z) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      test_e_coop_new = calc_mse(yhat_coop_new_test, test_y)
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
      train_e_coop = calc_mse(yhat_coop, train_y)
      err_train_coop[ii,j] = train_e_coop
      support_coop[ii, j] = sum(coop_fit$fit$thetax != 0) + sum(coop_fit$fit$thetaz != 0)
      
      yhat_coop_test = (test_Z%*%coop_fit$fit$thetax + coop_fit$fit$intercept_x) + 
        (test_X%*%coop_fit$fit$thetaz + coop_fit$fit$intercept_z)
      err_test_coop[ii,j] = calc_mse(yhat_coop_test, test_y)
      print("Iterative")
      print(err_test_coop[ii,j])
      
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
      new_train_coop[ii,j] = calc_mse(yhat_coop_new, train_y)
      
      new_cv_coop[ii, j] = min(full_fit$cvm)
      new_support_coop[ii, j] = full_fit$support[which.min(full_fit$cvm)]
      yhat_coop_new_test = cbind(test_Z,test_X) %*% full_fit$best_fit$beta + (full_fit$best_fit$a0*2)
      new_test_coop[ii,j] = calc_mse(yhat_coop_new_test, test_y)
      
      print("Full")
      print(new_test_coop[ii,j])
      cvm_min[j] = min(full_fit$cvm)
      test_MSE_min[j] = new_test_coop[ii,j]
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
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, train_y)
    
    new_cv_coop_no_pf[ii, j] = min(full_fit_no_pf$cvm)
    new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[which.min(full_fit_no_pf$cvm)]
    yhat_coop_new_test_no_pf = cbind(test_X, test_Z) %*% full_fit_no_pf$best_fit$beta + (full_fit_no_pf$best_fit$a0*2)
    test_e_coop_new_no_pf = calc_mse(yhat_coop_new_test_no_pf, test_y)
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

test_err = cbind(X_err_test_lasso, Z_err_test_lasso, err_test_lasso, err_fuse, 
                 coop_selected_by_cv_no_pf, coop_selected_by_cv,
                 new_test_coop_no_pf, new_test_coop)

support_df = cbind(X_lasso_support, Z_lasso_support,
                   support_fuse_late, lasso_support,
                   support_by_cv_no_pf, support_by_cv, 
                   new_support_coop_no_pf, 
                   new_support_coop)

p_name = paste("p = ", ncol(x_fil) + ncol(z_fil), "|")
n_name = paste("n = ", nrow(x))
plot_title = paste(p_name, n_name)
alpha_axis = as.character(round(alphalist,1))

pdf(file="./tcga_kirc_rna_mic.pdf", width=15, height=5)
par(mfrow=c(1,3), las=2, cex.main=2.1, cex.axis = 1.8, cex.lab=1.6, mar=c(12,3.2,3.2,0.9))
boxplot(as.vector(test_err[,1:6])~sort(rep(1:6, simN)), 
        names=c("Separate RNA", "Separate miRNA", 
                "Early", "Late",   
                "Coop", "Adap Coop"), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist))),
        xlab="", ylab="", main="Test MSE")

lasso_diff = cbind(test_err) - test_err[,3]
boxplot(as.vector(lasso_diff[,1:6])~sort(rep(1:6, simN)), 
        names=c("Separate RNA", "Separate miRNA", "Early",
                "Late",  "Coop", "Adap Coop"), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist))),
        xlab="", ylab="", main="Test MSE Difference with Early Fusion")
abline(h=0, lty=2)

boxplot(as.vector(support_df[,1:6])~
          sort(rep(1:6, simN)),
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist))),
        names=c("Separate RNA","Separate miRNA", "Early", "Late", "Coop", "Adap Coop"),
        xlab="", ylab="", main="Number of Features Selected")
dev.off()

saveRDS(test_err, paste0("tcga_test_err.rds"))
saveRDS(support_df, paste0("tcga_support.rds"))
saveRDS(alpha_by_cv_no_pf, paste0("alpha_no_pf.rds"))
saveRDS(alpha_by_cv, paste0("alpha_pf.rds"))

# counts = sapply(alphalist, function(x) sum(alpha_by_cv_no_pf == x))
# barplot(counts,
#         main="Selected Alpha Values by CV (Counts)",
#         xlab="", #ylab="Count",
#         names=as.character(round(alphalist,1)),
#         ylim=c(0,10), cex.lab=1.5)
# mtext("alpha", side=1, las=1, line=5, cex=1.6)
# 
# cl = #c("#ecffe8", "#ddf1d9", "#cfe2ca", "#c1d4bb", "#b3c7ac", "#a5b99e", "#97ab90",
#       # "#899e82", "#7c9174")
#   c(brewer.pal(length(alphalist), "RdGy"))
# cl[5] = "#edd3d3"
# plot(0,0, xlim = range(log(reg_path_all[[1]]$lambda)),
#      ylim=range(reg_path_all[[1]]$cvup,reg_path_all[[1]]$cvlo),
#      type = "n",
#      xlab="",
#      ylab="",
#      main="Regularization Curve: CV Error")
# 
# for(j in rev(1:length(alphalist))){
#     if (j != ind_best_full_data){
#       object = reg_path_all[[j]]
#       lines(x=log(object$lambda), y=object$cvm, 
#             ylim=range(object$cvup,object$cvlo),
#             pch=20, type="b", lty=1, lwd=0.5,
#             cex=0.6,
#             col=cl[j])} #alpha(cl[j],1.0))}
# }
# 
# object = reg_path_all[[ind_best_full_data]]
# lines(x=log(object$lambda), y=object$cvm, 
#       ylim=range(object$cvup,object$cvlo),
#       pch=18, type="b", lty=1, lwd=0.1,
#       cex=0.93,
#       col="#1c5c1c")
# axis(side=3,at=log(object$lambda),labels=paste(object$support),
#      tick=FALSE,line=-0.9, cex.axis=1.46, las=1)
# 
# mtext(expression(Log(lambda)), side=1, line=4.2, las=1, cex=1.5)
# abline(v=log(object$lambda.min), col="#1c5c1c", lty=2)
# legend(legend=as.character(round(alphalist,1)), 
#        x = "bottomleft",
#        col= c(cl[1:(ind_best_full_data-1)], "#1c5c1c", cl[ind_best_full_data:(length(alphalist))]),
#        pch= c(rep(20, ind_best_full_data-1), 18, rep(20,(length(alphalist)-ind_best_full_data))),
#        lty=rep(1,length(alphalist)), 
#        ncol=1)

# 
# pdf(file="path_tcga_kirc_rna_mic.pdf", width=5, height=5)
# par(mfrow=c(1,1), las=2, cex.main=1.06, cex.axis = 0.8, mar=c(5,3.9,3.9,0.9))
# cl = c(brewer.pal(length(alphalist), "RdBu"))
#   #brewer.pal(length(alphalist), "Paired")
# 
# plot(0,0, xlim = range(log(fit_null$lambda)),
#      ylim=range(fit_null$cvup,fit_null$cvlo),
#      type = "n",
#      xlab=expression(Log(lambda)),
#      ylab="CV Error",
#      main="Regularization Curve")
# 
# for(j in 1:length(alphalist)){
#   object = fit_no_pf[[j]]
#   lines(x=log(object$lambda), y=object$cvm, ylim=range(object$cvup,object$cvlo),
#         pch=20, type="b", lty=1, lwd=0.1,
#         cex=0.6,
#         col=cl[j])
# }
# 
# legend(legend=as.character(round(alphalist,1)), 
#        x = "topright",
#        col=cl,
#        pch=20,
#        lty=rep(1,length(alphalist)), 
#        ncol=1)
# dev.off()

##############exploratory##############

X_ind = which(log(apply(x_fil,2,var)) > 0.9)
Z_ind = which(log(apply(z_fil,2,var)) > -3)
X_df = x_fil[,X_ind]
Z_df = z_fil[,Z_ind]

preprocess_values_X = preProcess(X_df, method = c("center", "scale"))
X_std = predict(preprocess_values_X, X_df)

preprocess_values_Z = preProcess(Z_df, method = c("center", "scale"))
Z_std = predict(preprocess_values_Z, Z_df)

nfolds = 5
foldid = sample(rep_len(1:nfolds, dim(X_std)[1]))

alpha_selected = 3

pf_null = rep(1, ncol(X_std) + ncol(Z_std))
full_fit_no_pf = coop_cv_new(X_std,Z_std,y,
                             alpha=alpha_selected,foldid=foldid,
                             nfolds=max(foldid),
                             pf_values=pf_null)

x_fit = cv.glmnet(X_std, y, foldid=foldid)
z_fit = cv.glmnet(Z_std, y, foldid=foldid)

early_fit = cv.glmnet(cbind(X_std,Z_std),y,foldid=foldid)

coop_ind = which(coef(full_fit_no_pf$best_fit) != 0)
early_ind = which(coef(early_fit, s="lambda.min") != 0)

colnames(cbind(X_std,Z_std))[early_ind]
colnames(cbind(X_std,Z_std))[setdiff(coop_ind,early_ind)]



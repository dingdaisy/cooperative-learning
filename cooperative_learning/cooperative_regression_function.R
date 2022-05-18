library(glmnet)
library(randomForest)

# Cooperative Regression: Direct
coop_regression_full = function(x,z,y,alpha){
  n = length(y)
  xt = rbind(cbind(x, z),
             cbind(-sqrt(alpha)*x, sqrt(alpha)*z))
  yt = c(y, rep(0,n))
  g_fit = glmnet(xt, yt, standardize=F)
  return(g_fit)
}

#Stopping Criteria
calc_objective <- function(thetax,thetaz,x_intercept,z_intercept,alpha,x,z,y){
  res = sum((y-(x%*%thetax+x_intercept)-(z%*%thetaz+z_intercept))^2)/2 + 
    alpha*sum(((x%*%thetax+x_intercept)-(z%*%thetaz+z_intercept))^2)/2 
  return(res)
}

# Cooperative Regression: One-at-a-time
coop_regression_iter = function(x,z,y,alpha,
                                thetax=NULL, thetaz=NULL, 
                                foldid=NULL,
                                intercept_x=0, intercept_z=0,
                                max_iteration=3, thr=0.05){
  n = length(y)
  if(is.null(thetax)){thetax = rep(0,ncol(x))}
  if(is.null(thetaz)){thetaz = rep(0,ncol(z))}
  cv_error_x = NULL
  cv_error_z = NULL
  
  crit = NULL
  delta = 2 * thr
  crit_0 = calc_objective(thetax,thetaz,intercept_x,intercept_z,alpha,x,z,y)
  crit = c(crit, crit_0)
  n_iteration = 0
  
  lam_x = NULL
  lam_z = NULL
  
  while(delta > thr & n_iteration < max_iteration){
    old_crit = crit[length(crit)]
    n_iteration = n_iteration + 1
    
    y1=y/(1+alpha)- ((1-alpha)/(1+alpha))*(z%*%thetaz+intercept_z) 
    gfit1 = cv.glmnet(x, y1, standardize=F, foldid=foldid)
    thetax_temp=coef(gfit1,s=gfit1$lambda.min)[-1]
    intercept_x_temp=coef(gfit1,s=gfit1$lambda.min)[1]
    cv_error_x_temp=gfit1$cvm[which(gfit1$lambda == gfit1$lambda.min)]
    lam_x_temp = gfit1$lambda.min
    
    y2=y/(1+alpha)- ((1-alpha)/(1+alpha))*(x%*%thetax_temp+intercept_x_temp)
    gfit2 = cv.glmnet(z, y2, standardize=F, foldid=foldid)
    thetaz_temp=coef(gfit2,s=gfit2$lambda.min)[-1]
    intercept_z_temp=coef(gfit2,s=gfit2$lambda.min)[1]
    cv_error_z_temp=gfit2$cvm[which(gfit2$lambda == gfit2$lambda.min)]
    lam_z_temp = gfit2$lambda.min
    #print(cv_error_z_temp)
    
    current_obj = calc_objective(thetax_temp,thetaz_temp,
                                 intercept_x_temp,intercept_z_temp,
                                 alpha,x,z,y)
    
    if (current_obj < old_crit){
      thetax = thetax_temp
      thetaz = thetaz_temp
      intercept_z = intercept_z_temp
      intercept_x = intercept_x_temp
      cv_error_x = cv_error_x_temp
      cv_error_z = cv_error_z_temp
      lam_x = lam_x_temp
      lam_z = lam_z_temp
      crit = c(crit, current_obj)
    }
    
    new_crit = crit[length(crit)]
    delta = (old_crit-new_crit)/abs(old_crit)
  }
  
  return(list(thetax=thetax, thetaz=thetaz,
              intercept_x=intercept_x, intercept_z=intercept_z,
              cv_error_x=cv_error_x, cv_error_z=cv_error_z,
              lam_x=lam_x,lam_z=lam_z,
              alpha=alpha, crit=crit,
              n_iteration = n_iteration,
              n_thetax = sum(thetax!=0),
              n_thetaz = sum(thetax!=0)))
}

# Stopping Criteria: General
calc_objective_general <- function(X_fit,Z_fit,alpha,y){
  res = sum((y-X_fit-Z_fit)^2)/2 + alpha*sum((X_fit-Z_fit)^2)/2
  return(res)
}

# Cooperative Regression: General
coop_regression_general = function(x,z,y,alpha,
                            fitter = "gbm",#"gbm",
                            max_iteration=10, thr=0.05){
  n = length(y)
  X_fit = rep(0,n)
  Z_fit = rep(0,n)
  
  crit = NULL
  delta = 2 * thr
  crit_0 = calc_objective_general(X_fit,Z_fit,alpha,y)
  crit = c(crit, crit_0)
  n_iteration = 0
  
  while(delta > thr & n_iteration < max_iteration){
    old_crit = crit[length(crit)]
    n_iteration = n_iteration + 1
    
    y1=y/(1+alpha)- ((1-alpha)/(1+alpha))*Z_fit
    if (fitter == "rf") 
      gfit1_temp = randomForest(x,as.vector(y1))
    if (fitter == "gbm") 
      #gfit1_temp = cv.gbmfit(x,as.vector(y1),distribution="gaussian") #,n.trees=n.trees,shrinkage=shrinkage)
      gfit1_temp = gbm.fit(x,as.vector(y1),distribution="gaussian", interaction.depth=1,
                   shrinkage=0.05,n.trees=30, verbose=F)
    
    if (fitter == "rf")  
      X_fit_temp = predict(gfit1_temp, x)
    if (fitter == "gbm")  
      #X_fit_temp = predict.cv.gbmfit(gfit1_temp, x, s="ntrees.min")
      X_fit_temp = predict(gfit1_temp, x)
    
    y2=y/(1+alpha)- ((1-alpha)/(1+alpha))*X_fit_temp
    if (fitter == "rf") 
      gfit2_temp = randomForest(z,as.vector(y2))
    if (fitter == "gbm") 
      #gfit2_temp = cv.gbmfit(z,as.vector(y2),distribution="gaussian") #,n.trees=n.trees,shrinkage=shrinkage)
      gfit2_temp = gbm.fit(z,as.vector(y2),distribution="gaussian", interaction.depth=1,
                           shrinkage=0.05,n.trees=30, verbose=F)
      
    if (fitter == "rf")   
      Z_fit_temp = predict(gfit2_temp, z)
    if (fitter == "gbm")   
      #Z_fit_temp = predict.cv.gbmfit(gfit2_temp, z, s="ntrees.min")
      Z_fit_temp = predict(gfit2_temp, z)
    
    current_obj = calc_objective_general(X_fit_temp,Z_fit_temp,alpha,y)
    if (current_obj < old_crit){
      gfit1 = gfit1_temp
      gfit2 = gfit2_temp
      X_fit = X_fit_temp
      Z_fit = Z_fit_temp
      crit = c(crit, current_obj)
    }
    
    new_crit = crit[length(crit)]
    delta = (old_crit-new_crit)/abs(old_crit)
  }
  
  return(list(fX = gfit1, fZ = gfit2,
              alpha=alpha, crit=crit,
              n_iteration = n_iteration))
}

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

# Cooperative Regression: Adaptive, Determining Fitting Order
coop_regression_iter_order <- function(x,z,y,alpha,
                                       thetax=NULL, thetaz=NULL, 
                                       foldid=NULL,
                                       intercept_x=0, intercept_z=0,
                                       max_iteration=3, thr=0.05){
  
  xz = coop_regression_iter(x,z,y,alpha,
                            thetax=thetax, thetaz=thetaz, 
                            foldid=foldid,
                            intercept_x=intercept_x, intercept_z=intercept_z,
                            max_iteration=max_iteration, thr=thr)
  xz_error = xz$cv_error_x + xz$cv_error_z
  zx = coop_regression_iter(z,x,y,alpha,
                            thetax=thetax, thetaz=thetaz, 
                            foldid=foldid,
                            intercept_x=intercept_x, intercept_z=intercept_z,
                            max_iteration=max_iteration, thr=thr)
  zx_error = zx$cv_error_x + zx$cv_error_z
  
  if(xz_error < zx_error){
    fit = xz
    order = "xz"
  }else{
    fit = zx
    order = "zx"
  }
  return(list(fit=fit, order=order))
}

# Cooperative Regression: Cross Validation
coop_cv_new = function(x,z,y,alpha=0,foldid,nfolds=10,pf_values=NULL,n_iter=3){
  
  #get lambda sequence
  xt0 = rbind(cbind(x, z),
              cbind(-sqrt(alpha)*x, sqrt(alpha)*z))
  yt0 = c(y, rep(0, dim(x)[1]))
  coop0 = glmnet(xt0, yt0, standardize=F, penalty.factor = pf_values)
  #lasso0 = glmnet(cbind(x,z),y,standardize=F,lambda=coop0$lambda*2,penalty.factor = pf_values)
  #sum(abs(coop0$beta - lasso0$beta))
  lambda0 = coop0$lambda
  
  outlist = as.list(seq(nfolds))
  err_mat = matrix(NA, nfolds, length(lambda0))
  for (i in seq(nfolds)){
    which = (foldid == i)
    x_sub = x[!which, , drop = FALSE]
    z_sub = z[!which, , drop = FALSE]
    y_sub = y[!which]
    
    #centering inside seems necessary?
    x_sub = scale(x_sub, T, F)
    z_sub = scale(z_sub, T, F)
    
    #fit model
    xt = rbind(cbind(x_sub, z_sub),
               cbind(-sqrt(alpha)*x_sub, sqrt(alpha)*z_sub))
    N = dim(x_sub)[1]
    yt = c(y_sub, rep(0, N))
    #coop lambda = 1/2 lasso lambda
    #coop intercept = 1/2 lasso intercept
    coop = glmnet(xt, yt, standardize=F, lambda=lambda0, penalty.factor = pf_values) #relaxed TRUE??
    outlist[[i]] = coop
    #lasso = glmnet(cbind(x_sub,z_sub),y_sub,standardize=F, lambda=lambda0*2) #coop$lambda*2)
    
    #evaluate
    theta = coef(coop)[-1,] #matrix of coefficients, go through thetas, px100, column by column, find the support, LS estimates, new column of theta,
    #replace, prediction different, better fit with fewer features, by de-shrinking achieves better performance
    #coef, gamma 0, make sure the intercept is correct 
    
    #lasso_theta = coef(lasso)[-1,]
    #sum(abs(theta-lasso_theta))
    
    intercept = coef(coop)[1,] * 2
    pred_fold = cbind(x[which, , drop = FALSE], z[which, , drop = FALSE]) %*% theta + 
      rep(intercept, each=nrow(x[which, , drop = FALSE]))
    #lasso_pred = predict(lasso, cbind(x[which, , drop = FALSE], z[which, , drop = FALSE]))
    #sum(abs(pred_fold-lasso_pred))
    true_y = y[which]
    err_mse = (pred_fold - replicate(ncol(pred_fold), true_y))^2
    err_cv = apply(err_mse, 2, mean, na.rm = TRUE)
    err_mat[i, ] = err_cv
  }
  
  cvm = apply(err_mat, 2, mean, na.rm = TRUE)
  #cvsd = apply(err_mat, 2, sd)
  cvsd = sqrt(apply(scale(err_mat, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nfolds - 1))
  
  cvm_min_ind = which.min(cvm)
  #will remove repetitive fit later
  best_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda0[cvm_min_ind], penalty.factor = pf_values)
  full_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda0, penalty.factor = pf_values)
  n_nonzero = full_fit$df
  
  return(list(cvm=cvm, cvsd=cvsd, lambda=lambda0, lambda.min=lambda0[cvm_min_ind], best_fit = best_fit, 
              support=n_nonzero,cvup=cvm+cvsd, cvlow=cvm-cvsd, pf_values=pf_values))
}


# Cooperative Regression: Cross Validation with Both "min" and "1se" Fitting Mode
coop_cv_new_1se = function(x,z,y,alpha=0,foldid,nfolds=10,pf_values=NULL,n_iter=3,fit_mode="1se"){
  
  #get lambda sequence
  xt0 = rbind(cbind(x, z),
              cbind(-sqrt(alpha)*x, sqrt(alpha)*z))
  yt0 = c(y, rep(0, dim(x)[1]))
  coop0 = glmnet(xt0, yt0, standardize=F, penalty.factor = pf_values)
  #lasso0 = glmnet(cbind(x,z),y,standardize=F,lambda=coop0$lambda*2,penalty.factor = pf_values)
  #sum(abs(coop0$beta - lasso0$beta))
  lambda0 = coop0$lambda
  
  outlist = as.list(seq(nfolds))
  err_mat = matrix(NA, nfolds, length(lambda0))
  for (i in seq(nfolds)){
    which = (foldid == i)
    x_sub = x[!which, , drop = FALSE]
    z_sub = z[!which, , drop = FALSE]
    y_sub = y[!which]
    
    #centering inside is necessary
    x_sub = scale(x_sub, T, F)
    z_sub = scale(z_sub, T, F)
    
    #fit model
    xt = rbind(cbind(x_sub, z_sub),
               cbind(-sqrt(alpha)*x_sub, sqrt(alpha)*z_sub))
    N = dim(x_sub)[1]
    yt = c(y_sub, rep(0, N))
    #coop lambda = 1/2 lasso lambda
    #coop intercept = 1/2 lasso intercept
    coop = glmnet(xt, yt, standardize=F, lambda=lambda0, penalty.factor = pf_values) #relaxed TRUE??
    outlist[[i]] = coop
    #lasso = glmnet(cbind(x_sub,z_sub),y_sub,standardize=F, lambda=lambda0*2) #coop$lambda*2)
    
    #evaluate
    theta = coef(coop)[-1,] 
    
    #lasso_theta = coef(lasso)[-1,]
    #sum(abs(theta-lasso_theta))
    
    intercept = coef(coop)[1,] * 2
    pred_fold = cbind(x[which, , drop = FALSE], z[which, , drop = FALSE]) %*% theta + 
      rep(intercept, each=nrow(x[which, , drop = FALSE]))
    #lasso_pred = predict(lasso, cbind(x[which, , drop = FALSE], z[which, , drop = FALSE]))
    #sum(abs(pred_fold-lasso_pred))
    true_y = y[which]
    err_mse = (pred_fold - replicate(ncol(pred_fold), true_y))^2
    err_cv = apply(err_mse, 2, mean, na.rm = TRUE)
    err_mat[i, ] = err_cv
  }
  
  cvm = apply(err_mat, 2, mean, na.rm = TRUE)
  nn = apply(!is.na(err_mat), 2, sum, na.rm = T)
  cvse = sqrt(apply(err_mat, 2, var, na.rm = T)/nn)
  
  cvsd = apply(err_mat, 2, sd)
  #cvsd = sqrt(apply(scale(err_mat, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nfolds - 1))
  
  cvlo = cvm - cvse
  cvup = cvm + cvse
  cvm_min_ind = which.min(cvm)
  
  lambda.min = lambda0[cvm_min_ind]
  imin.1se = which(cvm < cvm[cvm_min_ind] + cvse[cvm_min_ind])[1]
  lambda.1se = lambda0[imin.1se]
  
  if (fit_mode == "1se"){
    best_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda.1se, penalty.factor = pf_values)
  } else if (fit_mode == "min"){
    best_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda.min, penalty.factor = pf_values)
  }
  full_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda0, penalty.factor = pf_values)
  n_nonzero = full_fit$df
  
  return(list(cvm=cvm, cvsd=cvse, lambda=lambda0, 
              lambda.min=lambda0[cvm_min_ind], lambda.1se = lambda.1se,
              best_fit = best_fit, 
              ind_min = cvm_min_ind, ind_1se = imin.1se,
              support=n_nonzero, cvup=cvup, cvlow=cvlo, pf_values=pf_values))
}

# Plot CV Curve for Cooperative Learning
plot_cv_coop = function(object, title=NULL, lamx=NULL, lamz=NULL){
  plot(x=log(object$lambda),y=object$cvm,ylim=range(object$cvup,object$cvlo),
       xlab=expression(Log(lambda)),ylab="MSE", type="p", pch=20,
       main=title, cex.axis=1, cex.main=0.9)
  segments(log(object$lambda), object$cvup, log(object$lambda), object$cvlow, col="grey")
  axis(side=3,at=log(object$lambda),labels=paste(object$support),
       tick=FALSE,line=0,cex.axis=0.5)
  abline(v=log(object$lambda.min),lty=3, col="deepskyblue3", lwd=1.2)
  abline(v=log(object$lambda.1se),lty=3, col="darkgreen", lwd=1.2)
  mtext(paste("lamx:",lamx," ","lamz:",lamz), side=1, las=1,cex=0.6)
}



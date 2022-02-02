source("cooperative_regression_function.R")

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

#coop cv implementation
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
    
    x_sub = scale(x_sub, T, F)
    z_sub = scale(z_sub, T, F)
    
    #fit model
    xt = rbind(cbind(x_sub, z_sub),
               cbind(-sqrt(alpha)*x_sub, sqrt(alpha)*z_sub))
    N = dim(x_sub)[1]
    yt = c(y_sub, rep(0, N))
    #coop lambda = 1/2 lasso lambda
    #coop intercept = 1/2 lasso intercept
    coop = glmnet(xt, yt, standardize=F, lambda=lambda0, penalty.factor = pf_values)
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

plot_cv_coop = function(object, title=NULL, lamx=NULL, lamz=NULL){
  lambda_iter_ind = which.min(abs(object$lambda-1))
  plot(x=log(object$lambda),y=object$cvm,ylim=range(object$cvup,object$cvlo),
       xlab=expression(Log(lambda)),ylab="MSE", type="p", pch=20,
       col=ifelse(object$lambda==object$lambda[lambda_iter_ind], "red", "blue"),
       main=title, cex.axis=1, cex.main=0.9)
  segments(log(object$lambda), object$cvup, log(object$lambda), object$cvlow, col="grey")
  axis(side=3,at=log(object$lambda),labels=paste(object$support),
       tick=FALSE,line=0,cex.axis=0.5)
  abline(v=log(object$lambda.min),lty=3)
  mtext(paste("lamx:",lamx," ","lamz:",lamz), side=1, las=1,cex=0.6)
}


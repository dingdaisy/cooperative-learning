library(glmnet)

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

#Cooperative Learning CV for 3 Data Views
coop_cv_multi = function(x1,x2,x3,y,alpha=0,foldid,nfolds=10,
                       pf_values=NULL,n_view=3,fit_mode='min'){
  
  #get lambda sequence
  xt0 = rbind(cbind(x1, x2, x3),
              cbind(-sqrt(alpha)*x1, sqrt(alpha)*x2, matrix(0,nrow(x3),ncol(x3))),
              cbind(-sqrt(alpha)*x1, matrix(0,nrow(x2),ncol(x2)), sqrt(alpha)*x3),
              cbind(matrix(0,nrow(x1),ncol(x1)), -sqrt(alpha)*x2, sqrt(alpha)*x3))
  n_pairs = (n_view * (n_view-1)) / 2
  yt0 = c(y, rep(0, dim(x2)[1]*n_pairs))
  coop0 = glmnet(xt0, yt0, standardize=F, penalty.factor = pf_values)
  lambda0 = coop0$lambda
  
  outlist = as.list(seq(nfolds))
  err_mat = matrix(NA, nfolds, length(lambda0))
  for (i in seq(nfolds)){
    which = (foldid == i)
    x1_sub = x1[!which, , drop = FALSE]
    x2_sub = x2[!which, , drop = FALSE]
    x3_sub = x3[!which, , drop = FALSE]
    y_sub = y[!which]
    
    #centering inside necessary
    x1_sub = scale(x1_sub, T, F)
    x2_sub = scale(x2_sub, T, F)
    x3_sub = scale(x3_sub, T, F)
    
    #fit model
    xt = rbind(cbind(x1_sub, x2_sub, x3_sub),
                cbind(-sqrt(alpha)*x1_sub, sqrt(alpha)*x2_sub, matrix(0,nrow(x3_sub),ncol(x3_sub))),
                cbind(-sqrt(alpha)*x1_sub, matrix(0,nrow(x2_sub),ncol(x2_sub)), sqrt(alpha)*x3_sub),
                cbind(matrix(0,nrow(x1_sub),ncol(x1_sub)), -sqrt(alpha)*x2_sub, sqrt(alpha)*x3_sub))
    
    N = dim(x1_sub)[1]
    yt = c(y_sub, rep(0, N*n_pairs))
    #coop lambda = 1/2 lasso lambda
    #coop intercept = 1/2 lasso intercept
    coop = glmnet(xt, yt, standardize=F, lambda=lambda0, penalty.factor = pf_values)
    outlist[[i]] = coop
    #lasso = glmnet(cbind(x_sub,z_sub),y_sub,standardize=F, lambda=lambda0*2) #coop$lambda*2)
    
    #evaluate
    theta = coef(coop)[-1,]
    intercept_adjustment_factor = (n_pairs + 1)
    intercept = coef(coop)[1,] * intercept_adjustment_factor
    
    pred_fold = cbind(x1[which, , drop = FALSE],x2[which, , drop = FALSE],x3[which, , drop = FALSE])%*%theta + 
      rep(intercept, each=nrow(x1[which, , drop = FALSE]))
      
    true_y = y[which]
    err_mse = (pred_fold - replicate(ncol(pred_fold), true_y))^2
    err_cv = apply(err_mse, 2, mean, na.rm = TRUE)
    err_mat[i, ] = err_cv
  }
  
  cvm = apply(err_mat, 2, mean, na.rm = TRUE)
  nn = apply(!is.na(err_mat), 2, sum, na.rm = T)
  cvse = sqrt(apply(err_mat, 2, var, na.rm = T)/nn)
  
  cvlo = cvm - cvse
  cvup = cvm + cvse
  cvm_min_ind = which.min(cvm)
  
  lambda.min = lambda0[cvm_min_ind]
  #print(lambda.min)
  imin.1se = which(cvm < cvm[cvm_min_ind] + cvse[cvm_min_ind])[1]
  lambda.1se = lambda0[imin.1se]
  #print(lambda.1se)
  
  if (fit_mode == "1se"){
    best_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda.1se, penalty.factor = pf_values)
  } else if (fit_mode == "min"){
    best_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda.min, penalty.factor = pf_values)
  }
  full_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda0, penalty.factor = pf_values)
  n_nonzero = full_fit$df
  
  #cvm = apply(err_mat, 2, mean, na.rm = TRUE)
  #cvsd = sqrt(apply(scale(err_mat, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nfolds - 1))
  
  return(list(cvm=cvm, cvsd=cvse, lambda=lambda0, 
              lambda.min=lambda0[cvm_min_ind], lambda.1se = lambda.1se,
              best_fit = best_fit, 
              best_fit_coef = best_fit$beta,
              best_fit_intercept = best_fit$a0 * intercept_adjustment_factor,
              ind_min = cvm_min_ind, ind_1se = imin.1se,
              support=n_nonzero, cvup=cvup, cvlow=cvlo, pf_values=pf_values))
}

#Plot CV Curve for Cooperative Learning
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


#Cooperative Regression: Iterative
calc_objective_multi <- function(thetax1,thetax2,thetax3,
                           x1_intercept,x2_intercept,x3_intercept,
                           alpha,x1,x2,x3,y){
  res = sum((y-(x1%*%thetax1+x1_intercept)-(x2%*%thetax2+x2_intercept)-(x3%*%thetax3+x3_intercept))^2)/2 + 
    alpha/2*(sum(((x1%*%thetax1+x1_intercept)-(x2%*%thetax2+x2_intercept))^2)+
               sum(((x1%*%thetax1+x1_intercept)-(x3%*%thetax3+x3_intercept))^2)+
               sum(((x3%*%thetax3+x3_intercept)-(x2%*%thetax2+x2_intercept))^2))
  return(res)
}

coop_regression_iter_multi = function(x1,x2,x3,y,alpha,
                                thetax1=NULL, thetax2=NULL, thetax3=NULL, 
                                foldid=NULL,
                                intercept_x1=0, intercept_x2=0,intercept_x3=0,
                                max_iteration=3, thr=0.05){
  n = length(y)
  if(is.null(thetax1)){thetax1 = rep(0,ncol(x1))}
  if(is.null(thetax2)){thetax2 = rep(0,ncol(x2))}
  if(is.null(thetax3)){thetax3 = rep(0,ncol(x3))}
  cv_error_x1 = NULL
  cv_error_x2 = NULL
  cv_error_x3 = NULL
  
  crit = NULL
  delta = 2 * thr
  crit_0 = calc_objective_multi(thetax1,thetax2,thetax3,
                                intercept_x1,intercept_x2,intercept_x3,
                                alpha,x1,x2,x3,y)
  crit = c(crit, crit_0)
  n_iteration = 0
  
  lam_x1 = NULL
  lam_x2 = NULL
  lam_x3 = NULL
  M = 3
  
  while(delta > thr & n_iteration < max_iteration){
    old_crit = crit[length(crit)]
    n_iteration = n_iteration + 1
    
    y1 = y/(1+(M-1)*alpha) - 
      ((1-alpha)/(1+(M-1)*alpha))*((x2%*%thetax2+intercept_x2) + (x3%*%thetax3+intercept_x3))
    gfit1 = cv.glmnet(x1, y1, standardize=F, foldid=foldid)
    
    thetax1_temp=coef(gfit1,s=gfit1$lambda.min)[-1]
    intercept_x1_temp=coef(gfit1,s=gfit1$lambda.min)[1]
    cv_error_x1_temp=gfit1$cvm[which(gfit1$lambda == gfit1$lambda.min)]
    lam_x1_temp = gfit1$lambda.min
    
    y2 = y/(1+(M-1)*alpha) - 
      ((1-alpha)/(1+(M-1)*alpha))*((x1%*%thetax1_temp+intercept_x1_temp) + (x3%*%thetax3+intercept_x3))
    gfit2 = cv.glmnet(x2, y2, standardize=F, foldid=foldid)
    
    thetax2_temp=coef(gfit2,s=gfit2$lambda.min)[-1]
    intercept_x2_temp=coef(gfit2,s=gfit2$lambda.min)[1]
    cv_error_x2_temp=gfit2$cvm[which(gfit2$lambda == gfit2$lambda.min)]
    lam_x2_temp = gfit2$lambda.min
    
    y3 = y/(1+(M-1)*alpha) - 
      ((1-alpha)/(1+(M-1)*alpha))*((x1%*%thetax1_temp+intercept_x1_temp) + (x2%*%thetax2_temp+intercept_x2_temp))
    gfit3 = cv.glmnet(x3, y3, standardize=F, foldid=foldid)
    
    thetax3_temp=coef(gfit3,s=gfit3$lambda.min)[-1]
    intercept_x3_temp=coef(gfit3,s=gfit3$lambda.min)[1]
    cv_error_x3_temp=gfit3$cvm[which(gfit3$lambda == gfit3$lambda.min)]
    lam_x3_temp = gfit3$lambda.min

    current_obj = calc_objective_multi(thetax1_temp,thetax2_temp,thetax3_temp,
                                 intercept_x1_temp,intercept_x2_temp,intercept_x3_temp,
                                 alpha,x1,x2,x3,y)
    
    if (current_obj < old_crit){
      thetax1 = thetax1_temp
      thetax2 = thetax2_temp
      thetax3 = thetax3_temp
      intercept_x1 = intercept_x1_temp
      intercept_x2 = intercept_x2_temp
      intercept_x3 = intercept_x3_temp
      cv_error_x1 = cv_error_x1_temp
      cv_error_x2 = cv_error_x2_temp
      cv_error_x3 = cv_error_x3_temp
      lam_x1 = lam_x1_temp
      lam_x2 = lam_x2_temp
      lam_x3 = lam_x3_temp
      crit = c(crit, current_obj)
    }
    
    new_crit = crit[length(crit)]
    delta = (old_crit-new_crit)/abs(old_crit)
  }
  
  return(list(thetax1=thetax1, thetax2=thetax2, thetax3=thetax3,
              intercept_x1=intercept_x1, intercept_x2=intercept_x2, intercept_x3=intercept_x3,
              cv_error_x1=cv_error_x1,  cv_error_x2=cv_error_x2,  cv_error_x3=cv_error_x3, 
              lam_x1=lam_x1, lam_x2=lam_x2, lam_x3=lam_x3,
              alpha=alpha, crit=crit,
              n_iteration = n_iteration,
              n_thetax1 = sum(thetax1!=0),
              n_thetax2 = sum(thetax2!=0),
              n_thetax3 = sum(thetax3!=0)))
}

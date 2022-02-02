library(glmnet)
library(randomForest)
#source("gbm_helper.R")

#Cooperative Regression: Full
coop_regression_full = function(x,z,y,alpha){
  n = length(y)
  xt = rbind(cbind(x, z),
             cbind(-sqrt(alpha)*x, sqrt(alpha)*z))
  yt = c(y, rep(0,n))
  g_fit = glmnet(xt, yt, standardize=F)
  return(g_fit)
}

#"Stopping Criteria"
calc_objective <- function(thetax,thetaz,x_intercept,z_intercept,alpha,x,z,y){
  res = sum((y-(x%*%thetax+x_intercept)-(z%*%thetaz+z_intercept))^2)/2 + 
    alpha*sum(((x%*%thetax+x_intercept)-(z%*%thetaz+z_intercept))^2)/2 
  return(res)
}

#Cooperative Regression: Iterative
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

#"Stopping Criteria: General"
calc_objective_general <- function(X_fit,Z_fit,alpha,y){
  res = sum((y-X_fit-Z_fit)^2)/2 + alpha*sum((X_fit-Z_fit)^2)/2
  return(res)
}

#Cooperative Regression: General
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



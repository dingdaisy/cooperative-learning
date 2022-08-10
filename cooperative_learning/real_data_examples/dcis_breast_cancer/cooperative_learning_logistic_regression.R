
#Calculate AUROC
myroc = function(ytest,rit){
  # ytest is test outcome, rit is predicted value (continous)
  sens=rep(0,100);spec=ppv=npv=sens
  qq=quantile(rit,probs=seq(0,1,len=100))
  for(i in 1:100){
    yhat=0*(rit<=qq[i])+1*(rit>qq[i])
    sens[i]=sum((ytest==1)&(yhat==1))/sum(ytest==1)
    spec[i]=sum((ytest==0)&(yhat==0))/sum(ytest==0)
    npv[i]=sum((ytest==0)&(yhat==0))/sum(yhat==0)
    ppv[i]=sum((ytest==1)&(yhat==1))/sum(yhat==1)
  }
  sens=c(1,sens); spec=c(0,spec)
  area=0
  for(i in 1:(length(sens)-1)){
    area=area+.5*(spec[i+1]-spec[i])*(sens[i+1]+sens[i])
  }
  
  a=area
  q1=a/(2-a)
  q2=2*a*a/(1+a)
  nn=sum(ytest==0)
  na=sum(ytest==1)
  se=sqrt( (a*(1-a)+(na-1)*(q1-a*a)+(nn-1)*(q2-a*a))/(na*nn))
  
  val=(spec+sens)/2
  cutp.max=qq[which.max(val)]
  return(list(sens=sens,spec=spec,ppv=ppv,npv=npv,area=area,se=se,cutp=qq,cutp.max=cutp.max))
}

#CV Measures
cv_measure = function(predmat, y, nfolds, foldid, type.measure="deviance"){
  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(y)
  raw_y = y 
  if (is.null(nc)) {
    y = factor(y, levels=0:1)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), , drop=FALSE]
  }
  
  predmat=1/(1+exp(-predmat))
  nlambda=ncol(predmat)
  nlams=rep(nlambda,nfolds)
  if (type.measure == "deviance"){
    predmat = pmin(pmax(predmat, prob_min), prob_max)
    lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
    ly = log(y)
    ly[y == 0] = 0
    ly = drop((y * ly) %*% c(1, 1))
    res = 2 * (ly - lp)
    N = nrow(y) - apply(is.na(predmat), 2, sum)
  }
  else if (type.measure == "auc") {
    res = matrix(NA, nfolds, nlambda)
    good = matrix(0, nfolds, nlambda)
    for (i in seq(nfolds)) {
      good[i, seq(nlams[i])] = 1
      which = foldid == i
      for (j in seq(nlams[i])) {
        res[i, j] = myroc(raw_y[which],predmat[which,j])$area 
      }
    }
    N = apply(good, 2, sum)
  }
  else if (type.measure == "mse"){
    res = (y[, 1] - (1 - predmat))^2 + (y[, 2] - predmat)^2
    N = nrow(y) - apply(is.na(predmat), 2, sum)
  }
  else if (type.measure == "class"){
    res = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5)
    N = nrow(y) - apply(is.na(predmat), 2, sum)
  }
  return(list(cv=res, N=N))
}

#Cooperative Logistic Regression
multiview.logistic = function(x, z, y, lambda, theta = 0, 
                              niter = 50, beta = rep(0, 1 + ncol(x) + ncol(z)),
                              meth="glmnet", verbose=F){
  
  # version that centers features + target
  # inputs lambda in standard parametrization
  my = mean(y)
  beta = c(log(my/(1-my)), rep(0, ncol(x) + ncol(z)))  #NOTE, with intercept
  n = length(y)
  
  nx = ncol(x)
  nz = ncol(z)
  
  xall = cbind(x, z)
  
  i = 0
  
  x_list=list(x=x, z=z)
  beta_all = c(beta[1], beta[-1])
  p1 <- ncol(x_list[[1L]])
  a0 <- beta_all[1L]
  betax <- beta_all[2L:(p1 + 1L)]
  betaz <- beta_all[-(1L:(p1 + 1L))]
  eta <- a0 + x %*% betax + z %*% betaz
  pr <- 1/(1 + exp(-eta))
  crit0 <- -sum(y * log(pr) + (1 - y) * log(1 - pr)) + 
    lambda * (sum(abs(betax)) + sum(abs(betaz))) + 
    (theta/2) * sum((x %*% betax - z %*% betaz)^2) 
  
  #crit0 = critf(x_list=list(x=x, z=z), y, lambda = lambda, beta_all = c(beta[1], beta[-1]),
  #              theta = theta)
  crit_imp = -10
  
  #for (i in 1:niter) {
  while (i < niter & crit_imp < 0){
    i = i + 1
    xx = xall
    eta = cbind(1, xx) %*% beta
    pr = 1/(1 + exp(-eta))
    w = pr * (1 - pr)
    zz = eta + (y - pr)/(pr*(1-pr))  #working response
    
    # center working response and features, using weights
    mzz = sum(w * zz)/sum(w)
    zzc = zz - mzz
    
    mx = rep(NA, ncol(xx))
    for (j in 1:ncol(xx)) {
      mx[j] = sum(w * xx[, j])/sum(w)
      xx[, j] = xx[, j] - mx[j]
    }
    
    features = xx
    target = zzc
    fac = n/sum(w)
    
    lambda2 = lambda*fac/n  #glmnet parametrization
    if (theta > 0) {
      features = rbind(features, cbind(-sqrt(theta) * xx[, 1:nx], 
                                       sqrt(theta) * xx[, -(1:nx)]))
      target = c(target, rep(0, n))
      w = c(w, rep(1,n))  #NOTE, weights for the agreement penalty
      fac = 2*n/sum(w)  #NOTE, 2 times sample size
      lambda2=lambda*fac/(2*n)  #NOTE, adjust lambda
    }
    
    if (meth=="cvx"){
      out = make_cvxr_problem_gaussian0(features, target, lambda2, w, intercept = FALSE, glmnet_lambda = TRUE)
      result = solve(out$prob, verbose = FALSE)
      beta = result$getValue(result$beta)
    }
    if (meth=="glmnet"){
      #fit = glmnet(features, target, lambda = lambda2, weights = w, standardize = F, intercept = F)
      fit = glmnet:::glmnet.fit(features, target, lambda = lambda2, weights = w/sum(w), intercept=F) #NOTE, intercept = F
      b = as.numeric(fit$beta)
      b0 = mzz - sum(mx * b) #mzz: working response mean; mx: feature mean
      beta = c(b0, b)
    }
    
    x_list = list(x=x, z=z)
    beta_all = c(beta[1], beta[-1])
    p1 <- ncol(x_list[[1L]])
    a0 <- beta_all[1L]
    betax <- beta_all[2L:(p1 + 1L)]
    betaz <- beta_all[-(1L:(p1 + 1L))]
    eta <- a0 + x %*% betax + z %*% betaz
    pr <- 1/(1 + exp(-eta))
    crit <- -sum(y * log(pr) + (1 - y) * log(1 - pr)) + 
      lambda * (sum(abs(betax)) + sum(abs(betaz))) + 
      (theta/2) * sum((x %*% betax - z %*% betaz)^2)
    
    #crit = critf(x_list=list(x=x, z=z), y, lambda = lambda, beta_all = c(beta[1], beta[-1]), theta = theta)
    crit_imp = (crit - crit0) / crit0
    
    #print(crit)
    #print(crit_imp)
    #print(lambda)
    #print(crit)
    #print(crit0)
    if (crit < crit0){
      crit0 = crit
    }
    
    if (verbose){
      print(i)
      cat(sprintf("mzz: %f, term: %f, b0: %f\n", mzz, sum(mx*b), b0))
      cat(c("crit=", round(crit, 5), "b= ", round(beta[1:6], 3)), fill = T)
    }
  }
  return(beta)
}


#Cooperative Logistic Regression with CV
coop_logistic_cv = function(x,z,y,alpha=0,foldid,nfolds=5,pf_values=NULL,
                            niter_LR=50,verbose=F,fit_mode='min',type.measure='deviance'){
  
  #get lambda sequence
  coop0 = glmnet(cbind(x, z), y, standardize=F, family="binomial", penalty.factor = pf_values)
  lambda0 = coop0$lambda * nrow(x)
  
  outlist = as.list(seq(nfolds))
  err_mat = matrix(NA, nfolds, length(lambda0))
  pred_all = matrix(NA, length(y), length(lambda0))
  
  for (i in seq(nfolds)){
    #print(i)
    which = (foldid == i)
    x_sub = x[!which, , drop = FALSE]
    z_sub = z[!which, , drop = FALSE]
    y_sub = y[!which]
    
    #centering inside necessary
    x_sub = scale(x_sub, T, F)
    z_sub = scale(z_sub, T, F)
    
    #fit model
    coop_beta = sapply(lambda0, FUN=function(lambda_each) multiview.logistic(x_sub, z_sub, 
                                                           y_sub, lambda_each, 
                                                           theta = alpha, 
                                                           niter = niter_LR, 
                                                           beta = rep(0,1+ncol(x_sub)+ncol(z_sub)),
                                                           meth="glmnet",verbose=verbose))
    theta = coop_beta[-1,]
    intercept = coop_beta[1,]
    
    pred_fold = cbind(x[which, , drop = FALSE], z[which, , drop = FALSE]) %*% theta + 
                      rep(intercept, each=nrow(x[which, , drop = FALSE]))
    pred_all[which,] = pred_fold
  }
  
  err_cv_all = cv_measure(pred_all, as.vector(y), nfolds=nfolds, 
                          foldid=foldid, type.measure=type.measure)
  cvm = apply(err_cv_all$cv, 2, mean, na.rm = TRUE)
  cvse = sqrt(apply(err_cv_all$cv, 2, var, na.rm = T)/ (err_cv_all$N-1))
  
  cvlo = cvm - cvse
  cvup = cvm + cvse
  
  if (type.measure == 'auc'){
    cvm_min_ind = which.max(cvm)
    imin.1se = which(cvm > cvm[cvm_min_ind] - cvse[cvm_min_ind])[1]
  } else{
    cvm_min_ind = which.min(cvm)
    imin.1se = which(cvm < cvm[cvm_min_ind] + cvse[cvm_min_ind])[1]
  }
  
  lambda.min = lambda0[cvm_min_ind]
  lambda.1se = lambda0[imin.1se]
  
  if (fit_mode == "min"){
    best_fit = multiview.logistic(x, z, y, lambda.min, theta = alpha, 
                niter = niter_LR, beta = rep(0, 1 + ncol(x) + ncol(z)),
                meth="glmnet")
  } else if (fit_mode == "1se"){
    best_fit = multiview.logistic(x, z, y, lambda.1se, theta = alpha, 
                                  niter = niter_LR, beta = rep(0, 1 + ncol(x) + ncol(z)),
                                  meth="glmnet")
  }
  #full_fit = glmnet(xt0, yt0, standardize=F, lambda = lambda0, penalty.factor = pf_values)
  #n_nonzero = full_fit$df
  
  return(list(cvm=cvm, cvsd=cvse, lambda=lambda0, 
              lambda.min=lambda.min, lambda.1se = lambda.1se,
              best_fit = best_fit, 
              ind_min = cvm_min_ind, ind_1se = imin.1se,
              #support=n_nonzero, 
              cvup=cvup, cvlow=cvlo, pf_values=pf_values))
}



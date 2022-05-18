library(ISLR2)
library(glmnet)
library(keras)
library(tensorflow)
library(sp)
library(dplyr)
library(gplots)
library(fields)
library(caret)
library(RColorBrewer)
library(ggplot2)

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))}

range_01 <- function(x){(x-min(x))/(max(x)-min(x))}

make_boundary = function(){
  a1 = runif(1,min=.5, max=1)
  b1 = runif(1,min=.5 ,max=1)
  a2 = runif(1,min=-1, max=-.5)
  b2 = runif(1,min=-1, max=-.5)
  vert = rbind(c(0,a1),c(b1,0),c(0,a2),c(b2,0),c(0,a1))
  return(vert)
}

make_image = function(n, isig, n_dim=32, max_intensity=6, score){
  #level of intensity corresponds to disease score
  score_01 = range_01(score)
  score_intensity = score_01 * max_intensity
  
  x = array(0,c(n,n_dim,n_dim,1))
  ind = array(NA,c(n,n_dim,n_dim))
  #isig = 1 * (sig > median(sig))
  res = matrix(NA,n_dim,n_dim)
  
  for (ii in 1:n){
    vert = make_boundary()
    for (j in sample(1:n_dim, size=n_dim)){
      #generate a random position in the image
      a = -1+2*(j-1)/n_dim
      for (k in sample(1:n_dim, size=n_dim)){
        b = -1+2*(k-1)/n_dim
        #if the point is in the polygon
        res[j,k] = point.in.polygon(a,b,vert[,1],vert[,2])
        if (res[j,k] == 1)
          #generate a random intensity for the point when it is in the polygon
          {x[ii,j,k,1] = runif(1)}
        }}
    
    #random matrix with n_dim x n_dim
    u = matrix(runif(n_dim*n_dim) < .01, n_dim, n_dim)
    for (j in 1:n_dim){
        for (k in 1:n_dim){
            #if a point is in the polygon + random matrix True + signal true, 
            #then the corresponding score_intensity
            if(res[j,k]==1 & u[j,k] & isig[ii]==1) x[ii,j,k,1] = score_intensity[ii]
        }
    }
    ind[ii,,] = u & res
    }
  return(list(x=x, ind=ind))
  }

#simulation settings
train_frac = 0.1
val_frac = 0.1

n = 2000
p_z = 100
p_r = 30
sigma_noise = 60
beta_strength = 1
level_intensity = 3.9 #4 is doing well too
u_std = 2 #sd of latent factors
s = 6 #strength of u

nsim = 10
niter = 2
nepoch = 5
#late_fusion_val_frac = 0.3
alphalist = c(0,.5,1,3,5,9,15,20)
set.seed(9)

res = matrix(NA,nsim,length(alphalist)+4)
alpha_selected = rep(NA, nsim)

for (ii in 1:nsim){
    z_all = matrix(rnorm(n*p_z), n, p_z)
    u_matrix = matrix(0, n, p_r)
    for (m in seq(p_r)){
      u = rnorm(n, sd = u_std)
      z_all[, m] = z_all[, m] + s*u
      u_matrix[, m] = s*u
    }

    #beta = rep(beta_strength, p_r)
    #signal_all = u_matrix %*% beta
    beta = c(rep(beta_strength, p_r), rep(0, p_z - p_r))
    signal_all = z_all %*% beta
    signal_all_noise = signal_all + sigma_noise * rnorm(n)
    snr = var(signal_all) / var(signal_all_noise-signal_all)
    
    prob_all = 1 / (1 + exp(-signal_all_noise)) 
    label_all = 1 * (runif(n) < prob_all) 

    x_all = make_image(n,isig=label_all,max_intensity=level_intensity,score=signal_all_noise)$x
    
    smp_size_train = floor(train_frac * nrow(x_all)) 
    smp_size_val = floor(val_frac * nrow(x_all))
    train_ind = sort(sample(seq_len(nrow(x_all)), size = smp_size_train))
    ind_no_train = setdiff(seq_len(nrow(x_all)), train_ind)
    val_ind = sort(sample(ind_no_train, size = smp_size_val))
    test_ind = setdiff(ind_no_train, val_ind)
    
    colnames(z_all) = seq(ncol(z_all))
    
    x_train = x_all[train_ind,,,,drop=F]
    x_val = x_all[val_ind,,,,drop=F]
    x_test = x_all[test_ind,,,,drop=F]
    z_train = z_all[train_ind,]
    z_val = z_all[val_ind,]
    z_test = z_all[test_ind,]
    
    label_train = label_all[train_ind]
    label_val = label_all[val_ind]
    label_test = label_all[test_ind]
    y_train = signal_all_noise[train_ind]
    y_val = signal_all_noise[val_ind]
    y_test = signal_all_noise[test_ind]
    
    preprocess_values_train = preProcess(z_train, method = c("center", "scale"))
    z_train = predict(preprocess_values_train, z_train)
    z_val = predict(preprocess_values_train, z_val)
    z_test = predict(preprocess_values_train, z_test)
    
    #look at 2 examples
    plotit=F
    #old_plot
    if(plotit){
      pdf("sample_images_vF.pdf", width=10.9, height=3)
      par(mfrow=c(1,3), las=2, cex.main=1.2, cex.axis = 1.03, mar=c(3,3,3,5))
      o1=which(label_all==0)[1]
      o2=which(label_all==1)[6]
      o3=which(label_all==1)[3]
      upper_limit = max(max(x_all[o3,,,]),max(x_all[o1,,,]),max(x_all[o2,,,]))
      
      image.plot(x_all[o1,,,], col = brewer.pal(9, "Greys"), zlim=c(0,upper_limit))
      mtext("Healthy Sample", side=3, line=0.6,las=1)
      
      image.plot(x_all[o2,,,], col = brewer.pal(9, "Greys"), zlim=c(0,upper_limit))
      #image.plot(x_all[o3,,,], col = gray.colors(12,rev=T), zlim=c(0,level_intensity))
      mtext("Disease Sample", side=3, line=0.6,las=1)
      
      image.plot(x_all[o3,,,], col = brewer.pal(9, "Greys"), zlim=c(0,upper_limit))
      mtext("Disease Sample", side=3, line=0.6,las=1)

      dev.off()
    }
    
    plotit_new=F
    if(plotit_new){
      pdf("sample_images_vF.pdf", width=6.9, height=3)
      
      par(mfrow=c(1,2), las=2, cex.main=1.2, cex.axis = 1.03, mar=c(3,3,3,3))
      p1=which(label_all==0)[1]
      p2=which(label_all==1)[3]
      upper_limit = max(max(x_all[p1,,,]),max(x_all[p2,,,]))
      
      image.plot(x_all[p1,,,], col = c("white", brewer.pal(9, "Blues")), zlim=c(0,upper_limit))
      mtext("Healthy Sample", side=3, line=0.6,las=1)
      
      image.plot(x_all[p2,,,], col = c("white",  brewer.pal(9, "Blues")), zlim=c(0,upper_limit))
      mtext("Disease Sample", side=3, line=0.6,las=1)
      
      dev.off()
    }

  # image data on its own---fit CNN
  model <- keras_model_sequential() %>%
     layer_conv_2d(filters = 32, kernel_size = c(3, 3),
        padding = "same", activation = "relu",
        input_shape = c(32, 32, 1)) %>%
     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
     layer_flatten() %>%
     layer_dense(units = 32, activation = "relu") %>%
     layer_dense(units = 1, activation = "linear")
  
  #summary(model)
  
  model %>% compile(loss = "mean_squared_error",
                   optimizer = optimizer_rmsprop(), metrics = c("mse"))
  history <- model %>% fit(x_train, label_train, epochs = nepoch)
  
  yhat_train = model %>% predict(x_train)
  yhat_val = model %>% predict(x_val)
  yhat_test = model %>% predict(x_test)
  
  yhatc_test = 1*(yhat_test>mean(yhat_test))
  res[ii,1] = mean(yhatc_test != label_test)
  print("Images")
  print(res[ii,1])

  #proteomics on its own
  gfit = cv.glmnet(z_train,label_train,standardize=F)
  yhatz_train = as.vector(predict(gfit,z_train,s="lambda.min"))
  yhatzc_train = 1*(yhatz0>mean(yhatz0))
  yhatz_val = as.vector(predict(gfit,z_val))
  yhatz_test = as.vector(predict(gfit,z_test))
  yhatzc_test = 1*(yhatz_test>mean(yhatz_test))
  #table(ytest,yhatzc)
  #err_lasso = calc_mse(yhatz_test,y_test)
  res[ii,2] = mean(yhatzc_test!=label_test)
  #res_score[ii,2] = err_lasso
  print("Genes")
  print(res[ii,2])
  
  #late fusion
  fuse_data = data.frame(y=label_val, X_pred=as.vector(yhat_val), Z_pred=as.vector(yhatz_val))
  fit_fuse = lm(y ~ X_pred + Z_pred, data=fuse_data)
  fuse_pred_test = predict(fit_fuse, data.frame(X_pred=as.vector(yhat_test), 
                                                Z_pred=as.vector(yhatz_test)))
  #fuse_pred_test = (yhat_test+yhatz_test)/2
  late_yhatc = 1*(fuse_pred_test>mean(fuse_pred_test))
  res[ii,3] = mean(late_yhatc!=label_test)
  #res_score[ii,3] = calc_mse(fuse_pred_test,y_test)
  print("Late")
  print(res[ii,3])
  
  #coop
  val_mis = rep(0, times = length(alphalist))
  test_mis = rep(0, times = length(alphalist))
  
  print("Coop")
  for (j in 1:length(alphalist)){
    alpha = alphalist[j]
    print(alpha)
    fz=rep(0,length(label_train))
    for (kk in 1:niter){
      r1 = label_train/(1+alpha) - (1-alpha)*fz/(1+alpha)
      gfit = cv.glmnet(z_train,r1)
      fx = predict(gfit,z_train,s="lambda.min")
      
      r2 = label_train/(1+alpha) - (1-alpha)*fx/(1+alpha)
      history = model %>% fit(x_train, r2, epochs = 5)
      fz = model %>% predict(x_train)
    }
    yhat_coop_val = predict(gfit,z_val) + model %>% predict(x_val)
    yhatc_coop_val = 1 * (yhat_coop_val>mean(yhat_coop_val))
    val_mis[j] = mean(yhatc_coop_val != label_val)
    
    yhat_coop = predict(gfit,z_test) + model %>% predict(x_test)
    yhatc_coop = 1 * (yhat_coop>mean(yhat_coop))
    test_mis[j] = mean(yhatc_coop != label_test)
    #table(ytest,yhatgc)
    res[ii,j+3+1] = mean(yhatc_coop != label_test)
    print(res[ii,j+3+1])
  }
  res[ii,4] = test_mis[which.min(val_mis)]
  alpha_selected[ii] = alphalist[which.min(val_mis)]
}

saveRDS(res, paste0("error_image_gene.rds"))
saveRDS(alpha_selected, paste0("alpha_image_gene.rds"))

out = rbind(colMeans(res), apply(res,2,sd)/sqrt(nsim))
dimnames(out) = list(c("mean","se"),
                     c("image","prot","late","coop",paste("coop",as.character(alphalist),sep="")))

alpha_axis = sapply(as.character(round(alphalist,1)),
         function(x) paste0("Coop (alpha=",x,")"))

#title still hard-coded right now
pdf(file="test_misclassification_images_genes.pdf", width=10, height=5)
par(mfrow=c(1,2), las=2, cex.main=1.2, cex.axis = 1.03, mar=c(7.8,3.3,3.2,0.9))
boxplot(c(res[,1:4])
        ~sort(rep(1:4, nsim)), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist)),
        #        "#ff8080",rep("#BCDEBC",length(alphalist))),
        names=c("Only Images","Only Proteomics", "Late", "Coop"), #,alpha_axis), 
        xlab="",ylab="", main="Test Misclassification Error (SNR = 1)")

images_diff = res - res[,3]
boxplot(c(images_diff[,1:4])
        ~sort(rep(1:4, nsim)), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist)),
        #        "#ff8080",rep("#BCDEBC",length(alphalist))),
        names=c("Only Images","Only Proteomics", "Late", "Coop"), #, alpha_axis), 
        xlab="",ylab="", main="Difference with Late Fusion (SNR = 1)")
abline(h=0, lty=2)
dev.off()


library(glmnet)
library(latex2exp)
library(RColorBrewer)
library(viridis)
library(scico)

set.seed(3333)
n=100
p=20
nsim=50
alphalist=c(.01,.05,.1,.25,.5,.75,1,10,100)

x=matrix(rnorm(n*p),n,p)
z=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)
z=scale(z,T,F)
b=rep(2,p)
sigma=6
mu=x%*%b
snr=var(mu)/sigma^2
cat(c("snr=",snr),fill=T)
y=mu+sigma*rnorm(n)
y=y-mean(y)

 coop00 = glmnet(cbind(x,z),y, standardize=F)

lambda00=coop00$lambda
lambda00=c(2*max(lambda00),lambda00)

nz0=norm0=matrix(NA,nsim,length(lambda00))
nz=norm=array(NA, c(nsim,length(alphalist),length(lambda00)))



for(ii in 1:nsim){
x=matrix(rnorm(n*p),n,p)
z=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)
z=scale(z,T,F)
y=rnorm(n)
y=y-mean(y)


 coop0= glmnet(cbind(x,z),y, standardize=F)
beta0=predict(coop0,s=lambda00,type="coef")

nz0[ii,]=colSums(beta0[-1,]!=0)
norm0[ii,]=colSums(abs(beta0[-1,]))

for(k in 1:length(alphalist)){
alpha=alphalist[k]
xt0 = rbind(cbind(x, z),
              cbind(-sqrt(alpha)*x, sqrt(alpha)*z))
  yt0 = c(y, rep(0, dim(x)[1]))

 coop = glmnet(xt0, yt0 , standardize=F)
 beta=predict(coop,s=lambda00, type="coef")
nz[ii,k,]=colSums(beta[-1,]!=0)
norm[ii,k,]=colSums(abs(beta[-1,]))
 }}

nz0m=colMeans(nz0)
nzm=apply(nz,c(2,3),mean)
norm0m=colMeans(norm0)
normm=apply(norm,c(2,3),mean)

pdf(file="df9.pdf",width=7.9,height=7)
par(mar=c(6,6,3,3))
color = c(scico(9,palette="bamako")[1:7],"#e1e01e", "#ffb21a")
plot(norm0m,nz0m,type="l", ylim=c(0,36),xlab="L1 norm",
     cex.lab = 1.5, cex.axis=1.2,
     ylab="Number of nonzero coefficients",lwd=3,xlim=c(0,2.7))

for(k in 1:length(alphalist)){
 #lines(normm[k,],nzm[k,],type="l",col=k+1)
 lines(normm[k,],nzm[k,],type="l",col=color[k], lwd=2)
 }
 legend('bottomright', legend=TeX(sprintf("$\\rho = %5.2f$", c(0,alphalist))),
        # lwd=c(2,rep(1,length(alphalist))),lty=1,col=1:(1+length(alphalist)),cex=.8)
         lwd=c(3,rep(2,length(alphalist))),
         lty=1,
         col=c('black',color[1:length(alphalist)]),
         cex=1.5,
         bty = "n")
dev.off()

library(RColorBrewer)
library(latex2exp)

general_filename = 'test_mse_3_pimp20_px500_pz500_sigma10_factorstr2_ustd1_factor2_SNR3_sim1.rds'

test_filename = general_filename
test_mse_df = readRDS(test_filename)

diff_filename = paste0("diff_", general_filename)
lasso_diff_df = readRDS(diff_filename)

support_filename = paste0("support", substring(general_filename,9))
support_df = readRDS(support_filename)

alpha_filename =  paste0("alpha", substring(general_filename,9))
alpha_df = readRDS(alpha_filename)

simN = 10
alphalist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)
alpha_axis = sapply(as.character(round(alphalist,1)),
         function(x) paste0("(alpha=",x,")"))
color = c(rep("white",4), "#B6D7A8", "#83B9E3")

plot_name = paste0("color_plot", substring(general_filename,9,nchar(general_filename)-3), "pdf")
pdf(file=plot_name, width=18, height=5)
par(mfrow=c(1,4), las=2, cex.main=2.33, cex.axis = 2.15, cex.lab=1.95, mar=c(10.3,6.6,6.1,0.9)) #, srt=45)
p1 = boxplot(unlist(cbind(test_mse_df[,1:5], test_mse_df[,15]))
        ~sort(rep(1:6, simN)), 
        #ylim=c(3000,5000),
        col=color,
        names=c("Separate X","Separate Z", "Early fusion", "Late fusion", "      Coop","Adap Coop"), 
        xlab="",ylab="", xaxt = "n")#, main="Test MSE")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Test MSE", line=1.8, cex.lab=2.75)

boxplot(unlist(cbind(lasso_diff_df[,1:5],lasso_diff_df[,15]))~
          sort(rep(1:6, simN)), 
       #ylim=c(-1000,1000),
       col=color,
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
        col=color,
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
        names=c(as.character(round(alphalist,1))),#as.character(round(alphalist,1))),
        ylim=c(0,10), cex.lab=2.0, 
        density=c(66,66),
        angle=c(36,36),
        col=c("#B6D7A8", "#83B9E3"))
legend("topleft", c("Coop","Adap Coop"), 
       cex=2.33, 
       pch=15,
       bty = "n",
       col=c("#B6D7A8", "#83B9E3"))
title(main=expression(bold(paste("Selected ", rho, " by CV (Counts)"))), line=2.2, cex.lab=2.75)
mtext(TeX(sprintf("$\\rho$")), side=1, las=1, line=6.3, cex=1.72)
dev.off()


library(RColorBrewer)
library(latex2exp)

general_filename = 'test_mse_123_pimp60_px1300_px2300_px3300_sigma25_factorstr2_factorx1_2_factorx2_2_factorx3_0_ustd1_factor2_SNR1_mode_comb_sim1.rds'

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
color = c(rep("white",5), "#B6D7A8", "#83B9E3")

plot_name = paste0("color_plot", substring(general_filename,9,nchar(general_filename)-3), "pdf")
pdf(file=plot_name, width=18, height=5)
par(mfrow=c(1,4), las=2, cex.main=2.33, cex.axis = 2.15, cex.lab=1.95, mar=c(10.3,6.6,6.1,0.9)) #, srt=45)
p1 = boxplot(unlist(cbind(test_mse_df[,1:7]))
        ~sort(rep(1:7, simN)), 
        ylim=c(700,1050),
        col=color,
        names=c("Separate X1","Separate X2", "Separate X3", "Early fusion", "Late fusion", "      Coop", "Adap Coop"), 
        xlab="",ylab="", xaxt = "n")#, main="Test MSE")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3])-9, p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Test MSE", line=1.8, cex.lab=2.75)

X1_diff_df = data.frame(test_mse_df - test_mse_df[,1])
boxplot(unlist(lasso_diff_df[,1:7])~ #cbind(X1_diff_df[,1:7])
          sort(rep(1:7, simN)), 
       ylim=c(-200,200),
       col=color,
       names=c("Separate X1","Separate X2", "Separate X3", "Early fusion", "Late fusion", "      Coop", "Adap Coop"), 
       xlab="",ylab="", main="", xaxt = "n")
mtext("Relative to early fusion", side=3, line=0.6, las=1, cex=1.6)
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3])-9, p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Test MSE Difference", line=3, cex.lab=2.75)
abline(h=0, lty=2)

boxplot(unlist(cbind(support_df[,1:7]))~
          sort(rep(1:7, simN)),
        #ylim=c(0,200),
        col=color,
        names=c("Separate X1","Separate X2", "Separate X3", "Early fusion", "Late fusion", "      Coop", "Adap Coop"), 
        xlab="", ylab="", main="", xaxt = "n")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.6, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3])-9, p1$names, srt = 45, xpd = T, cex=2.0)
title(main="Number of Features Selected", line=1.8, cex.lab=2.75)

counts = sapply(alphalist, function(x) sum(alpha_df$alpha_by_cv_no_pf == x))
counts_ada = sapply(alphalist, function(x) sum(alpha_df$alpha_by_cv == x))
barplot(rbind(counts,counts_ada),
        beside=TRUE,
        main="",
        #xlab="alpha", #ylab="Count",
        names=c(as.character(round(alphalist,1))),#as.character(round(alphalist,1))),
        ylim=c(0,10), cex.lab=2.0, 
        density=c(66,66),
        angle=c(36,36),
        col=c("#B6D7A8", "#83B9E3"))
#col=c("#a4aaa4", "black"))#"#89c589"))
legend("topleft", c("Coop","Adap Coop"), 
       cex=2.33, 
       #col=c("#a4aaa4", "black"), 
       #density=c(66,66),
       pch=15,
       #angle=c(36,36),
       bty = "n",
       col=c("#B6D7A8", "#83B9E3"))

title(main=expression(bold(paste("Selected ", rho, " by CV (Counts)"))), line=2.2, cex.lab=2.75)
mtext(TeX(sprintf("$\\rho$")), side=1, las=1, line=6.3, cex=1.72)
dev.off()


library(RColorBrewer)
library(latex2exp)

setwd('./../snr_6')

res = readRDS("error_image_gene.rds")
alpha_selected = readRDS("alpha_image_gene.rds")
nsim = 10

out = rbind(colMeans(res), apply(res,2,sd)/sqrt(nsim))
#dimnames(out) = list(c("mean","se"),
#                     c("image","prot","late","coop",paste("coop",as.character(alphalist),sep="")))
color = c(rep("white",3), "#B6D7A8")

#title still hard-coded right now
pdf(file="color_6_test_misclassification_images_genes.pdf", width=7.9, height=5)
par(mfrow=c(1,2), las=2, cex.main=1.39, cex.axis = 1.2, mar=c(6.9,3.3,4.2,1.5))
p1=boxplot(c(res[,1:4])
        ~sort(rep(1:4, nsim)), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist)),
        #        "#ff8080",rep("#BCDEBC",length(alphalist))),
        names=c("Only Images"," Only Omics", "  Late fusion", "           Coop"), #,alpha_axis), 
        col=color,
        xlab="",ylab="", main="Test Error", xaxt = "n")
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.39, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=1.3)

images_diff = res - res[,3]
p2 = boxplot(c(images_diff[,1:4])
        ~sort(rep(1:4, nsim)), 
        #col = c(rep("white",4),"#ff8080", rep("#BCDEBC",length(alphalist)),
        #        "#ff8080",rep("#BCDEBC",length(alphalist))),
        names=c("Only Images"," Only Omics", "  Late fusion", "          Coop"), #, alpha_axis), 
        col=color,
        xlab="",ylab="", main="", xaxt = "n")
title(main="Test Error Difference", line=2.2, cex.lab=2.75)
mtext("Relative to late fusion", side=3, line=0.5, las=1, cex=1.27)
tick = seq_along(p1$names)
axis(1, at = tick, labels = F)
text(tick-0.39, par("usr")[3]-0.2*(par("usr")[4]-par("usr")[3]), p1$names, srt = 45, xpd = T, cex=1.2)
abline(h=0, lty=2)
dev.off()


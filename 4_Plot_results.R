## R-Script - Plot results
## author: Javier Lopatin
## mail:javierlopatin@gmail.com
## Manuscript: Comparing Generalized Linear Models and random forest to model vascular plant species richness using LiDAR data in a natural forest in central Chile
## last changes: 12/11/2015

##############################################
### Importance of variables. Figure 2  ###
##############################################

library(lattice)

# load the table with the variable importances (gini index and hierarchical partitioning for RF and GLM respectively)
# the table must have 4 columns: 1. variable names, 2. the importance, 3. the layer (total, tree, shrub or herb), and 4. the model (GLM or RF)
imp<- read.table("importance.csv", header=T, sep=",", dec=".")
str(imp)

pdf(file = "Figures/Fig .pdf", width=10, height=5)
dotplot(factor(Variables, levels = rev(sort(unique(Variables))), ordered = T)~Importance | Layer, data=imp, cex=1.3, xlab="Importance (%)", groups=Model, layout = c(4,1),
        main="",index.cond=list(c(3,4,2,1)), box.width=2, origin=0, par.settings = simpleTheme(cex=1.5, pch=c(16,17), col=c("blue", "red")), auto.key=list(space="right", columns=1)) 
dev.off()


###################################################
### Distribution of model accuracies. Figure 3  ###
###################################################

library(beanplot)

# RF VS GLM
pdf(file = "Figures/RF-NB beanplot.pdf", width=10, height=4)
par(mfrow=c(1,3),lend = 1, mai = c(0.8, 0.6, 0.5, 0.1))
beanplot( unlist(r2.rf), unlist(r2.nb),  unlist(r2.rf.A), unlist(r2.nb.A),  unlist(r2.rf.AR), unlist(r2.nb.AR),  unlist(r2.rf.H), unlist(r2.nb.H), 
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", ll = 0, side = "b", log="",
          main = "Squared Pearson's correlation coefficient", names=c("Total", "Tree", "Shrub", "Herb"), ylab = expression(r^2), ylim = c(0,1), yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)
beanplot( unlist(Nrmse.rf), unlist(Nrmse.nb),  unlist(Nrmse.rf.A), unlist(Nrmse.nb.A),  unlist(Nrmse.rf.AR), unlist(Nrmse.nb.AR),  unlist(Nrmse.rf.H), unlist(Nrmse.nb.H),
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", ll = 0, side="b",log="" , ylim=c(10,40),
          main = "Normalized root mean square error", names=c("Total", "Tree", "Shrub", "Herb"), ylab = "nRMSE (%)", yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
beanplot( unlist(bias.rf), unlist(bias.nb),  unlist(bias.rf.A), unlist(bias.nb.A),  unlist(bias.rf.AR), unlist(bias.nb.AR),  unlist(bias.rf.H), unlist(bias.nb.H), 
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", ll = 0, side = "b", log="" ,
          main = "Bias", names=c("Total", "Tree", "Shrub", "Herb"), ylab = "Bias", yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)
legend("topleft", legend=c("RF", "GLM"), fill=c("black", "gray"), bty="n", cex=1.3)
dev.off()


####################################################################
### Scatter plots of observed versus predicted values. Figure 4  ###
####################################################################

par(oma = c(4.5, 1, 1, 1))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
par(mfrow=c(2,2),lend = 1, mai = c(0.5, 0.5, 0.5, 0.1))
##Total Richness
par(mar=c(5, 5, 3, 3))
Mytitle = "Total richness fitted model"
Pred1<- subset(rf, Obs < 32)$Pred
Obs1<-  subset(rf, Obs < 32)$Obs
Pred2<- subset(nb, Obs < 32)$Pred
Obs2<-  subset(nb, Obs < 32)$Obs
MyXlab <- expression(paste("Observed (N° spp)"))
MyYlab <- expression(paste( "Predicted (N° spp)"))
plot((Obs1 - 0.25), Pred1, col="black", xlim=c(0,40), ylim=c(0,40), xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=0.5, cex.lab=1.1, cex.axis=1.1, main = Mytitle, las= 1)
points((Obs2 + 0.25), Pred2, pch=17,col="gray", cex=0.5)
abline(0, 1, lty=1, lwd=2)
# RF
lm1 = lm(Pred1 ~ Obs1-1)
abline(lm1, lty=2, lwd=2)
# GLM
lm2 = lm(Pred2 ~ Obs2-1)
abline(lm2, lty=3, lwd=2)

## Tree 
par(mar=c(5, 5, 3, 3))
Mytitle = "Tree richness fitted model"
Pred1<- subset(rf.A, Obs < 16)$Pred
Obs1<-  subset(rf.A, Obs < 16)$Obs
Pred2<- subset(nb.A, Obs < 16)$Pred
Obs2<-  subset(nb.A, Obs < 16)$Obs
plot((Obs1 - 0.25), Pred1, xlim=c(0,20), ylim=c(0,20), xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=0.8, cex.lab=1.1, cex.axis=1.1, main = Mytitle, las= 1)
points((Obs2 + 0.25), Pred2, pch=17, col="grey", cex=0.8)
abline(0, 1, lty=1, lwd=2)
# RF
abline(lm1, lty=2, lwd=2)
# GLM
lm2 = lm(Pred2 ~ Obs2-1)
abline(lm2, lty=3, lwd=2)

## Shrub
par(mar=c(5, 5, 3, 3))
Mytitle = "Shrub richness fitted model"
Pred1<- rf.AR$Pred
Obs1<-  rf.AR$Obs
Pred2<- nb.AR$Pred
Obs2<-  nb.AR$Obs
plot((Obs1 - 0.25), Pred1, xlim=c(0,25), ylim=c(0,25), xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=0.8, cex.lab=1.1, cex.axis=1.1, main = Mytitle, las= 1)
points((Obs2 + 0.25), Pred2, pch=17, col="grey", cex=0.8)
abline(0, 1, lty=1, lwd=2)
# RF
lm1 = lm(Pred1 ~ Obs1-1)
abline(lm1, lty=2, lwd=2)
# GLM
lm2 = lm(Pred2 ~ Obs2-1)
abline(lm2, lty=3, lwd=2)

## Herb
par(mar=c(5, 5, 3, 3))
Mytitle = "Herb richness fitted model"
Pred1<- rf.H$Pred
Obs1<-  rf.H$Obs
Pred2<- nb.H$Pred
Obs2<-  nb.H$Obs
plot((Obs1 - 0.25), Pred1, xlim=c(0,20), ylim=c(0,20), xlab = MyXlab, ylab = MyYlab, pch=16, pty="s", cex=0.8, cex.lab=1.1, cex.axis=1.1, main = Mytitle, las= 1)
points((Obs2 + 0.25), Pred2, pch=17, col="grey", cex=0.8)
abline(0, 1, lty=1, lwd=2)
# RF
lm1 = lm(Pred1 ~ Obs1-1)
abline(lm1, lty=2, lwd=2)
# GLM
lm2 = lm(Pred2 ~ Obs2-1)
abline(lm2, lty=3, lwd=2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = T)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-0.08, -0.80,inset = 0 , xpd = T,  text.width=0.3,legend=c("RF","", "GLM",""), col=c("black","black","black", "gray"), pch=c(NA,16,NA,17), lty=c(2,NA,3,NA), lwd=c(2.5,2.5), horiz = F, bty = "n", cex=0.8) 


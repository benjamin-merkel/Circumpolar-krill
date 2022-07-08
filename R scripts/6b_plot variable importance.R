

es_vi <- readRDS("data/E.superba_circumpolar_GBM_100_estimated variable importance_Dec-Mar.rds")
ec_vi <- readRDS("data/E.crystallorophias_circumpolar_GBM_100_estimated variable importance_Dec-Mar.rds")
es_vi2<- apply(es_vi,2,mean)
ec_vi2<- apply(ec_vi,2,mean)


png(paste0("figures/both species variable importance.png"), res = 800, width=30, height = 17, units="cm")
opar <- par(mfrow=c(1,2),mar=c(4,8,1,1))

boxplot(es_vi[,order(es_vi2)], horizontal=T, las=1,
        col = grey(1), border = "skyblue3", lty = 1, lwd = 2,
        ylab = "", xlab = "variable importance")
points(es_vi2[order(es_vi2)],1:length(es_vi2),cex=1.5,pch=19,col="darkred")
mtext('a)', side=2, line= 5.5, at=14.5, cex=2, las = 1)


boxplot(ec_vi[,order(ec_vi2)], horizontal=T, las=1,
        col = grey(1), border = "skyblue3", lty = 1, lwd = 2,
        ylab = "", xlab = "variable importance")
points(ec_vi2[order(ec_vi2)],1:length(ec_vi2),cex=1.5,pch=19,col="darkred")
mtext('b)', side=2, line= 5.5, at=14.5, cex=2, las = 1)

par(opar)
dev.off()



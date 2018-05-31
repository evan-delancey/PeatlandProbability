

#Load libraries
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)
library(caret)
library(snow)
library(RPyGeo)
library(dismo)
library(gbm)
library(ggthemes)


#enter number of iteration of subsampling
iter <- 40

#set location of a temporary raster dump
#this can take up 100-300BG per run but is deleted after
rasterOptions(maxmemory = 1e+09,tmpdir = "J:/RtmpRasterDump")

tifs <- c("ARI.tif", "NDPOL.tif", "PC1.tif", "REIP.tif", "TPI.tif", "TWI.tif")
fls <- tifs


Validation <- read.csv("J:/LandCover/OrganicWetlandProbability/Tests/DFs/Validation.csv")
Validation <- Validation[,-c(4,8)]


#######################################################################################
#1)Use all varibles to get ideal learning rate in terms of accuracy, AUC, and deviance#
#######################################################################################


distance <- c("100", "250", "375", "500", "1000", "1500")
Npts <- c(91347, 14617, 6497, 3655, 915, 407)

Owetland <- raster("J:/LandCover/OrganicWetlandProbability/Training/Owetland.tif")

AUC <- matrix(nrow = iter, ncol = length(distance))
dev <- matrix(nrow = iter, ncol = length(distance))
accuracy <- vector()

for (i in 2:length(distance)){
	prediction <- matrix(nrow = length(Validation[,1]), ncol = iter)
	setwd(paste0("J:/LandCover/OrganicWetlandProbability/Training/ModelPts/Data", distance[i]))
	for (n in 1:iter){
	
		dat <- read.csv(paste0("d", n, ".csv"))
		dat <- dat[,-c(4,8)]
			
		#build model
		fit <- gbm.step(dat, 1:length(fls), length(fls)+1, family = "bernoulli", tree.complexity = 8,
				learning.rate = 0.005, bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
			
		#model stats
		AUC[n,i] <- fit$cv.statistics$discrimination.mean
		dev[n,i] <- (fit$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
		prediction[,n] <- predict.gbm(fit ,Validation[,1:length(fls)], n.trees = fit$gbm.call$best.trees, type = "response")
		print(n)
	}
	mean.prediction <- rowMeans(prediction)
	binary.prediction <- ifelse(mean.prediction > 0.5, 1, 0)
	evaluation <- cbind(Validation[,7], binary.prediction)
	correct <- evaluation[,1] - evaluation[,2]
	accuracy[i] <- length(correct[correct==0])/length(correct)
	
	
	print(paste0("done distance ", i, "out of ", length(distance)))
}

setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
write.csv(AUC, "AUC_Npts.csv", row.names=F)
write.csv(dev, "dev_Npts.csv", row.names=F)
write.csv(accuracy, "accuracy_Npts.csv", row.names=F)



#plot AUC
AUC.mean <- colMeans(AUC)
AUC.sd <- apply(AUC, 2, sd) 
AUC.nsd <- AUC.mean - AUC.sd
AUC.psd <- AUC.mean + AUC.sd
d <- data.frame (
	x = Npts,
	mean = AUC.mean,
	psd = AUC.psd,
	msd = AUC.nsd
) 

setwd("J:/LandCover/OrganicWetlandProbability/Tests/Figures")

tiff(paste0("AUC_Npts.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Number of training points") + ylab("AUROC") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot dev
dev.mean <- colMeans(dev)
dev.sd <- apply(dev, 2, sd) 
dev.nsd <- dev.mean - dev.sd
dev.psd <- dev.mean + dev.sd
d <- data.frame (
	x = Npts,
	mean = dev.mean,
	psd = dev.psd,
	msd = dev.nsd
) 

tiff(paste0("dev_Npts.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Number of training points") + ylab("Explained deviance") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot accuracy
d <- data.frame (
	x = Npts,
	mean = accuracy
) 

tiff(paste0("accuracy_Npts.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_line(data=d, aes(x=x, y=mean), size = 1.5, colour="forestgreen") + 
	xlab("Number of training points") + ylab("Accuracy") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


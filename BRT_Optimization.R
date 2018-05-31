

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

tifs <- c("ARI.tif", "NDPOL.tif", "PC1.tif", "PC2.tif", "REIP.tif", "TPI.tif", "TWI.tif", "VH.tif")
fls <- tifs


Validation <- read.csv("J:/LandCover/OrganicWetlandProbability/Tests/DFs/Validation.csv")


#######################################################################################
#1)Use all varibles to get ideal learning rate in terms of accuracy, AUC, and deviance#
#######################################################################################


LR <- seq(0.001, 0.012, by=0.001)
LR <- c(0.01, 0.011, 0.012)

Owetland <- raster("J:/LandCover/OrganicWetlandProbability/Training/Owetland.tif")

AUC.LR <- matrix(nrow = iter, ncol = length(LR))
dev.LR <- matrix(nrow = iter, ncol = length(LR))
accuracy.LR <- vector()

for (i in 1:length(LR)){
	prediction <- matrix(nrow = length(Validation[,1]), ncol = iter)
	setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
	for (n in 1:iter){
	
		dat <- read.csv(paste0("d", n, ".csv"))
			
		#build model
		fit <- gbm.step(dat, 1:length(fls), length(fls)+1, family = "bernoulli", tree.complexity = 6,
				learning.rate = LR[i], bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
			
		#model stats
		AUC.LR[n,i] <- fit$cv.statistics$discrimination.mean
		dev.LR[n,i] <- (fit$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
		prediction[,n] <- predict.gbm(fit ,Validation[,1:length(fls)], n.trees = fit$gbm.call$best.trees, type = "response")
	}
	mean.prediction <- rowMeans(prediction)
	binary.prediction <- ifelse(mean.prediction > 0.5, 1, 0)
	evaluation <- cbind(Validation[,9], binary.prediction)
	correct <- evaluation[,1] - evaluation[,2]
	accuracy.LR[n] <- length(correct[correct==0])/length(correct)
	
	
	print(paste0("done LR ", i, "out of ", length(LR)))
}

setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
write.csv(AUC.LR, "AUC_LR.csv", row.names=F)
write.csv(dev.LR, "dev_LR.csv", row.names=F)
write.csv(accuracy.LR, "accuracy_LR.csv", row.names=F)



#plot AUC
AUC.mean <- colMeans(AUC.LR)
AUC.sd <- apply(AUC.LR, 2, sd) 
AUC.nsd <- AUC.mean - AUC.sd
AUC.psd <- AUC.mean + AUC.sd
d <- data.frame (
	x = LR,
	mean = AUC.mean,
	psd = AUC.psd,
	msd = AUC.nsd
) 

setwd("J:/LandCover/OrganicWetlandProbability/Tests/Figures")

tiff(paste0("AUC_LR.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Learning rate") + ylab("AUROC") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot dev
dev.mean <- colMeans(dev.LR)
dev.sd <- apply(dev.LR, 2, sd) 
dev.nsd <- dev.mean - dev.sd
dev.psd <- dev.mean + dev.sd
d <- data.frame (
	x = LR,
	mean = dev.mean,
	psd = dev.psd,
	msd = dev.nsd
) 

tiff(paste0("dev_LR.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Learning rate") + ylab("Explained deviance") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot accuracy
d <- data.frame (
	x = LR,
	mean = accuracy,
) 

tiff(paste0("accuracy_LR.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Learning rate") + ylab("Accuracy") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()






#########################################################################################
#2)Use all varibles to get ideal tree complexity in terms of accuracy, AUC, and deviance#
#########################################################################################


TC <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)

Owetland <- raster("J:/LandCover/OrganicWetlandProbability/Training/Owetland.tif")

AUC.TC <- matrix(nrow = iter, ncol = length(TC))
dev.TC <- matrix(nrow = iter, ncol = length(TC))
accuracy.TC <- matrix(nrow = iter, ncol = length(TC))



for (i in 1:length(TC)){
	
	setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
	for (n in 1:iter){
	
		dat <- read.csv(paste0("d", n, ".csv"))
			
		#build model
		fit <- gbm.step(dat, 1:length(fls), length(fls)+1, family = "bernoulli", tree.complexity = TC[i],
				learning.rate = 0.0085, bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
			
		#model stats
		AUC.TC[n,i] <- fit$cv.statistics$discrimination.mean
		dev.TC[n,i] <- (fit$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
		prediction <- predict.gbm(fit ,Validation[,1:length(fls)], n.trees = fit$gbm.call$best.trees, type = "response")
		binary.prediction <- ifelse(prediction > 0.5, 1, 0)
		evaluation <- cbind(Validation[,9], binary.prediction)
		correct <- evaluation[,1] - evaluation[,2]
		accuracy.TC[n,i] <- length(correct[correct==0])/length(correct)
	
	}
	
	print(paste0("done TC ", i, "out of ", length(TC)))
}

setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
write.csv(AUC.TC, "AUC_TC.csv", row.names=F)
write.csv(dev.TC, "dev_TC.csv", row.names=F)
write.csv(accuracy.TC, "accuracy_TC.csv", row.names=F)



#plot AUC
AUC.mean <- colMeans(AUC.TC)
AUC.sd <- apply(AUC.TC, 2, sd) 
AUC.nsd <- AUC.mean - AUC.sd
AUC.psd <- AUC.mean + AUC.sd
d <- data.frame (
	x = TC,
	mean = AUC.mean,
	psd = AUC.psd,
	msd = AUC.nsd
) 

setwd("J:/LandCover/OrganicWetlandProbability/Tests/Figures")

tiff(paste0("AUC_TC.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Tree complexity") + ylab("AUROC") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot dev
dev.mean <- colMeans(dev.TC)
dev.sd <- apply(dev.TC, 2, sd) 
dev.nsd <- dev.mean - dev.sd
dev.psd <- dev.mean + dev.sd
d <- data.frame (
	x = TC,
	mean = dev.mean,
	psd = dev.psd,
	msd = dev.nsd
) 

tiff(paste0("dev_TC.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Tree complexity") + ylab("Explained deviance") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot accuracy
accuracy.mean <- colMeans(accuracy.TC)
accuracy.sd <- apply(accuracy.TC, 2, sd) 
accuracy.nsd <- accuracy.mean - accuracy.sd
accuracy.psd <- accuracy.mean + accuracy.sd
d <- data.frame (
	x = TC,
	mean = accuracy.mean,
	psd = accuracy.psd,
	msd = accuracy.nsd
) 

tiff(paste0("accuracy_TC.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Tree complexity") + ylab("Accuracy") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()







#########################################################################################
#3)Use all varibles to get ideal number of variables in terms of accuracy, AUC, and deviance#
#########################################################################################




VARs <- list(c(5,3), c(5,3,7), c(5,3,7,6), c(5,3,7,6,2), c(5,3,7,6,2,1), c(5,3,7,6,2,1,8), c(5,3,7,6,2,1,8,4))
VAR.num <- c(2,3,4,5,6,7,8)


AUC.VARs <- matrix(nrow = iter, ncol = length(VARs))
dev.VARs <- matrix(nrow = iter, ncol = length(VARs))
accuracy.VARs <- matrix(nrow = iter, ncol = length(VARs))

for (i in 1:length(VARs)){
	
	setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
	for (n in 1:iter){
	
		dat <- read.csv(paste0("d", n, ".csv"))
			
		#build model
		fit <- gbm.step(dat, VARs[[i]], length(fls)+1, family = "bernoulli", tree.complexity = 6,
				learning.rate = 0.005, bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
			
		#model stats
		AUC.VARs[n,i] <- fit$cv.statistics$discrimination.mean
		dev.VARs[n,i] <- (fit$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
		prediction <- predict.gbm(fit ,Validation[,VARs[[i]]], n.trees = fit$gbm.call$best.trees, type = "response")
		binary.prediction <- ifelse(prediction > 0.5, 1, 0)
		evaluation <- cbind(Validation[,9], binary.prediction)
		correct <- evaluation[,1] - evaluation[,2]
		accuracy.VARs[n,i] <- length(correct[correct==0])/length(correct)
	
	}
	
	print(paste0("done VARs ", i, "out of ", length(VARs)))
}

setwd("J:/LandCover/OrganicWetlandProbability/Tests/DFs")
write.csv(AUC.VARs, "AUC_VARs.csv", row.names=F)
write.csv(dev.VARs, "dev_VARs.csv", row.names=F)
write.csv(accuracy.VARs, "accuracy_VARs.csv", row.names=F)



#plot AUC
AUC.mean <- colMeans(AUC.VARs)
AUC.sd <- apply(AUC.VARs, 2, sd) 
AUC.nsd <- AUC.mean - AUC.sd
AUC.psd <- AUC.mean + AUC.sd
d <- data.frame (
	x = VAR.num,
	mean = AUC.mean,
	psd = AUC.psd,
	msd = AUC.nsd
) 

setwd("J:/LandCover/OrganicWetlandProbability/Tests/Figures")

tiff(paste0("AUC_VARs.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Number of variables") + ylab("AUROC") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot dev
dev.mean <- colMeans(dev.VARs)
dev.sd <- apply(dev.VARs, 2, sd) 
dev.nsd <- dev.mean - dev.sd
dev.psd <- dev.mean + dev.sd
d <- data.frame (
	x = VAR.num,
	mean = dev.mean,
	psd = dev.psd,
	msd = dev.nsd
) 

tiff(paste0("dev_VARs.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Number of variables") + ylab("Explained deviance") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()


#plot accuracy
accuracy.mean <- colMeans(accuracy.VARs)
accuracy.sd <- apply(accuracy.VARs, 2, sd) 
accuracy.nsd <- accuracy.mean - accuracy.sd
accuracy.psd <- accuracy.mean + accuracy.sd
d <- data.frame (
	x = VAR.num,
	mean = accuracy.mean,
	psd = accuracy.psd,
	msd = accuracy.nsd
) 

tiff(paste0("accuracy_VARs.tiff"), width = 1000, height = 700)
ggplot(d, aes(x=x, y=mean)) + 
	theme_minimal()+ 
	geom_ribbon(aes(ymin = msd, ymax = psd), fill="forestgreen", alpha=0.25) +
	geom_line(colour="forestgreen", size=1.9) +
	xlab("Number of variables") + ylab("Accuracy") +
	theme(axis.title.x = element_text(size=34), axis.title.y = element_text(size=34),axis.text = element_text(size=30))
dev.off()






	
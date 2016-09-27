Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)

# PBDB data
setwd(paste(Path_base, "fossilSalamanders/pbdb", sep=""))
fosRTT <- read.csv("raw_curve_data.csv", stringsAsFactors=F)

# Extant-only data
setwd(paste(Path_base, "bamm/extant_only", sep=""))
exTree <- read.tree("modernSal.tre")
exMCMC <- read.csv("ex_mcmc_out.txt", stringsAsFactors=F)
exEdata <- getEventData(exTree, "ex_event_data.txt", burnin=0.5)
maxT <- max(nodeHeights(exTree))

exPlot <- plot(exEdata, breaks="quantile")
addBAMMlegend(exPlot)
exRTT <- getRateThroughTimeMatrix(exEdata)

# fossilBAMM runs
### Post-generation addition of files
fbEdata <- list()
fbEdata_hi <- list()
fbRTT <- list()
fbRTT_hi <- list()
fmax <- c()
for (count in 1:10)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	fbTree <- read.tree("fossilTree.tre")
	fbEdata[[count]] <- getEventData(fbTree, "lo_event_data.txt", burnin=0.5)
	fbRTT[[count]] <- getRateThroughTimeMatrix(fbEdata[[count]])
	
	fbEdata_hi[[count]] <- getEventData(fbTree, "hi_event_data.txt", burnin=0.5)
	fbRTT_hi[[count]] <- getRateThroughTimeMatrix(fbEdata_hi[[count]])
	
	fmax[count] <- max(nodeHeights(fbTree))
}

Alpha <- 1
exCol <- rgb(217/255, 95/255, 2/255, Alpha)
hiCol <- rgb(27/255, 158/255, 119/255, Alpha)
loCol <- rgb(117/255, 112/255, 179/255, Alpha)

setwd(paste(Path_base, "figures", sep=""))
pdf("spex.pdf", height=5, width=10)
par(mfcol=c(1,2), mar=c(4,4,1,1), mgp=c(1.5,0.3,0), tck=-0.005, bty="n", las=1)
# Speciation plot
tDiff <- max(fmax) - maxT
plot(exRTT$times + tDiff, apply(exRTT$lambda,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="", ylab="", col=exCol, lwd=1.5)
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$lambda, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$lambda, 2, mean), col=loCol))
Xax <- c(-100, 0, 50, 100, 150, 204, 240, 270)
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))


# Extinction plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="", ylab="", col=exCol, lwd=1.5)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$mu, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
dev.off()

pdf("derivedPars.pdf", height=5, width=10)
par(mfcol=c(1,2), mar=c(4,4,1,1), mgp=c(1.5,0.3,0), tck=-0.005, bty="n", las=1)
# Net div plot
plot(exRTT$times + tDiff, apply(exRTT$lambda,2,mean) - apply(exRTT$mu,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="", ylab="", col=exCol, lwd=1.5)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$lambda, 2, mean) - apply(fbRTT_hi[[x]]$mu, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$lambda, 2, mean) - apply(fbRTT[[x]]$mu, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))

# Rel ext plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean) / apply(exRTT$lambda,2,mean), type='l', ylim=c(0,1), xlim=c(0,max(fmax)), axes=F, xlab="", ylab="", lwd=1.5, col=exCol)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean) / apply(fbRTT_hi[[x]]$lambda, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x],  apply(fbRTT[[x]]$mu, 2, mean) / apply(fbRTT[[x]]$lambda, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()
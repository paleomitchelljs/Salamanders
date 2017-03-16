Path_base <- "~/Documents/Salamanders/"

addFamily <- function(Taxa, Family, Cex=2, centering=T, Top=T)	{
	Names <- unique(c(tree$tip.label[grep(Taxa[1], tree$tip.label)], tree$tip.label[grep(Taxa[2], tree$tip.label)], tree$tip.label[grep(Taxa[3], tree$tip.label)]))
	allNames <- na.omit(tree$tip.label[getDescendants(tree, findMRCA(tree, tips=Names))])
	if (centering == T)		{
		Center <- allNames[round(mean(1:length(allNames)))]
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=-0.05, cex=Cex)
	}
	else	 if (centering == F){
		if (Top == T)	{
			Center <- rev(allNames)[1]
		}
		else	{
			Center <- allNames[1]
		}
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=-0.05, cex=Cex)
	}
}

library(phytools)
library(BAMMtools)
library(coda)

# PBDB data
setwd(paste(Path_base, "fossilSalamanders/pbdb", sep=""))
fosRTT <- read.csv("raw_curve_data.csv", stringsAsFactors=F)
fosNetDiv <- fosRTT$X3T_lambda - fosRTT$X3T_mu
fosRelExt <- fosRTT$X3T_mu/fosRTT$X3T_lambda

# Extant-only data
setwd(paste(Path_base, "bamm/extant_only", sep=""))
exTree <- ladderize(read.tree("modernSal.tre"), right=F)
exMCMC <- read.csv("ex_mcmc_out.txt", stringsAsFactors=F)
exEdata <- getEventData(exTree, "ex_event_data.txt", burnin=0.5)
maxT <- max(nodeHeights(exTree))

setwd(paste(Path_base, "figures", sep=""))
pdf("extant_sp.pdf", height=10, width=10)
addFamily <- function(Taxa, Family, Cex=2, centering=T, Top=T)	{
	Names <- unique(c(tree$tip.label[grep(Taxa[1], tree$tip.label)], tree$tip.label[grep(Taxa[2], tree$tip.label)], tree$tip.label[grep(Taxa[3], tree$tip.label)]))
	allNames <- na.omit(tree$tip.label[getDescendants(tree, findMRCA(tree, tips=Names))])
	if (centering == T)		{
		Center <- allNames[round(mean(1:length(allNames)))]
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=-0.05, cex=Cex)
	}
	else	 if (centering == F){
		if (Top == T)	{
			Center <- rev(allNames)[1]
		}
		else	{
			Center <- allNames[1]
		}
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=-0.05, cex=Cex)
	}
}
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
exPlot <- plot(exEdata, breaks="quantile")
addBAMMlegend(exPlot)

# Add family labels
tree <- as.phylo(exEdata)
#addFamily(c("Rhyacotriton", "Rhyacotriton", "Rhyacotriton"), "Rhyacotritonidae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Proteus", "Necturus", "Necturus"), "Proteidae", centering=F, Top=F)
dev.off()

setwd(paste(Path_base, "check_tree_files", sep=""))
pdf("bamm_ex_output.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,2))
exPlot <- plot(exEdata, breaks="quantile")
addBAMMlegend(exPlot)
tiplabels(text=exEdata$tip.label, cex=0.15, adj=-0.1, frame="none", bg=rgb(0,0,0,0))
dev.off()


# Extant-no pleth data
setwd(paste(Path_base, "bamm/extant_only", sep=""))
npTree <- ladderize(read.tree("no_pleth.tre"), right=F)
npMCMC <- read.csv("np_mcmc_out.txt", stringsAsFactors=F)
npEdata <- getEventData(npTree, "np_event_data.txt", burnin=0.5)
maxT <- max(nodeHeights(npTree))

setwd(paste(Path_base, "figures", sep=""))
pdf("np_sp.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
npPlot <- plot(npEdata, breaks="quantile")
addBAMMlegend(npPlot)

# Add family labels
tree <- as.phylo(npEdata)
addFamily(c("Rhyacotriton", "Rhyacotriton", "Rhyacotriton"), "Rhyacotritonidae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
addFamily(c("Proteus", "Necturus", "Necturus"), "Proteidae", centering=F, Top=F)
dev.off()

# Extant-only rates through time
exRTT <- getRateThroughTimeMatrix(exEdata)
npRTT <- getRateThroughTimeMatrix(npEdata)

# fossilBAMM runs
### Post-generation addition of files
fbEdata <- list()
fbEdata_all <- list()
fbRTT <- list()
fbRTT_all <- list()
fmax <- c()
fullmax <- c()
llES <- c()
nsES <- c()
Nshifts <- c()
meanNshifts <- c()
presRate <- c()
Nshifts_all <- c()
meanNshifts_all <- c()

for (count in 1:70)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	fbTree <- ladderize(read.tree("fossilTree.tre"), right=F)
	MCMC_dat <- read.csv("base_mcmc_out.txt")
	MCMC_dat <- MCMC_dat[floor(0.1*nrow(MCMC_dat)):nrow(MCMC_dat),]
	
	fbEdata[[count]] <- getEventData(fbTree, "base_event_data.txt", burnin=0.1, nsamples=200)
	fbRTT[[count]] <- getRateThroughTimeMatrix(fbEdata[[count]])
	llES[count] <- effectiveSize(MCMC_dat$logLik)
	nsES[count] <- effectiveSize(MCMC_dat$N_shifts)
	BFmat <- computeBayesFactors(MCMC_dat, 200, burnin=0)
	Nshifts[count] <- stepBF(BFmat, step.size=20, expectedNumberOfShifts=200, inputType="matrix")[1]
	meanNshifts[count] <- mean(MCMC_dat$N_shifts)
	presRate[count] <- mean(MCMC_dat$preservationRate)
	
	fullTree <- ladderize(read.tree("fullTree.tre"), right=F)
	fbEdata_all[[count]] <- getEventData(fullTree, "full_event_data.txt", burnin=0.1, nsamples=200)
	fbRTT_all[[count]] <- getRateThroughTimeMatrix(fbEdata_all[[count]])

	MCMC_dat <- read.csv("full_mcmc_out.txt")
	MCMC_dat <- MCMC_dat[floor(0.1*nrow(MCMC_dat)):nrow(MCMC_dat),]

	BFmat <- computeBayesFactors(MCMC_dat, 200, burnin=0)
	Nshifts_all[count] <- stepBF(BFmat, step.size=20, expectedNumberOfShifts=200, inputType="matrix")[1]
	meanNshifts_all[count] <- mean(MCMC_dat$N_shifts)	

}


setwd(paste(Path_base, "figures", sep=""))
pdf("fossil_topologies.pdf", height=8, width=8)
for (i in 1:length(fbEdata))	{
	plot(as.phylo(fbEdata[[i]]), cex=0.25, label.offset=0.02)
}
dev.off()

setwd(paste(Path_base, "figures", sep=""))
pdf("fossil_sp.pdf", height=10, width=10)
fosEdata <- fbEdata[[1]]
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
fosPlot <- plot(fosEdata, colorbreaks=exPlot$colorbreaks)
#addBAMMlegend(fosPlot, location="bottomleft")

# Add family labels
tree <- as.phylo(fosEdata)
addFamily(c("Batrachosauroides", "Opisthotriton", "Peratosauroides"), "Batrachosauroididae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
addFamily(c("Scapherpeton", "Piceoerpeton", "Piceoerpeton"), "Scapherpetonidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
dev.off()

setwd(paste(Path_base, "figures", sep=""))
pdf("fossil_sp_alt.pdf", height=10, width=10)
fosEdata <- fbEdata[[35]]
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
fosPlot <- plot(fosEdata, colorbreaks=exPlot$colorbreaks)
#addBAMMlegend(fosPlot, location="bottomleft")

# Add family labels
tree <- as.phylo(fosEdata)
addFamily(c("Batrachosauroides", "Opisthotriton", "Peratosauroides"), "Batrachosauroididae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
addFamily(c("Scapherpeton", "Piceoerpeton", "Piceoerpeton"), "Scapherpetonidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
dev.off()

Alpha <- 1
npCol <- rgb(166/255, 206/ 255, 227/255, Alpha) # extant w/o pleth
exCol <- rgb(31/255, 120/255, 180/255, Alpha) # extant w/ pleth
npfCol <- rgb(178/255, 223/255, 138/255, Alpha) # fossils w/o pleth
allCol <- rgb(51/255, 160/255, 44/255, Alpha) # fossils w/ pleth

LWD <- 2
Xax <- c(0, 50, 100, 150, 200, 250)
startT <- 225

setwd(paste(Path_base, "figures", sep=""))
pdf("spex.pdf", height=6.5, width=6.5)

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.2,0.3,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)

# Speciation plot
plot(max(exRTT$times) - exRTT$times, apply(exRTT$lambda,2,mean), type='n', ylim=c(0,0.1), xlim=c(startT, 0), axes=F, xlab="age (Ma)", ylab=expression(paste("speciation rate (", lambda, ")", sep="")), col=exCol, lwd=LWD)

silent <- sapply(1:length(fbRTT_all), function(x) lines(max(fbRTT_all[[x]]$times) - fbRTT_all[[x]]$times, apply(fbRTT_all[[x]]$lambda, 2, mean), col=allCol))
silent <- sapply(1:length(fbRTT), function(x) lines(max(fbRTT[[x]]$time) - fbRTT[[x]]$times, apply(fbRTT[[x]]$lambda, 2, mean), col=npfCol))

lines(max(npRTT$times) - npRTT$times, apply(npRTT$lambda, 2, mean), col=npCol, lwd=LWD)
lines(max(exRTT$times) - exRTT$times, apply(exRTT$lambda, 2, mean), col=exCol, lwd=LWD)

axis(1, at=Xax)
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
segments(66, 0, 66, 0.1, lty=3, lwd=1.2)

# Extinction plot
plot(max(exRTT$times) - exRTT$times, apply(exRTT$mu,2,mean), type='n', ylim=c(0,0.1), xlim=c(startT, 0), axes=F, xlab="age (Ma)", ylab=expression(paste("extinction rate (", mu, ")", sep="")), col=exCol, lwd=LWD)

silent <- sapply(1:length(fbRTT_all), function(x) lines(max(fbRTT_all[[x]]$times) - fbRTT_all[[x]]$times, apply(fbRTT_all[[x]]$mu, 2, mean), col=allCol))
silent <- sapply(1:length(fbRTT), function(x) lines(max(fbRTT[[x]]$time) - fbRTT[[x]]$times, apply(fbRTT[[x]]$mu, 2, mean), col=npfCol))

lines(max(npRTT$times) - npRTT$times, apply(npRTT$mu, 2, mean), col=npCol, lwd=LWD)
lines(max(exRTT$times) - exRTT$times, apply(exRTT$mu,2,mean), col=exCol, lwd=LWD)

axis(1, at=Xax)
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))

segments(66, 0, 66, 0.1, lty=3, lwd=1.2)
#dev.off()

#pdf("derivedPars.pdf", height=5, width=10)
#par(mfcol=c(1,2), mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
# Net div plot
plot(max(exRTT$times) - exRTT$times, apply(exRTT$lambda,2,mean) - apply(exRTT$mu,2,mean), type='n', ylim=c(0,0.1), xlim=c(startT, 0), axes=F, xlab="age (Ma)", ylab=expression(paste("net diversification (", lambda, "-", mu, ")", sep="")), col=exCol, lwd=LWD)

silent <- sapply(1:length(fbRTT_all), function(x) lines(max(fbRTT_all[[x]]$times) - fbRTT_all[[x]]$times, apply(fbRTT_all[[x]]$lambda, 2, mean) - apply(fbRTT_all[[x]]$mu, 2, mean), col=allCol))
silent <- sapply(1:length(fbRTT), function(x) lines(max(fbRTT[[x]]$times) - fbRTT[[x]]$times, apply(fbRTT[[x]]$lambda, 2, mean) - apply(fbRTT[[x]]$mu, 2, mean), col=npfCol))

lines(max(npRTT$times) - npRTT$times, apply(npRTT$lambda,2,mean) - apply(npRTT$mu, 2, mean), col=npCol, lwd=LWD)
lines(max(exRTT$times) - exRTT$times, apply(exRTT$lambda,2,mean) - apply(exRTT$mu,2,mean), col=exCol, lwd=LWD)

axis(1, at=Xax)
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
legend("topleft", bty="n", legend=c("extant", "extant-reduced", "fossils", "fossils-reduced"), text.col=c(exCol, npCol, allCol, npfCol), cex=1.5)

# Rel ext plot
plot(max(exRTT$times) - exRTT$times, apply(exRTT$mu,2,mean) / apply(exRTT$lambda,2,mean), type='n', ylim=c(0,1), xlim=c(startT, 0), axes=F, xlab="age (Ma)", ylab=expression(paste("relative extinction (", mu, "/", lambda, ")", sep="")), lwd=LWD, col=exCol)

silent <- sapply(1:length(fbRTT_all), function(x) lines(max(fbRTT_all[[x]]$times) - fbRTT_all[[x]]$times, apply(fbRTT_all[[x]]$mu, 2, mean) / apply(fbRTT_all[[x]]$lambda, 2, mean), col=allCol))
silent <- sapply(1:length(fbRTT), function(x) lines(max(fbRTT[[x]]$times) - fbRTT[[x]]$times,  apply(fbRTT[[x]]$mu, 2, mean) / apply(fbRTT[[x]]$lambda, 2, mean), col=npfCol))

lines(max(npRTT$times) - npRTT$times, apply(npRTT$mu, 2, mean) / apply(npRTT$lambda,2,mean), col=npCol, lwd=LWD)
lines(max(exRTT$times) - exRTT$times, apply(exRTT$mu,2,mean) / apply(exRTT$lambda,2,mean), lwd=LWD, col=exCol)

axis(1, at=Xax)
axis(2, at=c(-10, 0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()

sPlot <- plot(fbEdata_all[[1]], breaks="jenks", spex="s")
ePlot <- plot(fbEdata_all[[1]], breaks="jenks", spex="e")

setwd(paste(Path_base, "figures", sep=""))
pdf("spexComparison.pdf", height=10, width=10)
par(mfrow=c(2,2), mar=c(0,0,0,7), oma=c(0,0,0,0), pty="s")
x <- plot(exEdata, colorbreaks=sPlot$colorbreaks)
LOC <- c(25, 30, 210, 380)
addBAMMlegend(x, location=LOC)
tree <- as.phylo(exEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae", Cex=CEX)
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae", Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extant-only", side=3, line=0, xpd=NA, cex=1.5)
mtext("speciation rate", side=2, at=mean(c(LOC[3], LOC[4])), line=-2.5, xpd=NA, cex=1.5, las=0)
#segments(219, 177, 302, 475, xpd=NA, lty=2)
#segments(215, -14, 519, -14, xpd=NA, lty=2)

plot(fbEdata_all[[1]], colorbreaks=sPlot$colorbreaks)
tree <- as.phylo(fbEdata_all[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae", Cex=CEX)
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae", Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extant + fossils", side=3, line=0, xpd=NA, cex=1.5)

x <- plot(exEdata, colorbreaks=ePlot$colorbreaks, spex="e")
addBAMMlegend(x, location=LOC)
tree <- as.phylo(exEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae", Cex=CEX)
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae", Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extinction rate", side=2, line=-2.5, at=mean(c(LOC[3], LOC[4])), xpd=NA, cex=1.5)
#segments(219, 177, 302, 475, xpd=NA, lty=2)
#segments(215, -14, 519, -14, xpd=NA, lty=2)

plot(fbEdata_all[[1]], colorbreaks=ePlot$colorbreaks, spex="e")
tree <- as.phylo(fbEdata_all[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae", Cex=CEX)
addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae", Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
dev.off()

sPlot <- plot(fbEdata[[1]], breaks="jenks", spex="s")
ePlot <- plot(fbEdata[[1]], breaks="jenks", spex="e")

setwd(paste(Path_base, "figures", sep=""))
pdf("spexComparison_np.pdf", height=10, width=10)
par(mfrow=c(2,2), mar=c(0,0,0,7), oma=c(0,0,0,0), pty="s")
x <- plot(npEdata, colorbreaks=sPlot$colorbreaks)
LOC <- c(25, 30, 105, 160)
addBAMMlegend(x, location=LOC)
tree <- as.phylo(npEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extant-only", side=3, line=0, xpd=NA, cex=1.5)
mtext("speciation rate", side=2, at=mean(c(LOC[3], LOC[4])), line=-2.5, xpd=NA, cex=1.5, las=0)

plot(fbEdata[[1]], colorbreaks=sPlot$colorbreaks)
tree <- as.phylo(fbEdata[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extant + fossils", side=3, line=0, xpd=NA, cex=1.5)

x <- plot(npEdata, colorbreaks=ePlot$colorbreaks, spex="e")
addBAMMlegend(x, location=LOC)
tree <- as.phylo(npEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extinction rate", side=2, line=-2.5, at=mean(c(LOC[3], LOC[4])), xpd=NA, cex=1.5)

plot(fbEdata[[1]], colorbreaks=ePlot$colorbreaks, spex="e")
tree <- as.phylo(fbEdata[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
addFamily(c("Tylototriton", "Echinotriton", "Pleurodeles"), "Salamandridae", centering=T, Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
dev.off()


make.plot <- function(vec1, vec2, label=NULL, labline=3)	{
	R <- range(c(vec1, vec2))
	Ticks <- pretty(R)
	plot(1, 1, type="n", xlim=range(Ticks), ylim=range(Ticks), axes=F, xlab="", ylab="")
	points(vec1, vec2, pch=21, bg='gray70')
	axis(1, at=c(-100, Ticks))
	axis(2, at=c(-100, Ticks))
	segments(Ticks[1], Ticks[1], Ticks[length(Ticks)], Ticks[length(Ticks)], lty=2)
	if (!is.null(label))	{
		mtext(label, side=2, at=max(Ticks), line=labline, cex=1.25)
	}
}
options(scipen=6)


setwd(paste(Path_base, "figures", sep=""))
pdf("rateCors.pdf", height=6.5, width=9.75)
# Tip rates between fossil and extant-only species
par(mfrow=c(2,3), mar=c(4, 5, 0.5, 0.5), bty="n", las=1, mgp=c(2, 0.5, 0), tck=-0.01, cex.axis=0.9)
Liney <- 3
Linex <- 2

# all
edata <- fbEdata_all[[1]]
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label
exSpec <- exEdata$meanTipLambda
names(exSpec) <- exEdata$tip.label
exExt <- exEdata$meanTipMu
names(exExt) <- exEdata$tip.label

Inter <- intersect(exEdata$tip.label, edata$tip.label)

make.plot(exSpec[Inter], Spec[Inter], label="a)")
mtext("speciation rate - fossils", side=2, las=0, line=Liney)
mtext("speciation rate - extant", side=1, line=Linex)
make.plot(exExt[Inter], Ext[Inter], "b)")
mtext("extinction rate - fossils", side=2, las=0, line=Liney)
mtext("extinction rate - extant", side=1, line=Linex)
make.plot(exSpec[Inter] - exExt[Inter], Spec[Inter] - Ext[Inter], "c)")
mtext("net div. rate - fossils", side=2, las=0, line=Liney)
mtext("net div. rate - extant", side=1, line=Linex)

# No Pleths
edata <- fbEdata[[1]]
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label
exSpec <- npEdata$meanTipLambda
names(exSpec) <- npEdata$tip.label
exExt <- npEdata$meanTipMu
names(exExt) <- npEdata$tip.label

Inter <- intersect(npEdata$tip.label, edata$tip.label)

make.plot(exSpec[Inter], Spec[Inter], label="c)")
mtext("speciation rate - fossils", side=2, las=0, line=Liney)
mtext("speciation rate - extant", side=1, line=Linex)
make.plot(exExt[Inter], Ext[Inter], label="d)")
mtext("extinction rate - fossils", side=2, las=0, line=Liney)
mtext("extinction rate - extant", side=1, line=Linex)
make.plot(exSpec[Inter] - exExt[Inter], Spec[Inter] - Ext[Inter], "c)")
mtext("net div. rate - fossils", side=2, las=0, line=Liney)
mtext("net div. rate - extant", side=1, line=Linex)

dev.off()


#########################################################################################
# No Pleths
edata <- fbEdata_all[[1]]
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label
exSpec <- npEdata$meanTipLambda
names(exSpec) <- npEdata$tip.label
exExt <- npEdata$meanTipMu
names(exExt) <- npEdata$tip.label

# IUCN
IUCN <- read.csv(paste(Path_base, "datafiles/IUCN_list.csv", sep=""), stringsAsFactors=F)
IUCN_sp <- apply(IUCN, 1, function(x) paste(c(x[which(colnames(IUCN)=="Genus")], x[which(colnames(IUCN)=="Species")]), sep="_", collapse="_"))
Scores <- IUCN$Red.List.status
names(Scores) <- IUCN_sp

# Key: DD= data deficient, NE= not evaluated, 
# then: LC -> CD -> NT -> VU -> EN -> CR -> EW -> EX

source('~/Documents/BAMMtools_Functions/spindlePlot.R', chdir = TRUE)
par(mfrow=c(1, 2), mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(1, 0.5, 0), tck=-0.01)
COL <- 'gray70'
WF <- 1.5
Cats <- c("LC", "NT", "VU", "EN")

plot(meanSp[Cats], axes=F, xlim=c(0.5, 4.5), ylim=c(0, 0.2), xlab="", ylab="", type="n")
spindlePlot(Spec[names(which(Scores=="LC"))], at=1, col=COL, widthFactor=WF, breakPts=seq(from=0, to=1, by=0.01), Border=COL)
spindlePlot(Spec[names(which(Scores=="NT"))], at=2, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.01))
spindlePlot(Spec[names(which(Scores=="VU"))], at=3, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.01))
spindlePlot(Spec[names(which(Scores=="EN"))], at=4, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.01))

axis(1, at=c(1, 2, 3, 4), labels=Cats)
axis(2, at=c(0, 0.05, 0.1, 0.15, 0.2))

plot(meanEx[Cats], axes=F, xlim=c(0.5, 4.5), ylim=c(0, 0.08), xlab="", ylab="", type="n")
spindlePlot(Ext[names(which(Scores=="LC"))], at=1, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.005))
spindlePlot(Ext[names(which(Scores=="NT"))], at=2, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.005))
spindlePlot(Ext[names(which(Scores=="VU"))], at=3, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.005))
spindlePlot(Ext[names(which(Scores=="EN"))], at=4, col=COL, widthFactor=WF, breaks=20, Border=COL, breakPts=seq(from=0, to=1, by=0.005))

axis(1, at=c(1, 2, 3, 4), labels=Cats)
axis(2, at=c(0, 0.02, 0.04, 0.06, 0.08))

plot(meanEx[Cats], axes=F, xlim=c(0.5, 4.5), ylim=c(0, 0.1), xlab="", ylab="", type="n")
spindlePlot(Spec[names(which(Scores=="LC"))] - Ext[names(which(Scores=="LC"))], at=1, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Spec[names(which(Scores=="NT"))] - Ext[names(which(Scores=="NT"))], at=2, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Spec[names(which(Scores=="VU"))] - Ext[names(which(Scores=="VU"))], at=3, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Spec[names(which(Scores=="EN"))] - Ext[names(which(Scores=="EN"))], at=4, col=COL, widthFactor=WF, breaks=20, Border=COL)

axis(1, at=c(1, 2, 3, 4), labels=Cats)
axis(2, at=c(0, 0.02, 0.04, 0.06, 0.08, 0.1))

plot(meanEx[Cats], axes=F, xlim=c(0.5, 4.5), ylim=c(0, 1), xlab="", ylab="", type="n")
spindlePlot(Ext[names(which(Scores=="LC"))] / Spec[names(which(Scores=="LC"))], at=1, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Ext[names(which(Scores=="NT"))] / Spec[names(which(Scores=="NT"))], at=2, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Ext[names(which(Scores=="VU"))] / Spec[names(which(Scores=="VU"))], at=3, col=COL, widthFactor=WF, breaks=20, Border=COL)
spindlePlot(Ext[names(which(Scores=="EN"))] / Spec[names(which(Scores=="EN"))], at=4, col=COL, widthFactor=WF, breaks=20, Border=COL)
axis(1, at=c(1, 2, 3, 4), labels=Cats)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))


#########################################################################################
# Geographic range
Areas <- read.table(paste(Path_base, "datafiles/areas.txt", sep=""), row.names=1, stringsAsFactors=F)

edata <- fbEdata_all[[1]]
Inter <- intersect(rownames(Areas), edata$tip.label)
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label

par(mfrow=c(2,2))
Xax <- c(-10, 5.5, 7.5, 9.5, 11.5, 13.5)
Pch <- 16
Col <- rgb(0,0,0,0.25)
plot(log10(Areas[Inter,1]), Ext[Inter], axes=F, xlim=c(5.5,13.5), pch=Pch, col=Col, xlab="")
axis(1, at=Xax)
plot(log10(Areas[Inter,1]), Spec[Inter], axes=F, xlim=c(5.5,13.5), pch=Pch, col=Col, xlab="")
axis(1, at=Xax)
plot(log10(Areas[Inter,1]), Spec[Inter]-Ext[Inter], axes=F, xlim=c(5.5,13.5), pch=Pch, col=Col, xlab=expression(paste(log, "Area")))
axis(1, at=Xax)
plot(log10(Areas[Inter,1]), Ext[Inter]/Spec[Inter], axes=F, xlim=c(5.5,13.5), pch=Pch, col=Col, xlab=expression(paste(log, "Area")))
axis(1, at=Xax)

#########################################################################################
# Neot scores
Neot <- read.csv(paste(Path_base, "datafiles/neoteny/neoteny.csv", sep=""), stringsAsFactors=F)
rownames(Neot) <- apply(Neot, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
nScores <- Neot[,"neotenic"]
names(nScores) <- rownames(Neot)

EDAT <- fbEdata_all[[1]]
Sim <- make.simmap(as.phylo(EDAT), nScores[EDAT$tip.label], model="ARD", nsim=1e3, pi=c(0.95, 0.05, 0), Q="mcmc")

nodeVec <- seq(from=(length(EDAT$tip.label) + 1), to=(length(EDAT$tip.label) + EDAT$Nnode), by=1)
relRates <- sapply(nodeVec, nodeRate, edata=EDAT)
transP <- sapply(nodeVec, transProb, simmap=Sim)

Ylims <- c(0, 500)
Xlims <- c(-2, 2)

par(mfrow=c(2,3), mar=c(4,5,1,1))
# Spec
rateVec <- relRates[3,2:ncol(relRates)] / relRates[1,2:ncol(relRates)]
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["1",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=c(-1, 2), ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)

# Ext
rateVec <- relRates[4, 2:ncol(relRates)] / relRates[2, 2:ncol(relRates)]
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["1",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=Xlims, ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)

# Turn
rateVec <- (relRates[4, 2:ncol(relRates)] + relRates[3, 2:ncol(relRates)]) / (relRates[2, 2:ncol(relRates)] + relRates[1, 2:ncol(relRates)])
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["1",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=Xlims, ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)


# Spec
rateVec <- relRates[3,2:ncol(relRates)] / relRates[1,2:ncol(relRates)]
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["2",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=Xlims, ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)

# Ext
rateVec <- relRates[4, 2:ncol(relRates)] / relRates[2, 2:ncol(relRates)]
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["2",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=Xlims, ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)

# Turn
rateVec <- (relRates[4, 2:ncol(relRates)] + relRates[3, 2:ncol(relRates)]) / (relRates[2, 2:ncol(relRates)] + relRates[1, 2:ncol(relRates)])
subRateVec <- sample(rateVec, size=length(rateVec), replace=TRUE, prob=transP["2",2:ncol(transP)])
x <- hist(log(rateVec), main="", xlim=Xlims, ylim=Ylims, col=rgb(0, 0, 0, 0.25))
hist(log(subRateVec), main="", xlim=Xlims, ylim=Ylims, add=T, col=rgb(1,0,0,0.5), breaks=x$breaks)



setwd(paste(Path_base, "figures", sep=""))
pdf("dev_rates.pdf", height=6.5, width=6.5)
# Tip rates between fossil and extant-only species
par(mfrow=c(1,1), mar=c(4, 5, 0.5, 0.5), bty="n", las=1, mgp=c(2, 0.5, 0), tck=-0.01, cex.axis=0.9)

NetDiv <- fbEdata[[1]]$meanTipLambda - fbEdata[[1]]$meanTipMu
names(NetDiv) <- fbEdata[[1]]$tip.label
subScores <- nScores[names(NetDiv)]
plot(0, 0, type="n", bty='n', xlab='', ylab='', axes=F, xlim=c(-0.02, 0.1), ylim=c(0, 2e2))
polygon(density(NetDiv[names(subScores)[which(subScores==1)]]), col=rgb(1,0,0,0.5), border=NA)
polygon(density(NetDiv[names(subScores)[which(subScores==0)]]), col=rgb(0,0,1,0.5), border=NA)
polygon(density(NetDiv[names(subScores)[which(subScores==2)]]), col=rgb(0,1,0,0.5), border=NA)
axis(1, at=c(-100, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
axis(2, at=c(-100, 0, 50, 100, 150, 200))
mtext("net div.", 1, line=2, cex=1.5)
mtext("density", 2, line=3, las=0, cex=1.5)
legend("topleft", legend=c("metamorph.", "fac. neot.", "obl. neot."), text.col=c(rgb(0,0,1,0.9), rgb(1,0,0,0.9), rgb(0,1,0,0.9)), bty="n", cex=1.25)
dev.off()

Ext <- matrix(NA, nrow=70, ncol=5)
Spec <- matrix(NA, nrow=70, ncol=5)
Net <- matrix(NA, nrow=70, ncol=5)

for (i in 1:70)	{
	edat <- fbEdata_all[[i]]
	vec <- nScores[edat$tip.label]
	names(vec) <- edat$tip.label
	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="extinction", method="kruskal", logrates=FALSE)
	Ext[i,] <- c(unlist(Est$estimate)[c("-1", "0", "1", "2")], p=Est$p.value)

	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="speciation", method="kruskal", logrates=FALSE)
	Spec[i,] <- c(unlist(Est$estimate)[c("-1", "0", "1", "2")], p=Est$p.value)
	
	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="net", method="kruskal", logrates=FALSE)
	Net[i,] <- c(unlist(Est$estimate)[c("-1", "0", "1", "2")], p=Est$p.value)
	cat(i, "\n")
}

nfExt <- matrix(NA, nrow=70, ncol=4)
nfSpec <- matrix(NA, nrow=70, ncol=4)
nfNet <- matrix(NA, nrow=70, ncol=4)
for (i in 1:70)	{
	edat <- fbEdata_all[[i]]
	vec <- nScores[edat$tip.label]
	names(vec) <- edat$tip.label
	vec[which(vec==1)] <- 0
	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="extinction", method="kruskal", logrates=FALSE)
	nfExt[i,] <- c(unlist(Est$estimate)[c("-1", "0", "2")], p=Est$p.value)

	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="speciation", method="kruskal", logrates=FALSE)
	nfSpec[i,] <- c(unlist(Est$estimate)[c("-1", "0", "2")], p=Est$p.value)

	Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="net", method="kruskal", logrates=FALSE)
	nfNet[i,] <- c(unlist(Est$estimate)[c("-1", "0", "2")], p=Est$p.value)
	cat(i, "\n")
}


edat <- fbEdata[[1]]

vec <- nScores[edat$tip.label]
names(vec) <- edat$tip.label
vec[which(vec==1)] <- 2
Est <- traitDependentBAMM(ephy=edat, traits=vec, rate="n", method="m", logrates=FALSE)



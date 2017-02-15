Path_base <- "~/Documents/Salamanders/"

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
exTree <- read.tree("no_pleth.tre")
exMCMC <- read.csv("np_mcmc_out.txt", stringsAsFactors=F)
exEdata <- getEventData(exTree, "np_event_data.txt", burnin=0.5)
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

# Extant-only rates through time
exRTT <- getRateThroughTimeMatrix(exEdata)

# fossilBAMM runs
### Post-generation addition of files
fbEdata <- list()
fbEdata_hi <- list()
fbRTT <- list()
fbRTT_hi <- list()
fmax <- c()
fullmax <- c()
llES <- c()
nsES <- c()
Nshifts <- c()
meanNshifts <- c()
presRate <- c()

for (count in 1:70)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	fbTree <- read.tree("fossilTree.tre")
	MCMC_dat <- read.csv("base_mcmc_out.txt")
	MCMC_dat <- MCMC_dat[floor(0.1*nrow(MCMC_dat)):nrow(MCMC_dat),]
	
	fbEdata[[count]] <- getEventData(fbTree, "base_event_data.txt", burnin=0.1, nsamples=200)
	fbRTT[[count]] <- getRateThroughTimeMatrix(fbEdata[[count]])
	llES[count] <- effectiveSize(MCMC_dat$logLik)
	nsES[count] <- effectiveSize(MCMC_dat$N_shifts)
	BFmat <- computeBayesFactors(MCMC_dat, 200, burnin=0)
	Nshifts[count] <- stepBF(BFmat, step.size=20, expectedNumberOfShifts=20, inputType="matrix")
	meanNshifts[count] <- mean(MCMC_dat$N_shifts)
	presRate[count] <- mean(MCMC_dat$preservationRate)
	
	fullTree <- read.tree("fullTree.tre")
	fbEdata_hi[[count]] <- getEventData(fullTree, "full_event_data.txt", burnin=0.1, nsamples=200)
	fbRTT_hi[[count]] <- getRateThroughTimeMatrix(fbEdata_hi[[count]])
	
	fmax[count] <- max(nodeHeights(fbTree))
	fullmax[count] <- max(nodeHeights(fullTree))
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
#addFamily(c("Rhyacotriton", "Rhyacotriton", "Rhyacotriton"), "Rhyacotritonidae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
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
#addFamily(c("Rhyacotriton", "Rhyacotriton", "Rhyacotriton"), "Rhyacotritonidae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
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
exCol <- rgb(217/255, 95/255, 2/255, Alpha)
hiCol <- rgb(27/255, 158/255, 119/255, Alpha)
loCol <- rgb(117/255, 112/255, 179/255, Alpha)

setwd(paste(Path_base, "figures", sep=""))
pdf("spex.pdf", height=5, width=10)
par(mfcol=c(1,2), mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
# Speciation plot
tDiff <- max(fmax) - maxT
plot(exRTT$times + tDiff, apply(exRTT$lambda,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("speciation rate (", lambda, ")", sep="")), col=exCol, lwd=1.5)

silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$lambda, 2, mean), col=hiCol))

silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$lambda, 2, mean), col=loCol))
Xax <- c(-100, 0, 50, 100, 150, 204, 240, 270)
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
segments(max(fmax) - 66, 0, max(fmax) - 66, 0.1, lty=3, lwd=1.2)

# Extinction plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("extinction rate (", mu, ")", sep="")), col=exCol, lwd=1.5)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')

silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$mu, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))

legend("topright", bty="n", legend=c("extant", "fossils", "fossils-full"), lty=1, col=c(exCol, loCol, hiCol))
#legend("topright", bty="n", legend=c("extant", "fossils"), lty=1, col=c(exCol, loCol))
segments(max(fmax) - 66, 0, max(fmax) - 66, 0.1, lty=3, lwd=1.2)
dev.off()

pdf("derivedPars.pdf", height=5, width=10)
par(mfcol=c(1,2), mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
# Net div plot
plot(exRTT$times + tDiff, apply(exRTT$lambda,2,mean) - apply(exRTT$mu,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("net diversification (", lambda, "-", mu, ")", sep="")), col=exCol, lwd=1.5)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$lambda, 2, mean) - apply(fbRTT_hi[[x]]$mu, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$lambda, 2, mean) - apply(fbRTT[[x]]$mu, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
legend("topright", bty="n", legend=c("extant", "fossils", "fossils-full"), lty=1, col=c(exCol, loCol, hiCol))

# Rel ext plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean) / apply(exRTT$lambda,2,mean), type='l', ylim=c(0,1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("relative extinction (", mu, "/", lambda, ")", sep="")), lwd=1.5, col=exCol)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean) / apply(fbRTT_hi[[x]]$lambda, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x],  apply(fbRTT[[x]]$mu, 2, mean) / apply(fbRTT[[x]]$lambda, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()

sPlot <- plot(fbEdata[[1]], breaks="jenks", spex="s")
ePlot <- plot(fbEdata[[1]], breaks="jenks", spex="e")

setwd(paste(Path_base, "figures", sep=""))
pdf("spexComparison.pdf", height=10, width=10)
par(mfrow=c(2,2), mar=c(0,0,0,7), oma=c(0,5,2,0), pty="s")
plot(exEdata, colorbreaks=sPlot$colorbreaks)
addBAMMlegend(sPlot, location=c(-36,-30,152,313))
tree <- as.phylo(exEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
#addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
#addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae", Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extant-only", side=3, line=0, xpd=NA, cex=1.5)
mtext("speciation rate", side=2, line=3, xpd=NA, cex=1.5)
segments(219, 177, 302, 475, xpd=NA, lty=2)
segments(215, -14, 519, -14, xpd=NA, lty=2)

plot(fbEdata[[1]], colorbreaks=sPlot$colorbreaks)
tree <- as.phylo(fbEdata[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
#addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
#addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Paramesotriton", "Paramesotriton", "Pachytriton"), "Salamandridae", Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
mtext("extant + extinct", side=3, line=0, xpd=NA, cex=1.5)

#plot(fbEdata_hi[[1]], colorbreaks=sPlot$colorbreaks)
#tree <- as.phylo(fbEdata_hi[[1]])
#addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
#addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
#addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
#addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
#addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
#addFamily(c("Paramesotriton", "Paramesotriton", "Pachytriton"), "Salamandridae", Cex=CEX)
#addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)

plot(exEdata, colorbreaks=ePlot$colorbreaks, spex="e")
addBAMMlegend(ePlot, location=c(-36,-30,152,313))
tree <- as.phylo(exEdata)
CEX <- 1
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
#addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
#addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae", Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae", Cex=CEX)
mtext("extinction rate", side=2, line=3, xpd=NA, cex=1.5)
segments(219, 177, 302, 475, xpd=NA, lty=2)
segments(215, -14, 519, -14, xpd=NA, lty=2)

plot(fbEdata[[1]], colorbreaks=ePlot$colorbreaks, spex="e")
tree <- as.phylo(fbEdata[[1]])
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=T, Cex=CEX)
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T, Cex=CEX)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T, Cex=CEX)
#addFamily(c("Aneides", "Plethodon", "Ensatina"), "Plethodontinae")
#addFamily(c("Bolito", "Chiropterotriton", "Batrachoseps"), "Bolitoglossinae")
addFamily(c("Paramesotriton", "Paramesotriton", "Pachytriton"), "Salamandridae", Cex=CEX)
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae", Cex=CEX)
dev.off()



# Tip rates between fossil and extant-only species
edata <- fbEdata[[1]]
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label
exSpec <- exEdata$meanTipLambda
names(exSpec) <- exEdata$tip.label
exExt <- exEdata$meanTipMu
names(exExt) <- exEdata$tip.label

Inter <- intersect(exEdata$tip.label, edata$tip.label)

par(mfrow=c(1,2))
hist(exSpec[Inter], xlim=c(0,0.14), border=rgb(0,0,0,0), col=rgb(0,0,0,0.5), main="", ylim=c(0, 60))
hist(Spec[Inter], xlim=c(0,0.14), ylim=c(0,60), add=T, axes=F, col=rgb(1, 0, 0, 0.5), border=rgb(0, 0, 0, 0))

hist(exExt[Inter], xlim=c(0,0.1), border=rgb(0,0,0,0), col=rgb(0,0,0,0.5), main="", ylim=c(0, 60))
hist(Ext[Inter], xlim=c(0,0.1), add=T, axes=F, ylim=c(0,60), col=rgb(1,0,0,0.5), border=rgb(0,0,0,0))


#########################################################################################
# IUCN
IUCN <- read.csv(paste(Path_base, "datafiles/IUCN_list.csv", sep=""), stringsAsFactors=F)
IUCN_sp <- apply(IUCN, 1, function(x) paste(c(x[which(colnames(IUCN)=="Genus")], x[which(colnames(IUCN)=="Species")]), sep="_", collapse="_"))
Scores <- IUCN$Red.List.status
names(Scores) <- IUCN_sp

# Key: DD= data deficient, NE= not evaluated, 
# then: LC -> CD -> NT -> VU -> EN -> CR -> EW -> EX

#Inter <- intersect(names(Scores), edata$tip.label)

meanSp <- tapply(Spec[Inter], Scores[Inter], mean)
sdSp <- tapply(Spec[Inter], Scores[Inter], sd)
meanEx <- tapply(Ext[Inter], Scores[Inter], mean)
sdEx <- tapply(Ext[Inter], Scores[Inter], sd)

Cats <- c("LC", "NT", "VU", "EN")
plot(meanSp[Cats], axes=F)

plot(meanEx[Cats], axes=F, xlim=c(0,4), ylim=c(0,0.15))
silent <- sapply(1:length(Cats), function(x) segments(x, meanEx[Cats[x]]-sdEx[Cats[x]], x, meanEx[Cats[x]]+sdEx[Cats[x]]))


# Geographic range
Areas <- read.table(paste(Path_base, "datafiles/areas.txt", sep=""), row.names=1, stringsAsFactors=F)

edata <- fbEdata[[1]]
Inter <- intersect(rownames(Areas), edata$tip.label)
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label

# Let's check only Salamandridae
#newts <- edata$tip.label[c(grep("Salamandrina", edata$tip.label), grep("Taricha", edata$tip.label))]
#newt_tips <- edata$tip.label[getDescendants(as.phylo(edata), node=findMRCA(as.phylo(edata), newts))]
#newt_tips <- na.omit(newt_tips)
#Inter <- intersect(Inter, newt_tips)

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



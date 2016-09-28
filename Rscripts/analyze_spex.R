Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)

# PBDB data
setwd(paste(Path_base, "fossilSalamanders/pbdb", sep=""))
fosRTT <- read.csv("raw_curve_data.csv", stringsAsFactors=F)
fosNetDiv <- fosRTT$X3T_lambda - fosRTT$X3T_mu
fosRelExt <- fosRTT$X3T_mu/fosRTT$X3T_lambda

# Extant-only data
setwd(paste(Path_base, "bamm/extant_only", sep=""))
exTree <- read.tree("modernSal.tre")
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
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=F)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=F)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
addFamily(c("Proteus", "Necturus", "Necturus"), "Proteidae", centering=F, Top=F)
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
for (count in 1:10)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	fbTree <- read.tree("fossilTree.tre")
	fbEdata[[count]] <- getEventData(fbTree, "lo_event_data.txt", burnin=0.5)
	fbRTT[[count]] <- getRateThroughTimeMatrix(fbEdata[[count]])
	
	fbEdata_hi[[count]] <- getEventData(fbTree, "hi_event_data.txt", burnin=0.5)
	fbRTT_hi[[count]] <- getRateThroughTimeMatrix(fbEdata_hi[[count]])
	
	fmax[count] <- max(nodeHeights(fbTree))
}


setwd(paste(Path_base, "figures", sep=""))
pdf("fossil_sp.pdf", height=10, width=10)
fosEdata <- fbEdata[[1]]
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
fosPlot <- plot(fosEdata, colorbreaks=exPlot$colorbreaks)
#addBAMMlegend(fosPlot, location="bottomleft")

# Add family labels
tree <- as.phylo(fosEdata)
#addFamily(c("Rhyacotriton", "Rhyacotriton", "Rhyacotriton"), "Rhyacotritonidae")
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=F)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=F)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
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


# Extinction plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("extinction rate (", mu, ")", sep="")), col=exCol, lwd=1.5)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$mu, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
legend("topright", bty="n", legend=c("extant", "fossils", "fossils x 10"), lty=1, col=c(exCol, loCol, hiCol))
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
legend("topright", bty="n", legend=c("extant", "fossils", "fossils x 10"), lty=1, col=c(exCol, loCol, hiCol))

# Rel ext plot
plot(exRTT$times + tDiff, apply(exRTT$mu,2,mean) / apply(exRTT$lambda,2,mean), type='l', ylim=c(0,1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("relative extinction (", mu, "/", lambda, ")", sep="")), lwd=1.5, col=exCol)
#points(max(fmax) - fosRTT$Midpoint_Ma, fosRTT$X3T_mu, pch=16, type='b')
silent <- sapply(1:length(fbRTT_hi), function(x) lines(fbRTT_hi[[x]]$times + max(fmax) - fmax[x], apply(fbRTT_hi[[x]]$mu, 2, mean) / apply(fbRTT_hi[[x]]$lambda, 2, mean), col=hiCol))
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x],  apply(fbRTT[[x]]$mu, 2, mean) / apply(fbRTT[[x]]$lambda, 2, mean), col=loCol))
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.2, 0.4, 0.6, 0.8, 1))
dev.off()


# IUCN
IUCN <- read.csv(paste(Path_base, "datafiles/IUCN_list.csv", sep=""), stringsAsFactors=F)
IUCN_sp <- apply(IUCN, 1, function(x) paste(c(x[which(colnames(IUCN)=="Genus")], x[which(colnames(IUCN)=="Species")]), sep="_", collapse="_"))
Scores <- IUCN$Red.List.status
names(Scores) <- IUCN_sp

# Key: DD= data deficient, NE= not evaluated, 
# then: LC -> CD -> NT -> VU -> EN -> CR -> EW -> EX

edata <- fbEdata[[1]]
Inter <- intersect(names(Scores), edata$tip.label)
Ext <- edata$meanTipMu
names(Ext) <- edata$tip.label
Spec <- edata$meanTipLambda
names(Spec) <- edata$tip.label

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



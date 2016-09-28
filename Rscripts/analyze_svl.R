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

# Extant-only data
setwd(paste(Path_base, "bamm_morph/extant-only", sep=""))
exTree <- read.tree("extantTree.tre")
exMCMC <- read.csv("ex_svl_mcmc_out.txt", stringsAsFactors=F)
exEdata <- getEventData(exTree, "ex_svl_event_data.txt", burnin=0.5, type="trait")
nodeTrees <- read.tree("ex_svl_node_state.txt")
nodeMat <- sapply(nodeTrees, function(x) x$edge.length)
nodeVals <- apply(nodeMat[,101:1000], 1, mean)
nodeStates <- nodeTrees[[1]]
nodeStates$edge.length <- nodeVals
ex_nodeStates <- nodeStates
maxT <- max(nodeHeights(exTree))


setwd(paste(Path_base, "figures", sep=""))
pdf("extant_beta.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
exPlot <- plot(exEdata, breaks="quantile")
addBAMMlegend(exPlot)

# Add family labels
tree <- as.phylo(exEdata)
#tiplabels(text="Rhyacotritonidae", tip=grep("Rhyacotriton", tree$tip.label), bg='white', frame='none', xpd=NA, adj=-0.05, cex=0.5)
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae")
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae")
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae")
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Salamandra", "Neurergus", "Chioglossa"), "Salamandridae")
tiplabels(text="Cryptobranchidae", tip=grep("Andrias", tree$tip.label), bg='white', frame='none', xpd=NA, adj=-0.05, cex=1.5)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
addFamily(c("Proteus", "Necturus", "Necturus"), "Proteidae")
dev.off()

# Extant-only rates through time
exRTT <- getRateThroughTimeMatrix(exEdata)

# fossilBAMM runs
### Post-generation addition of files
fbEdata <- list()
fbRTT <- list()
fbNodes <- list()
fmax <- c()
for (count in 1:10)	{
	setwd(paste(Path_base, "bamm_morph/output-", count, sep=""))
	fbTree <- read.tree("fossilTree.tre")
	fbEdata[[count]] <- getEventData(fbTree, "svl_event_data.txt", burnin=0.5, type="trait")
	fbRTT[[count]] <- getRateThroughTimeMatrix(fbEdata[[count]])
	fmax[count] <- max(nodeHeights(fbTree))
	nodeTrees <- read.tree("svl_node_state.txt")
	nodeMat <- sapply(nodeTrees, function(x) x$edge.length)
	nodeVals <- apply(nodeMat[,101:1000], 1, mean)
	nodeStates <- nodeTrees[[1]]
	nodeStates$edge.length <- nodeVals
	fbNodes[[count]] <- nodeStates
}


setwd(paste(Path_base, "figures", sep=""))
pdf("fossil_beta.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
fosPlot <- plot(fbEdata[[1]], colorbreaks=exPlot$colorbreaks)

# Add family labels
tree <- as.phylo(fbEdata[[1]])
#tiplabels(text="Rhyacotritonidae", tip=grep("Rhyacotriton", tree$tip.label), bg='white', frame='none', xpd=NA, adj=-0.05, cex=0.5)
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae")
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae")
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=F, Top=T)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Taricha", "Taricha", "Taricha"), "Salamandridae")
tiplabels(text="Cryptobranchidae", tip=grep("Andrias", tree$tip.label), bg='white', frame='none', xpd=NA, adj=-0.05, cex=2)
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
addFamily(c("Proteus", "Necturus", "Necturus"), "Proteidae", centering=F, Top=F)
dev.off()


setwd(paste(Path_base, "figures", sep=""))
pdf("beta_curves.pdf", height=5, width=5)
Alpha <- 1
exCol <- rgb(217/255, 95/255, 2/255, Alpha)
hiCol <- rgb(27/255, 158/255, 119/255, Alpha)
loCol <- rgb(117/255, 112/255, 179/255, Alpha)

par(mfcol=c(1,1), mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
tDiff <- max(fmax) - maxT
plot(exRTT$times + tDiff, apply(exRTT$beta,2,mean), type='l', ylim=c(0,0.1), xlim=c(0,max(fmax)), axes=F, xlab="age (Ma)", ylab=expression(paste("svl rate (", beta, ")", sep="")), col=exCol, lwd=1.5)
silent <- sapply(1:length(fbRTT), function(x) lines(fbRTT[[x]]$times + max(fmax) - fmax[x], apply(fbRTT[[x]]$beta, 2, mean), col=loCol))
Xax <- c(-100, 0, 50, 100, 150, 204, 240, 270)
axis(1, at=Xax, labels=as.character(sapply(max(fmax)-Xax, round)))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))
legend("topright", bty="n", legend=c("extant", "fossils"), lty=1, col=c(exCol, loCol))
dev.off()



### LOOK AT NODES THROUGH TIME
sameNode <- function(node, fostree, extree)	{
	taxa <- extree$tip.label
	Tips <- fostree$tip.label[getDescendants(fostree, node)]	
	Tips <- Tips[Tips %in% taxa]
	exnode <- findMRCA(extree, Tips)
	return(c(node, exnode))
}
Alpha <- 0.75
exCol <- rgb(217/255, 95/255, 2/255, Alpha)
hiCol <- rgb(27/255, 158/255, 119/255, Alpha)
loCol <- rgb(117/255, 112/255, 179/255, Alpha)

eHeights <- nodeHeights(exTree)
fHeights <- nodeHeights(as.phylo(fbEdata[[1]]))

nodeNums <- seq(from=Ntip(as.phylo(fbNodes[[1]]))+1, to=nrow(fbNodes[[1]]$edge), by=1)
# some nodes in fossil tree have no extant descendants
nodeMat <- sapply(nodeNums, function(x) try(sameNode(x, fostree=as.phylo(fbNodes[[1]]), extree=ex_nodeStates)))
nodeMat <- nodeMat[sapply(nodeMat, class)!="try-error"]
nodeMat <- do.call(rbind, nodeMat)
nodeMat_save <- nodeMat
nodeMat <- nodeMat[-c(1,2,3),]

setwd(paste(Path_base, "figures", sep=""))
pdf("crownnodes.pdf", height=5, width=5)
par(mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
plot(1, 1, type="n", xlim=c(1, 6), ylim=c(1, 6), xlab="logSVL (extant-only)", ylab="logSVL (fossils+extant)", main="", axes=F)
for (x in unique(nodeMat[,2]))	{
	exNode <- x
	fosNode <- max(nodeMat[which(nodeMat[,2]==x),1])
	exEdge <- which(ex_nodeStates$edge[,2]==exNode)
	fosEdge <- which(fbNodes[[1]]$edge[,2]==fosNode)
	points(ex_nodeStates$edge.length[exEdge], fbNodes[[1]]$edge.length[fosEdge], pch=16, col='gray70',cex=1.5)
}
segments(1,1,5.8,5.8,lty=2)
axis(1, at=c(-10,1,2,3,4,5,6))
axis(2, at=c(-10,1,2,3,4,5,6))
dev.off()



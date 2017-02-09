###
source('~/Documents/BAMMtools_Functions/spindlePlot.R', chdir = TRUE)

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

# Neoteny
z <- list()
Ecor <- c()
Scor <- c()
Empty <- c()
counter <- 1
for (Tree in 1:70)	{
	cat(Tree, "\n")
	setwd(paste(Path_base, "bamm/output-", Tree, sep=""))
	tree <- read.tree("fossilTree.tre")
	load(file=paste(Path_base, "output/neotAnc/neotAnc", Tree, ".RData", sep=""))
	Heights <- nodeHeights(tree)
	rownames(Heights) <- tree$edge[,2]
	Liabs <- apply(ACE$liab[101:1001,], 2, mean)
	Alpha <- 1
	nCols <- c(rgb(254/255,232/255,200/255,Alpha), rgb(253/255,187/255,132/255,Alpha), rgb(227/255, 74/255, 	51/255, Alpha))
	names(nCols) <- c("0", "1", "2")
	Thresh <- apply(ACE$par[10:101,], 2, mean)
	nVec <- rep(1, length(Liabs))
	nVec[Liabs>=Thresh[2]] <- 2
	nVec[Liabs>=Thresh[3]] <- 3
	names(nVec) <- names(Liabs)
	Facul <- names(Liabs)[intersect(which(Liabs > 0), which(Liabs < Thresh[3]))]
	z[[Tree]] <- data.frame(tree=rep(Tree, length(Liabs)), names=names(Liabs), Liabs, nVec)

	# SpEx
	if (file.exists("base_event_data.txt"))	{
	edata <- getEventData(tree, "base_event_data.txt", burnin=0.5)
	mcmc <- read.csv("base_mcmc_out.txt", stringsAsFactors=F)

	Ext <- edata$meanTipMu
	names(Ext) <- edata$tip.label
	Spec <- edata$meanTipLambda
	names(Spec) <- edata$tip.label

	Int <- intersect(names(Ext), names(Liabs))
	Ecor[Tree] <- cor(Ext[Int], Liabs[Int])
	Scor[Tree] <- cor(Spec[Int], Liabs[Int])
	}
	else		{
		Empty[counter] <- Tree
		counter <- counter + 1
	}
}

setwd(paste(Path_base, "/manuscripts/neoteny/", sep=""))
pdf("salamanders.pdf", height=3, width=3)
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.2, 0.5, 0), tck=-0.005, las=1, cex.lab=1.5, bty="n", bg='white', fg='black', col.axis='black', col.lab='black')
plot(Ext[Int], Liabs[Int], pch=16, col=nCols[nVec], xlab="extinction rate", ylab="develomental lability", ylim=c(-15, 20), xlim=c(0, 0.1), axes=F)
axis(1, at=c(-10, -0.05, 0, 0.05, 0.1))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))
legend("topleft", legend=rev(c("Normal", "Fac. neoteny", "Obl. neoteny")), pch=16, col=rev(nCols), text.col=rev(nCols), bty="n", cex=0.75)
#segments(-0.02, 0, 0.2, 0, col=nCols[2], lwd=1.5)
#segments(-0.02, Thresh[3], 0.2, Thresh[3], col=nCols[3], lwd=1.5)
dev.off()

par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.2, 0.5, 0), tck=-0.005, las=1, cex.lab=1.5, bty="n", bg='white', fg='black', col.axis='black', col.lab='black')
SPAT <- 1
EXAT <- 2
LINE <- 0.5
BUFFER <- 0.5
plot(1, 1, xlim=c(SPAT-BUFFER,EXAT+BUFFER), ylim=c(-0.5, 0.5), type="n", ylab="", xlab="", axes=F)
spindlePlot(y=na.omit(Scor), at=SPAT, widthFactor=4, col='gray60')
spindlePlot(y=na.omit(Ecor), at=EXAT, widthFactor=4, col='gray60')
#axis(1, at=c(-10,1,3), labels=F)
axis(2, at=c(-0.5,-0.25,0,0.25,0.5))
mtext("speciation", side=1, line=LINE, at=SPAT)
mtext("extinction", side=1, line=LINE, at=EXAT)
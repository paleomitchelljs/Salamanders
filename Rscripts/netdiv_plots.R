#### Net div calculations
Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)
library(coda)
library(png)

source('~/Documents/BAMMtools_Functions/addScale.R', chdir = TRUE)


plotSil <- function(TipTax, xoffset, yoffset, Tree, Hpad=0, Vpad=0, plotPic=NULL, Download=TRUE, Address=NULL)	{
	require(png)
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (Download == TRUE)	{
		if (is.null(Address))	{
			stop("need to supply a web address if Download == TRUE")
		}
		picFile <- download.file(Address, dest=paste(TipTax, ".png", sep=""), mode="wb")
	}
	if (is.null(plotPic))	{
		plotPic <- readPNG(paste(TipTax, ".png", sep=""))
	}
	
	tip <- which(Tree$tip.label==TipTax)
	xloc <- lastPP$xx[tip] + Hpad
	yloc <- lastPP$yy[tip] + Vpad

	rasterImage(plotPic, xloc-xoffset, yloc-yoffset, xloc+xoffset, yloc+yoffset, xpd=NA)
}
addFamily <- function(Taxa, tree, Family="", Tip=TRUE, Cex=2, centering=T, Top=T, ncex=1.25, ntxt="")	{
	Names <- unique(c(tree$tip.label[grep(Taxa[1], tree$tip.label)], tree$tip.label[grep(Taxa[2], tree$tip.label)], tree$tip.label[grep(Taxa[3], tree$tip.label)]))
	Node <- findMRCA(as.phylo(tree), tips=Names)
	allNames <- na.omit(tree$tip.label[getDescendants(as.phylo(tree), Node)])

	if (Tip)	{
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
	else if (!Tip)	{
		Edge <- which(tree$edge[,2] == Node)
		Age <- tree$begin[Edge]
		fromRoot <- max(tree$end) - Age

		lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		subedge <- lastPP$edge[Edge, , drop = FALSE]
		xx <- lastPP$xx
		yy <- lastPP$yy
		YY <- yy[subedge[, 2]]
		XX <- xx[subedge[,1]] + (0.1 * xx[subedge[,1]])
		points(XX, YY, pch=21, bg='white', cex=ncex)
		points(XX, YY, pch=ntxt, col='black', cex=0.333*ncex)		
	}
}

makePlot <- function(x, y, Lims=c(-0.02, 0.1), Fit=T, cex=1.25, xLab="", yLab="", Linex=2, Liney=2.5, Label="")	{
	options(scipen=6)
	plot(1, 1, xlim=Lims, ylim=Lims, axes=F, xlab="", ylab="", type="n")
	axis(1, at=c(-100, pretty(Lims)))
	axis(2, at=c(-100, pretty(Lims)))
	segments(Lims[1], Lims[1], Lims[2], Lims[2], lty=2)
	if (Fit)	{
		Model <- lm(y~x)
		segments(0, Model$coefficients[1], 
				max(x, na.rm=T), Model$coefficients[2]*max(x, na.rm=T))
	}
	points(x, y, pch=21, bg='gray80', cex=cex)
#	legend("top", bty="n", 
#	legend=paste("r = ", round(cor(x, y, use="complete", method="pearson"), digits=2), sep=""))
	mtext(xLab, side=1, line=Linex)
	mtext(yLab, side=2, line=Liney, las=0)
	mtext(Label, side=2, line=3, at=Lims[2], cex=1.25, xpd=NA)
}


## Read in dataframe
totalMat <- read.csv(paste(Path_base, "datafiles/netdiv.csv", sep=""), stringsAsFactors=F)

## Parse dataframe
extant <- totalMat[totalMat$data=="extant",]
exAll <- extant[extant$pleth == 1, ]
npAll <- extant[extant$pleth == 0, ]

fossil <- totalMat[totalMat$data == "fossil",]
fossil_all <- fossil[fossil$pleth == 1, ]
fossil_np <- fossil[fossil$pleth == 0, ]

fnpStems <- tapply(fossil_np$gcr_stem, fossil_np$clade, mean)
fnpCrowns <- tapply(fossil_np$gcr_crown, fossil_np$clade, mean)

fStems <- tapply(fossil_all$gcr_stem, fossil_all$clade, mean)
fCrowns <- tapply(fossil_all$gcr_crown, fossil_all$clade, mean)

## Make plots
pdf(paste(Path_base, "figures/no_pleth_compPlot_hi.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fnpStems[npAll$clade], npAll$rstem_hi.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fnpCrowns[npAll$clade], npAll$rcrown_hi.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()

pdf(paste(Path_base, "figures/no_pleth_compPlot.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fnpStems[npAll$clade], npAll$rstem.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fnpCrowns[npAll$clade], npAll$rcrown.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()

pdf(paste(Path_base, "figures/no_pleth_compPlot_lo.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fnpStems[npAll$clade], npAll$rstem_lo.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fnpCrowns[npAll$clade], npAll$rcrown_lo.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()


pdf(paste(Path_base, "figures/all_compPlot_hi.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fStems[exAll$clade], exAll$rstem_hi.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fCrowns[exAll$clade], exAll$rcrown_hi.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()

pdf(paste(Path_base, "figures/all_compPlot.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fStems[exAll$clade], exAll$rstem.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fCrowns[exAll$clade], exAll$rcrown.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()

pdf(paste(Path_base, "figures/all_compPlot_lo.pdf", sep=""), height=4, width=8)
par(mfrow=c(1, 2), mar=c(2.5,4,0.5,0.5), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01, cex.axis=0.8)
makePlot(fStems[exAll$clade], exAll$rstem_lo.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A", Linex=1.5)
makePlot(fCrowns[exAll$clade], exAll$rcrown_lo.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="B", Linex=1.5)
dev.off()


pdf(paste(Path_base, "figures/netdiv_comp_all.pdf", sep=""), height=2.5, width=10)
par(mfrow=c(1, 4), mar=c(4,5,1,1), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01)
makePlot(fStems[exAll$clade], exAll$rstem.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A")
makePlot(fStems[exAll$clade], exAll$gcr_stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (BAMM)", Label="B")
makePlot(fCrowns[exAll$clade], exAll$rcrown.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="C")
makePlot(fCrowns[exAll$clade], exAll$gcr_crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (BAMM)", Label="D")
dev.off()

pdf(paste(Path_base, "figures/netdiv_comp_no_pleth.pdf", sep=""), height=2.5, width=10)
par(mfrow=c(1, 4), mar=c(4,5,1,1), bty="n", las=1, mgp=c(2,0.5,0), tck=-0.01)
makePlot(fnpStems[npAll$clade], npAll$rstem.stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (calc.)", Label="A")
makePlot(fnpStems[npAll$clade], npAll$gcr_stem, xLab="stem net div. (fossilBAMM)", yLab="stem net div. (BAMM)", Label="B")
makePlot(fnpCrowns[npAll$clade], npAll$rcrown.crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (calc.)", Label="C")
makePlot(fnpCrowns[npAll$clade], npAll$gcr_crown, xLab="crown net div. (fossilBAMM)", yLab="crown net div. (BAMM)", Label="D")
dev.off()

###########################################
exRTT <- getRateThroughTimeMatrix(exEdata)
npRTT <- getRateThroughTimeMatrix(npEdata)

fEdata <- subtreeBAMM(fbEdata[[1]], tips=npEdata$tip.label)
fnpRTT <- getRateThroughTimeMatrix(fbEdata[[1]])

rttPoly <- function(RTT, opt="lam", cPoly="gray70", cLine='black', Lty=2, Lwd=2, LoP=0.05, HiP=0.95)		{
	Times <- max(RTT$times) - RTT$times

	if (opt == "lam")	{
		Mean <- apply(RTT$lambda, 2, median)
		Lo <- apply(RTT$lambda, 2, quantile, probs=LoP)
		Hi <- apply(RTT$lambda, 2, quantile, probs=HiP)
	}
	if (opt == "mu")	{
		Mean <- apply(RTT$mu, 2, median)
		Lo <- apply(RTT$mu, 2, quantile, probs=LoP)
		Hi <- apply(RTT$mu, 2, quantile, probs=HiP)
	}
	if (opt == "net")	{
		Vec <- RTT$lambda - RTT$mu
		Mean <- apply(Vec, 2, median)
		Lo <- apply(Vec, 2, quantile, probs=LoP)
		Hi <- apply(Vec, 2, quantile, probs=HiP)
	}
	if (opt == "rel")	{
		Vec <- RTT$mu / RTT$lambda
		Mean <- apply(Vec, 2, mean)
		Lo <- apply(Vec, 2, quantile, probs=LoP)
		Hi <- apply(Vec, 2, quantile, probs=HiP)
	}
	if (opt == "turn")	{
		Vec <- RTT$mu + RTT$lambda
		Mean <- apply(Vec, 2, mean)
		Lo <- apply(Vec, 2, quantile, probs=LoP)
		Hi <- apply(Vec, 2, quantile, probs=HiP)
	}	
	Xvec <- c(Times[1], Times, Times[length(Times)], rev(Times))
	Yvec <- c(Mean[1], Lo, Mean[length(Mean)], rev(Hi))
	polygon(Xvec, Yvec, border=F, col=cPoly)
	lines(Times, Mean, col=cLine, lty=Lty, lwd=Lwd)
	
}

Cols <- c(rgb(166/255, 206/255, 227/255, 1), 
		  rgb(31/255, 120/255, 180/255, 1), 
		  rgb(178/255, 223/255, 138/255, 1), 
		  rgb(51/255, 160/255, 44/255, 1), 
		  rgb(251/255, 154/255, 153/255, 1))

Cols <- c('gray70', 'gray80', 'gray90',
		  'gray60', 'gray60')

pdf(paste(Path_base, "figures/ExtRate.pdf", sep=""), height=3.5, width=7)
layout(mat=matrix(c(1, 2, 3), byrow=T, nrow=1), widths=c(1,0.5,0.5), heights=c(2))
par(mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
plot(1, 1, xlim=c(225, 0), ylim=c(0, 0.06), type="n", axes=F, bty="n", xlab="", ylab="")
addScale(c(225, 200, 175, 150, 125, 100, 75, 50, 25, 0), units="Stage", YLim=c(-10, 0))
axis(1, at=seq(from=250, to=0, by=-25))
axis(2, at=c(-10, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06))

mtext("age (Ma)", side=1, line=2, cex=1.25)
mtext("extinction rate", side=2, line=3, las=0, cex=1.25)
mtext("A", side=2, at=0.06, line=3)
mtext("B", side=2, at=0.06, line=-21.5, xpd=NA)
legend("topleft", legend=c("extant & extinct", "extant"), lty=c(2, 3), bty="n", cex=1.25, lwd=2)

#polygon(x=c(250, 250, 201, 201), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))
#polygon(x=c(145, 145, 66, 66), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))

rttPoly(fnpRTT, opt="mu", cPoly=Cols[1], LoP=0.25, HiP=0.75)
rttPoly(npRTT, opt="mu", cPoly=Cols[5], Lty=3, LoP=0.25, HiP=0.75)

par(mar=c(0,0,0,1.5))
P <- plot(fEdata, breaksmethod="jenks", spex="e", lwd=2)
addBAMMlegend(P, location="topleft", cex.axis=0.9, labelDist=0.5, shortFrac=0.04)

addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=npEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=npEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=npEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=npEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=npEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Salamandrella"), tree=npEdata, Tip=F, ntxt="6", ncex=2)
addFamily(c("Cryptobranchus", "Andrias", "Andrias"), tree=npEdata, Tip=F, ntxt="7", ncex=2)
#addFamily(c("Siren", "Pseudobranchus", "Kababisha"), tree=fbEdata[[1]], Tip=F, ntxt="5", ncex=2)
#addFamily(c("Hynobiius", "Onychodactylus", "Iridotriton"), tree=fbEdata[[1]], Tip=F, ntxt="6", ncex=2)
#addFamily(c("Cryptobranchus", "Andrias", "Eoscapherpeton"), tree=fbEdata[[1]], Tip=F, ntxt="7", ncex=2)


par(mar=c(0,1.5,0,0))
plot(npEdata, colorbreaks=P$colorbreaks, spex="e", direction="leftwards", lwd=2)

addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=npEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=npEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=npEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=npEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=npEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Salamandrella"), tree=npEdata, Tip=F, ntxt="6", ncex=2)
addFamily(c("Cryptobranchus", "Andrias", "Andrias"), tree=npEdata, Tip=F, ntxt="7", ncex=2)
dev.off()


pdf(paste(Path_base, "figures/SpecRate.pdf", sep=""), height=3.5, width=7)
layout(mat=matrix(c(1, 2, 3), byrow=T, nrow=1), widths=c(1,0.5,0.5), heights=c(2))
par(mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
plot(1, 1, xlim=c(225, 0), ylim=c(0, 0.1), type="n", axes=F, bty="n", xlab="", ylab="")
addScale(c(225, 200, 175, 150, 125, 100, 75, 50, 25, 0), units="Stage", YLim=c(-10, 0))
axis(1, at=seq(from=250, to=0, by=-25))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))

mtext("age (Ma)", side=1, line=2, cex=1.25)
mtext("speciation rate", side=2, line=3, las=0, cex=1.25)
mtext("A", side=2, at=0.1, line=3)
mtext("B", side=2, at=0.1, line=-21.5, xpd=NA)
legend("topleft", legend=c("extant & extinct", "extant"), lty=c(2, 3), bty="n", cex=1.25, lwd=2)
#polygon(x=c(250, 250, 201, 201), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))
#polygon(x=c(145, 145, 66, 66), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))

rttPoly(fnpRTT, opt="lam", cPoly=Cols[1], LoP=0.25, HiP=0.75)
rttPoly(npRTT, opt="lam", cPoly=Cols[5], Lty=3, LoP=0.25, HiP=0.75)

par(mar=c(0,0,0,1.5))
P <- plot(fEdata, breaksmethod="jenks", spex="s", lwd=2)
addBAMMlegend(P, location="topleft", cex.axis=0.9, labelDist=0.5, shortFrac=0.04)

addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=npEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=npEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=npEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=npEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=npEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Salamandrella"), tree=npEdata, Tip=F, ntxt="6", ncex=2)
addFamily(c("Cryptobranchus", "Andrias", "Andrias"), tree=npEdata, Tip=F, ntxt="7", ncex=2)
#addFamily(c("Siren", "Pseudobranchus", "Kababisha"), tree=fbEdata[[1]], Tip=F, ntxt="5", ncex=2)
#addFamily(c("Hynobiius", "Onychodactylus", "Iridotriton"), tree=fbEdata[[1]], Tip=F, ntxt="6", ncex=2)
#addFamily(c("Cryptobranchus", "Andrias", "Eoscapherpeton"), tree=fbEdata[[1]], Tip=F, ntxt="7", ncex=2)

Image <- readPNG(paste(Path_base, "figures/images/newt2.png", sep=""))
plotSil("Notophthalmus_viridescens", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/ambystoma.png", sep=""))
plotSil("Ambystoma_tigrinum", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/siren.png", sep=""))
plotSil("Siren_lacertina", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/hellbender.png", sep=""))
plotSil("Cryptobranchus_alleganiensis", xoffset=20, yoffset=2, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

par(mar=c(0,1.5,0,0))
plot(npEdata, colorbreaks=P$colorbreaks, spex="s", direction="leftwards", lwd=2)

addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=npEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=npEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=npEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=npEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=npEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Salamandrella"), tree=npEdata, Tip=F, ntxt="6", ncex=2)
addFamily(c("Cryptobranchus", "Andrias", "Andrias"), tree=npEdata, Tip=F, ntxt="7", ncex=2)

dev.off()


###########################################
exRTT <- getRateThroughTimeMatrix(exEdata)

fRTT <- getRateThroughTimeMatrix(fbEdata_all[[1]])
fEdata <- subtreeBAMM(fbEdata_all[[1]], tips=exEdata$tip.label)

rttPoly <- function(RTT, opt="lam", cPoly="gray70", cLine='black', Lty=2, Lwd=2, LoP=0.05, HiP=0.95)		{
	Times <- max(RTT$times) - RTT$times

	if (opt == "lam")	{
		Mean <- apply(RTT$lambda, 2, median)
		Lo <- apply(RTT$lambda, 2, quantile, probs=LoP)
		Hi <- apply(RTT$lambda, 2, quantile, probs=HiP)
	}
	if (opt == "mu")	{
		Mean <- apply(RTT$mu, 2, median)
		Lo <- apply(RTT$mu, 2, quantile, probs=LoP)
		Hi <- apply(RTT$mu, 2, quantile, probs=HiP)
	}
	if (opt == "net")	{
		Vec <- RTT$lambda - RTT$mu
		Mean <- apply(Vec, 2, median)
		Lo <- apply(Vec, 2, quantile, probs=LoP)
		Hi <- apply(Vec, 2, quantile, probs=HiP)
	}
	if (opt == "rel")	{
		Vec <- RTT$mu / RTT$lambda
		Mean <- apply(Vec, 2, mean)
		Lo <- apply(Vec, 2, quantile, probs=LoP)
		Hi <- apply(Vec, 2, quantile, probs=HiP)
	}
	Xvec <- c(Times[1], Times, Times[length(Times)], rev(Times))
	Yvec <- c(Mean[1], Lo, Mean[length(Mean)], rev(Hi))
	polygon(Xvec, Yvec, border=F, col=cPoly)
	lines(Times, Mean, col=cLine, lty=Lty, lwd=Lwd)
	
}

Cols <- c(rgb(166/255, 206/255, 227/255, 1), 
		  rgb(31/255, 120/255, 180/255, 1), 
		  rgb(178/255, 223/255, 138/255, 1), 
		  rgb(51/255, 160/255, 44/255, 1), 
		  rgb(251/255, 154/255, 153/255, 1))

Cols <- c('gray70', 'gray80', 'gray90',
		  'gray60', 'gray60')

pdf(paste(Path_base, "figures/ExtRate_all.pdf", sep=""), height=3.5, width=7)
layout(mat=matrix(c(1, 2, 3), byrow=T, nrow=1), widths=c(1,0.5,0.5), heights=c(2))
par(mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
plot(1, 1, xlim=c(225, 0), ylim=c(0, 0.06), type="n", axes=F, bty="n", xlab="", ylab="")
addScale(c(225, 200, 175, 150, 125, 100, 75, 50, 25, 0), units="Stage", YLim=c(-10, 0))
axis(1, at=seq(from=250, to=0, by=-25))
axis(2, at=c(-10, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06))

mtext("age (Ma)", side=1, line=2, cex=1.25)
mtext("extinction rate", side=2, line=3, las=0, cex=1.25)
mtext("A", side=2, at=0.06, line=3)
mtext("B", side=2, at=0.06, line=-21.5, xpd=NA)
legend("topleft", legend=c("extant & extinct", "extant"), lty=c(2, 3), bty="n", cex=1.25, lwd=2)

#polygon(x=c(250, 250, 201, 201), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))
#polygon(x=c(145, 145, 66, 66), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))

rttPoly(fRTT, opt="mu", cPoly=Cols[1], LoP=0.25, HiP=0.75)
rttPoly(npRTT, opt="mu", cPoly=Cols[5], Lty=3, LoP=0.25, HiP=0.75)

par(mar=c(0,0,0,1.5))
P <- plot(fEdata, breaksmethod="jenks", spex="e", lwd=2)
addBAMMlegend(P, location="topleft", cex.axis=0.9, labelDist=0.5, shortFrac=0.04)

addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=exEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=exEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=exEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=exEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=exEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Andrias"), tree=exEdata, Tip=F, ntxt="6", ncex=2)


par(mar=c(0,1.5,0,0))
plot(exEdata, colorbreaks=P$colorbreaks, spex="e", direction="leftwards", lwd=2)

addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=exEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=exEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=exEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=exEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=exEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Andrias"), tree=exEdata, Tip=F, ntxt="6", ncex=2)
dev.off()


pdf(paste(Path_base, "figures/SpecRate_all.pdf", sep=""), height=3.5, width=7)
layout(mat=matrix(c(1, 2, 3), byrow=T, nrow=1), widths=c(1,0.5,0.5), heights=c(2))
par(mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
plot(1, 1, xlim=c(225, 0), ylim=c(0, 0.1), type="n", axes=F, bty="n", xlab="", ylab="")
addScale(c(225, 200, 175, 150, 125, 100, 75, 50, 25, 0), units="Stage", YLim=c(-10, 0))
axis(1, at=seq(from=250, to=0, by=-25))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1))

mtext("age (Ma)", side=1, line=2, cex=1.25)
mtext("speciation rate", side=2, line=3, las=0, cex=1.25)
mtext("A", side=2, at=0.1, line=3)
mtext("B", side=2, at=0.1, line=-21.5, xpd=NA)
legend("topleft", legend=c("extant & extinct", "extant"), lty=c(2, 3), bty="n", cex=1.25, lwd=2)
#polygon(x=c(250, 250, 201, 201), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))
#polygon(x=c(145, 145, 66, 66), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))

rttPoly(fRTT, opt="lam", cPoly=Cols[1], LoP=0.25, HiP=0.75)
rttPoly(npRTT, opt="lam", cPoly=Cols[5], Lty=3, LoP=0.25, HiP=0.75)

par(mar=c(0,0,0,1.5))
P <- plot(fEdata, breaksmethod="jenks", spex="s", lwd=2)
addBAMMlegend(P, location="topleft", cex.axis=0.9, labelDist=0.5, shortFrac=0.04)

addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=exEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=exEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=exEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=exEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=exEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Andrias"), tree=exEdata, Tip=F, ntxt="6", ncex=2)

Image <- readPNG(paste(Path_base, "figures/images/newt2.png", sep=""))
plotSil("Notophthalmus_viridescens", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/ambystoma.png", sep=""))
plotSil("Ambystoma_tigrinum", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/siren.png", sep=""))
plotSil("Siren_lacertina", xoffset=20, yoffset=5, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

Image <- readPNG(paste(Path_base, "figures/images/hellbender.png", sep=""))
plotSil("Cryptobranchus_alleganiensis", xoffset=20, yoffset=2, Tree=fbEdata[[1]], plotPic=Image, Vpad=0, Hpad=30, Download=F)

par(mar=c(0,1.5,0,0))
plot(exEdata, colorbreaks=P$colorbreaks, spex="s", direction="leftwards", lwd=2)

addFamily(c("Rhyacotriton", "Gyrinophilus", "Amphiuma"), tree=exEdata, Tip=F, ntxt="1", ncex=2)
addFamily(c("Necturus", "Necturus", "Proteus"), tree=exEdata, Tip=F, ntxt="2", ncex=2)
addFamily(c("Notophthalmus", "Neurergus", "Salamandrina"), tree=exEdata, Tip=F, ntxt="3", ncex=2)
addFamily(c("Ambystoma", "Dicamptodon", "Ambystoma"), tree=exEdata, Tip=F, ntxt="4", ncex=2)
addFamily(c("Siren", "Pseudobranchus", "Pseudobranchus"), tree=exEdata, Tip=F, ntxt="5", ncex=2)
addFamily(c("Hynobiius", "Onychodactylus", "Andrias"), tree=exEdata, Tip=F, ntxt="6", ncex=2)

dev.off()


pdf(paste(Path_base, "figures/netdiv_polar.pdf", sep=""), height=3.5, width=7)
LWD <- 2
par(mfrow=c(1, 2), mar=c(1, 1, 1, 1), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
x <- plot(fbEdata[[1]], method="polar", spex="netdiv", breaksmethod="jenks", lwd=LWD)
addBAMMlegend(x, side=4, location="topright", labelDist=0.25, cex.axis=0.75)
text(1.15, 1, "net div.", srt=90)
plot(npEdata, method="polar", spex="netdiv", colorbreaks=x$colorbreaks, lwd=LWD)
dev.off()

pdf(paste(Path_base, "figures/netdiv_polar_all.pdf", sep=""), height=3.5, width=7)
LWD <- 2
par(mfrow=c(1, 2), mar=c(1, 1, 1, 1), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
x <- plot(fbEdata_all[[1]], method="polar", spex="netdiv", breaksmethod="jenks", lwd=LWD)
addBAMMlegend(x, side=4, location="topright", labelDist=0.25, cex.axis=0.75)
text(1.15, 1, "net div.", srt=90)
plot(exEdata, method="polar", spex="netdiv", colorbreaks=x$colorbreaks, lwd=LWD)
dev.off()


par(mar=c(4, 5, 0.5, 0.5), las=1, mgp=c(2, 0.3, 0), tck=-0.01)
plot(1, 1, xlim=c(225, 0), ylim=c(0, 0.14), type="n", axes=F, bty="n", xlab="", ylab="")
addScale(c(225, 200, 175, 150, 125, 100, 75, 50, 25, 0), units="Stage", YLim=c(-10, 0))
axis(1, at=seq(from=250, to=0, by=-25))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14))

mtext("age (Ma)", side=1, line=2, cex=1.25)
mtext("turnover rate", side=2, line=3, las=0, cex=1.25)
mtext("A", side=2, at=0.14, line=3)
legend("topleft", legend=c("extant & extinct", "extant"), lty=c(2, 3), bty="n", cex=1.25, lwd=2)

#polygon(x=c(250, 250, 201, 201), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))
#polygon(x=c(145, 145, 66, 66), y=c(0, 0.06, 0.06, 0), border=NA, col=rgb(0, 0, 0, 0.05))

rttPoly(fnpRTT, opt="turn", cPoly=Cols[1], LoP=0.25, HiP=0.75)
rttPoly(npRTT, opt="turn", cPoly=Cols[5], Lty=3, LoP=0.25, HiP=0.75)

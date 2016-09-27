source('~/Documents/BAMMtools_Functions/placeFossil.R', chdir = TRUE)
library(phytools)
library(TreePar)

setwd("~/Documents/Salamanders/datafiles")
extant <- read.csv("salamanders.csv", stringsAsFactors=FALSE)
fossil <- read.csv("fossilSalamanders.csv", stringsAsFactors=FALSE)
salTree <- read.tree("fossilTree.tre")

meanSVL <- tapply(extant$SVL, extant$Species, mean, na.rm=T)
meanHW <- tapply(extant$Head.width, extant$Species, mean, na.rm=T)
f_meanSVL <- tapply(fossil$svl, fossilSp, mean)
f_meanHW <- tapply(fossil$head_width, fossilSp, mean)

par(mfrow=c(2,2), mar=c(3,4,0.5,0.5), las=1, bty="n")

plot(log(meanSVL), log(meanHW))
text(log(meanSVL), log(meanHW), names(meanSVL), pos=4, cex=0.5)
points(log(f_meanSVL), log(f_meanHW), pch=15, col='red')

extantSizeDist <- density(meanSVL, na.rm=T, adjust=4)
fossilSizeDist <- density(f_meanSVL, na.rm=T, adjust=0.5)

yDelim <- c(0, 0.06)
xDelim <- c(0, 500)

plot(1, 1, type="n", xlab="snout-vent length", ylab="density", xlim=xDelim, ylim=yDelim, axes=FALSE)
polygon(extantSizeDist, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0))
polygon(fossilSizeDist, col=rgb(1,0,0,0.5), border=rgb(1,0,0,0))
axis(1, at=c(-100, pretty(xDelim)))
axis(2, at=c(-100, pretty(yDelim)))

Missing <- c(254, 500, 359)
meanSVL[c("Cryptobranchus_alleganiensis", "Andrias_davidianus", "Andrias_japonicus")] <- Missing

Data <- c(meanSVL, f_meanSVL)
salTree2 <- drop.tip(salTree, setdiff(salTree$tip.label, names(Data)))
ancs <- fastAnc(salTree2, Data[salTree2$tip.label])
Heights <- nodeHeights(salTree2)
Vec <- c(Data[salTree2$tip.label], ancs)
plot(Heights[,2], Vec[salTree2$edge[,2]], axes=F, xlab="age (Ma)", ylab="SVL")
Times <- c(0, 50, 100, 150, 200, 250, 300)
Ages <- sapply(max(Heights) - Times, round) + 2
axis(1, at=Times, labels=Ages)
axis(2, at=c(-1000, 0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500))


BMfit <- fitContinuous(salTree2, Data[salTree2$tip.label])
EBfit <- fitContinuous(salTree2, Data[salTree2$tip.label], model="EB")
OUfit <- fitContinuous(salTree2, Data[salTree2$tip.label], model="OU")


phenogram(salTree2, Data[salTree2$tip.label], fsize=1)
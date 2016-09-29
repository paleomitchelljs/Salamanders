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

# Read in neoteny scores
setwd(paste(Path_base, "neoteny", sep=""))
neoteny <- read.csv("neoteny.csv", stringsAsFactors=F)
spNames <- apply(neoteny, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
neoteny <- neoteny[,3]
names(neoteny) <- spNames

Alpha <- 1
Cols <- c(rgb(254/255,232/255,200/255,Alpha), rgb(253/255,187/255,132/255,Alpha), rgb(227/255, 74/255, 51/255, Alpha))
names(Cols) <- c("0", "1", "2")

# Load trees
trees <- list()
for (count in 1:5)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	tree <- read.tree("fossilTree.tre")
	trees[[count]] <- tree
	Ngen <- 1e6
	ACE <- ancThresh(tree, neoteny[tree$tip.label], ngen=Ngen, control=list(sample=Ngen/1e3, burnin=Ngen*0.2, piecol=Cols))
	setwd(paste(Path_base, "output", sep=""))
	save(ACE, file=paste(Path_base, "output/neotAnc/neotAnc", count, ".RData", sep=""))
}

setwd(paste(Path_base, "bamm/output-1", sep=""))
tree <- read.tree("fossilTree.tre")
load(file=paste(Path_base, "output/neotAnc/neotAnc1.RData", sep=""))


setwd(paste(Path_base, "check_tree_files", sep=""))
pdf("neoteny_check.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
plot(tree, show.tip.label=T, label.offset=0.25, cex=0.15)
nodelabels(pie=ACE$ace, piecol=Cols, cex=0.25)
tiplabels(pch=16, col=Cols[as.character(neoteny[tree$tip.label])], cex=0.25)
dev.off()

setwd(paste(Path_base, "figures", sep=""))
pdf("neoteny.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,14))
plot(tree, show.tip.label=F)
nodelabels(pie=ACE$ace, piecol=Cols, cex=0.5)
tiplabels(pch=16, col=Cols[as.character(neoteny[tree$tip.label])], cex=0.5)
# Add families
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=F)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
legend("bottomleft", legend=c("Normal", "Fac. neoteny", "Obl. neoteny"), pch=16, col=Cols, pt.cex=2, cex=2, bty='n')
dev.off()


# Use liab and par to plot up Amphiumas, Sirens, Cryptobranchids and Ambystomatids
Burn <- floor(0.1 * nrow(ACE$liab))
Vals <-  apply(ACE$liab[Burn:nrow(ACE$liab),],2,mean)

trackNode <- function(tree, tip)		{
	tipNum <- which(tree$tip.label == tip)
	counter <- 1
	ancvec <- c()
	ancvec[counter] <- tree$edge[which(tree$edge[,2]==tipNum),1]
	while (ancvec[counter] > Ntip(tree)+2)	{
		ancvec[counter+1] <- tree$edge[which(tree$edge[,2]==ancvec[counter]),1]
		counter <- counter + 1
	}
	return(ancvec)
}

Heights <- nodeHeights(tree)
rownames(Heights) <- tree$edge[,2]
Thresh <- apply(ACE$par[10:101,],2,mean)

setwd(paste(Path_base, "figures", sep=""))
pdf("neotEvol.pdf", height=8, width=16)
# Sirenidae
par(mar=c(4,4,1,1), oma=c(0,0,0,3), las=1, mgp=c(2.2, 0.5, 0), tck=-0.005, bty="n", cex.lab=1.5, cex=1.5)
plot(0, 0, xlim=c(max(Heights), -50), ylim=c(-10, 10), xlab="time", ylab="liability", type="n", axes=F)

#Amphiuma
Names <- tree$tip.label[grep("Amphiuma", tree$tip.label)]
tipScores <- c()
counter <- 1
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
	tipScores[counter] <- zL[which(zH==max(zH))]
	counter <- counter+1
}
text(0, mean(tipScores), pos=4, labels="Amphiuma", xpd=NA, cex=1.5, font=3)

# Ambystoma
Names <- tree$tip.label[grep("Ambystoma", tree$tip.label)]
tipScores <- c()
counter <- 1
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
	tipScores[counter] <- zL[which(zH==max(zH))]
	counter <- counter+1
}
text(0, mean(tipScores), pos=4, labels="Ambystoma", xpd=NA, cex=1.5, font=3)

# Rhyacotriton
Names <- tree$tip.label[grep("Rhyacotriton", tree$tip.label)]
tipScores <- c()
counter <- 1
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
	tipScores[counter] <- zL[which(zH==max(zH))]
	counter <- counter+1
}
text(0, mean(tipScores), pos=4, labels="Rhyacotriton", xpd=NA, cex=1.5, font=3)

# Eurycea
#Names <- tree$tip.label[grep("Eurycea", tree$tip.label)]
#text(0, mean(tipScores), pos=4, labels="Eurycea", xpd=NA, cex=1.5, font=3)
#Names <- tree$tip.label[c(grep("Siren", tree$tip.label), grep("Pseudobranchus", tree$tip.label))]
#text(0, mean(tipScores), pos=4, labels="Sirenidae", xpd=NA, cex=1.5)

# Finish plot
axis(1, at=c(500, 250, 200, 150, 100, 50, 0))
axis(2, at=c(-100, -10, -5, 0, 5, 10))
segments(266, Thresh["1"], 0, Thresh["1"], col=Cols[3], lwd=1.5)
segments(266, Thresh["0"], 0, Thresh["0"], col=Cols[2], lwd=1.5)
points(270, -0.55, pch=21, bg=Cols[1])
text(270, -0.55, labels="Normal", pos=4)
points(270, 1.2, pch=21, bg=Cols[3])
text(270, 1.2, labels="Obl. neoteny", pos=4)
dev.off()




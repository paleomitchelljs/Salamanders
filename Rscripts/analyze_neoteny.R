Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)

# Read in neoteny scores
setwd(paste(Path_base, "neoteny", sep=""))
neoteny <- read.csv("neoteny.csv", stringsAsFactors=F)
spNames <- apply(neoteny, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
neoteny <- neoteny[,3]
names(neoteny) <- spNames

# Load trees
trees <- list()
for (count in 1:5)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	tree <- read.tree("fossilTree.tre")
	trees[[count]] <- tree
	ACE <- ancThresh(tree, neoteny[tree$tip.label], ngen=1e6, control=list(sample=1e4, burnin=2e5))
	setwd(paste(Path_base, "output", sep=""))
	save(ACE, file=paste(Path_base, "output/neotAnc/neotAnc", count, ".RData", sep=""))
}

tree <- trees[[1]]

Alpha <- 1
Cols <- c(rgb(254/255,232/255,200/255,Alpha), rgb(253/255,187/255,132/255,Alpha), rgb(227/255, 74/255, 51/255, Alpha))
names(Cols) <- c("0", "1", "2")

setwd(paste(Path_base, "figures", sep=""))
pdf("neoteny.pdf", height=10, width=10)
plot(tree, show.tip.label=F)
nodelabels(pie=ACE$ace, piecol=Cols, cex=0.25)
tiplabels(pch=16, col=Cols[as.character(neoteny[tree$tip.label])])
# Add families
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


par(mfcol=c(2,2))
# Amphiuma
plot(0, 0, xlim=c(max(Heights), 0), ylim=c(range(Vals)), xlab="time", ylab="liability", type="n")
Names <- tree$tip.label[grep("Amphiuma", tree$tip.label)]
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
}

# Sirenidae
plot(0, 0, xlim=c(max(Heights), 0), ylim=c(range(Vals)), xlab="time", ylab="liability", type="n")
Names <- tree$tip.label[c(grep("Siren", tree$tip.label), grep("Pseudobranchus", tree$tip.label))]
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
}

# Ambystoma
plot(0, 0, xlim=c(max(Heights), 0), ylim=c(range(Vals)), xlab="time", ylab="liability", type="n")
Names <- tree$tip.label[grep("Ambystoma", tree$tip.label)]
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
}

# Eurycea
plot(0, 0, xlim=c(max(Heights), 0), ylim=c(range(Vals)), xlab="time", ylab="liability", type="n")
Names <- tree$tip.label[grep("Eurycea", tree$tip.label)]
for (x in Names)		{
	Name <- x
	z <- trackNode(tree, Name)
	zH <- Heights[as.character(z),2]
	zL <- Vals[c(Name, as.character(z))]
	ID <- rep(1, length(zL))
	ID[which(zL>mean(ACE$par[,2]))] <- 2
	ID[which(zL>mean(ACE$par[,3]))] <- 3
	points(max(Heights)-c(max(Heights),zH), zL, pch=21, bg=Cols[ID], type='b')
}






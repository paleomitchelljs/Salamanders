Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)

# Load trees
trees <- list()
for (count in 1:100)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	trees[[count]] <- read.tree("fossilTree.tre")
}

# Read in neoteny scores
setwd(paste(Path_base, "neoteny", sep=""))
neoteny <- read.csv("neoteny.csv", stringsAsFactors=F)
spNames <- apply(neoteny, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
neoteny <- neoteny[,3]
names(neoteny) <- spNames


tree <- trees[[1]]
ACE <- ancThresh(tree, neoteny[tree$tip.label], ngen=1e6)

Alpha <- 1
Cols <- c(rgb(254/255,232/255,200/255,Alpha), rgb(253/255,187/255,132/255,Alpha), rgb(227/255, 74/255, 51/255, Alpha))
names(Cols) <- c("0", "1", "2")
plot(tree, show.tip.label=F)
nodelabels(pie=ACE$ace, piecol=Cols, cex=0.25)
tiplabels(pch=16, col=Cols[as.character(neoteny[tree$tip.label])])

# Use liab and par to plot up Amphiumas, Sirens, Cryptobranchids and Ambystomatids

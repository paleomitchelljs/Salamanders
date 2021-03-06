Path_base <- "/scratch/drabosky_flux/jonsmitc/Salamanders/"
library(phytools)

# Read in neoteny scores
setwd(paste(Path_base, "neoteny", sep=""))
neoteny <- read.csv("neoteny.csv", stringsAsFactors=F)
spNames <- apply(neoteny, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
neoteny <- neoteny[,3]
names(neoteny) <- spNames

# extant tree
tree <- read.tree(paste(Path_base, "bamm/extant_only/modernSal.tre", sep=""))

Alpha <- 1
Cols <- c(rgb(254/255, 237/255, 222/255, Alpha), rgb(253/255, 190/255, 133/255, Alpha), rgb(253/255, 141/255, 60/255, Alpha), rgb(217/255, 71/255, 1/255, Alpha))
names(Cols) <- c("-1", "0", "1", "2")
PrA <- matrix(rep(0.01666667, length(Cols)), nrow=1)
colnames(PrA) <- names(Cols)
rownames(PrA) <- Ntip(tree) + 1
PrA[1,"0"] <- 1 - (3 * PrA[1,"1"])
Ngen <- 1e8

setwd(paste(Path_base, "bamm/extant_only/", sep=""))
source(paste(Path_base, "Rscripts/ancThresh_write.R", sep=""), chdir = TRUE)

ACE <- ancThresh(tree, neoteny[tree$tip.label], ngen=Ngen, sequence=colnames(PrA), control=list(sample=Ngen/1e3, propthresh=0.15*max(nodeHeights(tree)), propliab=0.75*max(nodeHeights(tree)), burnin=Ngen*0.1, piecol=Cols, tipcol="input", pr.anc=PrA, plot=F))
save(ACE, file=paste(Path_base, "output/neotAnc/neotAnc_extant.RData", sep=""))


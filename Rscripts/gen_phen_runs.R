Path_base <- "~/Documents/Salamanders/"
setwd(paste(Path_base, "Rscripts", sep=""))
source('placeFossil.R', chdir = TRUE)
library(phytools)
library(TreePar)

# Read in the data
setwd(paste(Path_base, "datafiles", sep=""))
extant <- read.csv("salamanders.csv", stringsAsFactors=FALSE)
extant <- rbind(extant, c("Karsenia_koreana", "NA", 41.8, 0.15*41.8))	# From description in Nature, 2005 Min et al.

fossil <- read.csv("fossilSalamanders.csv", stringsAsFactors=FALSE)
tree <- read.tree("amphibians3309.tre")

fossilSp <- apply(fossil,1,function(x) paste(x[1], x[2], sep="_", collapse=""))
fossil2 <- cbind(fossil$svl, fossil$head_width, fossil$age_min, fossil$age_max)
rownames(fossil2) <- fossilSp

extantSVL <- cbind(extant[,1], log(as.numeric(extant$SVL)))
rownames(extantSVL) <- extant[,1]
colnames(extantSVL) <- c("species", "svl")
fossilSVL <- cbind(fossilSp, log(fossil$svl))
rownames(fossilSVL) <- fossilSp
colnames(fossilSVL) <- c("species", "svl")
morphData <- rbind(extantSVL, fossilSVL)

### Post-generation addition of files
for (count in 1:100)	{
	setwd(paste(Path_base, "bamm_morph/control_file/", sep=""))
	z <- read.table("control.txt", stringsAsFactors=F, row.names=1)
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	phy <- read.tree("fossilTree.tre")
	phy2 <- drop.tip(phy, setdiff(phy$tip.label, rownames(morphData)))
	tDat <- morphData[phy2$tip.label,]
	# remove NAs
	tDat2 <- tDat[!is.na(tDat[,2]),]
	phy3 <- drop.tip(phy2, setdiff(phy2$tip.label, rownames(tDat2)))
	
	dir.create(paste(Path_base, "bamm_morph/output-", count, sep=""))
	setwd(paste(Path_base, "bamm_morph/output-", count, sep=""))
	write.tree(phy3, "fossilTree.tre")
	write.table(tDat2, "traits.txt", quote=FALSE, row.names=F, col.names=F, sep="\t")
	priors <- setBAMMpriors(phy3, traits="traits.txt", outfile=NULL)
	
	z[names(priors)[2:4],2] <- priors[2:4]
	write.table(z, file="control.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
}


setwd(paste(Path_base, "bamm_morph/extant-only", sep=""))
phy <- read.tree("modernSal.tre")
phy2 <- drop.tip(phy, setdiff(phy$tip.label, rownames(morphData)))
tDat <- morphData[phy2$tip.label,]
# remove NAs
tDat2 <- tDat[!is.na(tDat[,2]),]
phy3 <- drop.tip(phy2, setdiff(phy2$tip.label, rownames(tDat2)))

write.tree(phy3, "extantTree.tre")
write.table(tDat2, "traits.txt", quote=FALSE, row.names=F, col.names=F, sep="\t")
	
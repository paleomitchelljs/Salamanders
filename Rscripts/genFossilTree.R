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

# Drop tips with no data
SalNode <- findMRCA(tree, tips=c("Plethodon_cinereus", "Cryptobranchus_alleganiensis", "Andrias_japonicus"))
GetSals <- na.omit(tree$tip.label[getDescendants(tree, SalNode)])

Drops <- setdiff(tree$tip.label, GetSals)
onlySal <- drop.tip(tree, Drops)
totalSal <- 701
write.tree(onlySal, "modernSal.tre")

salTimes <- getx(onlySal)
Frac <- Ntip(onlySal) / totalSal
tree_rates <- function(par)	{
	x <- LikConstant(par[1], par[2], Frac, salTimes)
	return(x)
}

# Need speciation, extinction and preservation rates to determine branch lengths for fossil tips
rates <- optim(c(0.01, 0.005), tree_rates)
psi_rate <- sum(fossil$numOcc) / sum(onlySal$edge.length)


### Plot showing probability of different fossil branch lengths
ds <- seq(from=0,to=100,by=0.1)
brlen_prob <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=psi_rate))
brlen_prob_nopres <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=0))
brlen_prob_tenpres <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=psi_rate*10))

setwd(paste(Path_base, "figures", sep=""))
pdf("brProb.pdf", height=5, width=5)
par(mar=c(4,4,1,1), mgp=c(2.5,0.3,0), tck=-0.005, las=1)
plot(ds, brlen_prob, type='l', ylim=c(0,0.08), xlim=c(0,100), axes=F, xlab=expression(paste("missing time (" %~~% "fossil branch length)", sep="")), ylab="probability")
lines(ds, brlen_prob_nopres, col='red', lty=3)
lines(ds, brlen_prob_tenpres, col='blue', lty=2)
axis(1, at=c(-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08))
legend("topright", legend=c(expression(paste(psi, " = 0")), expression(paste(psi, " = est.")), expression(paste(psi, " = 10X est."))), col=c('red', 'black', 'blue'), lty=c(3,1,2), bty='n')
dev.off()

# Keep one frog and caecillian b/c some fossil salamanders are stem to sal clade
Keeper <- c(grep("Rana", Drops)[1], grep("Ichthyophis", Drops)[1])
Frog <- Drops[Keeper]
Drops <- Drops[-Keeper]

# Make a tree that is salamanders + one frog
salTree <- drop.tip(tree, Drops)
saveTree <- salTree

# Generate trees
Stop <- FALSE
sim <- 1
while (Stop == FALSE)	{
	setwd(paste(Path_base, "/bamm/", sep=""))
	salTree <- saveTree
	
# Add each fossil iteratively, as some are sister to other fossils
for (count in 1:length(unique(fossilSp)))	{
	Row <- which(fossilSp == unique(fossilSp)[count])[1]
	Children <- fossil[Row,c("lchild", "rchild")]
	Fossil <- unique(fossilSp)[count]

	# Draw age from bounds
	#Age <- runif(1, min=fossil[Row, "age_min"], max=fossil[Row, "age_max"])
	
	# Set age to min age
	Age <- fossil[Row, "age_min"]
	
	if (Children$rchild == "")	{
		Children <- Children$lchild
		lTips <- salTree$tip.label[grep(Children[1], salTree$tip.label)]
		Tips <- lTips
	}
	else	{
		lTips <- salTree$tip.label[grep(Children[1], salTree$tip.label)]
		rTips <- salTree$tip.label[grep(Children[2], salTree$tip.label)]
		Tips <- unique(c(lTips, rTips))
	}
	
	newTree <- try(placeFossil(salTree, Age, Name=Fossil, taxa=Tips, grain=0.01, rates=c(rates$par, psi_rate), maxAge=310))

	if (class(newTree) == "try-error")	{
		print(Fossil)
	}
	if (class(newTree) == "phylo")	{
		plot(newTree, show.tip.label=F)
		salTree <- newTree
		cat(Fossil, ":", count, ": Tips =", Ntip(salTree), ": Nodes =", Nnode(salTree), "\n")
	}
}

	# Remove the cursed frog!
	salTree <- drop.tip(salTree, c(Frog, "Gerobatrachus_hottoni"))

	if (Ntip(salTree) >= 550)	{
		dir.create(paste("output-", sim, sep=""))
		setwd(paste("output-", sim, sep=""))
		write.tree(salTree, "fossilTree.tre")
	
		setwd(paste(Path_base, "bamm/control_file/", sep=""))
		z <- read.table("control.txt", stringsAsFactors=F, row.names=1)
		setwd(paste(Path_base, "bamm/output-", sim, sep=""))
		phy <- salTree
		priors <- setBAMMpriors(phy, outfile=NULL)
		z[names(priors)[2:4],2] <- priors[2:4]
		totalTime <- max(nodeHeights(phy))
		z["observationTime",2] <- totalTime
		z["numberOccurrences",2] <- sum(fossil$numOcc)
		write.table(z, file="control.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
		z["numberOccurrences",2] <- sum(fossil$numOcc)*10
		write.table(z, file="control_hi.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
		samp <- read.table(paste(Path_base, "bamm/control_file/sampling.txt", sep=""), sep="\t", row.names=1)
		samp <- samp[phy$tip.label,]
		samp <- cbind(rownames(samp), samp)
		colnames(samp) <- c("0.68367346938","","")
		write.table(samp, "sampling.txt", row.names=T, col.names=T, quote=FALSE)		
		sim <- sim + 1
	}
	if (sim >= 1)	{
		Stop <- TRUE
		break;
	}å
}



### Post-generation addition of files
#for (count in 1:100)	{
#	setwd(paste(Path_base, "bamm/control_file/", sep=""))
#	z <- read.table("control.txt", stringsAsFactors=F, row.names=1)
#	setwd(paste(Path_base, "bamm/output-", count, sep=""))
#	phy <- read.tree("fossilTree.tre")
#	priors <- setBAMMpriors(phy, outfile=NULL)
#	z[names(priors)[2:4],2] <- priors[2:4]
#	totalTime <- max(nodeHeights(phy))
#	z["observationTime",2] <- totalTime
#	z["numberOccurrences",2] <- sum(fossil$numOcc)
#	write.table(z, file="control.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
#	z["numberOccurrences",2] <- sum(fossil$numOcc)*10
#	write.table(z, file="control_hi.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
#	samp <- read.table("~/Documents/Salamanders/bamm/control_file/sampling.txt", sep="\t", row.names=1)
#	samp <- samp[phy$tip.label,]
#	samp <- cbind(rownames(samp), samp)
#	colnames(samp) <- c("0.68367346938","","")
#	write.table(samp, "sampling.txt", row.names=F, col.names=T, quote=FALSE)
#}


### Make diagnostic plots to look at missing/fossil data
setwd(paste(Path_base, "check_tree_files", sep=""))
missing <- setdiff(onlySal$tip.label, unique(extant$Species))
tipCols <- rep(1, Ntip(onlySal))
tipCols[onlySal$tip.label%in%missing] <- 2
plot(onlySal, tip.color=c('black', 'red')[tipCols], cex=0.25)
pdf("missingTree.pdf", height=10, width=10)
plot(onlySal, tip.color=c('black', 'red')[tipCols], cex=0.15)
dev.off()

pdf("fossilTree.pdf", height=10, width=10)
missing <- setdiff(salTree$tip.label, unique(fossilSp))
tipCols <- rep(1, Ntip(salTree))
tipCols[salTree$tip.label%in%missing] <- 2
plot(salTree, tip.color=c('red', 'black')[tipCols], cex=0.15)
dev.off()


### Make diagnostic plots to look at particular groups
newtNode <- findMRCA(salTree, tips=c("Palaeoproteus_klatti", "Batrachosauroides_dissimulans", "Opithisthotriton_kayi"))
newtTips <- salTree$tip.label[getDescendants(salTree, newtNode)]
newtTree <- drop.tip(salTree, setdiff(salTree$tip.label, newtTips))

newtNode <- findMRCA(onlySal, tips=c("Triturus_karelinii", "Salamandrina_perspicillata"))
newtTips <- onlySal$tip.label[getDescendants(onlySal, newtNode)]
newtTree <- drop.tip(onlySal, setdiff(onlySal$tip.label, newtTips))

newtNode <- findMRCA(onlySal, tips=c("Pseudobranchus_axanthus", "Siren_lacertina"))
newtTips <- onlySal$tip.label[getDescendants(onlySal, newtNode)]
newtTree <- drop.tip(onlySal, setdiff(onlySal$tip.label, newtTips))

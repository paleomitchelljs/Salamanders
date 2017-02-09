Path_base <- "~/Documents/Salamanders/"
setwd(paste(Path_base, "Rscripts", sep=""))
source('placeFossil.R', chdir = TRUE)
library(phytools)
library(TreePar)
library(BAMMtools)

# Read in the data
setwd(paste(Path_base, "datafiles", sep=""))
extant <- read.csv("salamanders.csv", stringsAsFactors=FALSE)
extant <- rbind(extant, c("Karsenia_koreana", "NA", 41.8, 0.15*41.8))	# From description in Nature, 2005 Min et al.
tree <- read.tree("amphibians3309.tre")

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
	Pars <- exp(par)
	x <- LikConstant(Pars[1], Pars[2], Frac, salTimes)
	return(x)
}

# Need speciation, extinction and preservation rates to determine branch lengths for fossil tips
rates <- optim(c(0.01, 0.005), tree_rates)
psi_rate <- 129 / sum(onlySal$edge.length) #sum(fossil$numOcc)

Save <- c("Plethodon_cinereus", "Ambystoma_tigrinum", "Amphiuma_means", "Andrias_japonicus", "Rhyacotriton_cascadae", "Necturus_maculosus", "Bolitoglossa_stuarti", "Dicamptodon_copei", "Siren_lacertina", "Hynobius_katoi", "Neurergus_kaiseri")
newtree <- drop.tip(onlySal, setdiff(onlySal$tip.label, Save))
newtree$tip.label <- c("Cryptobranchidae", "Hynobiidae", "Sirenidae", "Dicamptodontidae", "Ambystomatidae", "Salamandridae", "Proteidae", "Rhyacotritonidae", "Amphiumidae", "Plethodontinae", "Bolitoglossinae")
setwd(paste(Path_base, "figures", sep=""))
pdf("skeletontree.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,5))
plot(newtree, cex=2, font=1)
dev.off()
plethGens <- c("Plethodon", "Bolito", "Eurycea", "Aneides", "Pseudotriton", "Thorius", "Batrachoseps")
plethTips <- unlist(sapply(plethGens, function(x) onlySal$tip.label[grep(x, onlySal$tip.label)]))
#plethNode <- findMRCA(onlySal, plethTips)	#633
plethTipN <- na.omit(onlySal$tip.label[getDescendants(onlySal, node=633)])
plethTipN <- setdiff(plethTipN, "Gyrinophilus_porphyriticus")

hynGens <- c("Hynobius_kimurae", "Salamandrella_keyserlingii", "Onychodactylus_japonicus")
hynNode <- findMRCA(onlySal, hynGens)	# 171
hynTipN <- na.omit(onlySal$tip.label[getDescendants(onlySal, hynNode)])
hynTipN <- setdiff(hynTipN, "Salamandrella_keyserlingii")

### Plot showing probability of different fossil branch lengths
ds <- seq(from=0,to=100,by=0.1)
brlen_prob <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=psi_rate))
brlen_prob_nopres <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=0))
brlen_prob_tenpres <- sapply(ds, function(x) probDel(x, lam=rates$par[1], mu=rates$par[2], psi=psi_rate*10))

setwd(paste(Path_base, "figures", sep=""))
pdf("brProb.pdf", height=5, width=5)
par(mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, bty="n", las=1, cex.lab=1.5)
Alpha <- 1
exCol <- rgb(217/255, 95/255, 2/255, Alpha)
hiCol <- rgb(27/255, 158/255, 119/255, Alpha)
loCol <- rgb(117/255, 112/255, 179/255, Alpha)
plot(ds, brlen_prob, type='l', ylim=c(0,0.08), xlim=c(0,100), axes=F, xlab=expression(paste("missing time (" %~~% "fossil branch length)", sep="")), ylab="probability", col=loCol, lwd=3)
lines(ds, brlen_prob_nopres, col=exCol, lty=3, lwd=3)
lines(ds, brlen_prob_tenpres, col=hiCol, lty=2, lwd=3)
axis(1, at=c(-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
axis(2, at=c(-10, 0, 0.02, 0.04, 0.06, 0.08))
legend("right", legend=c(expression(paste(psi, " = 0")), expression(paste(psi, " = est.")), expression(paste(psi, " = 10X est."))), col=c(exCol, loCol, hiCol), lty=c(3,1,2), bty='n', lwd=3, cex=2)
dev.off()

# Keep one frog and caecillian b/c some fossil salamanders are stem to sal clade
Keeper <- c(grep("Rana", Drops)[1], grep("Ichthyophis", Drops)[1])
Frog <- Drops[Keeper]
Drops <- Drops[-Keeper]

# Make a tree that is salamanders + one frog
salTree <- drop.tip(tree, Drops)
salTree_uni <- salTree
saveTree <- salTree

# 7 13 20 60 68
Topologies <- c("fossilSalamanders.csv", "fossilSalamanders_BS+E.csv", "fossilSalamanders_BSE.csv", "fossilSalamanders_S+BE.csv", "fossilSalamanders_SBE.csv", "fossilSalamanders_SD+BHC.csv", "fossilSalamanders_SHC+BD.csv")
for (tcount in 1:length(Topologies))	{
	setwd(paste(Path_base, "datafiles/", sep=""))
	
	fossil <- read.csv(Topologies[tcount], stringsAsFactors=FALSE)
	fossilSp <- apply(fossil,1,function(x) paste(x[1], x[2], sep="_", collapse=""))
	fossil2 <- cbind(fossil$svl, fossil$head_width, fossil$age_min, fossil$age_max)
	rownames(fossil2) <- fossilSp

	# Generate trees
	Stop <- FALSE
	sim <- 1
	while (Stop == FALSE)	{
		simN <- sim + ( (tcount - 1) * 10 )
		Opt <- "bamm/"
	#	Opt <- "full/"
		setwd(paste(Path_base, Opt, sep=""))
	
		salTree <- saveTree
		salTree_uni <- saveTree
	
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
	
			newTree <- try(placeFossil(salTree, Age, Name=Fossil, taxa=Tips, grain=1e-3, rates=c(rates$par, psi_rate), maxAge=325))
			newTree_uni <- try(placeFossil(salTree_uni, Age, Name=Fossil, taxa=Tips, grain=1e-3, rates=NULL, maxAge=325))

			if (class(newTree) == "try-error")	{
				print(Fossil)
			}
			if (class(newTree) == "phylo")	{
		
				if (min(salTree$edge.length) > 0)	{
					salTree <- newTree
					cat(Fossil, ":", count, ": Tips =", Ntip(salTree), ": Nodes =", Nnode(salTree), "\n")
				}
				else	{
					minCount <- 0
					Stop2 <- FALSE
					while (Stop2 != TRUE)	{
						newTree <- try(placeFossil(salTree, Age, Name=Fossil, taxa=Tips, grain=1e-3, rates=c(rates$par, psi_rate), maxAge=325))
						if (class(newTree) == "phylo")	{
							if (min(newTree$edge.length) > 0)	{
								salTree <- newTree
								cat(Fossil, ":", count, ": Tips =", Ntip(salTree), ": Nodes =", Nnode(salTree), " || min used!\n")
								Stop2 <- TRUE
							}
						}
						minCount <- minCount + 1
						if (minCount > 10)	{
							Stop2 <- TRUE
							cat("Failed to slot in ", Fossil, "\n")
						}
					}
				}
			}
			if (class(newTree_uni) == "try-error")	{
				cat(Fossil, " Uniform error\n")
			}
			if (class(newTree_uni) == "phylo")	{
				salTree_uni <- newTree_uni
			}
		}

		if (Opt == "bamm/")	{
			# Remove the cursed frog! Remove crazy extants
			fullTree <- drop.tip(salTree, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", "Ichthyophis_bombayensis", "Rana_alticola"))
			fullTree_uni <- drop.tip(salTree_uni, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", "Ichthyophis_bombayensis", "Rana_alticola"))			
			salTree <- drop.tip(salTree, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", hynTipN, plethTipN, "Ichthyophis_bombayensis", "Rana_alticola"))
			salTree_uni <- drop.tip(salTree_uni, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", hynTipN, plethTipN, "Ichthyophis_bombayensis", "Rana_alticola"))
		}
		if (Opt == "full/")	{
			# Remove the cursed frog! Retain crazy extants
			salTree <- drop.tip(salTree, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", "Ichthyophis_bombayensis", "Rana_alticola"))
			salTree_uni <- drop.tip(salTree_uni, c(Frog, "Gerobatrachus_hottoni", "Celtedens_ibericus", "Ichthyophis_bombayensis", "Rana_alticola"))
		}
	
		if (Ntip(salTree) >= 170)	{
			dir.create(paste("output-", simN, sep=""))
			setwd(paste("output-", simN, sep=""))
			write.tree(salTree, "fossilTree.tre")
			write.tree(salTree_uni, "fossilTree_uni.tre")
			write.tree(fullTree, "fullTree.tre")
			write.tree(fullTree_uni, "fullTree_uni.tre")
			setwd(paste(Path_base, Opt, "control_file/", sep=""))
			z <- read.table("control.txt", stringsAsFactors=F, row.names=1)
			setwd(paste(Path_base, Opt, "output-", simN, sep=""))
			phy <- salTree
			priors <- setBAMMpriors(phy, outfile=NULL)
			z[names(priors)[2:4],2] <- priors[2:4]
			totalTime <- max(nodeHeights(phy))
			z["observationTime",2] <- totalTime
			z["numberOccurrences",2] <- sum(fossil$numOcc) + 26	# current count of unique pleth fossils
			write.table(z, file="control.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
			z["numberOccurrences",2] <- sum(fossil$numOcc)*10
			z["updateRateEventPosition",2] <- 0.8
			z["updateRateEventNumber",2] <- 1.2
			z["lambdaIsTimeVariablePrior",2] <- 1
			z["numberOfGenerations",2] <- 50000000
			z["mcmcWriteFreq",2] <- 25000
			z["eventDataWriteFreq",2] <- 25000
			z["minCladeSizeForShift",2] <- 1
			
			write.table(z, file="control_hi.txt", quote=FALSE, row.names=T, col.names=F, sep=" ")
			samp <- read.table(paste(Path_base, Opt, "control_file/sampling_nopleth.txt", sep=""), sep="\t", row.names=1)
			if (Opt == "full/")	{
				samp <- read.table(paste(Path_base, Opt, "control_file/sampling.txt", sep=""), sep="\t", row.names=1)
			}
			samp_uni <- samp
			Droptips <- setdiff(phy$tip.label, rownames(samp))
			if (length(Droptips) >= 1)	{
				cat("Dropping from tree: ", Droptips, "\n")
				phy <- drop.tip(phy, Droptips)
				write.tree(phy, "fossilTree.tre")
			}
			samp <- samp[phy$tip.label,]
			samp <- cbind(rownames(samp), samp)
			colnames(samp) <- c("0.1697575","","")
			if (Opt == "full/")	{
				colnames(samp) <- c("0.684593","","")
			}		
			write.table(samp, "sampling.txt", row.names=F, col.names=T, quote=FALSE)		
			Droptips <- setdiff(salTree_uni$tip.label, rownames(samp_uni))
			if (length(Droptips) >= 1)	{
				cat("Dropping from tree: ", Droptips, "\n")
				salTree_uni <- drop.tip(salTree_uni, Droptips)
				write.tree(salTree_uni, "fossilTree_uni.tre")
			}
			samp <- samp_uni[salTree_uni$tip.label,]
			samp <- cbind(rownames(samp), samp)
			colnames(samp) <- c("0.1697575","","")
			if (Opt == "full/")	{
				colnames(samp) <- c("0.684593","","")
			}	
			write.table(samp, "sampling_uni.txt", row.names=F, col.names=T, quote=FALSE)				
			sim <- sim + 1
		}
		if (sim >= 11)	{
			Stop <- TRUE
		}
	}
}

setwd(paste(Path_base, "datafiles", sep=""))
tree <- read.tree("fossilTree.tre")
Save <- c("Plethodon_cinereus", "Ambystoma_tigrinum", "Amphiuma_means", "Andrias_japonicus", "Rhyacotriton_cascadae", "Necturus_maculosus", "Bolitoglossa_stuarti", "Dicamptodon_copei", "Siren_lacertina", "Hynobius_katoi", "Neurergus_kaiseri", "Karaurus_sharovi", "Celtedens_ibericus", "Piceoerpeton_wilwoodense", "Peratosauroides_problematica")
newtree <- drop.tip(tree, setdiff(tree$tip.label, Save))
newtree$edge.length[28] <- newtree$edge.length[28] + 122.5 #push Albanerpetontidae to Pliocene
newtree$tip.label <- c("Cryptobranchidae", "Hynobiidae", "Batrachosauroididae", "Sirenidae", "Dicamptodontidae", "Ambystomatidae", "Salamandridae", "Proteidae", "Rhyacotritonidae", "Amphiumidae", "Plethodontinae", "Bolitoglossinae", "Scapherpetonidae", "Karauridae", "Albanerpetontidae")
setwd(paste(Path_base, "figures", sep=""))
pdf("skeletontree_fossils.pdf", height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,5))
plot(newtree, cex=2, font=1)
dev.off()


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


setwd(paste(Path_base, "datafiles", sep=""))
salTree <- read.tree("fossilTree.tre")
exTree <- read.tree("modernSal.tre")
missing <- setdiff(salTree$tip.label, exTree$tip.label)
setwd(paste(Path_base, "check_tree_files", sep=""))
pdf("fossilTree.pdf", height=10, width=10)
tipCols <- rep(1, Ntip(salTree))
tipCols[salTree$tip.label%in%missing] <- 2
#plot(salTree, tip.color=c('red', 'black')[tipCols], cex=0.15)
plot(salTree, show.tip.label=F)
tiplabels(pch=16, col=c(rgb(0,0,0,0), rgb(1,0,0,1))[tipCols], cex=1.1)
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

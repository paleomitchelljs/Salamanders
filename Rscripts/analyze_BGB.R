Path_base <- "~/Documents/Salamanders/"
Name <- "softFos"
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
setwd(paste(Path_base, "output", sep=""))
load(paste("sal_dec*_", Name, ".Rdata", sep=""))
load(paste("sal_dec*jv_", Name, ".Rdata", sep=""))

setwd(paste(Path_base, "datafiles", sep=""))
tree <- read.tree("fossilTree.tre")
Cols <- c('blue', "cyan", "chartreuse3", "goldenrod3", "red", "orchid1")

#################################################
pdf(paste(Path_base, "figures/sal_dec*_", Name, ".pdf", sep=""), height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,8))
plot_BGB(y_dec, analysis_titletxt="", plotwhat="pie", label.offset=0.25, tipcex=0.01, statecex=0.75, splitcex=0.75, titlecex=0, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=read.tree(z_dec$trfn), plotlegend=F)
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=F)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
dev.off()

pdf(paste(Path_base, "figures/sal_dec*jv_", Name, ".pdf",sep=""), height=10, width=10)
par(mar=c(0,0,0,0), oma=c(0,0,0,8))
plot_BGB(y_decjv, analysis_titletxt="", plotwhat="pie", label.offset=0.25, tipcex=0.01, statecex=0.5, splitcex=0.5, titlecex=0, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=read.tree(z_decjv$trfn), plotlegend=F)
addFamily(c("Ambystoma", "Ambystoma", "Ambystoma"), "Ambystomatidae", centering=F)
#addFamily(c("Dicamptodon", "Dicamptodon", "Dicamptodon"), "Dicamptodontidae")
addFamily(c("Siren", "Siren", "Pseudobranchus"), "Sirenidae", centering=T)
addFamily(c("Amphiuma", "Amphiuma", "Amphiuma"), "Amphiumidae", centering=T)
addFamily(c("Bolito", "Plethodon", "Desmognathus"), "Plethodontidae")
addFamily(c("Notophthalmus", "Neurergus", "Taricha"), "Salamandridae")
addFamily(c("Andrias", "Andrias", "Cryptobranchus"), "Cryptobranchidae")
addFamily(c("Hynobius", "Onychodactylus", "Liua"), "Hynobiidae")
#addFamily(c("Necturus", "Necturus", "Necturus"), "Proteidae", centering=T, Top=F)
if (Name == "hardFos")	{
	legend("bottomleft", legend=c("East NA", "West NA", "SA+CA", "Eur", "CAsia", "EAsia"), text.col=Cols, bty="n", cex=1.25)
}
dev.off()

#################################################
results_object <- y_decjv
BioGeoBEARS_run_object <- z_decjv
tipranges <- getranges_from_LagrangePHYLIP(BioGeoBEARS_run_object$geogfn)
areas <- getareas_from_tipranges_object(tipranges)
numareas <- length(areas)
max_range_size <- results_object$inputs$max_range_size
numstates <- numstates_from_numareas(numareas=length(areas), maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)
states_list_areaLetters <- areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=results_object$inputs$include_null_range)

Heights <- nodeHeights(tree)
occMat <- y_decjv$ML_marginal_prob_each_state_at_branch_top_AT_node
res <- 0.01
Inc <- seq(from=0, to=max(Heights), by=res)
Nregions <- 6

gMat <- matrix(0, nrow=length(Inc), ncol=numstates)
for (count in 1:length(Inc))	{
	# Edges beginning before Inc
	edgesBefore <- which(Heights[,1] <= Inc[count])

	# Edges ending after Inc
	edgesAfter <- which(Heights[,2] > Inc[count])

	# Edges that exist at Inc
	edgesAt <- intersect(edgesBefore, edgesAfter)

	Probs <- apply(occMat[edgesAt,], 2, function(x) x)
	gMat[count,] <- apply(Probs, 2, sum)
}

setwd(paste(Path_base, "fossilSalamanders/pbdb", sep=""))
pbdb_NA <- read.csv("raw_curve_dataNA.csv", stringsAsFactors=F)
pbdb_Eur <- read.csv("raw_curve_dataEur.csv", stringsAsFactors=F)
pbdb_Asia <- read.csv("raw_curve_dataAsia.csv", stringsAsFactors=F)

pdf(paste(Path_base, "figures/sal_dec*jv_LTT_", Name, ".pdf", sep=""), height=5, width=5)
par(las=1, bty="n", mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, cex.lab=1.5)
plot(1, 1, xlab="age (Ma)", ylab="log number of lineages", xlim=c(0, max(Heights)), ylim=c(-3, 5), type="n", axes=F)
Cols <- c('blue', "cyan", "chartreuse3", "goldenrod3", "red", "orchid1")
Times <- c(-100, 0, 63, 120, 199, 240, 265) + 0.8339
axis(1, at=Times, labels=sapply(Times, function(x) as.character(round(max(Heights)-x))))
axis(2, at=c(-100, -9, -7, -5, -3, -1, 1, 3, 5))
#axis(2, at=c(0.001, 0.2, 2, 20, 200), labels=c("n", "0.2", "2", "20", "200"))
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "A" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[1], lwd=2)
lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "W" %in% x)],1,sum)), col=Cols[2], lwd=2)
lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "S" %in% x)],1,sum)), col=Cols[3], lwd=2)
lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "Er" %in% x)],1,sum)), col=Cols[4], lwd=2)
lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "C" %in% x)],1,sum)), col=Cols[5], lwd=2)
lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "E" %in% x)],1,sum)), col=Cols[6], lwd=2)
if (Name == "hardFos")	{
	legend("topleft", legend=c("East NA", "West NA", "SA+CA", "Eur", "CAsia", "EAsia"), col=Cols, lty=1, lwd=1.5, bty="n", cex=0.75)
}
dev.off()


pdf(paste(Path_base, "figures/sal_dec*jv_LTT_", Name, "_focused_fos.pdf", sep=""), height=5, width=5)
par(las=1, bty="n", mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, cex.lab=1.5)
plot(1, 1, xlab="age (Ma)", ylab="log number of lineages", xlim=c(0, max(Heights)), ylim=c(-3, 5), type="n", axes=F)
Cols <- c('blue', "cyan", "chartreuse3", "goldenrod3", "red", "orchid1")
Times <- c(-100, 0, 63, 120, 199, 240, 265) + 0.8339
axis(1, at=Times, labels=sapply(Times, function(x) as.character(round(max(Heights)-x))))
axis(2, at=c(-100, -9, -7, -5, -3, -1, 1, 3, 5))
#axis(2, at=c(0.001, 0.2, 2, 20, 200), labels=c("n", "0.2", "2", "20", "200"))
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "A" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[1], lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "W" %in% x)],1,sum)), col='gray70', lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "S" %in% x)],1,sum)), col='gray70', lwd=2)
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "Er" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[4], lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "C" %in% x)],1,sum)), col='gray70', lwd=2)
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "E" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[6], lwd=2)
points(max(Heights)-pbdb_NA$Midpoint_Ma, log(pbdb_NA$X3T_rich), pch=16, col='blue', cex=1.5, lwd=1.5)
points(max(Heights)-pbdb_Eur$Midpoint_Ma, log(pbdb_Eur$X3T_rich), pch=16, col='goldenrod3', cex=1.5, lwd=1.5)
points(max(Heights)-pbdb_Asia$Midpoint_Ma, log(pbdb_Asia$X3T_rich), pch=16, col='orchid1', cex=1.5, lwd=1.5)
if (Name == "hardFos")	{
	legend("topleft", legend=c("East NA", "Eur", "EAsia"), col=Cols[c(1,4,6)], lty=1, lwd=1.5, bty="n", cex=0.75)
}
dev.off()


pdf(paste(Path_base, "figures/sal_dec*jv_LTT_", Name, "_focused_nofos.pdf", sep=""), height=5, width=5)
par(las=1, bty="n", mar=c(4,4,1,1), mgp=c(2.2,0.5,0), tck=-0.005, cex.lab=1.5)
plot(1, 1, xlab="age (Ma)", ylab="log number of lineages", xlim=c(0, max(Heights)), ylim=c(-3, 5), type="n", axes=F)
Cols <- c('blue', "cyan", "chartreuse3", "goldenrod3", "red", "orchid1")
Times <- c(-100, 0, 63, 120, 199, 240, 265) + 0.8339
axis(1, at=Times, labels=sapply(Times, function(x) as.character(round(max(Heights)-x))))
axis(2, at=c(-100, -9, -7, -5, -3, -1, 1, 3, 5))
#axis(2, at=c(0.001, 0.2, 2, 20, 200), labels=c("n", "0.2", "2", "20", "200"))
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "A" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[1], lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "W" %in% x)],1,sum)), col='gray70', lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "S" %in% x)],1,sum)), col='gray70', lwd=2)
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "Er" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[4], lwd=2)
#lines(Inc, log(apply(gMat[,sapply(states_list_areaLetters, function(x) "C" %in% x)],1,sum)), col='gray70', lwd=2)
Y <- log(apply(gMat[,sapply(states_list_areaLetters, function(x) "E" %in% x)],1,sum))
Y[Y < -3] <- NA
lines(Inc, Y, col=Cols[6], lwd=2)
if (Name == "hardFos")	{
	legend("topleft", legend=c("East NA", "Eur", "EAsia"), col=Cols[c(1,4,6)], lty=1, lwd=1.5, bty="n", cex=0.75)
}
dev.off()
#################################################
# Region of each node vs. proportion of all descendents still in that region
# Go through each node and find: what region has most extant tips? what is the prob the node existed in that region? what is the most likely region for node to have been in?


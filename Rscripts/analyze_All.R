### 
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
library(coda)

# Geographic range
Areas <- read.table(paste(Path_base, "datafiles/areas.txt", sep=""), row.names=1, stringsAsFactors=F)

# Read in extant tree and neotenic data
extantTree <- read.tree(paste(Path_base, "bamm/extant_only/no_pleth.tre", sep=""))
extantEdata <- getEventData(extantTree, paste(Path_base, "bamm/extant_only/np_event_data.txt", sep=""), burnin=0.1)
load(file=paste(Path_base, "output/neotAnc/neotAnc_extant.RData", sep=""))
exNeot <- ACE

# Read in fossil trees and neotenic data
corMat <- list()
specP <- c()
extP <- c()
for (count in 1:70)	{
	# Neoteny
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	tree <- read.tree("fossilTree.tre")
	load(file=paste(Path_base, "output/neotAnc/neotAnc", count, ".RData", sep=""))
	Heights <- nodeHeights(tree)
	rownames(Heights) <- tree$edge[,2]
	Liabs <- apply(ACE$liab[100:1001,], 2, mean)
	Thresh <- apply(ACE$par[100:1001,], 2, mean)

	# SpEx
	edata <- getEventData(tree, "base_event_data.txt", burnin=0.1, nsamples=500)
	mcmc <- read.csv("base_mcmc_out.txt", stringsAsFactors=F)

	Ext <- edata$meanTipMu
	names(Ext) <- edata$tip.label
	Spec <- edata$meanTipLambda
	names(Spec) <- edata$tip.label

	Int <- intersect(names(Ext), names(Liabs))
	if (length(Int) != length(Ext))	{
		cat("The number of neoteny values are not comparable to the number of tip rates!\n")
	}
	cat("Eff. Size N: ", effectiveSize(mcmc[floor(0.1*nrow(mcmc)):nrow(mcmc),"N_shifts"]), "\n")
	cat("Eff. Size LL: ", effectiveSize(mcmc[floor(0.1*nrow(mcmc)):nrow(mcmc),"logLik"]), "\n")	
	specP[count] <- traitDependentBAMM(edata, Liabs[edata$tip.label], reps=10000, method="s", rate="s")$p.value
	extP[count] <- traitDependentBAMM(edata, Liabs[edata$tip.label], reps=10000, method="s", rate="e")$p.value
	
	Inter <- intersect(names(Spec), rownames(Areas))
	corMat[[count]] <- cbind(Spec[Int], Ext[Int], Liabs[Int], rep(Thresh[3], length(Int)), rep(cor(Spec[Inter], log(Areas[Inter,1])), length(Int)), rep(cor(Ext[Inter], log(Areas[Inter,1])), length(Int)))
}


specR <- sapply(corMat, function(x) cor(x[,1], x[,3])^2)
extR <- sapply(corMat, function(x) cor(x[,2], x[,3])^2)

######################################################################################################################################
######################################################################################################################################
X <- sample(1:70, 1)
par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.2, 0.5, 0), tck=-0.005, las=1, cex.lab=1.5, bty="n", bg='white', fg='black', col.axis='black', col.lab='black')
Alpha <- 1
nCols <- c(rgb(254/255,232/255,200/255,Alpha), rgb(253/255,187/255,132/255,Alpha), rgb(227/255, 74/255, 51/255, Alpha))
names(nCols) <- c("0", "1", "2")

extantNeot <- load(file=paste(Path_base, "output/neotAnc/neotAnc_extant.RData", sep=""))
exLiabs <- apply(ACE$liab[10:101,], 2, mean)
exThresh <- apply(ACE$par[10:101,], 2, mean, na.rm=T)

nVec <- rep(1, length(exLiabs))
nVec[exLiabs>=exThresh[2]] <- 2
nVec[exLiabs>=exThresh[3]] <- 3

Lim <- pretty(c(0, max(extantEdata$meanTipLambda)))
plot(extantEdata$meanTipLambda, exLiabs[extantEdata$tip.label], pch=16, col=nCols[nVec], xlab="speciation rate", ylab="develomental lability", ylim=c(-15, 20), xlim=range(Lim), axes=F)
axis(1, at=c(-10, Lim))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))
legend("topright", legend=rev(c("Normal", "Fac. neoteny", "Obl. neoteny")), pch=16, col=rev(nCols), text.col=rev(nCols), bty="n", cex=0.75)

Lim <- pretty(c(0, max(extantEdata$meanTipMu)))
plot(extantEdata$meanTipMu, exLiabs[extantEdata$tip.label], pch=16, col=nCols[nVec], xlab="speciation rate", ylab="develomental lability", ylim=c(-15, 20), xlim=range(Lim), axes=F)
axis(1, at=c(-10, Lim))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))

Thresh <- c(0, corMat[[X]][1,4], Inf)
Liabs <- corMat[[X]][,3]
nVec <- rep(1, length(Liabs))
nVec[Liabs>=Thresh[1]] <- 2
nVec[Liabs>=Thresh[2]] <- 3

Lim <- pretty(c(0, max(corMat[[X]][,1])))
plot(corMat[[X]][,1], corMat[[X]][,3], pch=16, col=nCols[nVec], xlab="speciation rate", ylab="develomental lability", ylim=c(-15, 20), xlim=range(Lim), axes=F)
axis(1, at=c(-10, Lim))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))

Lim <- pretty(c(0, max(corMat[[X]][,2])))
plot(corMat[[X]][,2], corMat[[X]][,3], pch=16, col=nCols[nVec], xlab="extinction rate", ylab="develomental lability", ylim=c(-15, 20), xlim=range(Lim), axes=F)
axis(1, at=c(-10, Lim))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
# SVL 
#phendata <- getEventData(paste(Path_base, "bamm_morph/output-1/fossilTree.tre", sep=""), paste(Path_base, "bamm_morph/output-1/svl_event_data.txt", sep=""), burnin=0.5, type="trait")
#phenmcmc <- read.csv(paste(Path_base, "bamm_morph/output-1/svl_mcmc_out.txt", sep=""), stringsAsFactors=F)
#raw_svl <- read.table(paste(Path_base, "bamm_morph/output-1/traits.txt", sep=""), row.names=1, sep="\t")

# Geography
#setwd(paste(Path_base, "output", sep=""))
#Name <- "softFos"
#load(paste("sal_dec*jv_", Name, ".Rdata", sep=""))
#setwd(paste(Path_base, "datafiles", sep=""))
#results_object <- y_decjv
#tr_pruningwise <- reorder(tree, "pruningwise")
#tips <- 1:length(tr_pruningwise$tip.label)
#nodes <- (length(tr_pruningwise$tip.label)+1):(length(tr_pruningwise$tip.label)+tr_pruningwise$Nnode)
#relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
#relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
#bCols <- c('blue', "cyan", "chartreuse3", "goldenrod3", "red", "orchid1", rep('black', 57))

# Plots!
plot(tree, show.tip.label=F)
nodelabels(pie=relprobs_matrix_for_internal_states, piecol=bCols)

setwd(paste(Path_base, "figures", sep=""))
pdf("ratesVneot.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.2, 0.5, 0), tck=-0.005, las=1, cex.lab=1.5, bty="n")
plot(Spec, Liabs[names(Spec)], pch=21, bg=nCols[nVec], xlab="speciation rate", ylab="neoteny liability", ylim=c(-15, 20), xlim=c(0, 0.15), axes=F)
axis(1, at=c(-10, 0, 0.05, 0.1, 0.15))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))
segments(0.04, 0, 0.3, 0, col=nCols[2], lwd=1.5)
segments(0.04, Thresh[3], 0.3, Thresh[3], col=nCols[3], lwd=1.5)

plot(Ext, Liabs[names(Ext)], pch=21, bg=nCols[nVec], xlab="net diversification", ylab="", ylim=c(-15, 20), xlim=c(-0.05, 0.1), axes=F)
axis(1, at=c(-10, -0.05, 0, 0.05, 0.1))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15, 20))
legend("topright", legend=rev(c("Normal", "Fac. neo.", "Obl. neo.")), pch=21, pt.bg=rev(nCols), bty="n", cex=1.5)
segments(-0.02, 0, 0.2, 0, col=nCols[2], lwd=1.5)
segments(-0.02, Thresh[3], 0.2, Thresh[3], col=nCols[3], lwd=1.5)
dev.off()

setwd(paste(Path_base, "figures", sep=""))
pdf("neotVsvl.pdf", height=5, width=5)
par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.2, 0.5, 0), tck=-0.005, las=1, cex.lab=1.5, bty="n")
plot(raw_svl[,1],Liabs[rownames(raw_svl)], xlab="log SVL", ylab="neoteny liability", ylim=c(-15, 15), xlim=c(0.5, 6.5), axes=F, pch=21, bg=nCols[nVec[rownames(raw_svl)]])
axis(1, at=c(-10, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5))
axis(2, at=c(-50, -15, -10, -5, 0, 5, 10, 15))
segments(0.5, 0, 6.5, 0, col=nCols[2], lwd=1.5)
segments(0.5, Thresh[3], 6.5, Thresh[3], col=nCols[3], lwd=1.5)
dev.off()

Spec[grep("Gyrinophilus", names(Spec))], Liabs[grep("Gyrinophilus", names(Liabs))]

plot(Beta, Liabs[names(Beta)])

plot(raw_svl[,1], Liabs[rownames(raw_svl)])

plot(log10(Areas[Inter,1]), Spec[Inter])

plot(log10(Areas[Inter,1]), Liabs[Inter])


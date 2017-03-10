Path_base <- "~/Documents/Salamanders/grant_figure"
setwd(Path_base)

library(phytools)
library(BAMMtools)

addFamily <- function(Taxa, Family, tree, Cex=2, centering=T, Top=T, ADJ=-0.1)	{
	Names <- tree$tip.label[unique(unlist(sapply(Taxa, function(x) grep(x, tree$tip.label))))]
	if (length(Names) > 2)	{
		allNames <- na.omit(tree$tip.label[getDescendants(tree, findMRCA(tree, tips=Names))])
	}
	else if (length(Names <= 2))		{
		allNames <- Names
	}
	if (centering == T)		{
		Center <- allNames[round(mean(1:length(allNames)))]
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=ADJ, cex=Cex)
	}
	else if (centering == F){
		if (Top == T)	{
			Center <- rev(allNames)[1]
		}
		else	{
			Center <- allNames[1]
		}
		tiplabels(Family, tip=which(tree$tip.label == Center), bg='white', frame='none', xpd=NA, adj=ADJ, cex=Cex)
	}
}

NSamp <- 900
ex_tree <- ladderize(read.tree("extant.tre"))
ex_edata <- getEventData(ex_tree, "extant_event_data.txt", burnin=0.1, nsamples=NSamp)

f_tree <- ladderize(read.tree("fossils.tre"))
f_edata <- getEventData(f_tree, "fossils_event_data.txt", burnin=0.1, nsamples=NSamp)

SPEX <- "e"
CEX <- 0.75 	#Scale for the family labels
Scale <- plot(f_edata, breaksmethod="jenks", spex=SPEX)

pdf("salamanders_ext.pdf", height=4, width=8)
#dev.new(height=4, width=8)
par(mfrow=c(1, 2), mar=c(0, 0, 0, 2), oma=c(0, 4, 2, 0))
plot.bammdata(ex_edata, colorbreaks=Scale$colorbreaks, spex=SPEX)
addFamily(c("Ambystoma_jeffersonianum", "Ambystoma_maculatum", "Ambystoma_laterale"), "Amb.", as.phylo(ex_edata), centering=T, Cex=CEX)
addFamily("Pseudobranchus_axanthus", "Si.", as.phylo(ex_edata), centering=F, Cex=CEX)
addFamily("Amphiuma_means", "Amp.", as.phylo(ex_edata), centering=F, Cex=CEX)
addFamily("Cynops", "Sa.", as.phylo(ex_edata), Cex=CEX)
addFamily("Cryptobranchus_alleganiensis", "Cr.", as.phylo(ex_edata), Cex=CEX, centering=F)
Location <- c(-20, -10, 10, 65)
addBAMMlegend(Scale, side=2, location=Location, axisOffset=0.000002, cex.axis=0.75*CEX)
mtext("extant-only", side=3, line=0, cex=1.25 * CEX)
text(x=2.75 * Location[1], y=mean(Location[3:4]), "net div.", xpd=NA, srt=90)

plot(f_edata, colorbreaks=Scale$colorbreaks, spex=SPEX)
addFamily(c("Ambystoma_jeffersonianum", "Ambystoma_maculatum", "Ambystoma_laterale"), "Amb.", as.phylo(f_edata), centering=T, Cex=CEX)
addFamily(c("Pseudobranchus"), "Si.", as.phylo(f_edata), centering=T, Cex=CEX)
addFamily("Amphiuma", "Amp.", as.phylo(f_edata), centering=T, Cex=CEX)
addFamily("Cynops", "Sa.", as.phylo(f_edata), Cex=CEX)
addFamily(c("Andrias", "Cryptobranchus"), "Cr.", as.phylo(f_edata), Cex=CEX)
mtext("extant+extinct", side=3, line=0, cex=1.25 * CEX)
dev.off()
#### Net div calculations
Path_base <- "~/Documents/Salamanders/"

library(phytools)
library(BAMMtools)
library(coda)
source(paste(Path_base, 'Rscripts/gcr.R', sep=""), chdir = TRUE)

magND <- function(N, t_stem, t_crown, EpsS=0.9, EpsC=0.9)	{
	# Eqn 6
	R_stem <- (1 / t_stem) * log( N * (1 - EpsS) + EpsS)
	
	# Eqn 7
	SQR <- N * ((N * EpsC^2) - (8 * EpsC) + (2 * N * EpsC) + N)
	SQR <- sqrt(SQR)
	Component <- (0.5 * N * (1 - EpsC^2)) + (2 * EpsC) +
					(0.5 * (1 - EpsC) * SQR)
	R_crown <- (1 / t_crown) * (log(Component) - log(2))

	# Eqn 6, epsilon = 0.9
	R_stem_hi <- (1 / t_stem) * log( N * (1 - 0.9) + 0.9)
	
	# Eqn 7, epsilon = 0.9
	SQR <- N * ((N * 0.9^2) - (8 * 0.9) + (2 * N * 0.9) + N)
	SQR <- sqrt(SQR)
	Component <- (0.5 * N * (1 - 0.9^2)) + (2 * 0.9) +
					(0.5 * (1 - 0.9) * SQR)
	R_crown_hi <- (1 / t_crown) * (log(Component) - log(2))

	# Eqn 6, epsilon = 0
	R_stem_lo <- (1 / t_stem) * log( N * (1 - 0) + 0)
	
	# Eqn 7, epsilon = 0
	SQR <- N * ((N * 0^2) - (8 * 0) + (2 * N * 0) + N)
	SQR <- sqrt(SQR)
	Component <- (0.5 * N * (1 - 0^2)) + (2 * 0) +
					(0.5 * (1 - 0) * SQR)
	R_crown_lo <- (1 / t_crown) * (log(Component) - log(2))
	return(c(rstem=R_stem, rcrown=R_crown, 
			 rstem_hi=R_stem_hi, rcrown_hi=R_crown_hi, 
			 rstem_lo=R_stem_lo, rcrown_lo=R_crown_lo))
}

getAge <- function(tips, edata, Age=NULL)	{
	if (length(tips) > 1)	{
		maxT <- max(edata$end)
		Node <- findMRCA(as.phylo(edata), tips, type="node")
		Edge <- which(edata$edge[,2] == Node)
		crownAge <- maxT - edata$end[Edge]
		
		if (is.null(Age))	{
			stemAge <- maxT - edata$begin[Edge]
			gcr_stem <- gcr(edata, node=Node, drop.stem=F)
			gcr_ND <- mean(gcr_stem$lambda - gcr_stem$mu)
			EPS <- mean(gcr_stem$mu / gcr_stem$lambda)
		}
		if (!is.null(Age))	{
			stemAge <- maxT - Age
			
			# distance from each edge's start to the age
			EdgeVec <- abs(edata$begin - (maxT - Age))
			
			# the two edges that begin closest to the age
			ancEdge <- which(EdgeVec == min(EdgeVec))

			# nodes that descend from these edges
			descNodes <- edata$edge[ancEdge,2]
			if (length(descNodes) > 2)	{
				cat("Something went wrong with: ", Age, "\n")
			}
			# descendent tips
			descTips <- sapply(descNodes, function(x) edata$tip.label[getDesc(edata, x)$desc_set])
			
			# correct node
			Node <- descNodes[c(tips[1] %in% descTips[[1]], tips[1] %in% descTips[[2]])]
			gcr_stem <- gcr(edata, node=Node, drop.stem=F)
			gcr_ND <- mean(gcr_stem$lambda - gcr_stem$mu)
			EPS <- mean(gcr_stem$mu / gcr_stem$lambda)
		}
		gcr_crown <- gcr(edata, node=Node, drop.stem=T)
		gcr_NDc <- mean(gcr_crown$lambda - gcr_crown$mu)
		EPSc <- mean(gcr_crown$mu / gcr_crown$lambda)
		return(c(stem=stemAge, crown=crownAge, gcr_stem=gcr_ND, 
				gcr_crown=gcr_NDc, epsilon_stem=EPS, epsilon_crown=EPSc))
	}
	else	{
		return(c(stem=NA, crown=NA, gcr_stem=NA, gcr_crown=NA, 
				epsilon_stem=NA, epsilon_crown=NA))
	}
}


# Extant-only data
setwd(paste(Path_base, "bamm/extant_only", sep=""))

# Read in files
exTree <- ladderize(read.tree("modernSal.tre"), right=F)
exMCMC <- read.csv("ex_mcmc_out.txt", stringsAsFactors=F)
exEdata <- getEventData(exTree, "ex_event_data.txt", burnin=0.5)
exSampling <- read.table("extant_sampling.txt", 
				stringsAsFactors=F, header=F, sep="\t", skip=1)

# Generate species-clade-sampling matrix
exSpecies <- apply(exSampling, 1, function(x) unlist(strsplit(x, split=" ")))

# Number of species per family
exFams <- tapply(exSpecies[1,], exSpecies[2,], function(x) return(x))
exobsFamN <- sapply(exFams, length)
exFrac <- as.numeric(tapply(exSpecies[3,], exSpecies[2,], function(x) unique(x)))
exN <- exobsFamN / exFrac

# Extant-no pleth data
# Read in files
setwd(paste(Path_base, "bamm/extant_only", sep=""))
npTree <- ladderize(read.tree("no_pleth.tre"), right=F)
npMCMC <- read.csv("np_mcmc_out.txt", stringsAsFactors=F)
npEdata <- getEventData(npTree, "np_event_data.txt", burnin=0.5)
npSampling <- read.table("no_pleth_sampling.txt", 
				stringsAsFactors=F, header=F, sep="\t", skip=1)

# Generate species-clade matrix
npSpecies <- apply(npSampling, 1, function(x) unlist(strsplit(x, split=" ")))

# Number of species per family
npFams <- tapply(npSpecies[1,], npSpecies[2,], function(x) return(x))
npobsFamN <- sapply(npFams, length)
npFrac <- as.numeric(tapply(npSpecies[3,], npSpecies[2,], function(x) unique(x)))
npN <- npobsFamN / npFrac

# fossilBAMM runs
### Post-generation addition of files
fbEdata <- list()
fbEdata_all <- list()
llES <- c()
nsES <- c()

for (count in 1:70)	{
	setwd(paste(Path_base, "bamm/output-", count, sep=""))
	fbTree <- ladderize(read.tree("fossilTree.tre"), right=F)
	MCMC_dat <- read.csv("base_mcmc_out.txt")
	MCMC_dat <- MCMC_dat[floor(0.1*nrow(MCMC_dat)):nrow(MCMC_dat),]
	
	fbEdata[[count]] <- getEventData(fbTree, "base_event_data.txt", 
									burnin=0.1, nsamples=200)
	llES[count] <- effectiveSize(MCMC_dat$logLik)
	nsES[count] <- effectiveSize(MCMC_dat$N_shifts)
	
	fullTree <- ladderize(read.tree("fullTree.tre"), right=F)
	fbEdata_all[[count]] <- getEventData(fullTree, "full_event_data.txt", 
									burnin=0.1, nsamples=200)
}


### Calculate net diversification
# yes pleth, extant-only
exAges <- sapply(exFams, function(x) getAge(x, exEdata))
exData <- rbind(exAges, exN)
ex_ND_mag <- apply(exData, 2, function(x) 
			magND(N=x[7], t_stem=x[1], t_crown=x[2], EpsS=x[5], EpsC=x[6]))
exAll <- rbind(exData, ex_ND_mag)

# no pleth, extant-only
npAges <- sapply(npFams, function(x) getAge(x, npEdata))
npData <- rbind(npAges, npN)
np_ND_mag <- apply(npData, 2, function(x) 
			magND(N=x[7], t_stem=x[1], t_crown=x[2], EpsS=x[5], EpsC=x[6]))
npAll <- rbind(npData, np_ND_mag)

# yes pleth, fossil tree
fAll <- list()
for (count in 1:70)	{
	fEdata <- fbEdata_all[[count]]
	fAges <- sapply(1:length(exFams), function(x) getAge(exFams[[x]], fEdata, Age=exAges[1,x]))
	fData <- rbind(fAges, exN)
	f_ND_mag <- apply(fData, 2, function(x) 
			magND(N=x[7], t_stem=x[1], t_crown=x[2], EpsS=x[5], EpsC=x[6]))
	fAll[[count]] <- rbind(fData, f_ND_mag)
}
fStems <- sapply(fAll, function(x) x["gcr_stem",])
fCrown <- sapply(fAll, function(x) x["gcr_crown",])

# no pleth, fossil tree
fnpAll <- list()
for (count in 1:70)	{
	fnpEdata <- fbEdata[[count]]
	fnpAges <- sapply(1:length(npFams), function(x) getAge(npFams[[x]], fnpEdata, Age=npAges[1,x]))
	fnpData <- rbind(fnpAges, npN)
	fnp_ND_mag <- apply(fnpData, 2, function(x) 
			magND(N=x[7], t_stem=x[1], t_crown=x[2], EpsS=x[5], EpsC=x[6]))
	fnpAll[[count]] <- rbind(fnpData, fnp_ND_mag)
}
fnpStems <- sapply(fnpAll, function(x) x["gcr_stem",])
fnpCrown <- sapply(fnpAll, function(x) x["gcr_crown",])

# Compile data
reformatMat <- function(x, data_type="extant", pleth=0, treeN=1)	{
	x <- t(x)
	x <- data.frame(treeN=rep(treeN, nrow(x)), data=rep(data_type, nrow(x)), pleth=rep(pleth, nrow(x)), clade=rownames(x), x, row.names=NULL)
	return(x)
}
rexAll <- reformatMat(exAll, data_type="extant", pleth=1)
rnpAll <- reformatMat(npAll, data_type="extant", pleth=0)
rfAll <- sapply(1:length(fAll), function(x) reformatMat(fAll[[x]], data_type="fossil", pleth=1, treeN=x), simplify=F)
rfAll <- do.call(rbind, rfAll)
rfnpAll <- sapply(1:length(fnpAll), function(x) reformatMat(fnpAll[[x]], data_type="fossil", pleth=0, treeN=x), simplify=F)
rfnpAll <- do.call(rbind, rfnpAll)

allMat <- rbind(rexAll, rfAll)
npMat <- rbind(rnpAll, rfnpAll)
colnames(npMat) <- colnames(allMat)
totalMat <- rbind(allMat, npMat)

write.csv(totalMat, paste(Path_base, "datafiles/netdiv.csv", sep=""), quote=F, row.names=F)


#### Objective fit
exRates <- getCladeRates(exEdata)
exRates <- getCladeRates(exEdata)
 
### Make a PHYLIP style file
Path_base <- "~/Documents/Salamanders/"

setwd(paste(Path_base, "datafiles", sep=""))

geo1 <- read.csv("geography_morph.csv", stringsAsFactors=F)
geo2 <- read.csv("geography_tree.csv", stringsAsFactors=F)

fos <- read.csv("fossilSalamanders.csv", stringsAsFactors=F)[,c("genus", "species", "ENA", "WNA", "Neo", "Eur", "CASIA", "EASIA")]
fos_names <- apply(fos[,1:2], 1, paste, collapse="_")
fosMat <- cbind(fos_names, fos[,3:8])
colnames(fosMat) <- c("X", "ENA", "WNA", "NEO", "EUR", "CASIA", "EASIA")

combined <- rbind(geo1, geo2, fosMat)
combined[,1] <- gsub(" ", "_", combined[,1])

geoScore <- apply(combined[,2:ncol(combined)], 1, paste, collapse="")

Nregion <- 6
Ntaxa <- nrow(combined)

phylip <- cbind(combined[,1], geoScore)
colnames(phylip) <- c(Ntaxa, Nregion)

write.table(phylip, file="sal_geog.data", sep="\t", quote=F, row.names=F, col.names=T)

minorProb <- apply(combined[,2:7], 2, function(x) sapply(x, function(y) ifelse(y==0, 0.001, y)))
minorProb <- cbind(combined[,1], minorProb)
colnames(minorProb) <- c("X", "ENA", "WNA", "NEO", "EUR", "CASIA", "EASIA")

write.table(minorProb, file="sal_minor.txt", sep="\t", quote=F, row.names=F, col.names=T)
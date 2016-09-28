## Salamander diversity map
# Based on code from Pascal Title
Path_base <- "~/Documents/Salamanders/"
require(raster)
require(maptools)
require(rgeos)
require(rangeBuilder)
require(parallel)
require(geosphere)

# filename and location of IUCN salamanders shapefile
IUCNfile <- paste(Path_base, 'geography/CAUDATA/CAUDATA.shp', sep="")
extant <- read.csv(paste(Path_base, 'datafiles/salamanders.csv', sep=""), stringsAsFactors=F)
extantSp <- unique(extant$Species)
extantSp <- gsub("_", " ", extantSp)
fossil <- read.csv(paste(Path_base, 'datafiles/fossilSalamanders.csv', sep=""), stringsAsFactors=F)
fLat <- fossil$lat
fLong <- fossil$long

# load shapefile as SpatialPolygonsDataFrame
## We know it is unprojected long/lat, WGS84 datum
salamanders <- readShapeSpatial(IUCNfile, repair = TRUE, delete_null_obj = TRUE, force_ring = TRUE, proj4string = CRS('+proj=longlat +datum=WGS84'))

## Convert sp names to character
salamanders@data$binomial <- as.character(salamanders@data$binomial)
allsp <- unique(salamanders@data$binomial)

#nameVec <- intersect(extantSp, allsp)
#nameVec <- intersect(missing, allsp)
nameVec <- allsp
spList <- list()
for (i in 1:length(nameVec)) {
	ind <- which(salamanders@data$binomial == nameVec[i])
	tmp <- salamanders[ind,]	

	spList[[i]] <- tmp
	names(spList)[i] <- nameVec[i]
}

names(spList) <- gsub(' ', '_', names(spList))

# Get lat ranges for all sp:
latRange <- sapply(spList, function(x) x@bbox["y",])
latBins <- seq(from=-20, to=74, by=1)
Richness <- c()
latMid <- c()
for (count in 2:length(latBins))	{
	bin_min <- latBins[count-1]
	bin_max <- latBins[count]
	keep_min <- which(latRange[1,] <= bin_max & latRange[1,] >= bin_min)
	keep_max <- which(latRange[2,] <= bin_max & latRange[2,] >= bin_min)
	Richness[count-1] <- length(unique(c(keep_min, keep_max)))
	latMid[count-1] <- mean(c(bin_min, bin_max))
}

setwd(paste(Path_base, "figures", sep=""))
pdf("ldg.pdf", height=5, width=5)
par(mar=c(4,4,1,1), mgp=c(2.2,0.5,0), las=1, bty="n", cex.lab=1.5, tck=-0.005)
plot(latMid, Richness, type='b', col='black', pch=21, bg='gray70', lwd=1, cex=1.5, axes=F, xlim=c(-80, 80), ylim=c(0, 60), xlab="latitude", ylab="richness")
axis(1, at=c(-120, -80, -60, -40, -20, 0, 20, 40, 60, 80))
axis(2, at=c(-100, 0, 10, 20, 30, 40, 50, 60))
dev.off()


# -----------------------------------------------------------------------
# project these ranges to equal area projections
## For North America, we will use the North America Albers Equal Area
## http://spatialreference.org/ref/esri/102008/

EAproj <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"

for (i in 1:length(spList)) {
	spList[[i]] <- spTransform(spList[[i]], CRS(EAproj))
}

Areas <- sapply(spList, gArea)
Area_mat <- matrix(Areas, ncol=1, dimnames=list(names(Areas)))
write.table(Area_mat, file=paste(Path_base, "datafiles/areas.txt", sep=""), quote=F, col.names=F, sep="\t")
# -----------------------------------------------------------------------
listExtent <- getExtentOfList(spList)
template <- raster(xmn = listExtent$minLong, xmx = listExtent$maxLong, ymn = listExtent$minLat, ymx = listExtent$maxLat, res = c(50000, 50000), crs=proj4string(spList[[1]]))

# Plot salamander ranges & fossil occurrences
data(wrld_simpl)
wrld <- spTransform(wrld_simpl, CRS(proj4string(template)))

setwd(paste(Path_base, "geography/extant_maps", sep=""))
for (count in 1:length(spList))	{
	pdf(paste(names(spList)[count], ".pdf", sep=""), height=6, width=12)
	plot(wrld, col='gray80', lty=0, xlim=c(listExtent$minLong, listExtent$maxLong), ylim=c(listExtent$minLat, listExtent$maxLat))
	plot(spList[[count]], add=T, col=rgb(1,0,0,1), border='black')
	dev.off()
}

pdf("fullMap.pdf", height=6, width=12, encoding="MacRoman")
plot(wrld, col='gray80', lty=0, xlim=c(listExtent$minLong, listExtent$maxLong), ylim=c(listExtent$minLat, listExtent$maxLat))
silent <- sapply(spList, plot, add=T, col=rgb(1,0,0,0.1), border=rgb(0,0,0,0))

fosPts <- project(cbind(fLong,fLat), EAproj)
text(fosPts[,1], fosPts[,2], labels=expression("\u2020"), col='black', offset=0, cex=0.5)
dev.off()

pdf("fullMap_noFos.pdf", height=6, width=12, encoding="MacRoman")
plot(wrld, col='gray80', lty=0, xlim=c(listExtent$minLong, listExtent$maxLong), ylim=c(listExtent$minLat, listExtent$maxLat))
silent <- sapply(spList, plot, add=T, col=rgb(1,0,0,0.1), border=rgb(0,0,0,0))
dev.off()

# Pascal: this is where things stop functioned. richnessRaster() finds everything as NA
# rasterize
rasList <- list()
for (i in 1:length(spList)) {
	tmp <- rasterize(spList[[i]], template)
	values(tmp)[!is.na(values(tmp))] <- 1
	rasList[[i]] <- tmp
}
names(rasList) <- names(spList)

# We currently have a list of rasters, but as they are all of the same resolution, extent, etc, we can create a rasterStack
salStack <- stack(rasList)

# -----------------------------------------------------------------------
# We are now ready to generate a species richness raster
# Drop species with NA or recalculate at higher resolution. 
#salStack <- salStack[[setdiff(1:nlayers(salStack), c(1,2,3))]]

richness <- richnessRaster(ranges = salStack, nthreads = 4, resUnits="meters", resolution=50)

# Let's plot it:
ramp <- colorRampPalette(c('blue','yellow','red'))
plot(richness, col=ramp(100), legend = FALSE)
addRasterLegend(richness, location = c(-5e+06, -4.7e+06, -2e+06, 2e+06), ramp=c('blue','yellow','red'))
plot(wrld, add=TRUE)

# -----------------------------------------------------------------------






















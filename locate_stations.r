# File - locate_stations.r
# Version - 04.04.2013
# Author - Matthew Parkan
# Description - Compute inter-station distance matrix and create meta data file
# Data - ftp://ftp.ncdc.noaa.gov/pub/data/inventories/

#clear workspace
rm(list=ls())

#load librairies
library(rgdal)
library(sp)

# IMPORTANT!!! Define working directories
inputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\NCDC weather\\"
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\NCDC station metadata\\"

# IMPORTANT!!! Define inventory file name
# the input file should be placed in the input directory 
inputfile <- "ISH-HISTORY.CSV" #comma separated

# IMPORTANT!!! Define output file names
# the output file will be created in the work directory
outputfile1 <- "station_metadata"
outputfile2 <- "inter_station_distances"

#load data (comma separated)
file <- paste(inputpath,inputfile,sep="")
data <- read.csv(file, header=TRUE,colClasses=c(rep("character",7),rep("numeric",3),rep("character",2)))
data <- as.data.frame(data)

#Replace missing values by NA
data$WBAN[data$WBAN==99999] <- NA
data$LAT[data$LAT==-99999 | data$LAT==0] <- NA
data$LON[data$LON==-999999 | data$LON==0] <- NA
data$ELEV[data$ELEV==-99999] <- NA

# Dates
data$BEGIN <- as.Date(data$BEGIN, format="%Y%m%d")
data$END <- as.Date(data$END, format="%Y%m%d")

# Coordinates
data$LAT <- data$LAT/1000
data$LON <- data$LON/1000
data$ELEV <- data$ELEV/10

# Extract Greenland stations (GL) subset (only stations with known coordinates)
data.GL <- subset(data, data$CTRY=="GL" & !is.na(data$LON))
rownames(data.GL) <- NULL

# Add NSIDC Sea Ice Polar Stereographic North (EPSG:3411) coordinates in attributes
coords_ll = cbind(data.GL$LON, data.GL$LAT)
XY <- project(coords_ll, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs ")
XY <- round(XY,digits=0)
data.GL$XPS <- XY[,1]
data.GL$YPS <- XY[,2]
#coords_ps = cbind(data.GL$X, data.GL$Y)

# Extract Greenland stations (GL) subset with long time series
TSTART <- as.Date("19850101", format="%Y%m%d") #Date start
TSTOP <- as.Date("20121212", format="%Y%m%d") #Date stop
logind <- data.GL$BEGIN<=TSTART & data.GL$END>=TSTOP 
data.GL.long <- subset(data.GL, logind)
rownames(data.GL.long) <- NULL

# Export weather station metadata to .csv file
tablepath <- paste(outputpath,outputfile1,".csv",sep="") 
write.csv(data.GL.long, file=tablepath, row.names = FALSE)

# Export weather station metadata to .Rda file
save(data.GL.long,file=paste(outputpath,outputfile1,".Rda",sep=""))

# Export weather station metadata to .shp file
export.data <- data.GL.long[,c("USAF","WBAN","NAME","CTRY","FIPS","STATE","CALL","ELEV","BEGIN","END")]
spdf_ll = SpatialPointsDataFrame(subset(coords_ll,logind), proj4string=CRS("+init=epsg:4326"),export.data)
writeOGR(spdf_ll,paste(outputpath,"station_list_ll.shp"),"station_list_ll",driver="ESRI Shapefile",overwrite_layer=TRUE)

# Project coordinates to NSIDC Sea Ice Polar Stereographic North (EPSG:3411)
spdf_ps <- spTransform(spdf_ll, CRS("+init=epsg:3411"))
writeOGR(spdf_ps,paste(outputpath,"station_list_ps.shp"),"station_list_ps",driver="ESRI Shapefile",overwrite_layer=TRUE)

# Select station subset #settdiff, match
mystations <- c("043200","043390","043600","043900","042200","042020") #USAF number
mystations_logind <- data.GL.long$USAF %in% mystations
  
# Compute distances between all stations
dist_ll <- spDistsN1(spdf_ll, spdf_ll[1,], longlat=TRUE)
dist_ps <- spDistsN1(spdf_ps, spdf_ps[1,], longlat=FALSE)/1000 #as a comparison only, EPSG:3411 is conformal
dist <- cbind(dist_ll,dist_ps)

distmatrix_ll <- spDists(spdf_ll, spdf_ll, longlat=TRUE)
distmatrix_ps <- spDists(spdf_ps, spdf_ps, longlat=FALSE)/1000 #as a comparison only, EPSG:3411 is conformal

mydistmatrix_ll <- t(subset(distmatrix_ll, mystations_logind))
mydistmatrix_ps <- t(subset(distmatrix_ps, mystations_logind)) #as a comparison only, EPSG:3411 is conformal

orderind <- apply(mydistmatrix_ll,2,order)

sortednames <- matrix(rep(paste(spdf_ll$USAF," ",spdf_ll$NAME,sep="")),25,6)
sortedvalues <- round(mydistmatrix_ll, digits = 0)
sortednames[] <- sortednames[cbind(as.vector(orderind), as.vector(col(orderind)))]
sortedvalues[] <- sortedvalues[cbind(as.vector(orderind), as.vector(col(orderind)))]

#concatenate station names and distances
res <- matrix(paste(sortednames,", ~",sortedvalues, " km",sep=""),25,6) 
colnames(res) <- subset(spdf_ll$NAME, mystations_logind)

# Export distance matrix to .csv file
tablepath <- paste(outputpath,outputfile2,".csv",sep="") 
write.csv(res, file=tablepath, row.names = FALSE)

# Export distance matrix to .Rda file
save(res,file=paste(outputpath,outputfile2,".Rda",sep=""))


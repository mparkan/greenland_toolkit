# File - map_stations.r
# Version - 20.04.2013
# Author - Matthew Parkan
# Description - Plot a map of weather stations and buffers
# Data - ftp://ftp.ncdc.noaa.gov/pub/data/inventories/
#        http://www.naturalearthdata.com/downloads/

#clear workspace
rm(list=ls())

#load librairies
library(rgdal)
library(sp)

# IMPORTANT!!! Define path to Natural Earth data
coastlinepath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\Natural Earth Coastlines" 

# IMPORTANT!!! Define path to station metadata file (.Rda)
metadatapath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\NCDC station metadata\\station_metadata.Rda"

# IMPORTANT!!! Define path to folder for saving figures
figurepath <- "C:\\Users\\mat\\Google Drive\\Greenland\\figures\\"

# IMPORTANT!!! Define USAF number of desired stations (check station inventory)
myusaf <- c("043200","043390","043600","043900","042200","042020")

# IMPORTANT!!! Define buffer ranges (in km) around stations
circrange <- c(75,300)

# load Natural Earth coastline vector
ne_10m_coastline <- readOGR(dsn = coastlinepath, layer = "10m_coastline")

# project coastline to polar stereographic (EPSG:3411)
ne_10m_coastline_ps <- spTransform(ne_10m_coastline, CRS("+init=epsg:3411"))

# import .Rda station metadata
load(metadatapath)
numstat<-length(myusaf)
mystations <- subset(data.GL.long, data.GL.long$USAF %in% myusaf)
rownames(mystations)<-mystations$USAF
mystations<-mystations[myusaf,,drop=FALSE]
statcoords_ps <- cbind(mystations$XPS,mystations$YPS)
statcoords_ll <- cbind(mystations$LON,mystations$LAT) 
mystations_ps = SpatialPoints(statcoords_ps,proj4string=CRS("+init=epsg:3411"))

#define function to compute small circle coordinates (using spherical Earth approximation)
scircle <- function(clon,clat,r,id) {
  clon <- as.numeric(clon)
  clat <- as.numeric(clat)
  r <- as.numeric(r)
  id <- as.character(id)
  
  lonlatpoint=cbind(clon,clat)
  travelvector <- as.data.frame(cbind(direction = seq(0, 2*pi, by=2*pi/100), magnitude = r))
  Rearth <- 6372795
  Dd <- travelvector$magnitude / Rearth
  Cc <- travelvector$direction
  
  if (class(lonlatpoint) == "SpatialPoints") {
    lata <- coordinates(lonlatpoint)[1,2] * (pi/180)
    lona <- coordinates(lonlatpoint)[1,1] * (pi/180)
  }
  else {
    lata <- lonlatpoint[2] * (pi/180)
    lona <- lonlatpoint[1] * (pi/180)
  }
  latb <- asin(cos(Cc) * cos(lata) * sin(Dd) + sin(lata) * cos(Dd))
  dlon <- atan2(cos(Dd) - sin(lata) * sin(latb), sin(Cc) * sin(Dd) * cos(lata))
  lonb <- lona - dlon + pi/2
  
  lonb[lonb >  pi] <- lonb[lonb >  pi] - 2 * pi
  lonb[lonb < -pi] <- lonb[lonb < -pi] + 2 * pi
  latb <- latb * (180 / pi)
  lonb <- lonb * (180 / pi)
  c=cbind(longitude = lonb, latitude = latb)
  j = rbind(c, c[1, ])
  p= Polygons(list(Polygon(j)),id)
  return(p) 
}

#create circular buffers
radius=sort(rep(circrange,numstat))*1000
names=rep(myusaf,length(circrange))
bufnames=paste(myusaf,".",radius/1000,sep="")
bufinfo <- data.frame(statcoords_ll,radius,bufnames, stringsAsFactors=FALSE)
circles_ll <- SpatialPolygons(apply(bufinfo,1,function(x)scircle(x[1],x[2],x[3],x[4])),proj4string=CRS("+init=epsg:4326"))  
circles_ps <- spTransform(circles_ll, CRS("+init=epsg:3411"))


#plot map
##############################################################################

# Set plot boundaries (in data units)
xmin <- -1700000
xmax <- 1700000
ymin <- -4000000
ymax <- 100000

asratio = (ymax-ymin)/(xmax-xmin)


#print figure to .eps
postscript(paste(figurepath,"station_map.eps",sep=""), horizontal = FALSE, onefile = FALSE, paper = "special",
           width=4,height=4*asratio)

par(usr=c(xmin,xmax,ymin,ymax),
    mar=c(3,3,1,1), 
    oma=c(2,2,2,2),
    cex=0.5,
    col.axis="grey40",
    col.lab="grey40",
    mgp = c(2, 0.5, 0))

#plot basemap
plot(ne_10m_coastline_ps,xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     axes=FALSE,asp=1,xaxs="i",yaxs="i",lwd=0.2)
title(main = NULL, sub = NULL, xlab = "X [m]", ylab = "Y [m]", outer = FALSE)

#plot axis
box(lwd=0.5)
axis(side = 1, tck = -.01, lwd=0.5)
axis(side = 2, tck = -.01, lwd=0.5)

#plot selected stations
plot(mystations_ps,pch=21,cex=0.7,col="black",bg="tomato2",lwd=0.2,add=TRUE)

#plot circular buffers
plot(circles_ps,border="tomato2",lwd=0.3, lty=5,add=TRUE)

#plot station labels
text(statcoords_ps[,1], statcoords_ps[,2], 1:6, cex=1.5, adj=c(0.5,0.5), col="black")

dev.off()


#print figure to .pdf
pdf(paste(figurepath,"station_map.pdf",sep=""), onefile = FALSE, paper = "special",
    width=4,height=4*asratio)

par(usr=c(xmin,xmax,ymin,ymax),
    mar=c(3,3,1,1), 
    oma=c(2,2,2,2),
    cex=0.5,
    col.axis="grey40",
    col.lab="grey40",
    mgp = c(2, 0.5, 0))
#    cex.axis=0.6,


#plot basemap
plot(ne_10m_coastline_ps,xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     axes=FALSE,asp=1,xaxs="i",yaxs="i",lwd=0.2)
title(main = NULL, sub = NULL, xlab = "X [m]", ylab = "Y [m]", outer = FALSE)

#plot axis
box(lwd=0.5)
axis(side = 1, tck = -.01, lwd=0.5)
axis(side = 2, tck = -.01, lwd=0.5)


#plot selected stations
plot(mystations_ps,pch=21,cex=0.7,col="black",bg="tomato2",lwd=0.2,add=TRUE)
#text(statcoords_ps[,1], statcoords_ps[,2], 1:6, cex=1.5, pos=4, col="red") 

#plot circular buffers
plot(circles_ps,border="tomato2",lwd=0.3, lty=5,add=TRUE)

#plot station labels
text(statcoords_ps[,1], statcoords_ps[,2], 1:6, cex=1.5, adj=c(0.5,0.5), col="black")

dev.off()

#print figure to .png
png(paste(figurepath,"station_map.png",sep=""),width=2200,height=ceiling(2200*asratio),units = "px")

par(usr=c(xmin,xmax,ymin,ymax),
    mar=c(4,4,1,1), 
    oma=c(2,2,2,2),
    cex=2,
    col.axis="grey40",
    col.lab="grey40",
    mgp = c(2.5, 1, 0))

#plot basemap
plot(ne_10m_coastline_ps,xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     axes=FALSE,asp=1,xaxs="i",yaxs="i",lwd=0.2)
title(main = NULL, sub = NULL, xlab = "X [m]", ylab = "Y [m]", outer = FALSE)

#plot axis
box(lwd=0.5)
axis(side = 1, tck = -.01, lwd=0.5)
axis(side = 2, tck = -.01, lwd=0.5)

#plot selected stations
plot(mystations_ps,pch=21,cex=0.7,col="black",bg="tomato2",lwd=0.2,add=TRUE)

#plot circular buffers
plot(circles_ps,border="tomato2",lwd=0.3, lty=5,add=TRUE)

#plot station labels
text(statcoords_ps[,1], statcoords_ps[,2], 1:6, cex=1.5, adj=c(0.5,0.5), col="black")

dev.off()

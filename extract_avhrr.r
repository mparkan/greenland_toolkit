# File - extract_avhrr.r
# Version - 25.03.2013
# Author - Matthew Parkan
# Description - read OISST netcdf files and extract area statistics at specific locations
# Data - http://podaac.jpl.nasa.gov/dataset/NCDC-L4LRblend-GLOB-AVHRR_OI
#        http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html

#load librairies
library(rgdal)
library(ncdf)
library(sp)
library(raster) 
library(geosphere)
library(chron)
library(R.utils)

#clear workspace
rm(list=ls())

# IMPORTANT!!! Define working directories
#inputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\SST_AVHRR_OI"
inputpath <- "C:\\Users\\mat\\Desktop\\code\\DATA\\SST_AVHRR_OI\\data"
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\SST_AVHRR_OI\\"

# IMPORTANT!!! Define path to station metadata file (.Rda)
metadatapath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\NCDC station metadata\\station_metadata.Rda"

# IMPORTANT!!! Define path to quadrangle area file (.Rda)
quadareaspath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\SST_AVHRR_OI\\quadrangle_areas.Rda"

# IMPORTANT!!! Define global output file name (without extension)
outputfile <- "RS_STATISTICS"

# IMPORTANT!!! Define USAF number of desired stations (check station inventory)
myusaf <- c("043200","043390","043600","043900","042200","042020")

# IMPORTANT!!! Define ranges (in km) around stations to compute area statistics
circrange <- c(75,300)

#IMPORTANT!!! Specify number of lag days for each variable-range couple
SIC_lag_r75<-2
SST_lag_r75<-2
SIC_lag_r300<-2
SST_lag_r300<-2
lags<-data.frame(SIC_lag_r75,SST_lag_r75,SIC_lag_r300,SST_lag_r300)  

# IMPORTANT!!! Define USAF number of desired stations (check station inventory)
evalgeom=TRUE
computearea=FALSE  #WARNING may take several minutes to compute
extractvalues=FALSE  #WARNING may take several hours to extract


if (extractvalues == TRUE) {
  
  #list bz2 and netcdf files
  bzlist <- list.files(inputpath,recursive=TRUE,pattern = "\\.bz2$")
  bzlist <- gsub("/", "\\\\", bzlist)
  nclist <- gsub("[.]bz2$", "", bzlist)
  #nclist <- list.files(inputpath,recursive=TRUE,pattern = "\\.nc$")
  #nclist <- gsub("/", "\\\\", nclist)
  
  # import .Rda station metadata
  load(metadatapath)
  numstat<-length(myusaf)
  mystations <- subset(data.GL.long, data.GL.long$USAF %in% myusaf)
  rownames(mystations)<-mystations$USAF
  mystations<-mystations[myusaf,,drop=FALSE]
  statcoords_ps <- cbind(mystations$XPS,mystations$YPS)
  statcoords_ll <- cbind(mystations$LON,mystations$LAT) 
  
  #define function to rotate matrix 90° CCW 
  rot90 <- function(mat){apply(t(mat),2,rev)}
  
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
  
  
  if(evalgeom==TRUE){
    ncid <- open.ncdf(con=paste(inputpath,nclist[1],sep="\\"), write=FALSE)
    # create coordinates grids
    lat <- rev(get.var.ncdf(ncid, "lat")) #latitude
    lon <- get.var.ncdf(ncid, "lon") #longitude
    coords_ll <- as.matrix(expand.grid(lon,lat))
    coords_ps <-  project(coords_ll, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs ")
    grid_ps = SpatialPoints(coords_ps,proj4string=CRS("+init=epsg:3411"))
    grid_ll = SpatialPoints(coords_ll,proj4string=CRS("+init=epsg:4326"))
    height <- length(lat)
    width <- length(lon)
    
    mask <- rot90(get.var.ncdf(ncid, "mask")) #sea/land field composite mask
    mask_vec <- as.vector(t(mask))
    close.ncdf(ncid)
    
    #find all pixel (indices) that are inside the small circle
    radius=sort(rep(circrange,numstat))*1000
    names=rep(myusaf,length(circrange))
    bufnames=paste(myusaf,".",radius/1000,sep="")
    bufinfo <- data.frame(statcoords_ll,radius,bufnames, stringsAsFactors=FALSE)
    circles_ll <- SpatialPolygons(apply(bufinfo,1,function(x)scircle(x[1],x[2],x[3],x[4])),proj4string=CRS("+init=epsg:4326"))  
    circles_ps <- spTransform(circles_ll, CRS("+init=epsg:3411"))
  }
  
  #determine area of spherical quadrangles
  #WARNING may take several minutes to compute
  if(computearea==TRUE){
    tm1 <- system.time(
{
  res=0.25
  longitudes<-coords_ll[,1]
  latitudes<-coords_ll[,2]
  
  left<-longitudes-res/2
  right<-longitudes+res/2
  top<-latitudes+res/2
  bottom<-latitudes-res/2
  
  areaquad=cbind(left,left,right,right,top,bottom,bottom,top)
  area_vec<-apply(areaquad,1,function(x) areaPolygon(cbind(c(x[1],x[2],x[3],x[4]),c(x[5],x[6],x[7],x[8]))))
  sea_area_vec<-area_vec
  sea_area_vec[mask_vec==2]<-NA
  areas<-data.frame(area_vec,sea_area_vec)
  colnames(areas)<-c('ALL_AREAS','SEA_AREAS')
  
  #export to .Rda file
  save(areas,file=paste(inputpath,"quadrangle_areas.Rda",sep="\\"))
})} else{
  load(quadareaspath)
  sea_area_vec<-areas$SEA_AREAS
}
  
  #loop through netcdf file list and read variable values
  numfiles<-length(nclist)
  numcombi<-length(radius)
  ndfrows<-numfiles*numcombi
  start<-seq(1,ndfrows,by=numcombi)
  end<-seq(numcombi,ndfrows,by=numcombi)
  
  #initialize data frame
  final <- data.frame(AREA=numeric(ndfrows),
                      W_SIC=numeric(ndfrows),
                      W_SST=numeric(ndfrows),
                      W_ERR=numeric(ndfrows),
                      STATION=character(ndfrows),
                      RANGE=numeric(ndfrows),
                      DATE=as.Date(rep(0,ndfrows), origin = "1900-01-01"),
                      CHECK=character(ndfrows),
                      stringsAsFactors = FALSE)
  
  pb<-txtProgressBar(min = 0, max = numfiles, style = 3) #progress bar
  tm2 <- system.time({
    for(k in 1:3){ #numfiles
      
      # uncompress netcdf file
      bz2name=paste(inputpath,bzlist[k],sep="\\")
      ncdfname=gsub("[.]bz2$", "", bz2name)
      bunzip2(bz2name, ncdfname, overwrite=TRUE, remove=FALSE)
      
      # import OISST netcdf file
      ncid <- open.ncdf(con=paste(inputpath,nclist[k],sep="\\"), write=FALSE)
      # print(ncid) # print some basic information
      # vars = names(ncid$var)  #get a list of the data in the netcdf file
      
      # read data
      sea_ice_fraction <- rot90(get.var.ncdf(ncid, "sea_ice_fraction")) #sea ice area fraction
      #sort(unique(c(sea_ice_fraction)))
      sea_ice_fraction <- round(sea_ice_fraction,digits=2) 
      sea_ice_fraction[sea_ice_fraction == -1.28] <- 0
      sea_ice_fraction_vec <- as.vector(t(sea_ice_fraction))
      sea_ice_fraction_vec[mask_vec==2] <- NA
      
      analysed_sst <- rot90(get.var.ncdf(ncid, "analysed_sst")) #analysed sea surface temperature
      #sort(unique(c(analysed_sst)))
      analysed_sst <- round(analysed_sst,digits=2) 
      analysed_sst[analysed_sst==-54.53] <- NA
      analysed_sst <- analysed_sst-273.14999 #convert to celsius 
      analysed_sst_vec <- as.vector(t(analysed_sst))
      
      analysis_error <- rot90(get.var.ncdf(ncid, "analysis_error")) #estimated error standard deviation of analysed_sst
      #sort(unique(c(analysis_error)))
      analysis_error<- round(analysis_error,digits=2) 
      analysis_error[analysis_error==-327.68] <- NA
      analysis_error_vec <- as.vector(t(analysis_error))
      
      time <- get.var.ncdf(ncid, "time") #time
      tunits <- att.get.ncdf(ncid,"time","units") # time units
      
      #close ncdf file
      close.ncdf(ncid)
      
      # determine date
      tustr <- strsplit(tunits$value, " ") # parse tunits$value
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth=as.integer(unlist(tdstr)[2])
      tday=as.integer(unlist(tdstr)[3])
      tyear=as.integer(unlist(tdstr)[1])
      
      DATE <- as.Date(chron(time/86400,origin=c(tmonth, tday, tyear)), format="%d/%m/%Y")
      #DOY <- as.numeric(format(DATE, format = "%j")) 
      #YEAR = as.numeric(format(DATE, format = "%Y"))
      #MONTH = as.numeric(format(DATE, format = "%m"))
      #DAY = as.numeric(format(DATE, format = "%d"))
      
      # organise into spatial data frame
      rsdata<-data.frame(coords_ll,coords_ps,sea_area_vec,
                         sea_area_vec*sea_ice_fraction_vec,
                         sea_area_vec*analysed_sst_vec,
                         sea_area_vec*analysis_error_vec)
      colnames(rsdata) <- c('LONG','LAT','XPS','YPS','AREA','W_SIC','W_SST','W_ERR')
      rownames(rsdata)<-NULL
      spdf = SpatialPointsDataFrame(grid_ps,rsdata)
      # overlay polygons on grid (find points located in buffers)
      # values<-over(circles_ps,spdf,fn=function(x) sum(x,na.rm=TRUE),returnList = TRUE)
      # aggregate point values by polygons
      areastat<-as.data.frame(aggregate(spdf[,c('AREA','W_SIC','W_SST','W_ERR')], by=circles_ps, FUN = function(x) sum(x,na.rm=TRUE)))
      areastat[,c('W_SIC','W_SST','W_ERR')]<-areastat[,c('W_SIC','W_SST','W_ERR')]/areastat$AREA 
      areastat$STATION<-names
      areastat$RANGE<-radius
      areastat$DATE<-DATE
      areastat$CHECK<-rownames(areastat)
      final[start[k]:end[k],1:ncol(final)]<-areastat
      setTxtProgressBar(pb, k)
    }
  })
  close(pb)
  
  #round values
  final[,sapply(final,is.numeric)] <- round(final[,sapply(final,is.numeric)],digits=4)
  
  RS_STATISTICS<-final
  
  #export to .csv file
  tablepath <- paste(outputpath,outputfile,".csv",sep="") 
  write.csv(RS_STATISTICS, file=tablepath,row.names = FALSE)
  
  #export to .Rda file
  save(RS_STATISTICS,file=paste(outputpath,outputfile,".Rda",sep=""))
  
} else{
  load(paste(outputpath,outputfile,".Rda",sep=""))
}


#function to lag/lead, see http://www.r-bloggers.com/generating-a-laglead-variables/
shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  #stopifnot(is.numeric(x))
  
  if (length(shift_by)>1){
    return(sapply(shift_by,shift, x=x))}
    
    out<-NULL
    abs_shift_by=abs(shift_by)
  if (shift_by > 0 ){
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  }
  else if (shift_by < 0 ){
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  }
  else {out<-x}
  
  return(out)
}

#convert range column from m to km
RS_STATISTICS$RANGE <- RS_STATISTICS$RANGE/1000

for(k in 1:length(myusaf)){
  
  #subset by station
  RS_STATISTICS_SUB<-subset(RS_STATISTICS[,c("DATE","W_SIC","W_SST","RANGE","STATION")],RS_STATISTICS$STATION==myusaf[k])
  #split SIC and SST by range
  RS_split<-split(RS_STATISTICS_SUB[,c("W_SIC","W_SST")], RS_STATISTICS_SUB$RANGE, drop = FALSE)
  #split dates by range
  RS_dates<-split(RS_STATISTICS_SUB[,"DATE"], RS_STATISTICS_SUB$RANGE, drop = FALSE)
  #merge data frame columns
  merged<-do.call("cbind", RS_split)
  
  #lag
  lagn<-lapply(lags,function(x) rev(-x:0))
  lagged<-mapply(shift,merged,lagn,SIMPLIFY=FALSE)
  merged_lagged<-as.data.frame(do.call("cbind",lagged))
  
  #add columns names
  cnames<-unlist(mapply(rep, colnames(merged), lags+1))
  colnames(merged_lagged)<-paste(cnames,".M",abs(unlist(lagn)),sep="")
  
  #add date column
  merged_lagged$DATE <- RS_dates[[1]]
  
  RS_STATS<-merged_lagged
  
  #export to .csv file
  tablepath <- paste(outputpath,myusaf[k],"_",outputfile,".csv",sep="") 
  write.csv(RS_STATS, file=tablepath,row.names = FALSE)
  
  #export to .Rda file
  save(RS_STATS,file=paste(outputpath,myusaf[k],"_",outputfile,".Rda",sep=""))
  
}

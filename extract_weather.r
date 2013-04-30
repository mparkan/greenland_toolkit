# File - extract_weather.r
# Version - 31.03.2013
# Author - Matthew Parkan
# Description - Read NCDC weather data, create and export weather features to .csv and .Rda files 
# Data - http://www.ncdc.noaa.gov/cdo-web/

#clear workspace
rm(list=ls())

#IMPORTANT!!! Define working directories
inputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\NCDC weather\\"
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\NCDC weather\\"

# IMPORTANT!!! Define USAF number of desired stations (check station inventory)
#myusaf <- c("042030","042050","042280","042850","043300","043610")
# OR process all stations in input folder (comment the following two lines if you want to specify a subset of stations)
myusaf <- basename(list.dirs(inputpath,recursive=FALSE))
myusaf <- substr(myusaf, 1, 6) 

#IMPORTANT!!! Specify the minimum number of observations threshold 
# (features will not be created if the number of observations is below this threshold)
min_nobs <- 5000

#IMPORTANT!!! Specify the date range
daterange <- seq(as.Date("1978-01-01"),as.Date("2013-03-01"),by=1)

#IMPORTANT!!! Specify number of lag days for each variable
TEMP_MEAN_LAG <- 0
DEWPT_MEAN_LAG <- 2
SLP_MEAN_LAG <- 2
WIND_SPD_MEAN_LAG <- 2
WIND_DIR_RF_LAG <- 2
CEIL_HGT_RF_LAG <- 2

#histogram function
histcount <- function(input, bins){
  counts=hist(input,breaks=bins,plot=FALSE)[2];
  total=length(input);
  res <- lapply(counts, function(x) x/total)
  return(res) 
}

#function to lag/lead, see http://www.r-bloggers.com/generating-a-laglead-variables/
shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  #stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  out
}

#function to lag variables
lagvars<-function(variable,vardates,lagval,daterange){
  ORI <- with(variable, x[match(daterange,unique(vardates))])
  #rownames(ORI)<-NULL #########################################
  RES <- data.frame(daterange,shift(ORI,rev(-lagval:0)))
  varname <- deparse(substitute(variable))
  colnames(RES) <- c('DATE',paste(varname,'_M',as.character(0:lagval),sep=""))
  return(RES)
}

#function to find and unlist vector columns
#Cf. http://stackoverflow.com/questions/15930880/unlist-all-list-elements-in-a-dataframe
flat <- function(data) {
  temp1 <- sapply(data, is.list)
  temp2 <- do.call(
    cbind, lapply(data[temp1], function(x) 
      data.frame(do.call(rbind, x), check.names=FALSE)))
  cbind(data[!temp1], temp2)
}

# list csv files in the input folder
#dirlist <-list.dirs(inputpath, full.names = TRUE, recursive = TRUE)
#dirlist <- gsub("/", "\\\\", dirlist)
txtlist <- list.files(inputpath,recursive=TRUE,pattern = "dat\\.txt$")
station_id <- substr(txtlist, 1, 6) 
validind <- station_id %in% myusaf
txtlist <- subset(txtlist, validind==TRUE)
txtlist <- gsub("/", "\\\\", txtlist)
station_id <- subset(station_id, validind==TRUE)

numfiles <- length(txtlist)
nerrors=0
errors<-character()
pb<-txtProgressBar(min = 0, max = numfiles, style = 3) #progress bar
for(k in 1:numfiles){ 
  
  inputfile <- txtlist[k] # the input file should be placed in the input directory
  outputfile <- paste(station_id[k],"_WEATHER",sep="") # the output file will be created in the output directory
  
  #Import data
  file <- paste(inputpath,inputfile,sep="")
  column_cl <- c(rep("character",7),rep("numeric",2),"character",rep("numeric",4),rep("character",2),rep("numeric",2),"character",rep("numeric",8),"NULL")
  data <- read.csv(file, header=FALSE,skip=2,colClasses=column_cl)
  data <- as.data.frame(data)
  
  #Define column headers
  colnames(data) <- c( 'USAF', 'NCDC', 'YEARMODA','HRMN', 'OBS_I', 'OBS_TYPE','OBS_QCP',
                       'WIND_DIR','WIND_DIR_Q','WIND_OBS_I','WIND_SPD','WIND_SPD_Q',
                       'CEIL_HGT','CEIL_HGT_Q', 'CEIL_HGTQ_I', 'CAVOK',
                       'VISBY','VISBY_Q','VISBY_VAR','VISBY_VAR_Q',
                       'TEMP','TEMP_Q','DEWPT','DEWPT_Q','SLP','SLP_Q','RHX')
  
  #Subset date range of interest
  DATE <- as.Date(data$YEARMODA, format="%Y%m%d") #Date
  data <- subset(data, DATE>=daterange[1] & DATE<=tail(daterange,1))
  
  #Replace missing values by NA
  data$WIND_DIR[data$WIND_DIR==999] <- NA
  data$WIND_DIR[data$WIND_DIR_Q==2 | data$WIND_DIR_Q==3 | data$WIND_DIR_Q==6 | data$WIND_DIR_Q==7] <- NA
  data$WIND_SPD[data$WIND_SPD==9999 | data$WIND_SPD==999.9] <- NA
  data$WIND_SPD[data$WIND_SPD_Q==2 | data$WIND_SPD_Q==3 | data$WIND_SPD_Q==6 | data$WIND_SPD_Q==7] <- NA
  data$CEIL_HGT[data$CEIL_HGT==99999] <- NA
  data$CEIL_HGT[data$CEIL_HGT_Q==2 | data$CEIL_HGT_Q==3 | data$CEIL_HGT_Q==6 | data$CEIL_HGT_Q==7] <- NA
  data$VISBY[data$VISBY==999999] <- NA
  data$VISBY[data$VISBY_Q==2 | data$VISBY_Q==3 | data$VISBY_Q==6 | data$VISBY_Q==7] <- NA
  data$VISBY[data$VISBY==9] <- NA
  data$VISBY_VAR[data$VISBY_VAR==9] <- NA
  data$VISBY_VAR[data$VISBY_VAR_Q==2 | data$VISBY_VAR_Q==3 | data$VISBY_VAR_Q==6 | data$VISBY_VAR_Q==7] <- NA
  data$TEMP[data$TEMP==9999 | data$TEMP==999.9] <- NA
  data$TEMP[data$TEMP_Q==2 | data$TEMP_Q==3 | data$TEMP_Q==6 | data$TEMP_Q==7] <- NA
  data$DEWPT[data$DEWPT==9999 | data$DEWPT==999.9] <- NA
  data$DEWPT[data$DEWPT_Q==2 | data$DEWPT_Q==3 | data$DEWPT_Q==6 | data$DEWPT_Q==7] <- NA
  data$SLP[data$SLP==99999 | data$SLP==9999.9] <- NA
  data$SLP[data$SLP_Q==2 | data$SLP_Q==3 | data$SLP_Q==6 | data$SLP_Q==7] <- NA
  
  #Subset columns of interest
  data <- data[,c('USAF','YEARMODA','HRMN','WIND_DIR','WIND_SPD','CEIL_HGT','TEMP','SLP')]
  
  #delete rows with NA
  data <- data[!apply(data,1,function(y)any(is.na(y))),]
  rownames(data) <- NULL
  
  #Date and Time Code (UTC)
  YMDHM <- paste (data$YEARMODA,data$HRMN,sep =",")
  TIMESTAMP_UTC <- strptime(YMDHM, format="%Y%m%d,%H%M", tz="UTC")
  DATE <- as.Date(data$YEARMODA, format="%Y%m%d")
  DOY <- as.numeric(format(DATE, format = "%j")) 
  YEAR = as.numeric(format(DATE, format = "%Y"))
  MONTH = as.numeric(format(DATE, format = "%m"))
  DAY = as.numeric(format(DATE, format = "%d"))
  
  if (nrow(data) >= min_nobs) {
    
  #################################################################
  #TEMP: AIR-TEMPERATURE-OBSERVATION air temperature
  #################################################################
  TEMP_MEAN <- aggregate(data$TEMP,by=list(DATE),function(x) mean(x,na.rm=TRUE))
  TEMP_MEAN <- lagvars(TEMP_MEAN,unique(DATE),TEMP_MEAN_LAG,daterange)
  
  #################################################################
  #DEWPT: AIR-TEMPERATURE-OBSERVATION-DEWPOINT temperature
  #################################################################
  #DEWPT_MEAN <- aggregate(data$DEWPT,by=list(DATE),function(x) mean(x,na.rm=TRUE))
  #DEWPT_MEAN <- lagvars(DEWPT_MEAN,unique(DATE),DEWPT_MEAN_LAG,daterange)
  
  #################################################################
  #SLP: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure
  #################################################################
  SLP_MEAN <- aggregate(data$SLP,by=list(DATE),function(x) mean(x,na.rm=TRUE))
  SLP_MEAN <- lagvars(SLP_MEAN,unique(DATE),SLP_MEAN_LAG,daterange)
  
  #################################################################
  #SPD: WIND-OBSERVATION speed rate
  #################################################################
  WIND_SPD_MEAN <- aggregate(data$WIND_SPD,by=list(DATE),function(x) mean(x,na.rm=TRUE))
  WIND_SPD_MEAN <- lagvars(WIND_SPD_MEAN,unique(DATE),WIND_SPD_MEAN_LAG,daterange)
  
  #################################################################
  #WIND-OBSERVATION direction angle
  #################################################################
  #azimuth frequencies
  bins=seq(5,365,by=30)
  #hist(data$WIND_DIR,breaks=bins,freq=FALSE)
  
  #number of observations
  #WIND_DIR_NOBS <- aggregate(data$WIND_DIR,by=list(DATE),length)
  #colnames(WIND_DIR_NOBS) <- c('DATE','NOBS')
  
  #daily frequency count
  WIND_DIR_RF <- aggregate(data$WIND_DIR,by=list(DATE),function(x) histcount(x,bins))
  WIND_DIR_RF <- lagvars(WIND_DIR_RF,unique(DATE),WIND_DIR_RF_LAG,daterange)
  
  #################################################################
  #SKY-CONDITION-OBSERVATION ceiling height dimension
  #################################################################
  #ceiling height frequencies
  bins=c(0,500,seq(1000,8000,by=3500),22000)
  #bins=c(0,500,seq(1000,8000,by=2000),22000)
  #hist(data$CEIL_HGT,breaks=bins,freq=FALSE)
  
  #number of observations
  #CEIL_HGT_NOBS <- aggregate(data$CEIL_HGT,by=list(DATE),length)
  #hist(NOBS_CEIL_HGT$x,breaks=seq(0,24,by=1),freq=FALSE)
  
  #daily frequency count
  CEIL_HGT_RF <- aggregate(data$CEIL_HGT,by=list(DATE),function(x) histcount(x,bins))
  CEIL_HGT_RF <- lagvars(CEIL_HGT_RF,unique(DATE),CEIL_HGT_RF_LAG,daterange)
  
  #function to merge multiple data frames
  merge.all <- function(by, ...) {
    frames <- list(...)
    return (Reduce(function(x, y) {merge(x, y, by = by, all = TRUE)}, frames))
  }
  
  #weather <- merge.all(by = "DATE", TEMP_MEAN,DEWPT_MEAN,SLP_MEAN,WIND_SPD_MEAN,WIND_DIR_RF,CEIL_HGT_RF)
  weather <- merge.all(by = "DATE", TEMP_MEAN,SLP_MEAN,WIND_SPD_MEAN,WIND_DIR_RF,CEIL_HGT_RF)
  weather <- weather[!apply(weather,1,function(y)any(is.na(y))),]
  rownames(weather) <- NULL
  
  #find and unlist vector columns
  WEATHER <- flat(weather)
  
  #round values
  WEATHER[,sapply(WEATHER,is.numeric)] <- round(WEATHER[,sapply(WEATHER,is.numeric)],digits=3)
  
  #export to .csv file
  tablepath <- paste(outputpath,outputfile,".csv",sep="") 
  write.csv(WEATHER, file=tablepath,row.names = FALSE)
  
  #export to .Rda file
  save(WEATHER,file=paste(outputpath,outputfile,".Rda",sep=""))
  
  } else {
    nerrors=nerrors+1
    errors[nerrors]<-paste("Station ", myusaf[k], "has missing attributes or too few observations, unable to create features")
  }
  setTxtProgressBar(pb, k)
}
close(pb)

#write error log file
fileConn<-file(paste(outputpath,"error_log.txt",sep=""))
writeLines(errors, fileConn)
close(fileConn)


# Reference
##############################################################################

#WIND-OBSERVATION direction angle
#The angle, measured in a clockwise direction, between true north and the 
#direction from which the wind is blowing. 
#Default Value:999
#999: Missing.  If type code (below) = V, then 999 indicates variable wind direction.


#WIND-OBSERVATION direction quality code
#The code that denotes a quality status of a reported WIND-OBSERVATION direction
#Default Value:9
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#WIND-OBSERVATION type code
#The code that denotes the character of the WIND-OBSERVATION.
#Default Value:9
#Table of Values:
#A: Abridged Beaufort
#B: Beaufort
#C: Calm
#H: 5-Minute Average Speed
#N: Normal
#Q: Squall
#R: 60-Minute Average Speed
#T: 180 Minute Average Speed
#V: Variable


#SPD: WIND-OBSERVATION speed rate
#The rate of horizontal travel of air past a fixed point. 
#Unit:Meters per Second
#Default Value:9999


#WIND-OBSERVATION speed quality code
#The code that denotes a quality status of a reported WIND-OBSERVATION speed 
#rate.
#Default Value:9
#Table of Values: 
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#HGT: SKY-CONDITION-OBSERVATION ceiling height dimension
#The height above ground level (AGL) of the lowest cloud or obscuring phenomena 
#layer aloft with 5/8 or more summation total sky cover, which may be predominantly 
#opaque, or the vertical visibility into a surface-based obstruction.
#Unit:Meters
#Default Value:99999
#Table of Values: 
#22000: Unlimited


#Q: SKY-CONDITION-OBSERVATION ceiling quality code
#The code that denotes a quality status of a reported ceiling height dimension.
#Length:1
#Default Value:9
#Table of Values:
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check, data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#VISBY: VISIBILITY-OBSERVATION distance dimension
#The horizontal distance at which an object can be seen and identified. 
#Unit:Meters
#Table of Values:  
#Missing: 999999
#NOTE: Values greater than 160000 are entered as 160000


#Q: VISIBILITY-OBSERVATION distance quality code
#The code that denotes a quality status of a reported distance of a visibility 
#observation. 
#Length:1
#Default Value:9
#Table of Values:
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#I: VISIBILITY-OBSERVATION variability code
#The code that denotes whether or not the reported visibility is variable. 
#Length:1
#Default Value:9
#Table of Values: 
#N: Not variable
#V: Variable


#Q: VISIBILITY-OBSERVATION quality variability code
#The code that denotes a quality status of a reported VISIBILITY-OBSERVATION 
#variability code. 
#Length:1
#Default Value:9
#Table of Values: 
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#TEMP: AIR-TEMPERATURE-OBSERVATION air temperature
#The temperature of the air. 
#Length:5
#Scale:10
#Unit:Degrees Celsius
#Default Value:+9999


#Q: AIR-TEMPERATURE-OBSERVATION air temperature quality code
#The code that denotes a quality status of an AIR-TEMPERATURE-OBSERVATION. 
#Length:1
#Default Value:9
#Table of Values:
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#DEWPT: AIR-TEMPERATURE-OBSERVATION-DEWPOINT temperature
#The temperature to which a given parcel of air must be cooled at constant 
#pressure and water vapor content in order for saturation to occur. 
#Length:5
#Scale:10
#Unit:Degrees Celsius
#Default Value:+9999


#Q: AIR-TEMPERATURE-OBSERVATION-DEWPOINT quality code
#The code that denotes a quality status of the reported dew point temperature. 
#Length:1
#Default Value:9
#Table of Values: 
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check, data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present


#SLP: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure
#The air pressure relative to Mean Sea Level (MSL). 
#Length:5
#Scale:10
#Unit:Hectopascals
#Default Value:99999


#Q: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure quality code
#The code that denotes a quality status of the sea level pressure of an 
#ATMOSPHERIC-PRESSURE-OBSERVATION. 
#Length:1
#Default Value:9
#Table of Values: 
#0: Passed gross limits check
#1: Passed all quality control checks
#2: Suspect
#3: Erroneous
#4: Passed gross limits check , data originate from an NCDC data source
#5: Passed all quality control checks, data originate from an NCDC data source
#6: Suspect, data originate from an NCDC data source
#7: Erroneous, data originate from an NCDC data source
#9: Passed gross limits check if element is present

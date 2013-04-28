# File - extract_ephemeris.r
# Version - 19.03.2013
# Author - Matthew Parkan
# Description - computes daily daylight duration from JPL HORIZON ephemeris
# Data - http://ssd.jpl.nasa.gov/?horizons

#clear workspace
rm(list=ls())

#IMPORTANT!!! Define working directories
inputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\raw data\\JPL Horizons Ephemeris\\"
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\JPL Horizons Ephemeris\\"

# IMPORTANT!!! Define USAF number of desired stations (check station inventory)
myusaf <- c("043200","043390","043600","043900","042200","042020")
 
# Place the USAF number at the beginning of the file name (check reference below or station metadata file)
# "043200",NA,"DANMARKSHAVN","GL","GL","","BGDH",76.767,-18.667,12
# "043390",NA,"ITTOQQORTOORMIIT /S","GL","GL","","BGSC",70.483,-21.95,69
# "043600",NA,"TASIILAQ /AMMASSALI","GL","GL","","BGAM",65.6,-37.633,52
# "043900",NA,"PRINS CHRISTIAN SUN","GL","GL","","BGPC",60.05,-43.167,75
# "042200",NA,"AASIAAT /EGEDESMIND","GL","GL","","BGEM",68.7,-52.85,41
# "042020",NA,"PITUFFIK (THULE AB)","GL","GL","","BGTL",76.533,-68.75,59

# list csv files in the input folder
csvlist <- list.files(inputpath,recursive=FALSE,pattern = "\\.csv$")
station_id <- substr(csvlist, 1, 6) 
validind <- station_id %in% myusaf
csvlist <- subset(csvlist, validind==TRUE)
station_id <- subset(station_id, validind==TRUE)

numfiles <- length(csvlist)
#csvlist <- gsub("/", "\\\\", csvlist)

#time difference function
daydur <- function(input){
  if (length(input)==3){
    res <- round(as.numeric(difftime(input[3],input[1],units='hours')),digits = 2)
  } else {res <- round(24,digits = 2) }
  return(res) 
}

pb<-txtProgressBar(min = 0, max = numfiles, style = 3) #progress bar
for(k in 1:numfiles){ 
  
  outputfile <- paste(station_id[k],"_DAYLIGHT",sep="") # the output file will be created in the work directory
  
  #Import data
  file <- paste(inputpath,csvlist[k],sep="")
  data <- read.csv(file, header=FALSE,colClasses=c(rep("character",6)))
  data <- as.data.frame(data)
  
  #Define column headers
  colnames(data) <- c('YMDHM', 'IND', 'TYPE','M1','M2','M3')
  
  #Date and Time Code (UTC)
  TIMESTAMP_UTC <- strptime(data$YMDHM, format="%Y-%b-%d %H:%M", tz="UTC")
  DATE <- as.Date(TIMESTAMP_UTC, format="%Y-%b-%d")
  DOY <- as.numeric(format(TIMESTAMP_UTC, format = "%j")) 
  YEAR = as.numeric(format(TIMESTAMP_UTC, format = "%Y"))
  MONTH = as.numeric(format(TIMESTAMP_UTC, format = "%m"))
  DAY = as.numeric(format(TIMESTAMP_UTC, format = "%d"))
  
  #Create time vector (for dates with no daylight)
  DATEVECTOR <- seq(as.Date("1982-01-01"),as.Date("2013-01-01"),by=1)
  LIGHTVECTOR <- rep(0, length(DATEVECTOR))
  
  LONGNIGHTS <- data.frame(DATEVECTOR,LIGHTVECTOR)
  colnames(LONGNIGHTS) <- c('DATE', 'DAYLIGHT_HR') 
  
  #Aggregate ephemeris by date
  DAILY_RTS <- aggregate(TIMESTAMP_UTC,by=list(DATE),daydur)
  colnames(DAILY_RTS) <- c('DATE','DAYLIGHT_HR')
  
  #Find dates with no daylight
  missing_dates <- LONGNIGHTS$DATE %in% DAILY_RTS$DATE 
  NIGHTDATES <- subset(LONGNIGHTS, !missing_dates)
  rownames(NIGHTDATES) <- NULL
  
  #Bind and sort by date
  DAILY_RTS <- rbind(DAILY_RTS,NIGHTDATES)
  DAILY_RTS <- DAILY_RTS[order(DAILY_RTS$DATE),]
  rownames(DAILY_RTS) <- NULL

  EPHEMERIS<-DAILY_RTS
  
  #Export to .csv file
  tablepath <- paste(outputpath,outputfile,".csv",sep="") 
  write.csv(EPHEMERIS, file=tablepath,row.names = FALSE)
  
  #Export to .Rda file
  save(EPHEMERIS,file=paste(outputpath,outputfile,".Rda",sep=""))
  
  setTxtProgressBar(pb, k)
}
close(pb)
      
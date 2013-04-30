# File - merge_features.r
# Version - 17.04.2013
# Author - Matthew Parkan
# Description - merge ephemeris, weather and remote sensing features

#clear workspace
rm(list=ls())

#IMPORTANT!!! Define input directories
inputpath1 <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\NCDC weather\\" #path to folder containing weather feature files
inputpath2 <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\JPL Horizons Ephemeris\\" #path to folder containing weather feature files
inputpath3 <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\SST_AVHRR_OI\\" #path to folder containing remote sensing statistics files

#IMPORTANT!!! Define output directory
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\Features\\"

#IMPORTANT!!! Define USAF number of desired stations (check station inventory)
# "043200",NA,"DANMARKSHAVN","GL","GL","","BGDH",76.767,-18.667,12
# "043390",NA,"ITTOQQORTOORMIIT /S","GL","GL","","BGSC",70.483,-21.95,69
# "043600",NA,"TASIILAQ /AMMASSALI","GL","GL","","BGAM",65.6,-37.633,52
# "043900",NA,"PRINS CHRISTIAN SUN","GL","GL","","BGPC",60.05,-43.167,75
# "042200",NA,"AASIAAT /EGEDESMIND","GL","GL","","BGEM",68.7,-52.85,41
# "042020",NA,"PITUFFIK (THULE AB)","GL","GL","","BGTL",76.533,-68.75,59
myusaf <- c("043200","043390","043600","043900","042200","042020")
#myusaf <- c("042020")

#function to merge multiple data frames
merge.all <- function(by, ...) {
  frames <- list(...)
  return (Reduce(function(x, y) {merge(x, y, by = by, all = TRUE)}, frames))
}


numfiles<-length(myusaf)
pb<-txtProgressBar(min = 0, max = numfiles, style = 3) #progress bar
for(k in 1:numfiles){ 
  
  #load weather
  file1 <- list.files(inputpath1,recursive=FALSE,pattern = paste(myusaf[k],"\\w+\\.Rda$",sep=""))
  load(paste(inputpath1,file1,sep=""))
  F_WEATHER <- WEATHER
  
  #load ephemeris
  file2 <- list.files(inputpath2,recursive=FALSE,pattern = paste(myusaf[k],"\\w+\\.Rda$",sep=""))
  load(paste(inputpath2,file2,sep=""))
  F_EPHEMERIS <- EPHEMERIS
  
  #load remote sensing statistics
  file3 <- list.files(inputpath3,recursive=FALSE,pattern = paste(myusaf[k],"\\w+\\.Rda$",sep=""))
  load(paste(inputpath3,file3,sep=""))
  F_RS_STATS <- RS_STATS
  
  #merge features
  FEATURES <- merge.all(by = "DATE",F_WEATHER,F_RS_STATS,F_EPHEMERIS)
  FEATURES <- FEATURES[!apply(FEATURES,1,function(y)any(is.na(y))),]
  rownames(FEATURES) <- NULL
  
  #Export to .csv file
  tablepath <- paste(outputpath,myusaf[k],"_FEATURES",".csv",sep="") 
  write.csv(FEATURES, file=tablepath,row.names = FALSE)
  
  #Export to .Rda file
  save(FEATURES,file=paste(outputpath,myusaf[k],"_FEATURES",".Rda",sep=""))

  setTxtProgressBar(pb, k)
}
close(pb)

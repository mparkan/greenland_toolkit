# File - predict_self.r
# Version - 20.04.2013
# Author - Matthew Parkan
# Description - Support vector regression

#load librairies
library(e1071)

#clear workspace
rm(list=ls())

#IMPORTANT!!! Define input directory
inputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\Features\\"

#IMPORTANT!!! Define output directory
outputpath <- "C:\\Users\\mat\\Google Drive\\Greenland\\processed data\\"

#IMPORTANT!!! Define USAF number of desired stations (check feature availability before)
# "043200",NA,"DANMARKSHAVN","GL","GL","","BGDH",76.767,-18.667,12
# "043390",NA,"ITTOQQORTOORMIIT /S","GL","GL","","BGSC",70.483,-21.95,69
# "043600",NA,"TASIILAQ /AMMASSALI","GL","GL","","BGAM",65.6,-37.633,52
# "043900",NA,"PRINS CHRISTIAN SUN","GL","GL","","BGPC",60.05,-43.167,75
# "042200",NA,"AASIAAT /EGEDESMIND","GL","GL","","BGEM",68.7,-52.85,41
# "042020",NA,"PITUFFIK (THULE AB)","GL","GL","","BGTL",76.533,-68.75,59
#myusaf <- c("043200","043390","043600","043900","042200","042020")
myusaf <- c("042020")

#IMPORTANT!!! Define which month(s) should be predicted
months=seq(1,12,1)

#IMPORTANT!!! How many runs should be performed?
nruns=15

#IMPORTANT!!! Should support vector regression (epsilon) be performed? (TRUE= yes, FALSE=no)
svreg=TRUE

#IMPORTANT!!! Should linear regression be performed? (TRUE= yes, FALSE=no)
lreg=TRUE


runnum=sort(rep(1:nruns,length(months)))
months=rep(months,nruns)
nmonths=length(months)
nstat=length(myusaf)
step=0
pb<-txtProgressBar(min = 0, max = nmonths*nstat, style = 3) #progress bar

for(j in 1:nstat){
  
  #initialize performance data frame
  PERFORMANCE <- data.frame(RUN=numeric(nmonths),
                            MONTH=numeric(nmonths),
                            NTRAIN=numeric(nmonths),
                            NTEST=numeric(nmonths),
                            TRAINRATIO=numeric(nmonths),
                            LR.TEST.RMSE=numeric(nmonths),
                            LR.TEST.NRMSE=numeric(nmonths),
                            SVR.TEST.RMSE=numeric(nmonths),
                            SVR.TEST.NRMSE=numeric(nmonths))
  
  #initialize list of models
  MODELS <-list()
  
  #load features
  file <- list.files(inputpath,recursive=FALSE,pattern = paste(myusaf[j],"\\w+\\.Rda$",sep=""))
  load(paste(inputpath,file,sep=""))
  
  #split dataset by month
  #DOY <- as.numeric(format(FEATURES$DATE, format = "%j")) 
  #YEAR = as.numeric(format(FEATURES$DATE, format = "%Y"))
  MONTH = as.numeric(format(FEATURES$DATE, format = "%m"))
  monthly_features <- split(FEATURES, MONTH)
  
  for(k in 1:length(months)){
    
    #create training and test sets
    #############################################################
    FEATURES<-monthly_features[[months[k]]]
    
    #randomly sample observations to create training and test sets
    FEATURES$DATE <- NULL #remove date column
    nobs <- nrow(FEATURES)
    #index_train <- sample(nobs, ceiling(0.15*nobs))
    index_train <- sample(nobs, 230)
    index_test <- (1:nobs %in% index_train) == FALSE 
    train_set <- FEATURES[index_train, ]
    test_set <- FEATURES[index_test, ]
    
    #check for constant features and remove them if necessary
    train_const <- (apply(train_set, 2, sd))==0
    test_const <- (apply(test_set, 2, sd))==0
    const_features <- train_const | test_const
    train_set[,which(const_features)] <- list(NULL)
    test_set[,which(const_features)] <- list(NULL)
    
    
    #linear regression
    #############################################################
    
    if(lreg==TRUE){
      #create model
      model.lr <- lm(TEMP_MEAN_M0 ~ ., 
                     data = train_set)
      summary(model.lr)
      
      #predict
      pred.lr <- predict(model.lr, test_set, decision.values = TRUE)
      
      #performance
      LR.TEST.RMSE <- round(as.numeric(sqrt(crossprod(pred.lr-test_set$TEMP_MEAN_M0) / length(pred.lr))),digits=2)
      LR.TEST.NRMSE <- round(100*LR.TEST.RMSE/(max(FEATURES$TEMP_MEAN_M0)-min(FEATURES$TEMP_MEAN_M0)), digits=1)
    }   else {
      model.lr <- NA
      LR.TEST.RMSE <- NA
      LR.TEST.NRMSE <- NA
    }
    
    #support vector regression (epsilon)
    #############################################################
    
    if(svreg==TRUE){
      #tune free parameters (epsilon SVR)
      tm2 <- system.time({
        tobj <- tune.svm(TEMP_MEAN_M0 ~ ., 
                         data = train_set, 
                         type = "eps-regression", 
                         kernel = "linear",
                         cost = 2^seq(-10.5,-6,0.1)) #epsilon=seq(0.2,0.3,0.05)
      })
      summary(tobj)
      
      plot(tobj, xlab = "C", main="Parameter tuning")
      #plot(tobj,type = "contour", xlab = "C", ylab = expression(epsilon), main="Parameter tuning",color.palette=terrain.colors,nlevels=50 )
      #color.palette=heat.colors
      
      #plot(tobj, transform.x = log10,ylab = "C")
      bestC <- tobj$best.parameters[[1]]
      #bestEps <- tobj$best.parameters[[2]]
      
      #create model
      model.svr <- svm(TEMP_MEAN_M0 ~ ., 
                       data = train_set,
                       type = "eps-regression", 
                       kernel = "linear",
                       cost = bestC,#epsilon=bestEps,
                       cross = 10)
      summary(model.svr)
      
      #predict
      pred.svr <- predict(model.svr, test_set, decision.values = TRUE)
      
      #determine feature weights
      #weights <- abs(t(model.svr$coefs) %*% model.svr$SV)
      #weights.names <- colnames(weights)
      #names.sorted <- weights.names[order(weights)]
      #weigths.sorted <- weights[order(weights)]
      
      #performance
      SVR.TEST.RMSE <- round(as.numeric(sqrt(crossprod(pred.svr-test_set$TEMP_MEAN_M0) / length(pred.svr))),digits=2)
      SVR.TEST.NRMSE <- round(100*SVR.TEST.RMSE/(max(FEATURES$TEMP_MEAN_M0)-min(FEATURES$TEMP_MEAN_M0)),digits=1)  
    }  else {
      model.svr <- NA
      SVR.TEST.RMSE <- NA
      SVR.TEST.NRMSE <- NA
    }
    
    #save models
    MODELS[[k]] <- list(runnum[k],months[k],model.lr,model.svr)
    
    #save performance results
    PERFORMANCE$RUN[k]=runnum[k]
    PERFORMANCE$MONTH[k]=months[k]
    PERFORMANCE$NTRAIN[k]=nrow(train_set)
    PERFORMANCE$NTEST[k]=nrow(test_set)
    PERFORMANCE$TRAINRATIO[k]=round(PERFORMANCE$NTRAIN[k]/nobs, digits=2)
    PERFORMANCE$LR.TEST.RMSE[k]=LR.TEST.RMSE
    PERFORMANCE$LR.TEST.NRMSE[k]=LR.TEST.NRMSE
    PERFORMANCE$SVR.TEST.RMSE[k]=SVR.TEST.RMSE
    PERFORMANCE$SVR.TEST.NRMSE[k]=SVR.TEST.NRMSE
    
    step=step+1
    setTxtProgressBar(pb, step)
    
  }
  
  #export models to .Rda file
  save(MODELS,file=paste(outputpath,"Models\\S",myusaf[j],"_models_",nruns,"runs.Rda",sep=""))
  
  #export performance to .csv file
  tablepath <- paste(outputpath,"Predictions\\S",myusaf[j],"_performance_",nruns,"runs.csv",sep="") 
  write.csv(PERFORMANCE, file=tablepath,row.names = FALSE)
  
  #export performance to .Rda file
  save(PERFORMANCE,file=paste(outputpath,"Predictions\\S",myusaf[j],"_performance_",nruns,"runs.Rda",sep=""))
  
}

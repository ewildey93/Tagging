#load the packages you need
library(camtrapR)
library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(stringr)
#read in output from Timelapse, change to your file path
setwd("C:/Users/eliwi/OneDrive/Documents/NDOW-Lions/Tagging")
Timelapse <- read.csv(file = "C:/Users/eliwi/OneDrive/Documents/NDOW-Lions/Tagging/6F_DEL_230411-231013/6F_DEL_230411-231013.csv", na.strings = "")
#change name of DateTime column
colnames(Timelapse)[4] <- "DateTimeOriginal"

Timelapse$CAM_ID <- str_extract(Timelapse$RootFolder, "^.{2}")
Timelapse$MTN <- ifelse(grepl(pattern = "DEL", x = Timelapse$RootFolder) == TRUE, "DEL", "CLO")
######read temperature data off of images######################
#  slow step, could be problematic for bigger photo folders   #
###############################################################
#change path argument to path of photos folder
pics <- list.files(path="./6F_DEL_230411-231013/", pattern=".jpg$", full.names=T, recursive=T, ignore.case=T)
for (p in 1:length(pics)){
  IMAGE  <- magick::image_read(pics[p]) 
  #temperature
  IMAGE2  <- magick::image_crop(image = IMAGE, geometry = "273x60+958+1236") #crop image
  TempF <- tesseract::ocr_data(IMAGE2) # ocr the cropped image
  Timelapse[p,"TEMPERATURE"] <- TempF$word[1]
  
}
#check to see if temperature data was read correctly off images, should return: integer(0)
str_which(Timelapse$TEMPERATURE, "[:punct:]")
#remove F for fahrenheit from Temperature category and turn it into numeric data
Timelapse$TEMPERATURE <- as.numeric(gsub(pattern= "F", replacement= "",x=Timelapse$TEMPERATURE))

#for rows with SPP2 add a new row to dataframe for that detection
#this will show how many rows have SPP2 and therefore how many rows will be added
table(!is.na(Timelapse$SPP2))
#run the loop to add rows
loopnum <- nrow(Timelapse)
for(i in 1:loopnum){
if (!is.na(Timelapse$SPP2[i]) == TRUE){
  row <- Timelapse[i,]
  row[,"SPP"] <- Timelapse[i,"SPP2"]
  row[, "SPP2"] <- NA
  Timelapse[i,"SPP2"] <- NA
  Timelapse <- rbind(Timelapse, row)
  }
}
# do i need something else here to accounts for details in comments, or
#clear data from the rest of the row?

#add Date and time columns from datetime combined column
Timelapse$DATE <- as.Date(as.POSIXct(Timelapse$DateTimeOriginal, format="%Y-%m-%d %H:%M:%S"))
Timelapse$TIME <- format(as.POSIXct(Timelapse$DateTimeOriginal, format="%Y-%m-%d %H:%M:%S"), "%H:%M:%S")


#########################
#######function##########
# create function to assess temporal independence between records  for bouts ####
#should see assessTemporalIndependence appear in environment pane under functions once run#
#this is copied and modified from camtrapr:::assessTemporalIndependence function
#script will jump to line 199 when run, note this is just creating the function not running it on data
#I've only set this up to run for our specific set up only 1 camera at a site (CamerasIndependent=FALSE)

assessTemporalIndependence2 <- function(intable,
                                        deltaTimeComparedTo,
                                        columnOfInterest,     # species/individual column
                                        cameraCol,
                                        camerasIndependent,
                                        stationCol,
                                        minDeltaTime,
                                        eventSummaryColumn,
                                        eventSummaryFunction)
{
  # check if all Exif DateTimeOriginal tags were read correctly
  if(any(is.na(intable$DateTimeOriginal))){
    which.tmp <- which(is.na(intable$DateTimeOriginal))
    if(length(which.tmp) == nrow(intable)) stop("Could not read any Exif DateTimeOriginal tag at station: ", paste(unique(intable[which.tmp, stationCol])), " Consider checking for corrupted Exif metadata.")
    warning(paste("Could not read Exif DateTimeOriginal tag of", length(which.tmp),"image(s) at station", paste(unique(intable[which.tmp, stationCol]), collapse = ", "), ". Will omit them.\nConsider checking for corrupted Exif metadata. Or does your selected time zone have daylight saving time and the image(s) fall in the misisng hour at spring formward (cameras don't usually record DST)?. \n",
                  paste(file.path(intable[which.tmp, "Directory"],
                                  intable[which.tmp, "FileName"]), collapse = "\n")), call. = FALSE, immediate. = TRUE)
    intable <- intable[-which.tmp ,]
    rm(which.tmp)
  }
  
  # prepare to add time difference between observations columns
  intable <- data.frame(intable,
                        delta.time.secs  = NA,
                        delta.time.mins  = NA,
                        delta.time.hours = NA,
                        delta.time.days  = NA,
                        independent      = ifelse(minDeltaTime == 0, TRUE, NA),   # all independent if no temporal filtering
                        stringsAsFactors = FALSE,
                        check.names      = FALSE)        # to prevent ":" being converted to ".", e.g. in EXIF:Make
  
  # sort records by station, species, then time
  intable <- intable[order(intable[, stationCol], intable[, columnOfInterest], intable$DateTimeOriginal),] #, intable[, columnOfInterest]
  
  for(xy in 1:nrow(intable)){     # for every record
    
    
    which.columnOfInterest <- which(intable[, columnOfInterest]  == intable[xy, columnOfInterest])                     # same species/individual from either of first 2 columns
    which.stationCol       <- which(intable[, stationCol]        == intable[xy, stationCol])                # at same station
    which.independent      <- which(intable$independent          == TRUE)                                   # independent (first or only record of a species at a station)
    #which.earlier          <- which(intable$DateTimeOriginal     <  intable$DateTimeOriginal[xy])          # earlier than record xy (takes long)
    which.earlier          <- 1: (xy-1)                                                                  # earlier than record xy  (fast alternative, relies on table being sorted by date/time before anything else)
    if(camerasIndependent) {
      which.cameraCol      <- which(intable[, cameraCol]  == intable[xy, cameraCol])                        # at same camera
    }
    
    # set independent = TRUE and delta.time = 0 if it is the 1st/only  record of a species / individual
    
    if(camerasIndependent == TRUE){
      which.tmp <- Reduce(intersect, list(which.columnOfInterest, 
                                          which.stationCol, 
                                          which.cameraCol))
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which.tmp])){    # cameras at same station assessed independently
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    } else {
      #here
      #find common elements of multiple vectors, which index values are the same, subset to those
      which.tmp <- Reduce(intersect, list(which.columnOfInterest, #species
                                          which.stationCol))      #camera
      duplicateTime <- which(duplicated(intable$DateTimeOriginal) ==TRUE)
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which.tmp]) & !(xy %in% duplicateTime)) {
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      } 
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which.tmp]) & (xy %in% duplicateTime)){
        intable$independent[xy]       <- FALSE
        intable$delta.time.secs[xy]   <- 0
      }
    }
    
    # calculate time difference to previous records of same species at this station (if not the 1st/only record)
    if(is.na(intable$delta.time.secs[xy])) {
      
      if(deltaTimeComparedTo == "lastIndependentRecord"){
        
        if(camerasIndependent == TRUE){
          which_time2 <- Reduce(intersect, list(which.columnOfInterest, 
                                                which.stationCol,
                                                which.cameraCol,
                                                which.independent,
                                                which.earlier))
        } else {
          which_time2 <- Reduce(intersect, list(which.columnOfInterest, 
                                                which.stationCol,
                                                which.independent,
                                                which.earlier))
        }
      }  else {    # if(deltaTimeComparedTo == "lastRecord"){'
        if(camerasIndependent  == TRUE){
          which_time2 <- Reduce(intersect, list(which.columnOfInterest, 
                                                which.stationCol,
                                                which.cameraCol,
                                                which.earlier))
        } else {
          #here
          which_time2 <- Reduce(intersect, list(which.columnOfInterest, 
                                                which.stationCol,
                                                which.earlier))
        }
      }
      
      # time difference to last (independent) record
      diff_tmp <- min(na.omit(difftime(time1 = intable$DateTimeOriginal[xy],            # delta time to last independent record
                                       time2 = intable$DateTimeOriginal[which_time2],
                                       units = "secs")))
      
      # save delta time in seconds
      intable$delta.time.secs[xy] <-  diff_tmp
      if(intable$delta.time.secs[xy] >= (minDeltaTime * 60)){
        intable$independent[xy] <- TRUE
      } else {
        intable$independent[xy] <- FALSE
      }
      
    }   # end   if(intable$DateTimeOriginal[xy] == min(...)} else {...}
  }     # end for(xy in 1:nrow(intable))
  
  
  
  # keep only independent records
  outtable <- intable
  
  
  # compute delta time in hours and days
  outtable$delta.time.secs  <- round(outtable$delta.time.secs, digits = 0)
  outtable$delta.time.mins  <- round(outtable$delta.time.secs  / 60, digits = 1)
  outtable$delta.time.hours <- round(outtable$delta.time.mins  / 60, digits = 1)
  outtable$delta.time.days  <- round(outtable$delta.time.hours / 24, digits = 1)
  
  # remove "independent" column
  
  
  
  return(outtable)
}


##############################
####split data into bouts#####
##############################
#get only those photos where something was detected
Detections <- Timelapse[!is.na(Timelapse$SPP),]
#calculate time between bouts and "independence" of detections to determine bouts
outtable <- assessTemporalIndependence2(intable = Detections,deltaTimeComparedTo = "lastRecord",
                                         columnOfInterest = "SPP", camerasIndependent = FALSE,
                                         stationCol = "CAM_ID",minDeltaTime = 30)

#fill in bout column with sequential number
outtable$bout <- NA
for (i in 1:nrow(outtable)) {
  if (outtable$independent[i] == TRUE) {
    outtable$bout[i]= sum(outtable$independent[1:i] == TRUE)
  }else {
    outtable$bout[i] = outtable$bout[i-1]
  }
}

#######collapse data into bouts###########
#get first and last record for each bout since we need file and time start/end
head <- outtable%>%group_by(bout)%>%slice_head()
tail <- outtable%>%group_by(bout)%>%slice_tail()
tail2 <- tail[,c("File","TIME","bout")]
colnames(tail2)[1:2] <- c("FILE_END", "TIME_END")
#merge first and last record of each bout into same row
boutable <- merge(head, tail2)
#average temperature for each bout, based on Hannah's rule 
Temp <- as.numeric()
for (b in 1:length(unique(outtable$bout))){
  bouts <- outtable[outtable$bout == b,]
  print(b)
  diff <- max(bouts$TEMPERATURE) - min(bouts$TEMPERATURE)
  if (diff < 5) {
    Temp <- append(Temp,median(bouts$TEMPERATURE))
  } else {
    Temp <- append(Temp, mean(bouts$TEMPERATURE))
  }
}
#merge temperature data with bout data
boutable <- cbind(boutable, Temp )
#start to assemble dataframe in manner of Bout Template excel file, picking out columns we need
boutable2 <- boutable[,c("CAM_ID", "MTN", "DATE", "TIME", "TIME_END", "File", "FILE_END", 
                         "Temp", "SPP", "L_ANTLER", "R_ANTLER", "PHENO", "ENTER", "ENTER_2",
                         "NOTES")]
colnames(boutable2)[colnames(boutable2) == "TIME"] <- "TIME_START"
colnames(boutable2)[colnames(boutable2) == "File"] <- "FILE_START"

#add sex/age columns to dataframe
boutable3 <- data.frame(boutable2[,1:9], AD_M="", AD_F="", AD_U="", SUB_M="",
                        SUB_F="", SUB_U="", YOY_M="", YOY_F="", YOY_U="", 
                        U_U="", boutable2[,10:15])

#for filling in adult/sex columns if bout is only 1 photo
#for (b in 1:length(boutable3$bout)){
 # if (table(outtable$bout)[b] > 1) {boutable3[b, 10:19] <- ""} 
  #else{
    
#  }
#}
# age sex combo columns, has root folder, relative path columns

write.csv(boutable3, paste0(Timelapse$RootFolder, "bout.csv"))

#################################
###########scrap#################
#################################
#arguments to test assessTemporalIndependence2###
#######################
#intable <- Detections
#deltaTimeComparedTo = "lastRecord"
#columnOfInterest = "SPP"
#camerasIndependent = FALSE
#stationCol = "CAM_ID"
#minDeltaTime = 30



#check to see difference between lastIndependentRecord and lastRecord##
#################################################
#table(outtable$independent)
#table(outtable2$independent)
#which(outtable$independent ==TRUE)
#outtest <- as.data.frame(cbind(outtable$independent,outtable2$independent))
#table(outtest$V1 == outtest$V2)


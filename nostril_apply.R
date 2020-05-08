
setwd("C:/Users/vdemartsev/Google Drive/Thermal meerkats/Auto nose tracking/data")
library(Thermimage)
library(ggplot2)
library(stringr)
library(pracma)
library(KernSmooth)


######reading in nostril detection data and combining it with extracted color values######

#setting some cutoff values  
#lenght of sequential drame detection
frame_seq <- 100
# spike hight
temp_cutoff <- 1.5

#reading in data from loopy and imageJ

imageJ_output_dir <- "C:/Users/vdemartsev/Google Drive/Thermal meerkats/Auto nose tracking/data/ImageJ output"
loopy_output_dir <- "C:/Users/vdemartsev/Google Drive/Thermal meerkats/Auto nose tracking/data/Loopy output"

imageJ_files <- list.files (imageJ_output_dir, pattern = ".csv")
loopy_files <- list.files (loopy_output_dir,pattern = ".csv" )

#trim file names for matching

loopy_files_trimmed <- substring(loopy_files, 47, nchar(loopy_files)- 30)
loopy_files <- cbind(loopy_files, loopy_files_trimmed)
loopy_files[ ,2] <- paste(loopy_files[ ,2], "_out.csv", sep = "")                   
loopy_files <- data.frame(loopy_files, stringsAsFactors = F)
#make a mathing conversion table for video file names
proc_file_list <- data.frame(intersect(loopy_files[ ,2] ,imageJ_files ), stringsAsFactors = F)
proc_file_list <-  merge(loopy_files, proc_file_list, by.x = "loopy_files_trimmed", by.y = "intersect.loopy_files...2...imageJ_files.", all.x = F)
proc_file_list <- proc_file_list[!duplicated(proc_file_list$loopy_files_trimmed ), ]


#make empty data frame for collectina all input
temp_mes <- data.frame()

for (i in 1:nrow(proc_file_list)) {
  
  
  #read loopy output file
  loopy <- read.csv(paste(loopy_output_dir,"/", proc_file_list[i,2], sep=""))
  
  #correct slice count by adding  1 frame to each
  loopy$imageJframe <- loopy$frame_number+1
  
  #sorting rows by nostrilR and nostrilL
  loopy <-  loopy[order(loopy$name, decreasing = T),]
  
  #getting file name
  video_file_name <- proc_file_list[i,1]
  
  #read imageJ output file
  imageJ <- read.csv(paste(imageJ_output_dir, "/" ,  video_file_name, sep=""))
  
  #combine frame numbers with the correct nostril side
  imageJ <- cbind(imageJ, loopy$name, loopy$imageJframe, video_file_name)
  
  #collect all into one data frame
  temp_mes <- rbind(temp_mes, imageJ)
}

#####convert colour intensity to temperature ######

temp_mes$Median <- as.numeric(temp_mes$Median)
# Obtain the calibration constants from the csq files
# I used Exiftool (a command line utility to do this), but have also
# included this functionality in Thermimage 

PR1<-21546.2
PB<-1507
PF<-1
PO<-(-6331)
PR2<-0.01622949

# Set the object parameters required for the raw -> temperature
# conversion
# Emissivity, Object Distance, Reflected Temperature, Atmospheric Temperature,
# IR Window Temperature, IRT 
E<-0.95
OD<-1
RTemp<-10
ATemp<-10
IRWTemp<-10
IRT<-1
RH<-30

# We know that we only report on temperatures above 0C and below 40C, 
# so let's truncate the raw values to a much more reasonable range:
rawvalues<-10000:20000
tempvalues<-raw2temp(rawvalues, E, OD, RTemp, ATemp, IRWTemp=RTemp, IRT=1, RH, 
                     PR1, PB, PF, PO, PR2)

# add these to a data.frame
d<-data.frame(rawvalues, tempvalues)

lm1<-lm(tempvalues ~ stats::poly(rawvalues, 4), data=d)

# or the easier approach is to use the predict function
temp_mes$Median <- predict(lm1, newdata=data.frame(rawvalues=temp_mes$Median))



######quality control of data####

#filtering out short tracks

#make empty data frame to collect long tracks
long_tracks_all <- data.frame()

#standartize Median
#temp_mes$ZMedian <- scale(temp_mes$Median)


for (j in 1:nrow(proc_file_list))
{
  #selecting one video file output
  one <- subset(temp_mes, temp_mes$video_file_name == proc_file_list[j,1])
  
  x <- 1
  
  #adding column for marking frame seqiences
  one$seg <- NA
  
  #marking frame seqiences without gaps above 5 frames
  for (i in 1:(nrow(one)-1)) 
  {
    
    if (one[i+1, "Slice"] - one[ i, "Slice"] < 6 & one[i+1, "Slice"] - one[ i, "Slice"] > 0)
    {one[i, "seg"] <-x}
    else
    {one[i, "seg"] <- x
    x <- x+1}
  }
  
  
  #finding sequence lengts
  seg_lenghts <- data.frame(table(one$seg))
  
  
  #finding tracks over set number of frames lenght
  track_num <- which(seg_lenghts[,2] > frame_seq)
  
  #if no long frame sequence in the file, go to  next
  if (length(track_num) < 1)
  {next}
  
  #filtering the color values and filtering out jumps
  
  #choosing long tracks only
  
  long_tracks <- one[one$seg %in% which(seg_lenghts[ ,2] > frame_seq),]                       
  
  #converting frames to time in seconds
  long_tracks$time <- long_tracks$Slice/30
  
  
  #set the temp change cutoff. change this into replacing the outlied with local average !!!!
  z <- temp_cutoff
  
  #duplicating color value column to filter later 
  long_tracks$filt_median <- long_tracks$Median
  
  #getting rid of color jumps above cutoff
  for (i in track_num)
  { 
    seg_select <-  subset(long_tracks, long_tracks$seg == i)
    seg_select$dist <- NA
    #checking first and last frame in the sequence
    seg_select[1 , "dist"] <- abs(seg_select[2, "Median"]- seg_select[1, "Median"])
    seg_select[nrow(seg_select) , "dist"] <- abs(seg_select[nrow(seg_select)-1, "Median"]- seg_select[nrow(seg_select), "Median"])
    #checking the rest of the values
    for (y in 1:(nrow(seg_select)-1)) 
    {
      seg_select[y+1 , "dist"] <- abs(seg_select[y+1, "Median"]- seg_select[y, "Median"])
      seg_select[which(seg_select$dist > z), "filt_median"] <- NA
      
    }
    
    #collecting all filtered long tracks back together into one data frame
    long_tracks_all <- rbind(long_tracks_all, seg_select)
  }
}
#add video file number to each row
long_tracks_all$video_file_num <- substr(long_tracks_all$video_file_name, 16, 18)

#make unique segment ID

long_tracks_all$seqID <- paste(long_tracks_all$video_file_name, long_tracks_all$seg)


#finding segments with both nostrils tracked

filt_file_list <- unique(long_tracks_all$video_file_name )
#male empty df for filling with filtered segments <- 
cleaned_files_list <- matrix(NA, ncol = 3, nrow = 0)
for (i in filt_file_list) {
 
  
  select_file_nostrilL <- subset(long_tracks_all, long_tracks_all$video_file_name == i & long_tracks_all$`loopy$name` == "NostrilL")
  select_file_nostrilR <- subset(long_tracks_all, long_tracks_all$video_file_name == i & long_tracks_all$`loopy$name` == "NostrilR")
  # if only one nostril tracked skip to next file
  l_dur <- nrow(select_file_nostrilL)
  r_dur <- nrow(select_file_nostrilR)
  
  if (l_dur == 0 | r_dur == 0) {next}
  
  #merge both nostrils into one data frame
  both_nostrils <- merge(select_file_nostrilR, select_file_nostrilL[ , 3:10], by.x = "Slice", by.y = "Slice")
  
  #get long frame sequences 
  frame_seq <- unique(both_nostrils$seg.x)
  
  for (j in frame_seq) {
    
    
    #plot temp curves for visual inspection
    seq_plot <- ggplot(data = both_nostrils[which(both_nostrils$seg.x == j) ,], aes(x = Slice)) + 
      geom_line (aes(y= Median.x), color = "red") + geom_line(aes(y = Median.y), color = "blue") +
      ggtitle(paste(i,"\n", "Red _Right_nost_seg_", j,"\n", "Blue_Left_nost_seg_", both_nostrils[max(which(both_nostrils$seg.x == j)) , "seg.y"] ))
    print(seq_plot)
    
    #capture file and sequence name
    R_seq <- c(i, j, "NostrilR")
    L_seq <- c (i , both_nostrils[max(which(both_nostrils$seg.x == j)) , "seg.y"], "NostrilL" )
    cleaned_files_list <- rbind(cleaned_files_list, R_seq, L_seq)
  }
}

#####at this point visually inspect curves and manually select ones showing potentually clean breathing traces ######
#####make list of files and sequence numbers for further analysis 

#select only relevant frame sequences
row.names(cleaned_files_list) <- NULL

## output a file for manual selection
#write.csv(cleaned_files_list, "segments_for_filt.csv")

#this eelection is made manually from visual inspection of the plots
##list for 200 seq cutoff
#cleaned_files_list <- data.frame(cleaned_files_list[c(108,105,103,102,99,90,89,80,77,75,72,69,67,65,62,60,57,55,47,44,35,33,29,26,23,17,15,14,5,1) , ])
#cleaned_files_list$fs <- c(1650, 4180,4600,2550,1620,3200,2800,1780,4500,2200,1780,4750,3600,2800,4300,1980,900,1460,2900,1100,2100,4400,4850,2450,300,2750,1450,0,4450,2120)
#cleaned_files_list$ef <- c(1800, 4480,4750,2650,1850,3400,3000,1870,4700,2450,1950,4920,4000,3400,4420,2120,1280,1600,3020,1250,2400,4550,5000,2700,450,3100,1550,130,4650,2300)
##list for 150 seq cutoff  with 5 frame gaps
cleaned_files_list <- data.frame(cleaned_files_list[c(9,36,37,38,61,62,83,94,98,103,114,117,122,127,138,144,154,155,156,157,160,161,165,166,171,172,176,177,179,180,185,186,187,191,193,194,207,209,215,217,230,233,235,243,247,256,262,271,279,286) , ])
cleaned_files_list$fs <- c(4450,1450,2750,2800,2400,2400,4400,50,2500,3125,3050,90,200,1100,2900,2800,0,3500,3525,4750,150,650,500,500,1120,980,1980,4115,4290,4290,3600,2950,4800,1500,100,100,2250,3460,3120,4475,975,1870,2790,1150,3250,390,2550,1230,4175,1675
)
cleaned_files_list$ef <- c(4650,1700,3100,3100,2700,2700,4750,120,2750,3300,3250,220,400,1250,3150,3100,450,3680,3680,4900,270,800,750,750,1300,1080,2120,4180,4420,4420,4100,3450,4920,1920,260,300,2310,3510,3225,4700,1070,2000,3010,1900,3340,500,2635,1300,4440,1790
)
rows <- nrow(cleaned_files_list)
cleaned_files_list[rows+1 ,] <-  data.frame("IR_2019-06-19_0519-0001-5000_out.csv",	4,	"NostrilR",	1900	,2450, stringsAsFactors = F)
cleaned_files_list[rows+2 ,] <- data.frame("IR_2019-06-19_0519-0001-5000_out.csv",	4,	"NostrilR",	2500	,3200, stringsAsFactors = F)

cleaned_files_list$runID <- seq(1, nrow(cleaned_files_list), 1)
#filter al video data according to manul inspection

cleaned_files_list <- stacomirtools::killfactor(cleaned_files_list)
long_tracks_all <- stacomirtools::killfactor(long_tracks_all)

#long_tracks_filt <- long_tracks_all[long_tracks_all$video_file_name %in% cleaned_files_list$X1,]

long_tracks_filt_frames <- data.frame(stringsAsFactors = F)

for (i in 1:nrow(cleaned_files_list))
{
  file_selct <- subset(long_tracks_all, long_tracks_all$video_file_name == cleaned_files_list$X1[i])
  nostril_select <- subset(file_selct, file_selct$`loopy$name` == cleaned_files_list$X3[i])
  frame_slect <- nostril_select[which(nostril_select$Slice >= cleaned_files_list$fs[i] & nostril_select$Slice <= cleaned_files_list$ef[i]) ,]
  frame_slect$runID = cleaned_files_list[i, "runID"]
  long_tracks_filt_frames <- rbind(long_tracks_filt_frames, frame_slect)
 
}

long_tracks_filt_frames$seqID <- paste (long_tracks_filt_frames$video_file_name, long_tracks_filt_frames$runID)
#######sync audio lables with video frames#######


#load data file
call_times <- read.csv("thermal_audio_raw_timing.csv", stringsAsFactors = FALSE)
#separate file name column to DATE, GROUP and TRACK 
call_times$file <- substr(call_times$filename, 57, nchar(call_times$filename)-4)
call_times<- cbind(call_times, str_split_fixed(call_times$file, "_", 2))
call_times<- cbind(call_times, str_split_fixed(call_times$`2` , "\\\\", 2))
#call_times <- cbind(call_times, filenames)

#clean and rename columns
call_times_clean <- call_times[ , c(2, 3,4,5,12,14,16,17)]
colnames(call_times_clean)[6] <- "date"
colnames(call_times_clean)[7] <- "group"
colnames(call_times_clean)[8] <- "track"
call_times_clean$track <- gsub('L',"", call_times_clean$track)


#load records log file
logfile <- read.csv("thermal_logfile.csv", stringsAsFactors = F)
logfile <- logfile[, -c(11:14)]

#format data column
logfile$Date <- gsub('/',"", logfile$Date)
logfile$Date <- paste(substr(logfile$Date, 1 ,4), substr(logfile$Date , 7, 8), sep = "")

#merge call times with ID and video data 
merged_data <- merge(call_times_clean, logfile, by.x = c("date", "group", "track"),
                     by.y = c("Date", "Group", "Audio.file"), all.x = T)

#calculate sync tume with video
merged_data$frames <- as.numeric(merged_data$frames)
merged_data$frames[is.na(merged_data$frames)] <- 0
merged_data$syncT <- (merged_data$sec*30) 

#convert time from seconds into frames
merged_data$StartT <- merged_data$start.time*30 + merged_data$syncT + merged_data$frames
merged_data$EndT <- merged_data$end.time*30 + merged_data$syncT+ merged_data$frames
merged_data$AV.Sync.time.s. <-NULL
merged_data$Audio.labled <- NULL

#making all call lables lowercase
merged_data$label <- tolower(merged_data$label)


#### marking frames with calls ##### 
#list relevant video files
video_file_list <- as.numeric(unique(long_tracks_filt_frames$video_file_num))

#create empty frame to fill with data
AVsummary <- data.frame()

for (x in video_file_list)
{
  #print video file number
  print(x)
  #choose one audio file
  one_video_calls <-  subset(merged_data, merged_data$Video.file == x)
  #if no matching audio file go to next
  if( is.na(match(x, long_tracks_filt_frames$video_file_num)) == T)
  {next}
  #choose matching video filr
  one_video_colours <- subset(long_tracks_filt_frames, long_tracks_filt_frames$video_file_num == x)
  
  for (i in 1: nrow(one_video_calls))
  {
    #find matching frames for call times
    frame <- which(one_video_colours$Slice > one_video_calls$StartT[i] & one_video_colours$Slice < one_video_calls$EndT[i])
    #fill call type in corresponding frames
    one_video_colours [frame ,"call_type"] <- one_video_calls[i , "label"]
  }
  #collect all data in one big file
  AVsummary <- rbind(AVsummary,one_video_colours )
}

#AVsummary$seqID <- paste(AVsummary$video_file_name, AVsummary$seg, sep = "")

#add call indicatoe column and mark all existing focal call as 1

AVsummary$Calls <- 0
AVsummary[which(! is.na(AVsummary$call_type) & AVsummary$call_type != "nf" &   AVsummary$call_type != "nf sn"), "Calls"] <- 1
AVsummary[which(AVsummary$call_type == "nf") , "Calls"] <- 1/2  
AVsummary[which(AVsummary$call_type == "nf sn") , "Calls"] <- 1/2   


#plot curves against the call times and visually check for audio/video sync accuracy
upd_file_lst <- unique(AVsummary$seqID)

for (i in upd_file_lst)
{
  segment <- subset(AVsummary, AVsummary$seqID == i)
  plot(segment$time , segment$Median,type='l', lwd=1, col='red' , main = i)
  abline (v= segment[which(segment$Calls == 1), "time"], col = "blue")
  abline (v= segment[which(segment$Calls == 1/2), "time"], col = "green")
  
}

######digital filter data smoothing######


#library
library(lubridate)
library(readxl)
library(dplyr)
library(robfilter)
#creating empty matrix for end output
summ <- matrix(NA , 0, 28)
#making a vector for breathing rate check 
Breathrate_tot <- data.frame(stringsAsFactors = F)


for (k in 1:length(upd_file_lst)) 
{
 
  d <- subset(AVsummary, AVsummary$seqID == upd_file_lst[k])
  #check for segment lenght. Short segments are ignored
  if (nrow(d) < 80) {next}
 
   #filling NA with median values since some of the stages do not accept NAs
  CM <- mean(d[, "filt_median"], na.rm =T)
  for (b in 1:nrow(d))
  {
    
    d[b , "filt_median"]= ifelse(is.na(d[b, "filt_median"]), CM , d[b , "filt_median"])
  }
  
  #smoothing of raw data
  
  
   fs<-1/30 # camera frame rate is 30 frames/second.  fs=1/framerate
   d$time<-d$Slice/30 # create a time variable 
    
   #center temp values for avpiding filter edge effect
   d$m.center <- d$filt_median - mean(d$filt_median)
  
    library(signal) 
   bf <- butter(2, W=1/5, type="low") # low pass filter.  
   # The syntax is weird and poorly explained, but if you need to get rid of fast responses, 
   # set the W parameter appropriately.  
   
   filteredsignal <- filtfilt(bf, x= d$m.center)
   
   #adding filtered collumn to data and adding the mean of the original measurments
   d$filteredsignal<-  filteredsignal
   d$filt_med <- d$filteredsignal + mean(d$filt_median)
  
   #plot raw and processed data
   par(mfrow=c(2,1))
   plot(d$time, d$filt_med , type="l", lwd=2, main = paste(d$`loopy$name`[1],min(d$Slice),"-",max(d$Slice),"butterworth"))
   lines(d$filt_median~d$time, col="red")
 
   
  #### finding mins and maxs of the breathig curve ######
  
  dyad.dist <- d$filt_med #called this dyad.dist to match function
  
  plot.results <- F
  noise.thresh <- 0.25 #this is the starting max treshold value
 
  n.times <- length(dyad.dist)
 # setting a repeat counter for getting out of endless loops
   int <- 1
   
   
  repeat {
  
  
  #---find first minimum---
  #get first instance of a change in dyadic distance above noise.thresh, and determine whether it is going up or down
  i <- 1
  curr.min <- 1
  found = F
  while(found == F & i < length(dyad.dist)){
    dist.change <- dyad.dist[i] - dyad.dist[1]
    if(dyad.dist[i]<dyad.dist[curr.min]){
      curr.min <- i
    }
    if(abs(dist.change)>noise.thresh){
      found = T
      if(dist.change > 0){
        up.first <- T
      }
      else{
        up.first <- F
      }
    }
    else{
      i <- i + 1
    }
  }
  
  #if no first minimum found, return NULL
  if(!found){
    print('no first minimum found')
  }
  
  #find starting point (if it went up first, this is curr.min, otherwise, this is the first minimum)
  if(up.first){
    first <- curr.min
  } else{
    curr.min <- 1
    found <- F
    while(i < n.times & !found){
      dist.diff <- dyad.dist[i] - dyad.dist[curr.min]
      if(dist.diff < 0){
        curr.min <- i
      }
      else{
        if(dist.diff > noise.thresh){
          found <- T
          first <- curr.min
        }
      }
      i<-i+1
    }
    if(!found){ #if no starting point found, return NULL
      return(NULL)
    }
  }
  
  #----pull out maxes and mins in dyadic distance----
  min.max.min <- data.frame(t1 = NA, t2 = NA, t3 = NA)
  min.max.min$t1[1] <- first
  ref.idx <- first
  ref.val <- dyad.dist[first]
  data.idx <- 1
  going.up <- T
  for(i in first:n.times){
    curr.idx <- i
    curr.val <- dyad.dist[i]
    if(going.up){ #if going up
      if(curr.val >= ref.val){ #if the current value is greater than the reference value
        ref.idx <- curr.idx #replace the reference value with the current value
        ref.val <- dyad.dist[ref.idx]
      }
      if(abs(curr.val - ref.val) > noise.thresh){ #if the difference between current and reference exceeds the noise threshold
        min.max.min$t2[data.idx] <- ref.idx #record the local max
        going.up <- F #switch to going down
        ref.val <- curr.val
        ref.idx <- curr.idx
      } 
    }
    else{ #if going down
      if(curr.val <= ref.val){ #if the current value is less than the reference value
        ref.idx <- curr.idx #replace the reference value with the current value
        ref.val <- dyad.dist[ref.idx]
      }
      if(abs(curr.val - ref.val) > noise.thresh){ #if the difference between current and reference exceeds a noise threshold
        min.max.min$t3[data.idx] <- ref.idx #record the local min as the end of the current sequence
        data.idx <- data.idx + 1 #increment the data index
        min.max.min <- rbind(min.max.min,c(ref.idx,NA,NA)) #record the local min as the start of the next sequence
        going.up <- T #switch to going up
        ref.val <- curr.val
        ref.idx <- curr.idx
      }	
    }
  }
  
  
  #remove the last row if it is incomplete to get full sucles
  min.max.min.comp <- min.max.min[which(!is.na(min.max.min$t3)),]
  
  #calculate breathing rate based on the number of max peaks
  rate <-(d[max(min.max.min.comp$t3, na.rm = T) , "time"] - d[min(min.max.min.comp$t1, na.rm = T) , "time"])/(length(which(!is.na(min.max.min.comp[,2]))) - 1)
 
   #print the current peak detection threshold and the estimated breathing rate
  print(noise.thresh)
  print(rate)
  
  # conditions to end the repeat loop. 
  if (noise.thresh < 0.04) {break} #if noise.thresh gets bellow 0.04
  if (int == 51) {break}           # if loop goes over 50 iterations
  if (rate < 1.2)  {noise.thresh <- noise.thresh + 0.01} #if breathing rate is too high add to detection threshold
    else if (rate > 1.9) { noise.thresh <- noise.thresh - 0.01 ; int <- int+1} #if breathing rate is too low lower detection and add 1 to iterations
     
      else {break}
  }
  
   #secondary filtering of min peaks
   
   
   min_temp <- d[min.max.min$t1, "filt_med"]
   min_temp_1 <- d[min.max.min$t1, "Median"]
   
   f <- mean( min_temp_1)
   outliers <- which(min_temp > f+0.5)
   ifelse (length(outliers) > 0, new_peaks <- min.max.min[-which(min_temp > f+0.5) ,] , new_peaks <- min.max.min)
   
  
   #select max peaks
  makemax <-  new_peaks[,2] 
  
  #select min peaks
  makemin <-  new_peaks[ ,1]
  
 
  
  ####find plateau region for each peak####
  
  
  #pull min peak location
 peaks <- sort(unique(unlist(makemin)))
  
##here some manual correction conditions are listed. Should be commented out in the first run and edited after visual inspection of the 
 # peaks. At some point this could be made into an interqctive shinny.app interface. 
  if (k == 3) 
  { peaks <- peaks[-6]}
  else if (k == 5)
  {peaks <- peaks [-1]}
  else if (k== 6)
  {peaks <- peaks [-2]}
  else if (k== 7) 
  {peaks <- peaks [-5]}
  else if (k== 11) 
  {peaks <- peaks [-5]}
  else if (k == 12)
  {peaks <- peaks [-2]}
  else if (k== 13)
  {peaks <- peaks [-c(2,5)]}
  else if (k == 14)
  {peaks <- peaks [-1]}
  else if (k== 15)
  {peaks <- peaks [-4]}
  else if (k == 16)
  {peaks <- peaks [-4]}
  else if (k== 22)
  {peaks <- peaks [-4]}
  else if (k== 23)
  {peaks <- peaks [-5]}
  else if (k== 24)
  {peaks <- sort(c(peaks, 445))
  makemax <- c(makemax, 429, 477)}
  else if (k== 29)
  {peaks <- sort(c(peaks, 50))
  makemax <- c(makemax, 40 )}
  else if (k == 31)
  {peaks <- peaks [-2]}
  else if (k == 37)
  {peaks <- peaks [-c(1,2)]}
  else if (k== 38)
  {peaks <- peaks [-3]}
  else if (k== 45)
  {peaks <- peaks [-c(2,4)]}
  else if (k== 46)
  {peaks <- peaks [-1]
  peaks <- sort(c(peaks, 482))
  makemax <- c(makemax, 505)}
  else if (k ==52)
  {peaks <- c(40, 65, 110)
    makemax <- c(10, 75, 90)}
 ## end of manual peak correction bit 
  
  #make a vector of all peaks
  allpeaks <- sort(c(makemax, peaks))
  
  
  #make separate vectors for peaks from which to go forward or backward. 
  #this is done to prevent errors from peaks on the edges of the segment 
 for_peaks <-peaks[peaks <  (length(dyad.dist) - 10)]          
 rev_peaks <- peaks[peaks > 10]                           
                            
 #make empty vectors for plateau points
 plato.end <- numeric(0)
 plato.start <- numeric(0)
 

  #finding plateau start
 for(j in 1:length (for_peaks)) {
   #set the current minimum peak as a starting point and the next max peak(- 2 frames) as the end point
   
   slopes <- data.frame()
   local_stop <-  allpeaks[(which(allpeaks == for_peaks[j])+1)] - 2
   #if there is no next peak set the end of the segment (- 5 frames) as the end point 
   if (length(local_stop) == 0 || is.na(local_stop)) {local_stop = nrow(d)-5 }
  
   #start walking from start to end point making 5 frame regressions and extracting slope values
    for (i in for_peaks[j]:local_stop)
     
   {
     
     slopes[i , 1] <- lm(d$filt_med[i:(i+5)] ~ d$Slice[i:(i+5)])$coefficients[2]
     slopes[ i, 2] <- i
   }
   
   
   #set new start at the steepest slope point
   new_start <- which(slopes$V1 == max(slopes$V1, na.rm = T) )
   
   #from the new start point until the same local stop  make regressions on a 2 frame segments
   secondary_slopes <- data.frame()
   
   for (i in new_start:local_stop)
   {
     secondary_slopes[i,1] <- lm(d$filt_med[i:(i+2)] ~ d$Slice[i:(i+2)])$coefficients[2]
     secondary_slopes[i,2] <- i
   }
   
   
   #find a point where the slope is 1/3 of the steepest one or closest to zero (whatever comes first)
   plato_start <- which(abs(secondary_slopes$V1) == min(abs(secondary_slopes$V1), na.rm = T)  )
   plato_start_alt <- min(which(abs(secondary_slopes$V1) < max(abs(slopes$V1), na.rm = T)/3))
   start_point <- min(plato_start, plato_start_alt, na.rm = T)
   plato.start <- c(plato.start, start_point)
  
   
 }
 
 # finding plateau end
 #do the same but in the oposite direction to find plateau ends
 for(j in 1:length (rev_peaks)) {
   
   slopes <- data.frame()
   local_stop <-  allpeaks[(which(allpeaks == rev_peaks[j])-1)] +2 
   if (length(local_stop) == 0 || is.na(local_stop)) {local_stop =  5}
   for (i in rev_peaks[j]:local_stop)
     
     
     
   {
     
     slopes[i , 1] <- lm(d$filt_med[i:(i-5)] ~ d$Slice[i:(i-5)])$coefficients[2]
     slopes[ i, 2] <- i
   }
   
   
   new_start <- which(slopes$V1 == min(slopes$V1, na.rm = T) )
   
   
   secondary_slopes <- data.frame()
   
   for (i in new_start:local_stop)
   {
     secondary_slopes[i,1] <- lm(d$filt_med[i:(i-2)] ~ d$Slice[i:(i-2)])$coefficients[2]
     secondary_slopes[i,2] <- i
   }
   
   
   
   plato_end <- which(abs(secondary_slopes$V1) == min(abs(secondary_slopes$V1), na.rm = T)  )
   plato_end_alt <- max(which(abs(secondary_slopes$V1) < max(abs(slopes$V1), na.rm = T)/3))
   end_point <- max(plato_end, plato_end_alt, na.rm = T)
   plato.end <- c(plato.end, end_point)
   
   
 }
 
 if (k == 7) 
 { plato.end[1] <- 40}
 else if (k == 4) 
 { plato.start <- c(35, 105 ,141, 183, 245)}
 else if (k == 6)
 { plato.end[2] <- 110}
 else if (k==10)
 {plato.end[3] <- 85}
 else if (k==25)
 {plato.start[4] <- 89}
 else if (k==45)
 {plato.end <- plato.end[-1]}
 else if (k==46)
 {plato.start[8] <- 424}
 else if (k==47)
 {plato.start[3] <- 70}
  else if(k==52)
  {plato.end[2] <- 57 
    plato.start <- c( 7, 46 ,72)}
 
 #plot peaks on the filtered signal
 plot(dyad.dist~d$Slice,type='l', main = upd_file_lst[k]) #filtered signal
 
 abline (v = d$Slice[peaks], col = "red" ) #minimum
 abline(v = d$Slice[ plato.start], col = "gold") #start points
 abline(v = d$Slice[plato.end], col = "blue")    # end points
 
  #making empty matrices for populating with peak location and types
  t1n <- matrix(1, length(plato.end), 1)
  t3n <- matrix(3, length(plato.start), 1)
  #filling plateau start and end points
  t1 <- cbind(plato.end, t1n) 
  t3 <- cbind(plato.start, t3n) 
  
  
  # creating separate peaks column with all three transition type peak  locations
  
  plato.peaks <- matrix(NA, nrow(d), 2)
  #popplating column with t1 and t3 peaks
  for (i in 1:nrow(t1))
  {
    for (j in 1:2) {
      plato.peaks[t1[i], j] <- paste0(t1[i,j])
    }
  }
  
  for (i in 1:nrow(t3))
  {
    for (j in 1:2) { 
      plato.peaks[t3[i],j] <- paste0(t3[i,j])
    }
  }
  #adding  minimum (t2) peaks  to peaks column
  min.peaks <- peaks
  t2n <- matrix(2, length(min.peaks))
  t2 <- cbind(min.peaks, t2n) 
  for (i in 1:length(min.peaks))
  {
    {
      for (j in 1:2) {
        plato.peaks[t2[i], j] <- paste0(t2[i,j])
        
      }
    }
    
    
  }
  
  fpeaks <- as.numeric(na.omit(plato.peaks[ ,1]))
   
  
  ##making separate columns for change in intensity, frame distance and slope between peaks
  
  Intensity <- matrix(NA, nrow(d), 1)
  Distance <- matrix(NA, nrow(d), 1)
  Slope <- matrix(NA, nrow(d), 1)
  
  ##calculating all three peak parameters
  for (i in 1:length(fpeaks))
  {
    
    
    Intensity[fpeaks[i],1] <- abs(d[fpeaks[i], "Median"] - d[fpeaks[i+1], "Median"]) #intensity is calculated from the real data to avoid smoothing
    Distance[fpeaks[i], 1] <- abs(d[fpeaks[i], "time"] - d[fpeaks[i+1], "time"]) 
  }
  
  for (i in 1:(length(fpeaks)-1))
  {
    xs <- d[fpeaks[i]:fpeaks[i+1],"time"]
    ys<- d[fpeaks[i]:fpeaks[i+1], "Median"]
    
    
    Slope[fpeaks[i],1] <- lm(ys~xs)$coefficients[2]
  }
  
  #adding new columns to data
  
  d<- cbind(d, plato.peaks)
  colnames(d)[colnames(d)==1] <- "peak"
  colnames(d)[colnames(d)==2] <- "phase"
  d<- cbind(d, Intensity , Distance, Slope) 
  
  
  #making new columns to fill later
  addcolumns <- matrix(NA, nrow= nrow(d), ncol = 3)
  colnames(addcolumns) <- c("type", "NF", "Call.duration")
  d <- cbind(d, addcolumns)
  
  #clearing previous breaths, NF naming and call count (can be removed if input files are raw)
  d$NF = ""
  d$type = ""
  d$Call.duration = ""
  
  
  #make column of peak idx
  ll <- length(fpeaks)
  
  #adding the needed factor levels to breaths type column
  levels(d$type) <- c("c", "ac" , "bc" , "n" , "")
  
  #adding breathidx collumn
  d$breathidx <- NA
  #making phase column numeric
  d$phase <- as.numeric(d$phase)
  d$Calls <- as.numeric(d$Calls)
  d$peak <- as.numeric(as.character(d$peak))
  #setting the breaths counter
  r<-1
  
  # setting the first peak as breath #1
  d[fpeaks[1], "breathidx"] <- r
  
  #counting the breaths through the file. As the phase resets new count is added
  for (z in 2:(ll)) {
    if (d[fpeaks[z], "phase"] > d[fpeaks[z-1], "phase"]) {
      d[fpeaks[z], "breathidx"] <- r
    } else {
      r <- r+1
      d[fpeaks[z], "breathidx"] <- r
    }
    
  }
  
  
  #adding 0s insteadof NA to Calls collumn. Some later step fails wihout it
  d$Calls[is.na(d$Calls)] = 0
  
  #detecting and marking phases with calls in them
  for (i in 1:(ll-1)) {
    x <-  any(d[fpeaks[i]:fpeaks[i+1], "Calls"] == 1, na.rm = F)
    if  (x == T)
    {d[fpeaks[i], "type"] <- "c"
    
    } 
  }
  
  #detecting phases with non focal calls in them
  for (j in 1:(ll-1)) {
    y <-  any(d[fpeaks[j]:fpeaks[j+1], "Calls"] == 1/2, na.rm = F)
    if  (y == T)
    { d[fpeaks[j] , "NF"] <- "y"
    d[fpeaks[j+1] , "NF"] <- "y"
    
    }
  }
  
  #filling the rest of NF peaks as "n"
  for (peaknum in fpeaks) {
    if (d[peaknum ,"NF"] != "y") {
      
      d$NF[peaknum] = "n"
    }
    
  }
  IDX <- 0
  #finding and naming all peaks of "call" breaths
  for (peaknum in fpeaks) {
    if (d[peaknum, "type"] == "c") {
      IDX <- d$breathidx[peaknum]
      d[which(d$breathidx == IDX), "type"] = "c"
      
      #counting call duration in frames
      #making a table of peak index and breath counter of breaths with calls
      call.breaths <- (d[which(d$type == "c"), c("peak", "breathidx")])
     
       #summing up call duration in frames for each breath with calls
      for (breath in IDX) {
        breath.range <- range(call.breaths[which(call.breaths$breathidx == breath), 1])
       
         #removing NonFocal calls from the total call duration
        d2 = d
        d2[which(d2$Calls == 0.5), "Calls"] = 0
        d[which(d$breathidx == breath), "Call.duration"] = sum(d2[breath.range[1]:breath.range[2], "Calls"]) 
      }
    }
  }
  
  
  if (IDX != 0){
    
    #making a table of peak index and breath counter of breaths with calls
    call.breaths <- (d[which(d$type == "c"), c("peak", "breathidx")])
    
    
    #finding and naming all peaks of "after call" breaths
    for (peak1 in fpeaks) {
      if (d$type[peak1] == "c") {
        IDXAC <- d$breathidx[peak1]+1 
        d[which(d$breathidx == IDXAC), "type"] = "ac"
      }
    }
    #finding and naming all peaks of "before call" breaths
    for (peak2 in fpeaks) {
      if (d$type[peak2] == "c") {
        IDXBC <- d$breathidx[peak2]-1
        d[which(d$breathidx == IDXBC), "type"] = "bc"
      }
      
    }
    
  }
  #filling the rest of the unassigned peaks with "n"
  for (peaknum in fpeaks) {
    if (d$type[peaknum] == "") {
      
      d$type[peaknum] = "n"
    }
    #selecting rows with actual peak values
    filt <- d[d$type !="",]
    
    #Calculating mean breathing rate for benchmark of peak detection
    BRate_file <-data.frame(  (filt[max(which(filt$phase == 3)) , "time"] -   filt[min(which(filt$phase == 3)) , "time"])/(length( which(filt$phase == 3)) -1), filt$seqID[1])
   
    
  }
  
  
  #filling up data summary file
  summ <- rbind(summ, filt)
  #filling us breathing rates
  Breathrate_tot <- rbind(Breathrate_tot, BRate_file)
  
}

colnames(Breathrate_tot) <- c("rate", "file")

#Getting column names for the summary files 
columns <- colnames(filt)
colnames(summ) <- columns


print(mean(Breathrate_tot[ , 1]))
Breathrate_tot[ , 1]
#write files
write.csv(summ, "peaksummary4stats.csv"  )


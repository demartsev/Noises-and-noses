##### This is a script for summarizing, analyzing and visualizing the respiration results
##### it uses the full anonymised data-set and NOT the sample output from the minimal 
##### data processing script


library(DescTools)
library(ggplot2)
library(permute)
library(coin)
library(tidyverse)
library(fitdistrplus)
library(lmtest)
library(lme4)
library(jtools)
library(multcomp)
library(dplyr)

## read in the data csv ##
all_peaks <- data.frame(read.csv("/cloud/project/data_outputs/data_summary.csv", stringsAsFactors = F , sep = ";"))


# create function to calculate mean respiration phase proportions (and remove any NA values present)
meannona<-function(x) return(mean(x, na.rm=T))

breathsummary <- all_peaks %>% dplyr::select(type, NF, Distance, phase) %>%
  dplyr::filter(type=="n", NF=="y") %>%
  group_by(type, phase) %>%
  summarise_all(meannona) 
breathsummary
ti<-breathsummary$Distance[1]
te<-breathsummary$Distance[2]
tepause<-sum(breathsummary$Distance)
ti_te<-ti/te
ti_te  #inhalation/exhalation proportion
ti_ttot<-ti/(ti+te+tepause)
ti_ttot


#get  unique individual identifier
all_peaks$uniqueID <- paste(all_peaks$group, all_peaks$ID, sep = "")

#reorder factors for a better visualization
all_peaks$type <- factor(all_peaks$type, c("bc", "c", "ac", "n"))

#add new column for better Non Focal call separation
all_peaks$NF_ord <- NA
seqID <- unique((all_peaks$seqID))
all_peaks_n <- data.frame()

#mark all respiration phases that follow Non Focal calls
for(i in seqID){
  
  test <- all_peaks[which(all_peaks$seqID == i),]
  
  for (i in 1:(nrow(test) - 1))
  { ifelse(test[i,"NF"] == "y", test[ i+1, "NF_ord"] <- 1, test[i+1, "NF_ord"] <- 0)}
  all_peaks_n <- rbind(all_peaks_n, test)
}

#mark phases which have or directly follow Non Focal calls
all_peaks_n$NF_main <- 0
rows <- c(which(all_peaks_n$NF == "y") , which(all_peaks_n$NF_ord == 1 ))
all_peaks_n[rows, "NF_main"] <- 1

#rename the data frame for compatibility and to avoid mix-up with subsequent analysis
all_peaks <- all_peaks_n

#select only individuals that produced calls to allow permutations
table(all_peaks$uniqueID, all_peaks$type) #check the table and select individuals who called
callers <-  c(     115,            12,     17 ,            23 ,           25,            310,             311,           313,             314,           38 )

#filter data from caller only
all_peaks_callers <- all_peaks[which(all_peaks$uniqueID %in% callers) ,]

all_peaks_callers$uniqueID <- as.character(all_peaks_callers$uniqueID)

#add 0.01 to avoid zero values
all_peaks_callers$Intensity_1 <- all_peaks_callers$Intensity +0.01

#separate the file to breathing phases for easier handling
phase_1 <- all_peaks_callers[which(all_peaks_callers$phase == 1) ,] #inspiration
phase_2 <- all_peaks_callers[which(all_peaks_callers$phase == 2) ,] #expiration
phase_3 <- all_peaks_callers[which(all_peaks_callers$phase == 3) ,] #expiratory pause


resp_phases = c("inspiration", "expiration",  "expiratory pause")

###permutation tests are done here #####

perm_rounds <- 10000 #set the number of iterations here!!!

#run permutation tests with individual ID as a blocking parameter. 
for (type in c("bc", "c", "ac")) { #setting the types to compare with Non-call
  
  #function to subtract means   
  meanDif <- function(x, grp) {
    mean(x[grp == type], na.rm = T) - mean(x[grp == "n"], na.rm = T)
  }
  #make a vector of tested metrics
  for (metric in c("Intensity", "Distance", "Slope"))
    
  {
    par(mfrow=c(1,3))
    
    for (j in 1:3) {
      DIntensity <- numeric(length = perm_rounds) #make object for collecting the permuted values
      N <- nrow(all_peaks_callers[which(all_peaks_callers$phase == j),]) #number of samples
      
      set.seed(42)
      
      ##set the ID as the blocking variable. Permutations are only done within individuals 
      CTRL <- how(blocks  = all_peaks_callers[which(all_peaks_callers$phase == j), "uniqueID"]) 
      
      #shuffle the type values against each one of the metrics in all 3 phases
      
      for(i in seq_len(length(DIntensity) - 1)) {
        perm <- shuffle(N,  control = CTRL)
        h <-  all_peaks_callers[which(all_peaks_callers$phase == j), metric]
        DIntensity[i] <- with(all_peaks_callers[which(all_peaks_callers$phase == j) ,], meanDif(h, type[perm]))
      }
      #calculate the mean difference from real data
      DIntensity[perm_rounds] <- with(all_peaks_callers[which(all_peaks_callers$phase == j) ,], meanDif(h, type))
      
      #calculating pseudo-P value
      Dbig <- sum(abs(DIntensity) >= abs(DIntensity[perm_rounds]))
      P <- Dbig / length(DIntensity)
      
      #plotting hist
      hist(DIntensity, main = paste("P =", P, "Phase-",resp_phases[j]),
           xlab = paste("mean", metric ,"difference (", type, "- nc)" ))
      rug(DIntensity[perm_rounds], col = "red", lwd = 2) #marking the real data value
    }
    
  }
}


#make box_plots to summarize the results of calling phases
par(mfrow=c(1,3), mai = c(0.1, 0.1, 0.1, 0.1), oma=c(5,5,5,5))

#slope box plots
boxplot(Slope~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "bc" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab=" ", 
        xlab = "",
        xaxt = "n")
legend("topright",c("Pre_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

title("Respiration curve slope", line = 1, outer = TRUE, cex=3, font=2)
mtext("Respiration phase", line = -50, outer = T, cex = 1, font = 1)

boxplot(Slope~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "c" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab = "",
        xlab="Breathing phase",
        xaxt = "n", yaxt="n")
legend("topright",c("Call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

boxplot(Slope~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "c" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        xlab = "", ylab = "",
        xaxt = "n",  yaxt="n")
legend("topright",c("Post_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))


#amplitude box plots
boxplot(Intensity~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "bc" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab="Amplitude of change", 
        xlab="" ,
        xaxt = "n")
legend("topright",c("Pre_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

title("Respiration phase amplitude", line = 1, outer = TRUE, cex=3, font=2)
mtext("Respiration phase", line = -50, outer = T, cex = 1, font = 1)

boxplot(Intensity~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "c" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab=" ", 
        xlab=" " ,
        xaxt = "n", yaxt="n")
legend("topright",c("Call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

boxplot(Intensity~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "ac" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab=" ", 
        xlab=" " ,
        xaxt = "n", yaxt="n")
legend("topright",c("Post_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))


#duration box plots
boxplot(Distance~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "bc" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab=" ", 
        xlab=" " ,
        xaxt = "n")
legend("topright",c("Pre_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

title("Respiration phase duration", line = 1, outer = TRUE, cex=3, font=2)
mtext("Respiration phase", line = -50, outer = T, cex = 1, font = 1)

boxplot(Distance~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "c" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab=" ", 
        xlab=" " ,
        xaxt = "n", yaxt="n")
legend("topright",c("Call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

boxplot(Distance~type*phase, data=droplevels(all_peaks_callers[which(all_peaks_callers$type %in% c("n", "ac" )) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab="Phase duration", 
        xlab=" " ,
        xaxt = "n", yaxt="n")
legend("topright",c("Post_call","Quiet"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))


####looking at the non_call data subset#####
all_peaks_n <- all_peaks_n[which(all_peaks$type == "n"), ]

#run permutation tests with individual ID as a blocking parameter. 

#function to subtract means   
meanDif <- function(x, grp) {
  mean(x[grp == 1], na.rm = T) - mean(x[grp == 0], na.rm = T)
}
#make a vector of tested metrics
for (metric in c("Intensity", "Distance", "Slope"))
  
{
  par(mfrow=c(1,3))
  
  for (j in 1:3) {
    DIntensity <- numeric(length = perm_rounds) #make object for collecting the permuted values
    N <- nrow(all_peaks_n[which(all_peaks_n$phase == j),]) #number of samples
    
    set.seed(42)
    ##set the ID as the blocking variable. Permutations are only done within individuals 
    CTRL <- how(blocks  = all_peaks_n[which(all_peaks_n$phase == j), "uniqueID"]) 
    
    #shuffle the type values against each one of the metrics in all 3 phases
    for(i in seq_len(length(DIntensity) - 1)) {
      perm <- shuffle(N,  control = CTRL)
      h <-  all_peaks_n[which(all_peaks_n$phase == j), metric]
      DIntensity[i] <- with(all_peaks_n[which(all_peaks_n$phase == j) ,], meanDif(h, NF_main[perm]))
    }
    #calculate the mean difference from real data
    DIntensity[perm_rounds] <- with(all_peaks_n[which(all_peaks_n$phase == j) ,], meanDif(h, NF_main))
    
    #calculating P value
    Dbig <- sum(abs(DIntensity) >= abs(DIntensity[perm_rounds]))
    P <- Dbig / length(DIntensity)
    
    #plotting hist
    
    hist(DIntensity, main = paste("P =", P, "Phase -",resp_phases[j]),
         xlab = paste("mean", metric ,"difference (NF calls - No_Calls)" ))
    rug(DIntensity[perm_rounds], col = "red", lwd = 2)
  }
  
}



#make box_plots to summarize the results of non focal effect on respiration phases
par(mfrow=c(1,3), mai = c(0.5, 0.2, 0.5, 0.1), oma=c(2,2,2,2))

boxplot(Intensity~NF_main*phase, data=droplevels(all_peaks_n[which(all_peaks_n$type %in% c("n")) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab="Phase duration", 
        main="Respiration phase amplitude", xlab="" ,
        xaxt = "n")
legend("topright",c("No calls","Non focal calls"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

boxplot(Distance~NF_main*phase, data=droplevels(all_peaks_n[which(all_peaks_n$type %in% c("n")) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab="Phase duration", 
        main="Respiration phase duration", xlab="Respiration phase" ,
        xaxt = "n")
legend("topright",c("No calls","Non focal calls"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))

boxplot(Slope~NF_main*phase, data=droplevels(all_peaks_n[which(all_peaks_n$type %in% c("n")) , ]), notch=T,
        col=(c("gray","floralwhite")),
        ylab="Phase duration", 
        main="Respiration curve slope", xlab="" ,
        xaxt = "n")
legend("topright",c("No calls","Non focal calls"), fill=c("gray","floralwhite"))
axis(1, at = c(1.5, 3.5, 5.5) ,labels = c("inspiration", "expiration",  "expiratory pause"   ))





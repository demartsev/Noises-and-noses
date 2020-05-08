library(stringr)

#set workinf directiry
setwd("D:/VideoLibrary/Loopy_output")

#list all CSV files in WD
file_list <- list.files(pattern = "csv$")

for (i  in 1:length(file_list)) {
  

#read one file
loopy_output <- read.csv(file_list[i], stringsAsFactors = F)
#set the name to match video file
Loopy_name <- file_list[i]
Loopy_name <- str_sub(Loopy_name, start = 47, end = - 25)
Loopy_name <- sub("_[^_]+$", "", Loopy_name)

#subset relevant columns
loopy_output1 <- loopy_output[ , c("frame_number","name","x" ,"y")]

#separateing left nostril

left <- subset(loopy_output, loopy_output$name == "nostrilL")
left <- left[ , c(1,3:4)]
colnames(left)[1] <- "f"
left$f <- left$f+1
write.csv(left, paste(Loopy_name,"_L.csv",sep = "" ))
#separateing left nostril
right <- subset(loopy_output, loopy_output$name == "nostrilR")
right <- right[ , c(1,3:4)]
colnames(right)[1] <- "f"
right$f <- right$f+1
write.csv(right, paste(Loopy_name,"_R.csv",sep = "" ))

both <- rbind(right, left)
write.csv(both, paste(Loopy_name,".csv",sep = "" ))
}

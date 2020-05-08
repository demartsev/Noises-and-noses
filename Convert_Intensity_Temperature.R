
# How to convert the intensity (colour) value exported from ImageJ macro
# analysis of FLIR csq -> avi conversion

# Need to obtain the calibration constants from the csq files
# I used Exiftool (a command line utility to do this), but have also
# included this functionality in Thermimage (package for R)

PR1<-21546.2
PB<-1507
PF<-1
PO<-(-6331)
PR2<-0.01622949

# Need to set the object parameters required for the raw -> temperature
# conversion
# Emissivity, Object Distance, Reflected Temperature, Atmospheric Temperature,
# IR Window Temperature, IRT (these 2 are technically not required since
# we did not have a window, but are included for clarity)

E<-0.95
OD<-1
RTemp<-10
ATemp<-10
IRWTemp<-10
IRT<-1
RH<-30

library(Thermimage) # access to the raw2temp algorithm

rawvalues<-7000:65535
tempvalues<-raw2temp(rawvalues, E, OD, RTemp, ATemp, IRWTemp=RTemp, IRT=1, RH, 
         PR1, PB, PF, PO, PR2)
# This plot is for the full range of the camera. 
plot(tempvalues ~ rawvalues, type="l")

# We know that we only report on temperatures above 0C and below 40C, 
# so let's truncate the raw values to a much more reasonable range:
rawvalues<-10000:20000
tempvalues<-raw2temp(rawvalues, E, OD, RTemp, ATemp, IRWTemp=RTemp, IRT=1, RH, 
                     PR1, PB, PF, PO, PR2)
plot(tempvalues ~ rawvalues, type="l")

# add these to a data.frame
d<-data.frame(rawvalues, tempvalues)

lm1<-lm(tempvalues ~ stats::poly(rawvalues, 4), data=d)
plot(tempvalues ~ rawvalues, data=d)
lines(fitted(lm1) ~ rawvalues, col="red")

summary(lm1)
# R squared is 1, so we have a perfect fit with a 4th order polynomial
summary(lm1)$r.squared

# formula coeficients: 
cbind(coef(lm1))

# So this is awkward, but this is the empirical equation we can use:
paste0(coef(lm1)[1], " + ", 
       coef(lm1)[2], "*rawvalues ",
       coef(lm1)[3], "*rawvalues^2 ",
       coef(lm1)[4], "*rawvalues^3 ",
       coef(lm1)[5], "*rawvalues^4 ")

# or the easier approach is to use the predict function
predict(lm1, newdata=data.frame(rawvalues=13000))

# Might be easier to wrap the conversion into a function:
convert_rawtotemp<-function(rawvalues=rawvalues, modelfit=lm1){
  newdata<-data.frame(rawvalues)
  predictedtemperature<-predict(lm1, newdata)
  names(predictedtemperature)<-NULL
  return(predictedtemperature)
}

# then use the function on any of the rawvalues (what imagej calls Intensity)
# obtained from ImageJ

convert_rawtotemp(13000, lm1)
convert_rawtotemp(rawvalues, lm1)

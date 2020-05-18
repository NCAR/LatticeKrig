setwd('~/Dropbox/Home/Projects/BRACEData')
remove( list=ls())
library( ncdf4)
library( fields)

fileNames<- system( 'ls ~/Dropbox/Data/BRACEModelOutput/*85*TS.200601-208012.nc', intern =TRUE)
N<- length( fileNames)

# get common info from the first ensemble ncdf file
handle<- nc_open(fileNames[1]  )

lat <- ncvar_get( handle, 'lat')
lon <- ncvar_get( handle, 'lon')

Longitude.0360<- lon 
Longitude.trans <- ifelse(Longitude.0360<=180, Longitude.0360, Longitude.0360-360 ) 
index.sort <- order(Longitude.trans)
Longitude <- Longitude.trans[index.sort]
Latitude  <- lat
# boulder location 
# 40.0150° N, 105.2705° W
kBoulder<- 139
jBoulder<- 60

nx<- length( lon)
ny<- length( lat)

NTime <- length( date) 
date  <- ncvar_get( handle, 'date')
month<- substr( as.character( date), 5,6)
month<- as.numeric( month)
# assume that the month reported is actually the following month for the mean. 
# e.g. 2 is really January, 12 is November
month <- month - 1
month[ month == 0] <- 12
year<- as.numeric(substr( as.character( date), 1,4))
nt<-  sum( month==7)
year<- year[ month==7]
# Gaussian weights
# gw  <- ncvar_get( handle, 'gw')
#gw<- gw/ sum( gw)
# Gaussian weights for finding global surface average.
#gwField<- matrix( gw, 288, 192, byrow=TRUE )
# gwField<- gwField/ sum( gwField)
# these are Stacey's weights
cos.weights <- cos(lat*pi/180)

#testField<- matrix( cos.weights, 288, 192, byrow=TRUE )
#testField<- testField/sum( testField)
#test.for.zero( testField, gwField )

gwField<- matrix( cos.weights, 288, 192, byrow=TRUE )
gwField<- gwField/ sum( gwField)

seasonalAverage<- function( x){
  filter(x,rep(1/3,3), sides=2)}

globalAverage<- function( field){
  sum(gwField*field)
}

#JJASlope<-array(NA, c( nx, ny, N) )
#JJAGlobal<- matrix(NA, sum( month==7),N )


NYears<- sum( month==7)
tempX<- tempY<- matrix( NA,  NYears,N)
for( member in 1:30){
  cat("member", member, fill=TRUE)
handle<- nc_open(fileNames[member]  )
# convert to degrees C
TS<- ncvar_get( handle, 'TS')- 273.15
# reshuffle to center on prime meridian
TS<- TS[ index.sort,,]
globalTemp<- apply( TS, c(3), globalAverage)
globalTempSeasonal<- seasonalAverage( globalTemp)
tempX[,member]<- globalTempSeasonal[month==7]

tempY[,member]<- (seasonalAverage(TS[jBoulder,kBoulder,]))[month==7]
}

#JJASlope<- lsfit( tempX, tempY)$coef[2]
save( year, 
      tempX, tempY,
      jBoulder, kBoulder,
      file="slopeBoulder.rda" )  

  
  

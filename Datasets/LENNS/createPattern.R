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

JJASlope<-array(NA, c( nx, ny, N) )
JJAGlobal<- matrix(NA, sum( month==7),N )

for( member in 1:N ){
  cat( ' ', fill=TRUE)
  cat( ' Working on  member ', member, fill=TRUE)
  handle<- nc_open(fileNames[member]  )
  TS<- ncvar_get( handle, 'TS')
# convert to degrees C
  TS<- TS - 273.15
  globalTemp<- apply( TS, c(3), globalAverage)
  globalTempSeasonal<- seasonalAverage( globalTemp)
# plot( 1:100,globalTemp[1:100] ); lines(1:100, globalTempSeasonal[1:100] )
  tempX<- globalTempSeasonal[month==7]
  JJAGlobal[,member]<- tempX
  for( j in 1:nx){
    if( j%%50==0){ cat(j, ' ') }
    for( k in 1:ny){
      tempY<- (seasonalAverage(TS[j,k,]))[month==7]
      JJASlope[j,k,member]<- lsfit( tempX, tempY)$coef[2]
   }
}

}
# reshuffle to center on prime meridian 
JJASlope<- JJASlope[ index.sort,,]

JJASlopeMean<- apply(JJASlope, c(1,2), mean )
JJASlopeSd<- apply(JJASlope, c(1,2), sd )
#
#This may seem weird -- include this source code as part of the output object. 
#
sourceCode<- scan('createPattern.R', what='a', sep='\n')
# save everything of value
timeStamp<- date()
save(sourceCode,timeStamp,  JJASlope,
     JJASlopeMean, JJASlopeSd, JJAGlobal, Longitude,
     Latitude,
          file='JJAPatternScalingSlope.rda' )



load("JJAPatternScalingSlope.rda")
# stacey's version
load('~/Dropbox/Home/Projects/HPC4StatsBRACE/data/BRACEUfields.rda')

quartz()
set.panel( 1,3)
image.plot( Longitude,Latitude,JJASlopeMean,zlim=c( 0,5.5))
map( 'world', add=TRUE, col='magenta')
image.plot( Umean, zlim=c( 0,5.5))
image.plot( 100*(Umean- JJASlopeMean)/Umean, zlim=c( -10,10) )
title("percent relative error")
 
quartz()
set.panel(3,1)
image.plot( Longitude,Latitude,JJASlope[,,1],zlim=c( 0,5.5))
image.plot( Longitude,Latitude,Ufield[,,1] + Umean,zlim=c( 0,5.5))

image.plot( Longitude,Latitude,JJASlope[,,1] - (Ufield[,,1] + Umean) )

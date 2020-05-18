setwd('~/Dropbox/Home/Projects/BRACEData')
library( ncdf4)
library( fields)

##############################
########## setup
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
gw  <- ncvar_get( handle, 'gw')
gw<- gw/ sum( gw)
# Gaussian weights for finding global surface average.
gwField<- matrix( gw, 288, 192, byrow=TRUE )
gwField<- gwField/ sum( gwField)
# these are Stacey's weights
cos.weights <- cos(lat*pi/180)
cos.weights<- cos.weights/ sum(cos.weights)
testField<- matrix( cos.weights, 288, 192, byrow=TRUE )
testField<- testField/sum( testField)
test.for.zero( testField, gwField )

gwField<- testField

seasonalAverage<- function( x){
  filter(x,rep(1/3,3), sides=2)}

globalAverage<- function( field){
  sum(gwField*field)
}

patternScaling<- function( Tmp, X){
  return(lsfit( X,Tmp)$coef)
}
######################################
######## end setup

load( "StaceyWork/tas.rcp85.20062080.1.JJA.Rdata")
j<-100
k<- 100
staceyTest<- tas.k.JJA.sort[j,k,]

member<- 1
handle<- nc_open(fileNames[member]  )

TS<- ncvar_get( handle, 'TS')
TS<- TS - 273.15
TS<- TS[index.sort,,]
set.panel(2,1)
matplot(Longitude, TS[,1:10,1],type="l")
matplot(Longitude, TS[,182:192,1],type="l")



tempY<- (seasonalAverage(TS[j,k,]))[month==7]

test.for.zero(tempY, staceyTest )

tempX<- (seasonalAverage(apply( TS, 3, globalAverage)) )[month==7]

lsfit( tempX, tempY)$coef[2]
Ufield[j,k,1] + Umean[j,k]

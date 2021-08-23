LKinfo.CO2<- LKrigSetup(CO2$lon.lat, NC=150,nlevel=1, lambda=.2,
                        a.wght=4.1, alpha=1, 
                        LKGeometry="LKRing",
                        normalize=FALSE, 
                        choleskyMemory = list(nnzR=8e6) )
print(LKinfo.CO2)                                          
system.time( 
  obj1<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo.CO2)
)

system.time( 
  obj2<- LKrig( CO2$lon.lat,CO2$y, wX= obj1$wX,
                LKinfo=LKinfo.CO2, NtrA= 0)
)
system.time( look<- LKrig.sim(CO2$lon.lat, 
                              LKinfo = LKinfo.CO2 ))


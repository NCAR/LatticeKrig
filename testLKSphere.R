
data(CO2)
# estimate lambda ( should be around 0.003)
# NOTE: lambda values will tend to be sensitive to the model choice
dType0<- "GreatCircle"
attr(dType0,"Radius")<- .99
LKinfo0<- LKrigSetup( CO2$lon.lat, startingLevel=2 ,nlevel=2,
                      a.wght=1.1, alpha=c(1,.25),
                      LKGeometry="LKSphere",
                      distance.type= dType0) 
obj0B<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo0)
surface( obj0B, col=terrain.colors(256))
world( add=TRUE, col="magenta")

# use chordal distance 
dtype<- "Chordal"
attr(dtype,"Radius")<- 1.0

LKinfo1<- LKrigSetup( CO2$lon.lat, startingLevel=2 ,nlevel=2,
                      a.wght=1.1, alpha=c(1,.25),
                      LKGeometry="LKSphere", distance.type=dtype) 

obj0C<-  LatticeKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo1)
surface( obj0C, col=terrain.colors(256))
world( add=TRUE, col="magenta")



x1<- cbind( 5,-5)
look<- LKrig.basis(x1, LKinfo0)
look<- LKrig.basis(x1, LKinfo1)



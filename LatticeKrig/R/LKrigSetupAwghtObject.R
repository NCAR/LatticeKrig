LKrigSetupAwghtObject<- function( object){
# object is an LKinfo object
# should already have the lattice node information. 
  nLevel<- object$nlevel
  a.wght<- list()
  for( k in 1: nLevel ){
    latticeGrid<-  make.surface.grid((object$latticeInfo$grid )[[k]])
# if there is only one object use this for the a.wghts at
# every level
    aTemp <- predict( object$a.wghtObject, latticeGrid )
# check for missing values -- most likely outside range of spatial domain.    
    if( any( is.na(aTemp)) ){
      print(cbind( latticeGrid, aTemp))
      stop("Some a.wghts set to NAs from predict")
    }
# NOTE to check this for a 2-d geometry
# image.plot( as.surface( latticeGrid, aTemp))
# coerce to a matrix, if aTemp is a vector this will be a one column 
# matrix
    aTemp<- as.matrix( aTemp )
    a.wght<- c( a.wght, list( aTemp)  )
  }
  return(a.wght)
}
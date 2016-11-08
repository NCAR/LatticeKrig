omega2Awght<- function( omega, LKinfo){
  # for a rectangle this should be:
  #  Awght <- 4 +  exp( omega)^2
  xDimension<- dim(LKinfo$x)[2]
  Awght<- LKinfo$floorAwght + exp( omega * xDimension )
  return( Awght )
}

Awght2Omega<- function( Awght, LKinfo){
  # for a rectangle this should be:
  #  Awght <- 4 +  exp( omega)^2
  # omega  <- log(  Awght -4)/2
  if( Awght <= LKinfo$floorAwght){
    stop(paste("Awght is less than or equal to " , LKinfo$floorAwght) )
  }
  xDimension<- dim(LKinfo$x)[2]
  omega<-  log(Awght - LKinfo$floorAwght)/xDimension
  return( omega )
}
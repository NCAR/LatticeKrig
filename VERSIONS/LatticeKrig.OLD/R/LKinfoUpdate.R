
LKinfoUpdate<- function( LKinfo, ...){
   LKinfoCall<- as.list(LKinfo$call)
   argList<- list( ...)
   LKinfoCall[1] <- NULL
  
   for( argName in names( argList)){
 #    print(  LKinfoCall[ argName])
  #   print(  argList[[ argName]])
  #    LKinfo[ argName] <- NULL
      LKinfoCall[[ argName]] <- argList[[ argName]]}
#   print( LKinfoCall)
   do.call( "LKrig.setup" , LKinfoCall )
 }

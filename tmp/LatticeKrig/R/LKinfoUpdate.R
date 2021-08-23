# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

LKinfoUpdate<- function( LKinfo, ...){
  argList<- list( ...)
  # Is there anything to do?
  nArgs<- length( argList)
  if(  nArgs==0  ){
    return( LKinfo)
  }
  argNames<- names( argList)
# if only lambda is being changed then just substitute this value  
  if( (argNames[1]=="lambda")& (nArgs==1) ){
    LKinfo$lambda<- argList$lambda
    return( LKinfo)
  }
# Other arguments may require updating the LKinfo object
  LKinfoNames<- names( LKinfo)
  LKinfoSetupArgsNames<- names( LKinfo$setupArgs)
  # find where (or if) the new args match in LKinfo or setupArgs
  # NOTE: match("foo", NULL) gives an NA which is what we want!
  ind1<- match( argNames, LKinfoNames)
  ind2<- match( argNames, LKinfoSetupArgsNames )
  bad<- is.na(ind1)& is.na( ind2)
  if( any( bad)) {
    stop(
      cat(    "No match in current LKinfo for these new argument(s) ",
               argNames[bad]
         )
      )
  }
  setupArgsNew<- LKinfo$setupArgs
  LKinfoNew<- LKinfo
  # wipe out setupArgs in new LKinfo
  # otherwise they will appear twice in final list
  LKinfoNew$setupArgs<- NULL
  #
  # overwrite or add new values for extra ... arguments. 
  if( any( !is.na(ind1) ) ){
  for( argN in LKinfoNames[ind1]){
    LKinfoNew[[ argN]] <- argList[[ argN ]]
  }
  }
 #
  if( any( !is.na(ind2) ) ){
  for( argN in LKinfoSetupArgsNames[ind2]){
    setupArgsNew[[argN]] <- argList[[ argN]]
  }
  }
  # now call LKrigsetup again to update the LKinfo object with new values. 
  # this will also check that the arguments are consistent.
  LKinfoNew <- do.call( "LKrigSetup" , c( LKinfoNew, setupArgsNew) )
  # maybe reorder the components of this new list so in same order as
  # LKinfo (see example in help file for LKrigSetup)
  return(LKinfoNew)
}


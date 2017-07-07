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
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

summary.LKrig <- function(object, digits = 4, stripAwght=TRUE, ...) {
  x<- object
  obj<- list()
  LKinfo <- x$LKinfo
 
  if (is.matrix(x$residuals)) {
    n <- nrow(x$residuals)
    NData <- ncol(x$residuals)
  }
  else {
    n <- length(x$residuals)
    NData <- 1
  }
  c1 <- "Number of Observations:"
  c2 <- n
  if (NData > 1) {
    c1 <- c(c1, "Number of data sets fit:")
    c2 <- c(c2, NData)
  }
  c1 <- c(c1, "Number of parameters in the fixed component")
  c2 <- c(c2, x$nt)
  if (x$nZ > 0) {
    c1 <- c(c1, "Number of covariates")
    c2 <- c(c2, x$nZ)
  }
  if (!is.na(x$eff.df)) {
    c1 <- c(c1, " Effective degrees of freedom (EDF)")
    c2 <- c(c2, signif(x$eff.df, digits))
    c1 <- c(c1, "   Standard Error of EDF estimate: ")
    c2 <- c(c2, signif(x$trA.SE, digits))
  }
  c1 <- c(c1, "Smoothing parameter (lambda)")
  c2 <- c(c2, signif(x$lambda, digits))
  
  c1 <- c(c1, "MLE sigma ")
  c2 <- c(c2, signif(x$sigma.MLE.FULL, digits))
  
  c1 <- c(c1, "MLE rho")
  c2 <- c(c2, signif(x$rho.MLE.FULL, digits))

  
  c1 <- c(c1, "Total number of basis functions")
  c2 <- c(c2,  LKinfo$latticeInfo$m)
  
  c1 <- c(c1, "Multiresolution levels")
  c2 <- c(c2,  LKinfo$nlevel)
  
  c1<- c(c1,"log Profile Likelihood")
  c2<- c( c2, signif(x$lnProfileLike.FULL,10))
  
  c1<- c(c1,"log  Likelihood (if applicable)")
  c2<- c( c2, x$lnLike.FULL)
  
  c1 <- c(c1, "Nonzero entries in Ridge regression matrix")
  c2 <- c(c2, x$nonzero.entries)
  
  M<- length( c1)
  summary <-  data.frame(c1, c2)
  names( summary)<- c("", "")
    
  obj$call<-x$call 
  obj$inverseModel <- x$inverseModel
  obj$parameters<- summary
  obj$timingLKrig <- x$timingLKrig
  # fastNormalization matrices will add to the size of LKinfo
  # just reset to message
  if( stripAwght){
  attr(LKinfo$a.wght, 
  which = "fastNormDecomp") <- "Omitted by summary.LKrig function"
  }
  obj$LKinfo<- LKinfo
  obj$MLE<- x$MLE
  return(obj)
}




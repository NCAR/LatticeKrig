# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

print.LKrig <- function(x, digits = 4, ...) {
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
    if (NData == 1) {
        c1 <- c(c1, "MLE sigma ")
        c2 <- c(c2, signif(x$shat.MLE, digits))
        c1 <- c(c1, "MLE rho")
        c2 <- c(c2, signif(x$rho.MLE, digits))
    }
    c1 <- c(c1, "Nonzero entries in Ridge regression matrix")
    c2 <- c(c2, x$nonzero.entries)
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat(" ", fill = TRUE)
#    
    if( is.null( x$fixedFunction)){  
        cat("No fixed part of model", fill = TRUE)
    }
    else{
        if( x$fixedFunction=="LKrigDefaultFixedFunction"){
            cat("Fixed part of model is a polynomial of degree", x$fixedFunctionArgs$m - 1, "(m-1)", fill=TRUE)
        }  
        else{  
          cat("Fixed part of model uses the function:", x$fixedFunction, fill = TRUE)
          cat("with the argument list:", fill = TRUE)
          print( x$fixedFunctionArgs)
        }
    }
    cat("Radial basis function used: ", LKinfo$RadialBasisFunction, 
        fill = TRUE)
    cat(" ", fill = TRUE)
    
    cat(LKinfo$nlevel, "level(s)", LKinfo$m, " basis functions", 
        fill = TRUE)
    for (k in 1:LKinfo$nlevel) {
        cat("Lattice level", k, "is ", LKinfo$mx[k], "X", LKinfo$my[k], 
            "spacing", LKinfo$delta[k], fill = TRUE)
    }
    cat("Total number of basis functions: ", x$m, "  with overlap of ", 
        LKinfo$overlap, fill = TRUE)
    cat(" ", fill = TRUE)
#
    cat("Type of distance metric used: ", LKinfo$distance.type, 
        fill = TRUE)
#
    cat(" ", fill = TRUE)
    if (length(LKinfo$alpha[[1]]) == 1) {
        cat("Value(s) for weighting (alpha): ", unlist(LKinfo$alpha), 
            fill = TRUE)
    }
    else {
        cat("alpha values passed as a vector for each level", 
            fill = TRUE)
    }
#    
    cat(" ", fill = TRUE)
    if (length(LKinfo$a.wght[[1]]) == 1) {
        a.wght <- unlist(LKinfo$a.wght)
        cat("Value(s) for lattice dependence (a): ", a.wght, 
            fill = TRUE)
    }
    else {
        cat("Value(s) for weighting in GMRF (a.wght): ", unlist(LKinfo$alpha), 
            fill = TRUE)
    }
#    
    cat(" ", fill = TRUE)
    if (LKinfo$normalize) {
        cat("Basis functions normalized so marginal process variance is stationary", 
            fill = TRUE)
    }
    
    cat("Number of  lattice buffer points added in the x and y coordinates:",
        LKinfo$NC.buffer.x,",",LKinfo$NC.buffer.y, fill = TRUE)
    
    invisible(x)
}


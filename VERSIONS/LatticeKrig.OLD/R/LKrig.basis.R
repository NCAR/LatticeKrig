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

LKrig.basis <- function(x1, LKinfo, verbose = FALSE, 
    spam.format = TRUE) {
    grid.info <- LKinfo$grid.info
    nlevel <- LKinfo$nlevel
    overlap <- LKinfo$overlap
    normalize <- LKinfo$normalize
    scale.basis <- LKinfo$scale.basis
    distance.type <- LKinfo$distance.type
    V <- LKinfo$V
    # We will transform (scale and rotate) x matrix of locations by   x%*% t(solve(V))
    #
    # For the radial basis functions this
    # will imply that the Euclidean distance norm used between two vectors X1 and X2
    # is  t(X1-X2) solve(A) (X1-X2)  with A = (V %*%t(V))
    # Here's the deal on the linear transformation V:
    # It should be equivalent to just transforming the locations outside of LatticeKrig.
    # Where ever the observation locations are used
    # they are first transformed with t(solve(V))).
    # Surprisingly this needs to happen in one place below and in LKrig.setup to
    # determine the grid centers.
    #
    # The RBF centers and the delta scale factor, however, assumed to be
    # in the transformed scale and so a  transformation is not needed.
    # see LKrig.setup for how the centers are determined.
    # The concept is that if LKrig was called with the transformed locations
    # ( x%*%t( solve(V)) instead of x
    # and V set to diag( c(1,1) ) one would get the same results.
    # as passing a non identity V.
    #
    # accumulate matrix column by column in PHI
    PHI <- NULL
    for (l in 1:nlevel) {
        # Loop over levels and evaluate basis functions in that level.
        # Note that all the center information based on the regualr grids is
        # taken from the LKinfo object
        delta <- LKinfo$delta[l]
        grid.list <- LKinfo$grid[[l]]
        centers <- make.surface.grid(grid.list)
        #  set the range of basis functions, they are assumed to be zero outside
        #  the radius basis.delta
        basis.delta <- delta * overlap
        #
        x1Transformed<- x1 %*% t(solve(V))
        PHItemp <- Radial.basis(x1Transformed, centers, 
            basis.delta, max.points = NULL, mean.neighbor = 50, 
            just.distance = FALSE, RadialBasisFunction = get(LKinfo$RadialBasisFunction), 
            distance.type = LKinfo$distance.type)
normtime<-  system.time(
        if (normalize) {
            if( LKinfo$fastNormalization){
               wght <- LKrig.normalize.basis.fast(l, LKinfo, x1Transformed)
             }
            else{
               tempB <- LKrig.MRF.precision(LKinfo$mx[l], LKinfo$my[l], 
                     a.wght = (LKinfo$a.wght)[[l]], edge = LKinfo$edge, 
                     distance.type = LKinfo$distance.type)
               tempB <- LKrig.spind2spam(tempB)
               wght <- LKrig.quadraticform(  t(tempB) %*% tempB, PHItemp)
             }
# now normalize the basis functions by the weights            
# awkwardness of handling the 1 row case separately.
            if (nrow(x1) > 1) {
                PHItemp <- diag.spam(1/sqrt(wght)) %*% PHItemp
            }
            else {
                PHItemp@entries <- PHItemp@entries/sqrt(wght)
            }
        }
 )
# cat( "Level",l, normtime, fill=TRUE)
        # accumulate new level of basis function.
        PHI <- cbind(PHI, PHItemp)
    }
    # include a spatially varying multiplication of process.
    # wght are in scale of inverse marginal variance of process
    if (scale.basis) {
        wght <- c(predict(LKinfo$rho.object, x1))
        if (nrow(x1) > 1) {
            PHI <- diag.spam(sqrt(wght)) %*% PHI
        }
        else {
            PHI@entries <- PHI@entries * sqrt(wght)
        }
    }
    # attach  LKinfo list to the matrix to help identify how the basis functions
    # are organized.
    attr(PHI, which = "LKinfo") <- LKinfo
    return(PHI)
}

LKrig.quadraticform <- function(Q, PHI) {
    #   finds     the quadratic forms    PHI_j^T Q.inverse PHI_j  where PHI_j is the jth row of
    #   PHI.
    #   The loop is to avoid using memory for the entire problem if more than 2000 elements.
    nrow <- dim(PHI)[1]
    if (nrow > 1) {
        BLOCKSIZE <- 2000
        wght <- rep(NA, nrow)
        counter <- 1
        while (counter < nrow) {
            ind <- counter:min((counter + BLOCKSIZE), nrow)
            A <- forwardsolve(l = chol(Q), transpose = TRUE, 
                x = t(PHI[ind, ]), upper.tri = TRUE)
            wght[ind] <- c(colSums(A^2))
            counter <- counter + BLOCKSIZE
        }
    }
    else {
        # Unfortunately the case when there is one row needs to be handled separately.
        A <- forwardsolve(l = chol(Q), transpose = TRUE, x = t(PHI), upper.tri = TRUE)
        wght <- sum(A^2)
    }
    return(wght)
}





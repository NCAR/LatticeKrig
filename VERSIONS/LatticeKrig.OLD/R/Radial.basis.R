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

Radial.basis <- function(x1, centers, delta, max.points = NULL, 
    mean.neighbor = 50, just.distance = FALSE, RadialBasisFunction = "WendlandFunction", 
    distance.type = "Euclidean") {
    
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(centers)
    if (is.null(max.points)) {
        Nmax <- n1 * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
    # if delta is a scalar repeat the length of columns
    # then check for length n2
    if (length(delta) == 1) {
        delta <- rep(delta, n2)
    }
    if (length(delta) != n2) {
        stop("delta is not the length of n2")
    }
    # evalute RBF basis functions at the x1 locations with
    # nodes given by centers. Returned is a sparse matrix in
    # a simpler format than spam  (row/column indices and value)
    # returned values ind and rd below in 'out' object.
    # If the model is on the sphere or cylinder, convert to 3-d
    # coordinates to find chordal distances.
    # NOTE: longitude coordinates assumed to be in degrees.
    #
    # check distance choices
    if (is.na(match(distance.type, c("Euclidean", "cylinder")))) {
        stop("distance type is not supported (or is misspelled!).")
    }
    #
    if (distance.type == "Euclidean") {
        out <- .Fortran("dfind2d", x1 = as.double(x1), n1 = as.integer(n1), 
            center = as.double(centers), n2 = as.integer(n2), 
            delta2 = as.double(delta^2), ind = as.integer(rep(0, 
                Nmax * 2)), rd = as.double(rep(-1, Nmax)), Nmax = as.integer(Nmax), 
            iflag = as.integer(1), PACKAGE = "LatticeKrig")
    }
    if (distance.type == "cylinder") {
        # convert to radians
        R <- 360/(2 * pi)
        out <- .Fortran("dfind3d", x1 = as.double(LKrig.cyl(x1, 
            R)), n1 = as.integer(n1), center = as.double(LKrig.cyl(centers, 
            R)), n2 = as.integer(n2), delta2 = as.double(delta^2), 
            ind = as.integer(rep(0, Nmax * 2)), rd = as.double(rep(-1, 
                Nmax)), Nmax = as.integer(Nmax), iflag = as.integer(1), 
            PACKAGE = "LatticeKrig")
    }
    # negative iflag means one has run out of space
    if (out$iflag == -1) {
        cat("temp space set at", Nmax, fill = TRUE)
        stop("Ran out of space, increase max.points")
    }
    # trim down to a sparse matrix object where the elements that have
    # nonvalues (there are Nmax of these)
    N <- out$Nmax
#    print(N)
#    return(  list(ind = matrix(out$ind, Nmax, 2)[1:N, ], ra = out$rd[1:N], da = c(n1, n2)))
#     out@entries <- do.call(RadialBasisFunction, list(d = out@entries/delta[out@colindices]))
    # in spind format:
    out<- list(ind = matrix(out$ind, Nmax, 2)[1:N, ], ra = out$rd[1:N], da = c(n1, n2))
    if (just.distance) {
        return(out)
    }
    out$ra <- do.call(RadialBasisFunction, list(d = out$ra/delta[ out$ind[,2] ]) )
    out <- LKrig.spind2spam(out)
 # evaluate distance  with RBF ---  usually Wendland2.2
    return(out)
}

LKrig.cyl <- function(U, R = 1) {
    # x[,1] assumed to be in degrees scale to radians
    theta <- 2 * pi * U[, 1]/360
    # distance units are in degrees for x and should be scaled (using V in LKinfo object)
    # to be comparable for y
    cbind(R * cos(theta), R * sin(theta), U[, 2])
}

WendlandFunction <- function(d) {
    # dimension =2 k =2
    if (any(d < 0)) 
        stop("d must be nonnegative")
    return(((1 - d)^6 * (35 * d^2 + 18 * d + 3))/3 * (d < 1))
}

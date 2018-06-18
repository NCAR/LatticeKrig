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

LKrig.coef <- function(Mc, wPHI, wT.matrix, wy, lambda, 
    weights) {
    if (length(lambda) > 1) {
        stop("lambda must be a scalar")
    }
    A <- forwardsolve(Mc, transpose = TRUE, t(wPHI) %*% wT.matrix, 
        upper.tri = TRUE)
    A <- backsolve(Mc, A)
    A <- t(wT.matrix) %*% (wT.matrix - wPHI %*% A)/lambda
    #   A is  (T^t M^{-1} T)
    b <- forwardsolve(Mc, transpose = TRUE, t(wPHI) %*% wy, upper.tri = TRUE)
    b <- backsolve(Mc, b)
    b <- t(wT.matrix) %*% (wy - wPHI %*% b)/lambda
    # b is   (T^t M^{-1} y)
    # Save the intermediate matrix   (T^t M^{-1} T) ^{-1}
    # this the GLS covariance matrix of estimated coefficients
    # should be small for this to be efficient code -- the default is 3X3
    Omega <- solve(A)
    # GLS estimates
    d.coef <- Omega %*% b
    # coefficients of basis functions.
    c.coef <- forwardsolve(Mc, transpose = TRUE, t(wPHI) %*% 
        (wy - wT.matrix %*% d.coef), upper.tri = TRUE)
    c.coef <- backsolve(Mc, c.coef)
    c.mKrig <- sqrt(weights) * (wy - wPHI %*% c.coef - wT.matrix %*% 
        d.coef)/lambda
    
    return(list(c.coef = c.coef, d.coef = d.coef, Omega = Omega, 
        c.mKrig = c.mKrig))
}


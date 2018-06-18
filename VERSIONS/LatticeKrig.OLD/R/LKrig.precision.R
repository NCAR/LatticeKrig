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

LKrig.precision <- function(LKinfo, return.B = FALSE, level.index = NA,
                                   verbose=FALSE) {
    mx <- LKinfo$mx
    my <- LKinfo$my
    
    grid.info <- LKinfo$grid.info
    L <- LKinfo$nlevel
    if( any(unlist(LKinfo$a.wght)<4)){
        stop("a.wght less than 4")}
    if (L != length(my)) {
        stop("number of levels and mx and my are not consistent")
    }
    offset <- LKinfo$offset
    distance.type = LKinfo$distance.type
    
    # some checks on arguments
    if ((length(LKinfo$alpha) != L) | (length(LKinfo$a.wght) != 
        L)) {
        stop("number of levels for alpha or a.wght not right length")
    }
    # ind holds non-zero indices and ra holds the values
    ind <- NULL
    ra <- NULL
    da <- rep(0, 2)
    if (is.na(level.index)) {
        # loop over levels
        for (j in 1:L) {
            # evaluate the H matrix at level j.
            # each row of this matrix has an 'a.wght[j]'  on diagonal and
            # -1  for the nearest neighbor basis functions.
            tempB <- LKrig.MRF.precision(mx[j], my[j], a.wght = (LKinfo$a.wght)[[j]], 
                stationary = LKinfo$stationary, edge = LKinfo$edge, 
                distance.type = distance.type)
            # multiply this block by 1/ sqrt(diag( alpha[[j]]))
            alpha.level <- (LKinfo$alpha)[[j]]
            if (length(alpha.level) == 1) {
                tempra <- 1/sqrt(alpha.level[1]) * tempB$ra
            }
            else {
                rowindices <- tempB$ind[, 1]
                tempra <- 1/sqrt(alpha.level[rowindices]) * tempB$ra
            }
            # accumulate the new block
            # for the indices that are not zero
            ra <- c(ra, tempra)
            ind <- rbind(ind, tempB$ind + offset[j])
            # increment the dimensions
            da[1] <- da[1] + tempB$da[1]
            da[2] <- da[2] + tempB$da[2]
        }
        # dimensions of the full matrix
        # should be da after loop
        # check this against indices in LKinfo
        #
        if ((da[1] != LKinfo$offset[L + 1]) | (da[2] != LKinfo$offset[L + 
            1])) {
            stop("Mismatch of dimension with size in LKinfo")
        }
        # convert to spam format:
        temp <- LKrig.spind2spam(list(ind = ind, ra = ra, da = da))
    }
    else {
        # case to find precision at just one arbitrary level
        tempB <- LKrig.MRF.precision(mx[level.index], my[level.index], 
            a.wght = (LKinfo$a.wght)[[level.index]], stationary = LKinfo$stationary, 
            edge = LKinfo$edge, distance.type = distance.type)
        alpha.level <- (LKinfo$alpha)[[level.index]]
        if (length(alpha.level) == 1) {
            tempra <- 1/sqrt(alpha.level[1]) * tempB$ra
        }
        else {
            rowindices <- tempB$ind[, 1]
            tempra <- 1/sqrt(alpha.level[rowindices]) * tempB$ra
        }
        temp <- LKrig.spind2spam(list(ind = tempB$ind, ra = tempra, 
            da = tempB$da))
    }
    if (return.B) {
        return(temp)
    }
    else {
        # find precision matrix Q = t(B)%*%B and return
        return(t(temp) %*% (temp))
    }
}


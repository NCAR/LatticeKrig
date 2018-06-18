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

which.max.matrix <- function(z) {
    if (!is.matrix(z)) {
        stop("Not a matrix")
    }
    m <- nrow(z)
    n <- ncol(z)
    # take care of NAs
    ind <- which.max(z)
    iy <- trunc((ind - 1)/m) + 1
    ix <- ind - (iy - 1) * m
    return(cbind(ix, iy))
}


which.max.image <- function(obj) {
    ind.z <- which.max.matrix(obj$z)
    return(list(x = obj$x[ind.z[, 1]], y = obj$y[ind.z[, 2]], 
        z = obj$z[ind.z], ind = ind.z))
}


LKrig.rowshift <- function(A, shift.row, shift.col) {
    mx <- dim(A)[1]
    my <- dim(A)[2]
    if (abs(shift.row) > mx) {
        stop("shift exceeds matrix dimension")
    }
    if (shift.row < 0) {
        i.target <- 1:(mx + shift.row)
        i.source <- (-shift.row + 1):mx
    }
    else {
        i.target <- (shift.row + 1):mx
        i.source <- (1:(mx - shift.row))
    }
    
    B <- matrix(NA, mx, my)
    B[i.target, ] <- A[i.source, ]
    return(B)
}

LKrig.rowshift.periodic <- function(A, shift.row) {
    mx <- dim(A)[1]
    if (abs(shift.row) > mx) {
        stop("shift exceeds matrix dimension")
    }
    if (shift.row < 0) {
        i.source <- c((-shift.row + 1):mx, 1:(-shift.row))
    }
    else {
        i.source <- c(mx - (shift.row:1) + 1, 1:(mx - shift.row))
    }
    return(A[i.source, ])
}

LKrig.shift.matrix <- function(A, shift.row = 0, shift.col = 0, 
    periodic = c(FALSE, FALSE)) {
    if (shift.row != 0) {
        if (!periodic[1]) {
            A <- LKrig.rowshift(A, shift.row = shift.row)
        }
        else {
            A <- LKrig.rowshift.periodic(A, shift.row = shift.row)
        }
    }
    if (shift.col != 0) {
        if (!periodic[2]) {
            A <- t(LKrig.rowshift(t(A), shift.row = shift.col))
        }
        else {
            A <- t(LKrig.rowshift.periodic(t(A), shift.row = shift.col))
        }
    }
    return(A)
}


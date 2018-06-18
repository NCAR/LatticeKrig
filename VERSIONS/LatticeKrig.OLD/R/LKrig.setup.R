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

LKrig.setup <- function(x = NULL, NC = NULL, NC.buffer = 5, 
    nlevel=NULL, grid.info = NULL, lambda = NA, sigma = NA, rho = NA, 
    alpha = NA, nu = NULL, a.wght = NA, overlap = 2.5, normalize = TRUE, 
    normalize.level = NULL, edge = FALSE, rho.object = NULL, 
    RadialBasisFunction = "WendlandFunction", distance.type = "Euclidean", 
    V = diag(c(1, 1)), verbose = FALSE) {
#
# this serves as an interface from more simply described parameters and fills in defaults
# and other information that make the computational code cleaner.
#  
# update grid.info if it has not been specified and do some checks on the distance type.
    grid.info<- LKrig.make.grid.info(grid.info, NC, x, V, distance.type, verbose)
    # find the grid centers
    out.grid <- LKrig.make.centers(grid.info, nlevel, NC.buffer, distance.type)
    mx<- out.grid$mx
    my<- out.grid$my
    
#
# figure out form of model and form of lattice links,  reformat a.wght as a list over levels 
    out.a.wght<- LKrig.make.a.wght( a.wght, nlevel,mx,my)
    fastNormalization<- out.a.wght$fastNormalization
    a.wght<- out.a.wght$a.wght
# now fix up alpha    
    out.alpha<- LKrig.make.alpha( alpha, nu, nlevel)
    alpha<- out.alpha$alpha
# Check some details about scaling the basis functions and how they are
# normalized
    scale.basis <- !is.null(rho.object)
    if (scale.basis & !normalize) {
        stop("Can not scale an unnormalized basis")
    }
    if (is.null(normalize.level)) {
        normalize.level = rep(normalize, nlevel)
    }
# if just the center a.wght and it is constant for each level, setup some matrices for fast normalization
    if( fastNormalization){
        NormalizeList<- LKrig.make.Normalization( mx,my, a.wght)
    }
    else{
        NormalizeList<-NULL
    }
# set lambda if sigma and rho are passed.
    if (is.na(lambda[1])) {
        lambda <- sigma^2/rho
    }
#Note call argument in return allows the function to be re evaluated
    out <- list(nlevel = nlevel, grid.info = grid.info, NC = NC, 
        NC.buffer = NC.buffer, delta = out.grid$delta.save, mx = mx, 
        my = my, m = out.grid$m, m.domain = out.grid$m.domain, 
        NC.buffer.x = out.grid$NC.buffer.x, NC.buffer.y = out.grid$NC.buffer.y, 
        offset = out.grid$offset, grid = out.grid$grid, overlap = overlap, 
        alpha = alpha, nu = nu, a.wght = a.wght, stationary = out.a.wght$stationary, 
        lambda = lambda, sigma = sigma, rho = rho, normalize = normalize, 
        normalize.level = normalize.level, edge = edge, scalar.alpha = out.alpha$scalar.alpha, 
        scale.basis = scale.basis, rho.object = rho.object, RadialBasisFunction = RadialBasisFunction, 
        distance.type = distance.type, V = V,
        first.order= out.a.wght$first.order,
        fastNormalization=fastNormalization, NormalizeList=NormalizeList,
        call=match.call() )
    class(out) <- "LKinfo"
    return(out)
}


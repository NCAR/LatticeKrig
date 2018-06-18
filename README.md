# LatticeKrig: Multiresolution Kriging Based on Markov Random Fields

 
This repository contains some supplemental material in this top level and see the subdirectory **LatticeKrig** for the "head" of the standard R package. 
The most current package on CRAN  is listed here as *LatticeKrig_VERSION.tar.gz* . At the time of writing the version is *7.0*.
To create a possibly new version from this repository download the **LatticeKrig** subdirectory and in UNIX
```
 R CMD build --force LatticeKrig
```
To install this version for your system try
``` 
R CMD INSTALL LatticeKrig
```

or from the tar ball
```
R CMD INSTALL LatticeKrig_VERSION.tar.gz
```
where VERSION are the correct version numbers in this file.

# Package description
  Methods for the interpolation of large spatial
  datasets. This package follows a "fixed rank Kriging" approach using
  a large number of basis functions and provides spatial estimates
  that are comparable to standard families of covariance functions.
  Using a large number of basis functions allows for estimates that
  can come close to interpolating the observations (a spatial model
  with a small nugget variance.)  Moreover, the covariance model for
  this method can approximate the Matern covariance family but also
  allows for a multi-resolution model and supports efficient
  computation of the profile likelihood for estimating covariance
  parameters. This is accomplished through compactly supported basis
  functions and a Markov random field model for the basis
  coefficients. These features lead to sparse matrices for the
  computations and this package makes of the R spam package for this.
  An extension of this version over previous ones ( < 5.4 ) is the
  support for different geometries besides a rectangular domain. The
  Markov random field approach combined with a basis function
  representation makes the implementation of different geometries
  simple where only a few specific functions need to be added with
  most of the computation and evaluation done by generic routines that
  have been tuned to be efficient.  One benefit of this package's
  model/approach is the facility to do unconditional and conditional
  simulation of the field for large numbers of arbitrary points. There
  is also the flexibility for estimating non-stationary covariances
  and also the case when the observations are a linear combination
  (e.g. an integral) of the spatial process.  Included are generic
  methods for prediction, standard errors for prediction, plotting of
  the estimated surface and conditional and unconditional simulation.

-


The reference  DOI **10.5065/D6W957CT**is linked to the specific package version 5.5: [Versions/LatticeKrig_5.5.tar.gz](LatticeKrig_5.5.tar.gz)

MD5 check sum: **ebefee1efcf3b6da395d7c9c540fee4a**
  
For the most recent, distributed version of LatticeKrig  use a CRAN mirror site such as [R studio CRAN mirror site](http://cran.rstudio.com/) to download and install the package. Use ````citation("LatticeKrig")```` in R to generate a citation for this package with the current version number. 


Package contact:  Doug Nychka **nychka@ucar.edu**







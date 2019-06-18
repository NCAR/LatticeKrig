\name{registeredFORTRAN}
\alias{findnorm}
\alias{lkdist}
\alias{lkdistcomp}
\alias{lkdistgrid}
\alias{lkdistgridcomp}
\docType{data}
\title{
Internal FORTRAN routines for working with grids and finding distances. 
}
\description{
These are objects of class  \code{FortranRoutine} and also 
\code{\link{NativeSymbolInfo}}
They provide information for compiledfunctions called with \code{.Call},
or \code{.Fortran}.
Ordinarily one would not need to consult these and they are used to make
the search among dynamically loaded libraries ( in particular the fields library
) have less ambiguity and
also be faster. These are created when the package/library is loaded 
are have their definitions from the compliation of  \code{init.c} in the
package source (src) directory. 
}
\format{
  The format is a list with components:
  \describe{
  \item{name}{The (registration ?) name of the C function.}
  \item{address}{See \link{NativeSymbolInfo}. }
  \item{dll}{Dynamically linked library information.}
  \item{numParameters}{Number of calling arguments in function.}
  }
  }

\details{
Registered routines are
\describe{
 \item{findnorm}{Finds marginal variance quickly for the case of a rectangular lattice,
 a specific Wendland basis, and constant a.wght at a given level. Called from 
 \code{LKRectangleFastNormalization}}
 \item{lkdist}{Euclidean distance between two sets of coordinates but restricted to a
 maximum distance.Called from \code{LKDist}.}
 \item{lkdistcomp}{Same as lkdist but returns the component distances for each
 coordinate.Called by \code{LKDistComponents}}
 \item{lkdistgrid}{Finds all distances in x1 that are within delta of a set of grid
 points. Results are returned in the row, column,  value sparse matrix format.called from \code{LKDistGrid}}
 \item{lkdistgridcomp}{Same as lkdistgrid but the component distances are found for
 each coordinate.\code{LKDistGridComponents}}
 }
 See \code{package_native_routine_registration_skeleton} for the 
utility used to create these data objects.
%%  ~~ If necessary, more details than the __description__ above ~~
}

\references{
For background on registering C, C++ and Fortran functions see 5.4 of
Writing R Extensions. 
See \url{http://r.789695.n4.nabble.com/Registration-of-native-routines-td4728874.html} for additional dicsussion of code registration. 
}
\examples{
print(lkdistgridcomp)

}
\keyword{datasets}

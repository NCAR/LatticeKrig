\name{IcosahedronGrid}
\alias{IcosahedronGrid}
\alias{IcosahedronFaces}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Icosahedral multi-resolution grids
%%  ~~function to do ... ~~
}
\description{
Creates a multi-resolution grid based on subdividing triangles starting with faces of an 
icosahedron.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IcosahedronGrid(K)

IcosahedronFaces(K)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{ Number of levels.}
 
%%     ~~Describe \code{K} here~~
}
\details{
Creates a nearly regular grid by taking the first
level as the 12 points from a regular icosahedron. The subsequent levels generate a finer set
of points by subdividing each triangular face into 4 new triangles. The three new mid points
from the subdivision are added to the previous nodes to give the new level of resolution. The
triangles tend to be roughly  equilateral and so the nodes will tend to be roughly equally
space but there is some variation in distances among nearest neighbors. 

To depict the faces and nodes in a snazzy way use the \code{rgl} package and the following code 
\preformatted{
library( rgl)
# show level 3
Level<- 3
    SGrid <- IcosahedronFaces(4)
    Tri <- SGrid$Face[[Level]]
    L <- dim(Tri)[3]
    plot3d(rbind(c(0, 0, 0)), xlim = c(-1, 1), ylim = c(-1, 1), 
        zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "", 
        zlab = "", type = "p", box = TRUE)
    for (k in 1:L) \{
        U <- Tri[, , k]
        rgl.triangles(U, col = "grey80")
    \}
    plot3d(SGrid$nodes[[Level]], col = "green4", type = "s", 
        radius = 0.03, add = TRUE)
}

		
%%%	\figure{s1.pdf}{First Level} \figure{s2.png}{Second Level} \figure{s3.png}{Third Level}
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\strong{IcosahedronGrid} A list with \code{K} components each component is a three column matrix
giving the direction cosines for each grid point. 

\strong{IcosahedronFaces} A list with components:
\describe{
\item{MultiGrid}{The same list returned by IcosahedronGrid. }
\item{Faces}{A list with \code{K-1} components. Each components are the faces at a given level
represented as a three dimensional array (3X3XN with N the number of faces and a given level). The
array indices are vertices of triangle, coordinates and faces within a resolution level.   e.g. to
extract the 10th face (out of 80) for  the 2nd level: 
\preformatted{
look<- IcosahedronFaces(3)$Faces
triangle <- (look[[2]])[,,10]
print (triangle)
           [,1]       [,2]      [,3]
[1,] -0.5257311 -0.8506508  0.000000
[2,] -0.8090170 -0.5000000  0.309017
[3,] -0.8090170 -0.5000000 -0.309017
rowSums( triangle^2)
[1] 1 1 1
}

\code{triangle} will be the 10th face for the second level where
the columns are the 3d
coordinates of the direction cosines and the rows index the three
vertices. 
}

}

}
\seealso{ \link{toSphere}, \link{directionCosines}, \link{projectionSphere} }
\author{
Doug Nychka and Zachary Thomas
%%  ~~who you are~~
}


\examples{
# second level in lon lat coordinates 
look<- IcosahedronGrid(3)
lonlat<- toSphere( look[[3]])
plot( lonlat, xlab="lon", ylab="lat")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }% use one of  RShowDoc("KEYWORDS")

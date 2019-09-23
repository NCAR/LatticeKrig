LKTomographyPoints <- function(latticeInfo) {
  levelSizes <- t(latticeInfo$mx)
  dimension <- as.integer(dim(levelSizes)[1])
  nLevels <- as.integer(dim(levelSizes)[2])
  nPoints <- as.integer(sum(apply(levelSizes, 2, prod)))
  levelSizes <- as.integer(levelSizes)
  points <- as.double(matrix(0, nrow = dimension, ncol = nPoints))
  coordinates <- as.double(unlist(latticeInfo$grid))
  if (dimension == 3) {
    out <- .Fortran("LKTomPoints3D", coorinates=coordinates, levelSizes=levelSizes,
                    nLevels=nLevels, points=points, nPoints=nPoints, PACKAGE="LatticeKrigInverseProblem")
  } else if (dimension == 2) {
    out <- .Fortran("LKTomPoints2D", coorinates=coordinates, levelSizes=levelSizes,
                    nLevels=nLevels, points=points, nPoints=nPoints, PACKAGE="LatticeKrigInverseProblem")
  } else {
    stop("LKTomography only works for 2 or 3 dimensional problems")
  }
  points <- matrix(out$points, nrow = dimension)
  return(points)
}

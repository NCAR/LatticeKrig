LKBasisFunctionIntegral <- function(normDists, completions, ranges, rangeReps, points, weights, lineIntegralMode = 3L) {
  dists <- as.double(normDists)
  completions <- matrix(as.double(completions), nrow=2)
  N = length(dists)
  
  #this loop finds the range (basis function radius) corresponding to the point of each entry
  rangeVec <- rep(0, length(dists))
  for(i in 1:length(rangeReps)) {
    points <- points - rangeReps[i]
    rangeVec[points <= 0 & rangeVec == 0] <- (ranges[i] * weights[[i]])
  }
  
  if (lineIntegralMode == 1L) {
    if (!exists("lineIntegralGrid")) {
      data(LineIntegralSplines)
    }
    lineIntegralGridLength <- dim(lineIntegralGrid)[2]
    completionGridLength <- dim(completionGrid)[2]
    lineIntegralDelta <- 1/lineIntegralGridLength
    completionDelta <- 1/completionGridLength
    lineIntegrals <- .Fortran("CubicInterp", grid = lineIntegralGrid, nGrid = lineIntegralGridLength, delta = lineIntegralDelta,
                              points = dists, nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
    starts <- .Fortran("CubicInterp", grid = completionGrid, nGrid = completionGridLength, delta = completionDelta,
                       points = completions[1,], nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
    stops <- .Fortran("CubicInterp", grid = completionGrid, nGrid = completionGridLength, delta = completionDelta,
                      points = completions[2,], nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
    return (lineIntegrals * (stops - starts) * rangeVec)
  } else {
    quadLevels = as.integer(2^(lineIntegralMode - 2))
    lineIntegrals <- .Fortran("WendlandGaussQuad", offset = dists, completions = as.double(completions), nEntries = N,
                              nLevels = quadLevels, entries = rep(-1.8, N), PACKAGE = "LatticeKrigInverseProblem")
    return(lineIntegrals$entries * rangeVec)
  }
}

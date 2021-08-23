LKBasisFunctionIntegralSphere <- function(normDists, completions, ranges, rangeReps, points) {
  dists <- as.double(normDists)
  completions <- matrix(as.double(completions), nrow=2)
  ranges <- as.double(ranges)
  rangeReps <- as.integer(rangeReps)
  N = length(dists)
  
  if(!exists("lineIntegralGrid")) {
    data(LineIntegralSplines)
  }
  if(!exists("sphereRangeGrid")) {
    data(SphereErrorCorrectionSplines)
  }

  
  #this loop finds the range (basis function radius) corresponding to the point of each entry
  rangeVec <- rep(0, length(dists))
  for(i in 1:length(rangeReps)) {
    points <- points - rangeReps[i]
    rangeVec[points <= 0 & rangeVec == 0] <- ranges[i]
  }
  
  lineIntegralGridLength <- dim(lineIntegralGrid)[2]
  completionGridLength <- dim(completionGrid)[2]
  sphereRangeGridLength <- dim(sphereRangeGrid)[2]
  sphereDistGridLength <- dim(sphereDistGrid)[2]
  #evaluating all the splines
  lineIntegralDelta <- 1/lineIntegralGridLength
  completionDelta <- 1/completionGridLength
  sphereRangeDelta <- (pi/2)/sphereRangeGridLength
  sphereDistDelta <- 1/sphereDistGridLength
  lineIntegrals <- .Fortran("CubicInterp", grid = lineIntegralGrid, nGrid = lineIntegralGridLength, delta = lineIntegralDelta,
                  points = dists, nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
  starts <- .Fortran("CubicInterp", grid = completionGrid, nGrid = completionGridLength, delta = completionDelta,
                  points = completions[1,], nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
  stops <- .Fortran("CubicInterp", grid = completionGrid, nGrid = completionGridLength, delta = completionDelta,
                     points = completions[2,], nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
  sphereCorrectionR <- .Fortran("CubicInterp", grid = sphereRangeGrid, nGrid = sphereRangeGridLength, delta = sphereRangeDelta,
                  points = rangeVec, nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
  sphereCorrectionD <- .Fortran("CubicInterp", grid = sphereDistGrid, nGrid = sphereDistGridLength, delta = sphereDistDelta,
                  points = dists, nPoints = N, output = rep(0, N), PACKAGE = "LatticeKrigInverseProblem")$output
  
  return((lineIntegrals - sphereCorrectionR * sphereCorrectionD) * (stops - starts) * rangeVec)
}

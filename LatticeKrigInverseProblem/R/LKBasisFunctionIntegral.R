LKBasisFunctionIntegral <- function(normDists) {
  dists <- as.double(normDists)
  
  if(!exists("lineIntegralGrid")) {
    data(LineIntegralSplines)
  }
  nGrid = dim(lineIntegralGrid)[2]
  lineIntegralDelta <- 1/nGrid;
  out <- .Fortran("CubicInterp", grid = lineIntegralGrid, nGrid = nGrid, delta = lineIntegralDelta,
                  points = dists, nPoints = length(dists), output = rep(0, length(dists)), PACKAGE = "LatticeKrigInverseProblem")
  return(out$output)
}

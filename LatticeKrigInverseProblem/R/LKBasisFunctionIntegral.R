LKBasisFunctionIntegral <- function(normDists) {
  dists <- as.double(normDists)
  
  if(!exists("lineIntegralGrid")) {
    data(CubicLineIntegralGrid)
  }
  
  lineIntegralDelta <- 1/dim(lineIntegralGrid)[2]
  out <- .Fortran("CubicInterp", grid = lineIntegralGrid, nGrid = length(lineIntegralGrid), delta = lineIntegralDelta,
                  points = dists, nPoints = length(dists), output = rep(0, length(dists)), PACKAGE = "LatticeKrigInverseProblem")
  return(out$output)
}

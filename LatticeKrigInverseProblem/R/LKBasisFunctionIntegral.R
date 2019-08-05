LKBasisFunctionIntegral <- function(normDists) {
  dists <- as.double(normDists)
  
  # #loading in the integral grid, downloading it from github if needed
  # if(!exists("lineIntegralGrid")) {
  #   if(!file.exists("LineIntegralGrid.rda")) {
  #     tryCatch({
  #       download.file("https://github.com/NCAR/LatticeKrig/raw/master/TomographyIntegralTable/LineIntegralGrid.rda", "LineIntegralGrid.rda")
  #       print("Retrieved table from GitHub")
  #     }, error = function(e) {
  #       stop("Integral table not found in working directory and downloading failed (https://github.com/NCAR/LatticeKrig/raw/master/TomographyIntegralTable/LineIntegralGrid.rda)")
  #     })
  #   }
  #   load("LineIntegralGrid.rda", envir = globalenv())
  #   #print("Loaded table successfully")
  # }
  
  #out <- .Fortran("interp", grid = lineIntegralGrid, nGrid = length(lineIntegralGrid), delta = lineIntegralDelta,
  #         points = dists, nPoints = length(dists), output = rep(0, length(dists)), PACKAGE = "LatticeKrig")
  
  if(!exists("lineIntegralGrid")) {
    data(CubicLineIntegralGrid)
  }
  lineIntegralDelta <- 1/dim(lineIntegralGrid)[2]
  out <- .Fortran("CubicInterp", grid = lineIntegralGrid, nGrid = length(lineIntegralGrid), delta = lineIntegralDelta,
                  points = dists, nPoints = length(dists), output = rep(0, length(dists)), PACKAGE = "LatticeKrigInverseProblem")
  return(out$output)
}

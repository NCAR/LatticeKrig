LKBasisFunctionIntegral <- function(normDists) {
  dims <- dim(normDists)
  dists <- matrix(normDists, nrow = 1)
  
  #loading in the integral grid, downloading it from github if needed
  if(!exists("lineIntegralGrid")) {
    if(!file.exists("LineIntegralGrid.rda")) {
      tryCatch({
        download.file("https://github.com/NCAR/LatticeKrig/raw/master/TomographyIntegralTable/LineIntegralGrid.rda", "LineIntegralGrid.rda")
        print("Retrieved table from GitHub")
      }, error = function(e) {
        stop("Integral table not found in working directory and downloading failed (https://github.com/NCAR/LatticeKrig/raw/master/TomographyIntegralTable/LineIntegralGrid.rda)")
      })
    }
    load("LineIntegralGrid.rda")
    print("Loaded table successfully")
  }
  
  out <- .Fortran("interp", grid = lineIntegralGrid, nGrid = length(lineIntegralGrid), delta = lineIntegralDelta,
           points = dists, nPoints = length(dists), output = rep(0, length(dists)), PACKAGE = "LatticeKrig")
  
  return(out$output)
}
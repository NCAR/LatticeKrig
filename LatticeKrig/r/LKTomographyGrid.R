LKTomographyGrid <- function(lines, points, ranges, rangeReps) {
  lines <- as.matrix(lines)
  points <- as.matrix(points)
  ranges <- as.matrix(ranges)
  dims <- as.integer(dim(points)[1])
  nPoints <- as.integer(dim(points)[2])
  nLines <- as.integer(dim(lines)[2])
  nRanges <- as.integer(length(ranges))
  rangeReps <- as.integer(rangeReps)
  lines <- as.double(lines)
  points <- as.double(points)
  ranges <- as.double(ranges)
  
  start <- Sys.time()
  output <- .Fortran("LKTomGridCount", dim=dims, points=points, nPoints = nPoints,
                     lines=lines, nLines = nLines, ranges=ranges, rangeReps=rangeReps,
                     nRanges=nRanges, nEntries = 0L, PACKAGE = "LatticeKrig")
  stop <- Sys.time()
  delta <- stop-start
  cat("Counting done in", as.double(delta), "seconds\n")
  nEntries <- output$nEntries
  ind <- as.integer(matrix(0, nrow=2, ncol = nEntries))
  entries <- as.double(rep(0, nEntries))
  start <- Sys.time()
  output <- .Fortran("LKTomGrid", dim=dims, points=points, nPoints = nPoints,
                     lines=lines, nLines = nLines, ranges=ranges, rangeReps=rangeReps,
                     nRanges=nRanges, ind=ind, entries=entries, nEntries = nEntries, PACKAGE = "LatticeKrig")
  stop <- Sys.time()
  delta <- stop-start
  cat("Storage done in", as.double(delta), "seconds\n")
  ra = output$entries
  ind = matrix(output$ind, nrow = 2)
  return(list(ind = t(ind), da = c(nLines, nPoints), ra = ra))
}
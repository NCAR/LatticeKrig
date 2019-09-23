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
  
  output <- .Fortran("LKTomGridCount", dim=dims, points=points, nPoints = nPoints,
                     lines=lines, nLines = nLines, ranges=ranges, rangeReps=rangeReps,
                     nRanges=nRanges, nEntries = 0L, PACKAGE = "LatticeKrigInverseProblem")
  
  nEntries <- output$nEntries
  ind <- as.integer(matrix(0, nrow=2, ncol = nEntries))
  entries <- as.double(rep(0, nEntries))
  
  output <- .Fortran("LKTomGrid", dim=dims, points=points, nPoints = nPoints,
                     lines=lines, nLines = nLines, ranges=ranges, rangeReps=rangeReps,
                     nRanges=nRanges, ind=ind, entries=entries, nEntries = nEntries, PACKAGE = "LatticeKrigInverseProblem")
  
  ra = output$entries
  ind = t(matrix(output$ind, nrow = 2))
  ord <- order(ind[,1], ind[,2])
  ind <- ind[ord,]
  ra <- ra[ord]
  return(list(ind = ind, da = c(nLines, nPoints), ra = ra))
}

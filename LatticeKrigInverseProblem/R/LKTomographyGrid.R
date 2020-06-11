LKTomographyGrid <- function(lines, points, ranges, rangeReps, weights, lineIntegralMode = 3L) {
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
  completions <- as.double(matrix(0, nrow=2, ncol=nEntries))
  
  output <- .Fortran("LKTomGrid", dim=dims, points=points, nPoints = nPoints,
                     lines=lines, nLines = nLines, ranges=ranges, rangeReps=rangeReps,
                     nRanges=nRanges, ind=ind, entries=entries, completions=completions,
                     nEntries = nEntries, PACKAGE = "LatticeKrigInverseProblem")
  
  com <- matrix(output$completions, nrow = 2)
  ind <- t(matrix(output$ind, nrow = 2))
  ra <- LKBasisFunctionIntegral(output$entries, com, ranges, rangeReps, ind[,2], weights, lineIntegralMode)
  ord <- order(ind[,1], ind[,2])
  ind <- ind[ord,]
  ra <- ra[ord]
  return(list(ind = ind, da = c(nLines, nPoints), ra = ra))
}

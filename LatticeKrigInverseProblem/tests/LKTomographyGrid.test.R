options(echo=FALSE, width=999)
suppressWarnings(suppressMessages(library(LatticeKrigInverseProblem)))
lines = matrix(c(0, 0, 1, 1, 0, 1, 1, 0), nrow = 4)
points = matrix(c(0, 0.5, 0.5, 0, 2, 1), nrow = 2)
ranges = 1
rangeReps = 3
out <- LKTomographyGrid(lines, points, ranges, rangeReps)
spind2full(out)
LKBasisFunctionIntegral(out$ra)
LKBasisFunctionIntegral(seq(0,0.95,,6))

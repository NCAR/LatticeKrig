LKTomography <- function(lines, obs, points, LKinfo) {
  ranges <- LKinfo$basisInfo$overlap * LKinfo$latticeInfo$delta
  rangeReps <- LKinfo$latticeInfo$mLevel
  points <- LKTomographyPoints(LKinfo$latticeInfo$grid)
  
  tomMatrix <- LKTomographyGrid(lines, points, ranges, rangeReps)
  tomMatrix$ra <- LKBasisFunctionIntegral(tomMatrix$ra)
  idx = 1L
  for(i in 1:length(ranges)) {
    for(j in 1:rangeReps[i]) {
      tomMatrix$ra[idx] = tomMatrix$ra[idx] * ranges[i]
      idx = idx + 1
    }
  }
  kFit <- LatticeKrig(x=points, y=obs, X = tomMatrix)
  
}
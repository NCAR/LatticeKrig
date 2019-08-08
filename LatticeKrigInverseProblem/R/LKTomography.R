LKTomography <- function(lines, obs, LKinfo) {
  ranges <- LKinfo$basisInfo$overlap * LKinfo$latticeInfo$delta
  rangeReps <- LKinfo$latticeInfo$mLevel
  points <- LKTomographyPoints(LKinfo$latticeInfo)
  
  tomMatrix <- LKTomographyGrid(lines, points, ranges, rangeReps)
  tomMatrix$ra <- LKBasisFunctionIntegral(tomMatrix$ra)
  idx = 1L
  for(i in 1:length(ranges)) {
    for(j in 1:rangeReps[i]) {
      tomMatrix$ra[idx] = tomMatrix$ra[idx] * ranges[i]
      idx = idx + 1
    }
  }
  N <- length(obs)
  kFit <- LKrig(x=points, y=obs, U=cbind(rep(1, N), runif(N)), X=spind2spam(tomMatrix), LKinfo=LKinfo, lambda=0.05)
  return(kFit)
}
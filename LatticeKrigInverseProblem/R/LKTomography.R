LKTomography <- function(lines, obs, LKinfo) {
  ranges <- LKinfo$basisInfo$overlap * LKinfo$latticeInfo$delta
  rangeReps <- LKinfo$latticeInfo$mLevel
  points <- LKTomographyPoints(LKinfo$latticeInfo)
  dimension <- dim(points)[1]
  
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
  #TODO create a more reasonable U matrix; ideally everything but the first column should be 0, so there is no fixed effect, but I don't know how to
  #do that in LatticeKrig without creating a singular matrix
  kFit <- LKrig(x=t(points), y=obs, U=cbind(rep(1, N), matrix(runif(N*dimension), nrow = N)), X=spind2spam(tomMatrix), LKinfo=LKinfo, lambda=0.05)
  return(kFit)
}
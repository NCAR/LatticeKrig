LKTomography <- function(lines, obs, LKinfo) {
  ranges <- LKinfo$basisInfo$overlap * LKinfo$latticeInfo$delta
  rangeReps <- LKinfo$latticeInfo$mLevel
  points <- LKTomographyPoints(LKinfo$latticeInfo)
  dimension <- dim(points)[1]
  if (class(LKinfo)[2] == "LKSphere") {
    tomMatrix <- LKTomographyGridSphere(lines, points, ranges, rangeReps)
  } else {
    tomMatrix <- LKTomographyGrid(lines, points, ranges, rangeReps, LKinfo$alpha)
  }
  N <- length(obs)
  kFit <- LKrig(x=t(points), y=obs, X=spind2spam(tomMatrix), LKinfo=LKinfo, lambda=0.05)
  return(kFit)
}
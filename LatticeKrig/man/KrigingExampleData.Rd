\name{KrigingExampleData}
\alias{KrigingExampleData}
\docType{data}
\title{Synthetic data for kriging examples}
\description{
This is a sample of 50 (x,y) ordered pairs: the x values are in sorted order,
drawn from the interval [0,1]. The smallest x value is exactly 0, and the largest is exactly 1.
The y values are the function y=9*x*(1-x)^3 evaluated at each x value, with a small measurement error (normal with mean 0, standard deviation 0.01) added on.
}
\usage{data(KrigingExampleData)}
\format{The data in \code{KrigingExampleData} is listed in two vectors, \code{x} and \code{y}.}
\examples{
\dontrun{
data(KrigingExampleData)
plot(x,y)
kFit <- LatticeKrig(x, y)
xGrid <- seq(0,1,0.001)
lines(xGrid, predict(kFit, xGrid), col='red')
}
}

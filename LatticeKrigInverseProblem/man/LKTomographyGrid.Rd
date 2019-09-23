\name{LKTomographyGrid}
\alias{LKTomographyGrid}
\title{Compute a sparse matrix of line-point distances less than a threshold}

\description{
This function computes the matrix of line point distances from a list of points and a list of lines.
There is also a range (or list of ranges) specified: distances from a line to a point that are outside
of that point's range are not recorded, and distances below that range are divided by that range to be
in the interval [0,1]. The output of this function is recorded in the spind format, meaning we record
the row, column, and value of every nonzero matrix entry. The lines and points are formatted in columns
as ordered pairs; lines are specified by a column of two points on top of each other. For example, 100
points in 2D space would take up a 2x100 matrix, and 250 lines in 3D space would take up a 6x250 matrix.
The ranges are specified for each point; to save space, one can specify the number of times to repeat a
range entry. For example, ranges=c(1, 0.5, 0.25), rangeReps = c(10, 50, 100) would mean the first 10
points have range 1, the next 50 points have range 0.5, and the next 100 points have range 0.25.
}

\usage{
LKTomographyGrid(lines, points, ranges, rangeReps)
}

\arguments{
\item{lines}{The list of lines, specified in a matrix where each column has two sets of points.}
\item{points}{The list of points, specified in a matrix where each column contains the coordinates of a point.}
\item{ranges}{The ranges of each group of points, in order.}
\item{rangeReps}{The number of times to repeat each range.}
}

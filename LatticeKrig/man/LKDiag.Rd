\name{LKDiag}
\alias{LKDiag}
\title{
Create a matrix with given entries on the given diagonals
}

\description{
This function builds a matrix in which each entry on a given diagonal has the same value. 
The matrix can be any size, not necessarily square, and the user can specify arbitrary diagonals to fill. 
If \code{ncol} is not set, the output will be a square matrix. If \code{diags} is not set, the entries 
will be placed as close as possible to the main diagonal; e.g. passing in one value for \code{entries} 
will produce a diagonal matrix, and passing in a vector of 3 values will produce a tridiagonal matrix. 
If an even number of entries are provided, they will be placed symmetrically around the main diagonal, 
leaving the main diagonal as zeroes.
}

\usage{
LKDiag(entries, nrow, diags = NULL, ncol = nrow,
full = FALSE)
}

\arguments{
\item{entries}{The values to put on each diagonal.}
\item{nrow}{The number of rows in the matrix; if \code{ncol} isn't set, the matrix will be square with
\code{nrow} rows and columns.}
\item{diags}{The diagonals to put the entries on. 0 corresponds to the main diagonal, positive values 
go above the main diagonal, and negatives go below.}
\item{ncol}{The number of columns in the matrix, if it's different from the number of rows.}
\item{full}{If TRUE return a dense matrix in the usual R matrix format. If FALSE return a sparse
matrix in spam format. }

}

\value{
\strong{LKDiag} returns a (dense) matrix with the given diagonals.
}

\author{
Matt Iverson
}

\examples{
#produces a 5x5 identity matrix
LKDiag(1, 5)

#produces a 6x6 tridiagonal matrix
LKDiag(c(1, -2, 1), 6)
LKDiag(c(1, -2, 1), 6, diags=(-1:1))

#produces a 5x5 matrix with 3 on the main diagonal and the corners
n <- 5
LKDiag(3, n, diags=c(0, n-1, -(n-1)))

# each of the following produces a 4x6 matrix with a 2 in the diagonal 
# below the main diagonal and 5 in the diagonal  above it
LKDiag(c(2, 5), nrow = 4, ncol = 6)
LKDiag(c(2, 5), nrow = 4, ncol = 6, diags = c(-1, 1))
LKDiag(c(5, 2), nrow = 4, ncol = 6, diags = c(1, -1))
}

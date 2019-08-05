\name{LKBasisFunctionIntegral}
\alias{LKBasisFunctionIntegral}
\title{Computes the corresponding line integral of a vector of offsets on the 2D Wendland function}

\description{
This function computes the line integral across a line passing a distance $d$ from the Wendland function
in two dimensions (or any higher number of dimensions, since it's radially symmetric). Symbolically, if
f(x) is the Wendland function, LKBasisFunctionIntegral(d) evaluates to 
$\int_{-\infty}^{\infty} f(x^2 + y^2) dx$. This assumes the Wendland function has radius 1; to use this
function for a different radius, divide the distance between the center of the function and the line by
the radius, compute the line integral with this function, and multiply that value by the radius. 
}

\usage{LKBasisFunctionIntegral(normDists)}

\arguments{
\item{normDists}{The distances from the basis function to evaluate the line integrals at, normalized by the basis function radius as described above.}
}

\author{Matt Iverson}



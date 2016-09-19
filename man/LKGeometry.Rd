% # LatticeKrig  is a package for analysis of spatial data written for
% # the R software environment .
% # Copyright (C) 2016
% # University Corporation for Atmospheric Research (UCAR)
% # Contact: Douglas Nychka, nychka@ucar.edu,
% # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
% #
% # This program is free software; you can redistribute it and/or modify
% # it under the terms of the GNU General Public License as published by
% # the Free Software Foundation; either version 2 of the License, or
% # (at your option) any later version.
% # This program is distributed in the hope that it will be useful,
% # but WITHOUT ANY WARRANTY; without even the implied warranty of
% # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% # GNU General Public License for more details.

\name{LKGeometry}
\alias{LKGeometry}
\alias{LKBox}
\alias{LKRectangle}
\alias{LKInterval}
\title{Tailoring the LatticeKrig model to a specific geometry.
}
\description{To implement a spatial model for specific geometry, e.g. 
	a rectangular spatial domain, one needs to specify a few functions
	that set reasonable default choices, create the lattice, define the Markov
	random field and
	define the basis functions. One coudl also include some specific functions
	that are particularly
	efficient for some of the computations for that geometry. This help section
	gives an overview of
	what is needed to introduce a new geometry. }
	
	\section{ A quick start}{
		At the outset if you are an example person (and an impatient one) take a look at the functions
		created for the 1-d LatticeKrig model (\code{LKInterval}).
		(In the source for the package this is the file ModelLKInterval.R in the
		R subdirectory.) The specfic functions are:
		\describe{
     	\item{LKrigSetupLattice.LKInterval}{Creates a 1-d lattice.}
		\item{LKrigSAR.LKInterval}{Creates a sparse 1-d spatial autoregressive matrix.}
		\item{LKrigLatticeCenters.LKInterval}{Returns the spatial locations of the lattice points.}
	}
	Nothing else is needed to coax the LatticeKrig computation to fit a 1-d spatial model. There are more details and flexibility in specifying the model, but the ones listed above are generic in that they make sense for any LatticeKrig model and are not specific to the 1-d case.  
The \code{.LKInterval} tag on each of these functions indicates that these are S3 versions of a "method". The actual source code in LatticeKrig just relies on calling the generic 
method, e.g. LKrigSetupLattice, and the geometry, added as a class indicates, which function is actually used. Checking what each of these functions does and then 
seeing how they are called in \code{LKrigSetup} and \code{LKrig.basis} will give a simple introduction to how the code is structured. 
}

\section{Use of the LKInfo object}{
The basic idea behind this package is to consolidate all the
details of the spatial model including the geometry in the \code{LKinfo} object.
The way to create this object is by the function \code{SetupLKinfo} and we believe that this function has enough flexibility to avoid modifying this object outside of this function call. Typically for deliberate spatial modeling the first step will be creating the LKinfo object; however, the  quick and more "black box" function LatticeKrig is designed to take a minimal amount of information and then 
calls \code{SetupLKinfo} internally to fill out the LKinfo object. 
In fact, LatticeKrig( x,y) with just \code{x} as locations and \code{y} as observations will fit a spatial model using reasonable defaults  and the returned object from this fit will include the full LKinfo object that it used.  

A subtle and perhaps a potential disadvantage with the coding is that the LKinfo object is built in steps within the LKrigSetup function. In particular the method that setups the lattice is called before information is added about the alpha and
a.wght parameters. The source code LKrigSetup has been kept concise and it should be clear how the LKinfo is formed. Conceptually, the sequence of steps are:
\describe{
     	\item{The calling arguments}{Initialize the LKinfo object with the
     		 all the passed arguments.}
     	\item{Call: setDefaultsLKinfo}{If necessary add some components that fill in 
     		specific default parameters for a given method.}
        \item{Create basis information}{The component \code{ LKinfo\$basisInfo} holds the specific components for the basis functions. Currently there is not much flexibility here are as it is just derived from the calling arguments.}
     	\item{Call: LKrigSetupLattice}{Setup the lattice information and put in the 
     		component \code{LKinfo\$latticeInfo} }
         \item{Call: LKrigSetupAlpha}{Setup the format for the alpha parameters based on 
         	what was passed in and out this in the component \code{ LKinfo\$alpha}}
         \item{Call: LKrigSetupAwght}{Setup the format for the alpha parameters based on 
         	what was passed in and out this in the component \code{ LKinfo\$a.wght}}
         \item{Check using: LKinfoCheck}{Check for the required components in the object.}	
 }

 The LKinfo object is a list with  
class \code{LKinfo} and comes with a print method to make it easier to  read and also has the class of the geometry, e.g \code{LKInterval} for the 1-d model. The functions that are required for a geometry all take the LKinfo object as their first argument and the S3 object dispatching then calls the correct version. 
Put more simply, if  \code{LKinfoExample} has class \code{LKInterval} then the code
\code{LKrigLatticeCenters(LKinfoExample, Level= 1)} will result in computing the 
lattice centers at the first level and using the 1-d geometry.
In this way much of this package can use "generic" functions for the 
computational steps without revealing the extra complications and 
exceptions for particular models and geometries. 
The structure of the \code{LKinfo} object is not completely rigid. Certain components of the list are required and these are checked  by \code{LKinfoCheck}. Other components of this list can be geometry and problem dependent and this is a useful way to pass information to the computations. In general when faced with what should be included in LKinfo follow the style in \code{LKInterval} but also include
additional components that are needed to define the model.  Some examples are mentioned below. 
}

\section{Required methods}{
	For any model there are three functions required for a new geometry. By a new geometry we mean any variation from the regular grids and beyond  a
second order SAR. As a first step choose a name for the new geometry, e.g. 
MyNewIdea, and this should also be the name passed as the argument LKGeometry in the 
\code{LKrigSetup} function. Three new functions need to be supplied and that end in .MyNewIdea and the user will recognize these as S3 methods
being added to LatticeKrig.

The specific functions you will have to add to LatticeKrig package are:

\strong{LKrigSetupLattice.MyNewIdea}(LKinfo) This function will take the information in the LKinfo list and create any indexing and do other computations to create the specific lattice of spatial points. The  return value is a list with several basic
components that are required.  For example \code{m} is the total number of basis functions in the model.  See help( LKrigSetupLattice) for a description of these required arguments. In addition one can add other arguments that help in coding the following 
two functions. It possible that this function could be avoided by doing all the
setup whenever the SAR and  the centers functions are used but we suspect this would result in more complicated code and be less modular. 
The 1-d implementation in the source file ModelInterval.R is a  good example to get an idea of the whole process.  

\strong{LKrigSAR.MyNewIdea}(LKinfo, Level) Returns a sparse spatial autoregressive matrix at a specific level of the multiresolution and in \code{spind} sparse matrix format. See \code{help( LKrigSAR)}. Note that the lattice information from the previous function is in the component \code{latticeInfo} and so this function needs to 
access any lattice information in this way from the LKinfo object. e.g. 
\code{ LKinfo\$latticeInfo\$m } are the total number of basis functions.
				
\strong{LKrigLatticeCenters.MyNewIdea}(LKinfo, Level)  Returns the spatial locations of the lattice points at a specific level of the  multiresolution. This can either be as a matrix with columns indexing
dimension and rows indexing points or another class of object. These locations (or object) are subsequently passed to the second argument of the distance function \code{LKrigDistance} to create the basis functions.

For initially creating a new geometry one can just rely on this 
function returning the lattice locations. Once the whole geometry 
works one can consider returning a more specialized object and pairing these with more efficient distance finding algorithms. 
For example if the 
lattice points are on a regular, equally spaced grid the nearest neighbor locations can be found quickly just by simple arithmetic and this is much faster than a search. This is implemented for the
LKInterval, LKBox, and LKRectangle geometries by creating the class \code{gridList} for the grid of lattice centers and passing this object directly into the second argument of the distance function
LKrigDistance. See also \code{help(LKrigLatticeCenters)} and 
\code{help(LKrigDistance)}
}
		
\section{Optional Methods}{
	To make the modeling more flexible and take advantage of 
faster computation some additional methods can be added but are not required. Any new version should add your geometry name to the end.
e.g.  \code{setDefaultsLKinfo} should have the 
new name:  \code{setDefaultsLKinfo.MyNewIdea} and have the default arguments listed in the generic function (\code{setDefaultsLKinfo}).

\strong{setDefaultsLKinfo} This function takes the skeleton LKinfo list based on the arguments passed to LKrigSetup and can modify and add to them. The returned value is now taken as the initial LKinfo object. The default method is that no changes are made. But adding versions to this method can be very useful for  specific cases (e.g. see \code{setDefaultsLKinfo}) without complicating the main function with lots of conditional statements. 

\strong{LKrigSetupAlpha} These methods convert the alpha  information passed into 
LKrigSetup to a list format. A default is provided if you do not supply a specific version.  See help(LKrigSetupAlpha).

\strong{LKrigSetupAwght} These methods convert the a.wght information passed into 
LKrigSetup to a list format.

\strong{LKrigNormalizeBasisFast} This is a method to implement any fast methods for 
normalizing the basis functions. There is no default method because it is assumed this is geometry dependent. See \code{\link{LKrigNormalizeBasisFast.LKRectangle}}
as an example of how this is done. Note that there is already a "slow" method, 
\code{\link{LKrigNormalizeBasis}} involving a Cholesky solve of the precision matrix for each point where the basis function is evaluated. 

\strong{LKinfoCheck} This function checks the LKinfo object for necessary components and does some rudimentary comparisons to make sure arrays are the right size. It is always good to add more checks!

}

\author{
	Doug  Nychka
}


\seealso{
    \code{\link{LKrigSetup}},
	\code{\link{LKrigSAR}},
	\code{\link{LKrigLatticeCenters}}
}


\keyword{ spatial }
 
	


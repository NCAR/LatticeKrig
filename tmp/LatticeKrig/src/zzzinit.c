#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

/*
          subroutine   findNorm(mx,my,offset,Ux,Dx,Uy,Dy,nLocations, xLocations, weights, Z)
LKDist.R: #  subroutine               lkdist( x1, n1, x2, n2, dim,  delta2, ind,    rd,  Nmax, iflag)   
LKDistComponents.R: #  subroutine lkdistComp( x1, n1, x2, n2, dim,   delta, ind,    rd,  Nmax, iflag)
LKDistGrid.R:    # subroutine lkdistgrid(     x1, n1, nGrid, nDim,   delta, irow, jcol,   ra, Nmax, iflag)
LKDistGridComponents.R:#   subroutine lkdistgridcomp( 
                                              x1, n1, nGrid, nDim,   delta, irow, jcol,   ra, Nmax, iflag)
subroutine LKDiag(entries, nEntries, diags, nRow, nCol, matrix)
*/

/* .Fortran calls */
extern void F77_NAME(findnorm)(      void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,  void *);
extern void F77_NAME(lkdist)(        void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistcomp)(    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgrid)(    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgridcomp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdiag)(        void *, void *, void *, void *, void *, void *);


SEXP LinePointDistMat(SEXP dim_, SEXP points_, SEXP nPoints_, SEXP lines_, SEXP nLines_, SEXP ranges_, SEXP rangeReps_, SEXP nRanges_);

static const R_FortranMethodDef FortranEntries[] = {
    {"findnorm",       (DL_FUNC) &F77_NAME(findnorm),       11},
    {"lkdist",         (DL_FUNC) &F77_NAME(lkdist),         10},
    {"lkdistcomp",     (DL_FUNC) &F77_NAME(lkdistcomp),     10},
    {"lkdistgrid",     (DL_FUNC) &F77_NAME(lkdistgrid),     10},
    {"lkdistgridcomp", (DL_FUNC) &F77_NAME(lkdistgridcomp), 10},
    {"lkdiag",         (DL_FUNC) &F77_NAME(lkdiag),          6},
    {NULL, NULL, 0}
};

void R_init_LatticeKrig(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

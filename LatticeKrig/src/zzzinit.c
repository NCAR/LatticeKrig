#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(findnorm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdist)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistcomp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgrid)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgridcomp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdiag)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(interp)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgridcount)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgrid)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);


SEXP LinePointDistMat(SEXP dim_, SEXP points_, SEXP nPoints_, SEXP lines_, SEXP nLines_, SEXP ranges_, SEXP rangeReps_, SEXP nRanges_);

static const R_FortranMethodDef FortranEntries[] = {
    {"findnorm",       (DL_FUNC) &F77_NAME(findnorm),       11},
    {"lkdist",         (DL_FUNC) &F77_NAME(lkdist),         10},
    {"lkdistcomp",     (DL_FUNC) &F77_NAME(lkdistcomp),     10},
    {"lkdistgrid",     (DL_FUNC) &F77_NAME(lkdistgrid),     11},
    {"lkdistgridcomp", (DL_FUNC) &F77_NAME(lkdistgridcomp), 11},
    {"lkdiag",         (DL_FUNC) &F77_NAME(lkdiag),          6},
    {"interp",         (DL_FUNC) &F77_NAME(interp),          6},
    {"lktomgridcount", (DL_FUNC) &F77_NAME(lktomgridcount),  9},
    {"lktomgrid",      (DL_FUNC) &F77_NAME(lktomgrid),      11},
    {NULL, NULL, 0}
};

void R_init_LatticeKrig(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

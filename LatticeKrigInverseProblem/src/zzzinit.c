#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cubicinterp)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgridcount)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgrid)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgridspherecount)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgridsphere)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lktomgridspheremock)(void *, void *);
extern void F77_NAME(lktompoints3d)(void *, void *, void *, void *, void *);
extern void F77_NAME(lktompoints2d)(void *, void *, void *, void *, void *);


SEXP LinePointDistMat(SEXP dim_, SEXP points_, SEXP nPoints_, SEXP lines_, SEXP nLines_, SEXP ranges_, SEXP rangeReps_, SEXP nRanges_);

static const R_FortranMethodDef FortranEntries[] = {
    {"cubicinterp",          (DL_FUNC) &F77_NAME(cubicinterp),          6},
    {"lktomgridcount",       (DL_FUNC) &F77_NAME(lktomgridcount),       9},
    {"lktomgrid",            (DL_FUNC) &F77_NAME(lktomgrid),           11},
    {"lktomgridspherecount", (DL_FUNC) &F77_NAME(lktomgridspherecount), 8},
    {"lktomgridsphere",      (DL_FUNC) &F77_NAME(lktomgridsphere),     11},
    {"lktomgridspheremock",  (DL_FUNC) &F77_NAME(lktomgridspheremock),  2},
    {"lktompoints3d",        (DL_FUNC) &F77_NAME(lktompoints3d),        5},
    {"lktompoints2d",        (DL_FUNC) &F77_NAME(lktompoints2d),        5},
    {NULL, NULL, 0}
};

void R_init_LatticeKrig(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

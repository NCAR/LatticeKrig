#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(findnorm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdist)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistcomp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgrid)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lkdistgridcomp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"findnorm",       (DL_FUNC) &F77_NAME(findnorm),       11},
    {"lkdist",         (DL_FUNC) &F77_NAME(lkdist),         10},
    {"lkdistcomp",     (DL_FUNC) &F77_NAME(lkdistcomp),     10},
    {"lkdistgrid",     (DL_FUNC) &F77_NAME(lkdistgrid),     11},
    {"lkdistgridcomp", (DL_FUNC) &F77_NAME(lkdistgridcomp), 11},
    {NULL, NULL, 0}
};

void R_init_LatticeKrig(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

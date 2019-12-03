#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void checkopenmp(void *);
extern void Tdemix(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"checkopenmp", (DL_FUNC) &checkopenmp,  1},
    {"Tdemix",      (DL_FUNC) &Tdemix,      24},
    {NULL, NULL, 0}
};

void R_init_DeMixT(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

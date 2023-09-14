#include <R_ext/RS.h>                                                                                                                                                                         
#include <stdlib.h> // for NULL                                                                                                                                                               
#include <R_ext/Rdynload.h>                                                                                                                                                                   

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(r_dmzip)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dmzip_init_param)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dmzipt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dmzipt_init_param)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"r_dmzip",             (DL_FUNC) &F77_NAME(r_dmzip),             17},
  {"r_dmzip_init_param",  (DL_FUNC) &F77_NAME(r_dmzip_init_param),  12},
  {"r_dmzipt",            (DL_FUNC) &F77_NAME(r_dmzipt),            15},
  {"r_dmzipt_init_param", (DL_FUNC) &F77_NAME(r_dmzipt_init_param), 10},
  {NULL, NULL, 0}
};

void R_init_crimCV(DllInfo *dll)
{
     R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
     R_useDynamicSymbols(dll, FALSE);
}

#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void SurvBeranKernel(void *, void *,void *,void *,void *,void *,void *,void *,void *,
                     void *);

extern void SurvBeranNNE(void *, void *,void *,void *,void *,void *,void *,void *,void *,
                     void *);

extern void LLWeightsKernel(void *, void *,void *,void *,void *,void *);

extern void NWWeightsKernel(void *, void *,void *,void *,void *,void *);



static const R_CMethodDef CEntries[] = {
    {"SurvBeranKernel", (DL_FUNC) &SurvBeranKernel, 10},
    {"SurvBeranNNE", (DL_FUNC) &SurvBeranNNE, 10},
    {"LLWeightsKernel", (DL_FUNC) &LLWeightsKernel, 6},
    {"NWWeightsKernel", (DL_FUNC) &NWWeightsKernel, 6},
    {NULL, NULL, 0}
};

void R_init_CondTimeROC(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
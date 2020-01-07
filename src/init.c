#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
/* .Call calls */
extern SEXP _DeMixT_Alpha_search_2D(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Alpha_search_2D_sigma(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Alpha_search_MuT_2D(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Alpha_search_Pi_2D(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Alpha_search_SigmaT_2D(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1f0MuT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1f0Pi_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1f0SigmaT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1Loglikelihood_log_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1MuT_Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1MuTSigmaT_Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1Pi_Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D1SigmaT_Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0MuT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0MuTSigmaT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0Pi_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0PiMuT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0PiSigmaT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2f0SigmaT_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2Loglikelihood_log_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2Loglikelihood_unit_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2MuT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2MuT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2MuTSigmaT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2MuTSigmaT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2Pi_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2Pi_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2PiMuT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2PiMuT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2PiSigmaT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2PiSigmaT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2SigmaT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_D2SigmaT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DMuT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DMuT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DPi_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DPi_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DSigmaT_inner_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_DSigmaT_outer_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_f0_func_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_ft_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Gfunc_2D_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenLine_search_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenLine_search_MuT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenLine_search_Pi_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenLine_search_SigmaT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenSection_Loglikelihood_MuT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenSection_Loglikelihood_Pi_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_GoldenSection_Loglikelihood_SigmaT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Gt(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Gt_vec(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_inner_trapez_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_innerCPP_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_log_divide(SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_2D_L1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_ft_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_MuT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_Pi_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Loglikelihood_SigmaT_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_mint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_pf_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_pmin_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_SoftThreshold(SEXP, SEXP);
extern SEXP _DeMixT_SoftThreshold_vec(SEXP, SEXP);
extern SEXP _DeMixT_tf_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_tmin_y(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_tmin_y2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_Unit1(SEXP, SEXP, SEXP);
extern SEXP _DeMixT_x_update_2D(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DeMixT_x_update_inv_2D(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_DeMixT_Alpha_search_2D",                       (DL_FUNC) &_DeMixT_Alpha_search_2D,                       4},
    {"_DeMixT_Alpha_search_2D_sigma",                 (DL_FUNC) &_DeMixT_Alpha_search_2D_sigma,                 4},
    {"_DeMixT_Alpha_search_MuT_2D",                   (DL_FUNC) &_DeMixT_Alpha_search_MuT_2D,                   3},
    {"_DeMixT_Alpha_search_Pi_2D",                    (DL_FUNC) &_DeMixT_Alpha_search_Pi_2D,                    3},
    {"_DeMixT_Alpha_search_SigmaT_2D",                (DL_FUNC) &_DeMixT_Alpha_search_SigmaT_2D,                3},
    {"_DeMixT_D1f0MuT_func_2D",                       (DL_FUNC) &_DeMixT_D1f0MuT_func_2D,                       6},
    {"_DeMixT_D1f0Pi_func_2D",                        (DL_FUNC) &_DeMixT_D1f0Pi_func_2D,                        6},
    {"_DeMixT_D1f0SigmaT_func_2D",                    (DL_FUNC) &_DeMixT_D1f0SigmaT_func_2D,                    6},
    {"_DeMixT_D1Loglikelihood_2D",                    (DL_FUNC) &_DeMixT_D1Loglikelihood_2D,                    6},
    {"_DeMixT_D1Loglikelihood_log_2D",                (DL_FUNC) &_DeMixT_D1Loglikelihood_log_2D,                6},
    {"_DeMixT_D1MuT_Loglikelihood_2D",                (DL_FUNC) &_DeMixT_D1MuT_Loglikelihood_2D,                6},
    {"_DeMixT_D1MuTSigmaT_Loglikelihood_2D",          (DL_FUNC) &_DeMixT_D1MuTSigmaT_Loglikelihood_2D,          6},
    {"_DeMixT_D1Pi_Loglikelihood_2D",                 (DL_FUNC) &_DeMixT_D1Pi_Loglikelihood_2D,                 6},
    {"_DeMixT_D1SigmaT_Loglikelihood_2D",             (DL_FUNC) &_DeMixT_D1SigmaT_Loglikelihood_2D,             6},
    {"_DeMixT_D2f0MuT_func_2D",                       (DL_FUNC) &_DeMixT_D2f0MuT_func_2D,                       6},
    {"_DeMixT_D2f0MuTSigmaT_func_2D",                 (DL_FUNC) &_DeMixT_D2f0MuTSigmaT_func_2D,                 6},
    {"_DeMixT_D2f0Pi_func_2D",                        (DL_FUNC) &_DeMixT_D2f0Pi_func_2D,                        6},
    {"_DeMixT_D2f0PiMuT_func_2D",                     (DL_FUNC) &_DeMixT_D2f0PiMuT_func_2D,                     6},
    {"_DeMixT_D2f0PiSigmaT_func_2D",                  (DL_FUNC) &_DeMixT_D2f0PiSigmaT_func_2D,                  6},
    {"_DeMixT_D2f0SigmaT_func_2D",                    (DL_FUNC) &_DeMixT_D2f0SigmaT_func_2D,                    6},
    {"_DeMixT_D2Loglikelihood_2D",                    (DL_FUNC) &_DeMixT_D2Loglikelihood_2D,                    6},
    {"_DeMixT_D2Loglikelihood_log_2D",                (DL_FUNC) &_DeMixT_D2Loglikelihood_log_2D,                6},
    {"_DeMixT_D2Loglikelihood_unit_2D",               (DL_FUNC) &_DeMixT_D2Loglikelihood_unit_2D,               6},
    {"_DeMixT_D2MuT_inner_2D",                        (DL_FUNC) &_DeMixT_D2MuT_inner_2D,                        7},
    {"_DeMixT_D2MuT_outer_2D",                        (DL_FUNC) &_DeMixT_D2MuT_outer_2D,                        6},
    {"_DeMixT_D2MuTSigmaT_inner_2D",                  (DL_FUNC) &_DeMixT_D2MuTSigmaT_inner_2D,                  7},
    {"_DeMixT_D2MuTSigmaT_outer_2D",                  (DL_FUNC) &_DeMixT_D2MuTSigmaT_outer_2D,                  6},
    {"_DeMixT_D2Pi_inner_2D",                         (DL_FUNC) &_DeMixT_D2Pi_inner_2D,                         8},
    {"_DeMixT_D2Pi_outer_2D",                         (DL_FUNC) &_DeMixT_D2Pi_outer_2D,                         6},
    {"_DeMixT_D2PiMuT_inner_2D",                      (DL_FUNC) &_DeMixT_D2PiMuT_inner_2D,                      7},
    {"_DeMixT_D2PiMuT_outer_2D",                      (DL_FUNC) &_DeMixT_D2PiMuT_outer_2D,                      6},
    {"_DeMixT_D2PiSigmaT_inner_2D",                   (DL_FUNC) &_DeMixT_D2PiSigmaT_inner_2D,                   7},
    {"_DeMixT_D2PiSigmaT_outer_2D",                   (DL_FUNC) &_DeMixT_D2PiSigmaT_outer_2D,                   6},
    {"_DeMixT_D2SigmaT_inner_2D",                     (DL_FUNC) &_DeMixT_D2SigmaT_inner_2D,                     7},
    {"_DeMixT_D2SigmaT_outer_2D",                     (DL_FUNC) &_DeMixT_D2SigmaT_outer_2D,                     6},
    {"_DeMixT_DMuT_inner_2D",                         (DL_FUNC) &_DeMixT_DMuT_inner_2D,                         7},
    {"_DeMixT_DMuT_outer_2D",                         (DL_FUNC) &_DeMixT_DMuT_outer_2D,                         6},
    {"_DeMixT_DPi_inner_2D",                          (DL_FUNC) &_DeMixT_DPi_inner_2D,                          7},
    {"_DeMixT_DPi_outer_2D",                          (DL_FUNC) &_DeMixT_DPi_outer_2D,                          6},
    {"_DeMixT_DSigmaT_inner_2D",                      (DL_FUNC) &_DeMixT_DSigmaT_inner_2D,                      7},
    {"_DeMixT_DSigmaT_outer_2D",                      (DL_FUNC) &_DeMixT_DSigmaT_outer_2D,                      6},
    {"_DeMixT_f0_func_2D",                            (DL_FUNC) &_DeMixT_f0_func_2D,                            6},
    {"_DeMixT_ft_y",                                  (DL_FUNC) &_DeMixT_ft_y,                                  7},
    {"_DeMixT_Gfunc_2D_C",                            (DL_FUNC) &_DeMixT_Gfunc_2D_C,                            8},
    {"_DeMixT_GoldenLine_search_2D",                  (DL_FUNC) &_DeMixT_GoldenLine_search_2D,                  9},
    {"_DeMixT_GoldenLine_search_MuT_2D",              (DL_FUNC) &_DeMixT_GoldenLine_search_MuT_2D,              9},
    {"_DeMixT_GoldenLine_search_Pi_2D",               (DL_FUNC) &_DeMixT_GoldenLine_search_Pi_2D,               9},
    {"_DeMixT_GoldenLine_search_SigmaT_2D",           (DL_FUNC) &_DeMixT_GoldenLine_search_SigmaT_2D,           9},
    {"_DeMixT_GoldenSection_Loglikelihood_MuT_2D",    (DL_FUNC) &_DeMixT_GoldenSection_Loglikelihood_MuT_2D,    6},
    {"_DeMixT_GoldenSection_Loglikelihood_Pi_2D",     (DL_FUNC) &_DeMixT_GoldenSection_Loglikelihood_Pi_2D,     6},
    {"_DeMixT_GoldenSection_Loglikelihood_SigmaT_2D", (DL_FUNC) &_DeMixT_GoldenSection_Loglikelihood_SigmaT_2D, 6},
    {"_DeMixT_Gt",                                    (DL_FUNC) &_DeMixT_Gt,                                    3},
    {"_DeMixT_Gt_vec",                                (DL_FUNC) &_DeMixT_Gt_vec,                                3},
    {"_DeMixT_inner_trapez_2D",                       (DL_FUNC) &_DeMixT_inner_trapez_2D,                       6},
    {"_DeMixT_innerCPP_2D",                           (DL_FUNC) &_DeMixT_innerCPP_2D,                           7},
    {"_DeMixT_log_divide",                            (DL_FUNC) &_DeMixT_log_divide,                            2},
    {"_DeMixT_Loglikelihood_2D",                      (DL_FUNC) &_DeMixT_Loglikelihood_2D,                      6},
    {"_DeMixT_Loglikelihood_2D_L1",                   (DL_FUNC) &_DeMixT_Loglikelihood_2D_L1,                   7},
    {"_DeMixT_Loglikelihood_ft_y",                    (DL_FUNC) &_DeMixT_Loglikelihood_ft_y,                    6},
    {"_DeMixT_Loglikelihood_MuT_2D",                  (DL_FUNC) &_DeMixT_Loglikelihood_MuT_2D,                  7},
    {"_DeMixT_Loglikelihood_Pi_2D",                   (DL_FUNC) &_DeMixT_Loglikelihood_Pi_2D,                   7},
    {"_DeMixT_Loglikelihood_SigmaT_2D",               (DL_FUNC) &_DeMixT_Loglikelihood_SigmaT_2D,               7},
    {"_DeMixT_mint",                                  (DL_FUNC) &_DeMixT_mint,                                  6},
    {"_DeMixT_pf_y",                                  (DL_FUNC) &_DeMixT_pf_y,                                  7},
    {"_DeMixT_pmin_y",                                (DL_FUNC) &_DeMixT_pmin_y,                                9},
    {"_DeMixT_SoftThreshold",                         (DL_FUNC) &_DeMixT_SoftThreshold,                         2},
    {"_DeMixT_SoftThreshold_vec",                     (DL_FUNC) &_DeMixT_SoftThreshold_vec,                     2},
    {"_DeMixT_tf_y",                                  (DL_FUNC) &_DeMixT_tf_y,                                  7},
    {"_DeMixT_tmin_y",                                (DL_FUNC) &_DeMixT_tmin_y,                                8},
    {"_DeMixT_tmin_y2",                               (DL_FUNC) &_DeMixT_tmin_y2,                               9},
    {"_DeMixT_Unit1",                                 (DL_FUNC) &_DeMixT_Unit1,                                 3},
    {"_DeMixT_x_update_2D",                           (DL_FUNC) &_DeMixT_x_update_2D,                           5},
    {"_DeMixT_x_update_inv_2D",                       (DL_FUNC) &_DeMixT_x_update_inv_2D,                       3},
    {NULL, NULL, 0}
};


/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void checkopenmp(void *);
extern void Tdemix(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"checkopenmp", (DL_FUNC) &checkopenmp,  1},
    {"Tdemix",      (DL_FUNC) &Tdemix,      27},
    {NULL, NULL, 0}
};

void R_init_DeMixT(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

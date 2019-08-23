#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORNNLO_H
#define KFACTORNNLO_H
//needed for qqbar channel
double vegas_DY_NNLO_qqbar_LP(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_zdep(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_zdepcorr(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_delta(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_power(double *k, size_t dim, void *params);
//fit functions
double vegas_DY_NNLO_qqbar_LP_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_zdep_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_zdepcorr_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_delta_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qqbar_power_fit(double *k, size_t dim, void *params);

double DY_NNLO_qqbar_delta();
double DY_NNLO_qqbar_plus(double x);
double DY_NNLO_qqbar_NS(double x);
double DY_NNLO_qqbar_NS_expansion(double x, int power);
double DY_NNLO_BB_full(double x);
double DY_NNLO_BB_expansion(double x, int power);
double DY_NNLO_BC_full(double x);
double DY_NNLO_BC_expansion(double x, int power);

//needed for qq, qbarqbar
double vegas_DY_NNLO_qq_full(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qq_power(double *k, size_t dim, void *params);
// fitted
double vegas_DY_NNLO_qq_full_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qq_power_fit(double *k, size_t dim, void *params);

double DY_NNLO_CC_full(double x);
double DY_NNLO_CC_expansion(double x, int power);
double DY_NNLO_CD_full(double x);
double DY_NNLO_CD_expansion(double x, int power);
double DY_NNLO_CE_full(double x);
double DY_NNLO_CE_expansion(double x, int power);
double DY_NNLO_CF_full(double x);
double DY_NNLO_CF_expansion(double x, int power);


//needed for gg
double vegas_DY_NNLO_gg_full(double *k, size_t dim, void *params);
double vegas_DY_NNLO_gg_power(double *k, size_t dim, void *params);
// fitted
double vegas_DY_NNLO_gg_full_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_gg_power_fit(double *k, size_t dim, void *params);

double DY_NNLO_gg_full(double x);
double DY_NNLO_gg_expansion(double x, int power);

//needed for qg
double vegas_DY_NNLO_qg_full(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qg_power(double *k, size_t dim, void *params);
// fitted
double vegas_DY_NNLO_qg_full_fit(double *k, size_t dim, void *params);
double vegas_DY_NNLO_qg_power_fit(double *k, size_t dim, void *params);

double DY_NNLO_qg_full(double x);
double DY_NNLO_qg_expansion(double x, int power);

#endif

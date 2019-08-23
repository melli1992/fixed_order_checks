#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTOR_DY
#define KFACTOR_DY

//LO part
double DY_LO_factor();
double vegas_DY_LO(double *k, size_t dim, void *params);
//fit functions
double vegas_DY_LO_fit(double *k, size_t dim, void *params);


//NLO part

//needed for gg
double vegas_DY_NLO_qqbar_LP(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_zdep(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_zdepcorr(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_delta(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_power(double *k, size_t dim, void *params);
//fit functions
double vegas_DY_NLO_qqbar_LP_fit(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_zdep_fit(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_zdepcorr_fit(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_delta_fit(double *k, size_t dim, void *params);
double vegas_DY_NLO_qqbar_power_fit(double *k, size_t dim, void *params);

double DY_NLO_qqbar_delta();
double DY_NLO_qqbar_reg(double x);
double DY_NLO_qqbar_plus(double x);
double DY_NLO_qqbar_expansion(double x, int power);

//needed for qg
double vegas_DY_NLO_qg_full(double *k, size_t dim, void *params);
double vegas_DY_NLO_qg_power(double *k, size_t dim, void *params);
// fitted
double vegas_DY_NLO_qg_full_fit(double *k, size_t dim, void *params);
double vegas_DY_NLO_qg_power_fit(double *k, size_t dim, void *params);

double DY_NLO_qg_full(double x);
double DY_NLO_qg_expansion(double x, int power);

#endif

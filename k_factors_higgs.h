#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORHIGGS_H
#define KFACTORHIGGS_H

//LO part
double higgs_LO_factor();
std::complex<double> AQ(double x);
double vegas_higgs_LO(double *k, size_t dim, void *params);
//fit functions
double vegas_higgs_LO_fit(double *k, size_t dim, void *params);


//NLO part

//needed for gg
double vegas_higgs_NLO_gg_LP(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_zdepcorr(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_delta(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_power(double *k, size_t dim, void *params);
//fit functions
double vegas_higgs_NLO_gg_LP_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_zdepcorr_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_delta_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_power_fit(double *k, size_t dim, void *params);

double higgs_NLO_gg_delta();
double higgs_NLO_gg_reg(double x);
double higgs_NLO_gg_plus(double x);
double higgs_NLO_gg_expansion(double x, int power);

//needed for qg
double vegas_higgs_NLO_qg_full(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qg_power(double *k, size_t dim, void *params);
// fitted
double vegas_higgs_NLO_qg_full_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qg_power_fit(double *k, size_t dim, void *params);

double higgs_NLO_qg_full(double x);
double higgs_NLO_qg_expansion(double x, int power);


//needed for qqbar
double vegas_higgs_NLO_qqbar_full(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qqbar_power(double *k, size_t dim, void *params);
// fitted
double vegas_higgs_NLO_qqbar_full_fit(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qqbar_power_fit(double *k, size_t dim, void *params);

double higgs_NLO_qqbar_full(double x);
double higgs_NLO_qqbar_expansion(double x, int power);

#endif

#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORHIGGS_H
#define KFACTORHIGGS_H


double higgs_LO_factor();
double vegas_higgs_LO(double *k, size_t dim, void *params);

//needed for gg
double vegas_higgs_NLO_gg_LP(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_LP_corr(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_full(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_delta(double *k, size_t dim, void *params);
double vegas_higgs_NLO_gg_power(double *k, size_t dim, void *params);
double higgs_NLO_gg_LP(double x);
double higgs_NLO_gg_delta();
double higgs_NLO_gg_full(double x);
double higgs_NLO_gg_expansion(double x, int power);

//needed for qg
double vegas_higgs_NLO_qg_full(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qg_power(double *k, size_t dim, void *params);
double higgs_NLO_qg_full(double x);
double higgs_NLO_qg_expansion(double x, int power);


//needed for qqbar
double vegas_higgs_NLO_qqbar_full(double *k, size_t dim, void *params);
double vegas_higgs_NLO_qqbar_power(double *k, size_t dim, void *params);
double higgs_NLO_qqbar_full(double x);
double higgs_NLO_qqbar_expansion(double x, int power);

#endif

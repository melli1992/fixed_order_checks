#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORHIGGSNNLO_H
#define KFACTORHIGGSNNLO_H

//needed for gg
double vegas_higgs_NNLO_gg_LP(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_zdepcorr(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_delta(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_power(double *k, size_t dim, void *params);
//fitted functions
double vegas_higgs_NNLO_gg_LP_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_zdepcorr_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_delta_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_gg_power_fit(double *k, size_t dim, void *params);

double higgs_NNLO_gg_delta();
double higgs_NNLO_gg_reg(double x);
double higgs_NNLO_gg_plus(double x);
double higgs_NNLO_gg_expansion(double x, int power);

// qg (+qbar g)
double vegas_higgs_NNLO_qg_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qg_power(double *k, size_t dim, void *params);
//fitted
double vegas_higgs_NNLO_qg_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qg_power_fit(double *k, size_t dim, void *params);

double higgs_NNLO_qg_reg(double x);
double higgs_NNLO_qg_expansion(double x, int power);

// qq (+qbar qbar)
double vegas_higgs_NNLO_qq_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qq_power(double *k, size_t dim, void *params);
//fitted
double vegas_higgs_NNLO_qq_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qq_power_fit(double *k, size_t dim, void *params);

double higgs_NNLO_qq_reg(double x);
double higgs_NNLO_qq_expansion(double x, int power);

// qq' (+qbar qbar')
double vegas_higgs_NNLO_qqp_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qqp_power(double *k, size_t dim, void *params);
//fitted
double vegas_higgs_NNLO_qqp_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qqp_power_fit(double *k, size_t dim, void *params);

double higgs_NNLO_qqp_reg(double x);
double higgs_NNLO_qqp_expansion(double x, int power);

// qqbar
double vegas_higgs_NNLO_qqbar_zdep(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qqbar_power(double *k, size_t dim, void *params);
//fitted
double vegas_higgs_NNLO_qqbar_zdep_fit(double *k, size_t dim, void *params);
double vegas_higgs_NNLO_qqbar_power_fit(double *k, size_t dim, void *params);

double higgs_NNLO_qqbar_reg(double x);
double higgs_NNLO_qqbar_expansion(double x, int power);


#endif

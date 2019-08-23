#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

#ifndef MELLIN_H
#define MELLIN_H

// test functions
double vegas_fofx2_full(double *k, size_t dim, void *params);
double vegas_fofx2_deriv(double *k, size_t dim, void *params);
double vegas_fofx2_defor(double *k, size_t dim, void *params);
double vegas_fofx2_Nspace(double *k, size_t dim, void *params);

// the LO versions of DY(check of the implementations)
double vegas_DY_true(double *k, size_t dim, void *params);
double vegas_DY_fit(double *k, size_t dim, void *params);
double vegas_DY_mellin(double *k, size_t dim, void *params);
double vegas_DY_deriv(double *k, size_t dim, void *params);
double vegas_DY_defor(double *k, size_t dim, void *params);
// the LO versions of higgs 
double vegas_higgs_true(double *k, size_t dim, void *params);
double vegas_higgs_fit(double *k, size_t dim, void *params);
double vegas_higgs_mellin(double *k, size_t dim, void *params);

// the resummed functions
double vegas_resum_DY(double *k, size_t dim, void *params);
double vegas_resum_DY_expanded_NLO(double *k, size_t dim, void *params);
double vegas_resum_DY_expanded_NNLO(double *k, size_t dim, void *params);
double vegas_resum_higgs(double *k, size_t dim, void *params);
double vegas_resum_higgs_expanded_NLO(double *k, size_t dim, void *params);
double vegas_resum_higgs_expanded_NNLO(double *k, size_t dim, void *params);

#endif

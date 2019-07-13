#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

#ifndef MELLIN_H
#define MELLIN_H

double vegas_fofx2_full(double *k, size_t dim, void *params);
double vegas_fofx2_deriv(double *k, size_t dim, void *params);
double vegas_fofx2_defor(double *k, size_t dim, void *params);
double vegas_fofx2_Nspace(double *k, size_t dim, void *params);

double vegas_sigma0_nomel(double *k, size_t dim, void *params);
double vegas_sigma0_deriv(double *k, size_t dim, void *params);
double vegas_sigma0_defor(double *k, size_t dim, void *params);
double vegas_test_nomel(double *k, size_t dim, void *params);
double vegas_test_deriv(double *k, size_t dim, void *params);
double vegas_test_defor(double *k, size_t dim, void *params);
double vegas_resum_defor(double *k, size_t dim, void *params);
std::complex<double> LP_LL_q(std::complex<double>N);
#endif

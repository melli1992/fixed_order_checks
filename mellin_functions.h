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

// the LO version (check of the implementation)
double vegas_sigma0_true(double *k, size_t dim, void *params);
double vegas_sigma0_nomel(double *k, size_t dim, void *params);
double vegas_sigma0_exact(double *k, size_t dim, void *params);
double vegas_sigma0_deriv(double *k, size_t dim, void *params);
double vegas_sigma0_defor(double *k, size_t dim, void *params);

// the resummed functions
double vegas_resum_defor(double *k, size_t dim, void *params);

// the resummation coefficients
std::complex<double> LP_LL_q(std::complex<double>N);
std::complex<double> h0q(std::complex<double>lambda);
std::complex<double> h0g(std::complex<double>lambda);
std::complex<double> h0qNLP(std::complex<double> N, std::complex<double>lambda);
std::complex<double> h0gNLP(std::complex<double> N, std::complex<double>lambda);
std::complex<double> h1q(std::complex<double>lambda);
std::complex<double> h1g(std::complex<double>lambda);
#endif

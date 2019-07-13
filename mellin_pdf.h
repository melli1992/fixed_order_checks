#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef MELLINPDF_H
#define MELLINPDF_H



std::complex<double> mellin_pdf_sum_qqbar_charge_weighted(double x1, double x2, std::complex<double> N);
double deriv_pdf(int i, double x, double eps = 1.0E-5);
std::complex<double> mellin_test_pdf(double x1, double x2, std::complex<double> N);
std::complex<double> fit_sum_qqbar_charge_weighted(std::complex<double> x1, std::complex<double> x2);
std::complex<double> normal_sum_qqbar_charge_weighted(double x1,double x2);
double test_deriv_pdf(double x);
std::complex<double> xfit_pdfs(int i, std::complex<double> x);
std::complex<double> fit_pdfs(int i, std::complex<double> x);
double vegas_coefficients_0(double *k, size_t dim, void *params);
double vegas_coefficients_n(double *k, size_t dim, void *params);
double a0coeff(int i,double x);
double ancoeff(int i, int n, double x);
#endif

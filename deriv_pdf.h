#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef DERIVPDF_H
#define DERIVPDF_H

double pdf_sum_qq_charge_weighted(double x, double tau_over_z);
double pdf_sum_qg_charge_weighted(double x, double tau_over_z);
double derivative_qq_pdf(double x, double z, double tau, double eps= 0.0001);
double derivative_qg_pdf(double x, double z, double tau, double eps= 0.0001);
double derivative_qq_pdf_jac(double x, double z, double tau, double eps= 0.0001);
double derivative_qg_pdf_jac(double x, double z, double tau, double eps= 0.0001);
double vegas_sum_pdf_weigthed(double *k, size_t dim, void *params);
double vegas_pdf_up_minus_upbar(double *k, size_t dim, void *params);
double vegas_pdf_mom_consv(double *k, size_t dim, void *params);
#endif

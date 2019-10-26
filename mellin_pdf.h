#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef MELLINPDF_H
#define MELLINPDF_H



std::complex<double> mellin_pdf_sum_qqbar_charge_weighted(double x1, double x2, std::complex<double> N);
double deriv_xpdf(int i, double x, double eps = 1.0E-5);
double deriv_pdf(int i, double x, double eps = 1.0E-5);
std::complex<double> mellin_test_pdf(double x1, double x2, std::complex<double> N);
//weighted fitted pdfs
std::complex<double> fit_sum_qqbar_charge_weighted(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qqbar_charge_unweighted(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qq_charge_weighted_double(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qq_charge_weighted_double_vivj(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qq_charge_weighted_single(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qq_charge_weighted_single_vivi(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qg_charge_weighted(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_gg_charge_weighted(std::complex<double> x1, std::complex<double> x2);
//unweigthed fitted pdfs (checked with the real pdfs)
std::complex<double> fit_sum_gg(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qqbar(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qg(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qq(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qqNI(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_mellin_pdf_sum_gg(std::complex<double> Nint);
std::complex<double> fit_mellin_pdf_sum_qqbar_charge_weighted(std::complex<double> Nint);
std::complex<double> fit_mellin_pdf_sum_qqbarUP(std::complex<double> Nint);
std::complex<double> fit_mellin_pdf_sum_qqbarDOWN(std::complex<double> Nint);
//unweigthed fitted pdfs UP and DOWN seperately
std::complex<double> fit_sum_qqbarUP(std::complex<double> x1, std::complex<double> x2);
std::complex<double> fit_sum_qqbarDOWN(std::complex<double> x1, std::complex<double> x2);

std::complex<double> xfit_pdfs(int i, std::complex<double> x);
std::complex<double> Dxfit_pdfs(int i, std::complex<double> x);
std::complex<double> NDxfit_pdfs(int i, std::complex<double> x, double eps = 1.0E-5);
std::complex<double> fit_pdfs(int i, std::complex<double> x);
std::complex<double> Dfit_pdfs(int i, std::complex<double> x);
std::complex<double> NDfit_pdfs(int i, std::complex<double> x, double eps = 1.0E-5);
std::complex<double> xfit_Nspace_pdfs(int i, std::complex<double> N);
#endif

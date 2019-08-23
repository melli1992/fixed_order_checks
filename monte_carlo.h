#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include "parameters.h"
#include "k_factors_dy.h"

#ifndef MONTE_H
#define MONTE_H

struct results{double res; double err;};
struct functionint{gsl_monte_function G; std::vector<double> xl; std::vector<double> xu;};
void display_results(std::string title, double &result, double &error);
functionint init_vegas_coefficients(int n);
functionint init_vegas_mellin(std::string process = "DY", std::string order = "resum");
functionint init_vegas_DY(std::string order, std::string power="full", std::string process = "gg", int power_number = 0, bool fitted_pdfs = false);
functionint init_vegas_higgs(std::string order, std::string power="full", std::string process = "gg", int power_number = 0, bool fitted_pdfs = false);
functionint init_vegas_pf(std::string order, std::string power="full", std::string process = "qqbar", int power_number = 0, bool integrated = false);
results call_vegas(functionint integrand, lumni_params params, bool verbal = false,  bool high_prec = false);
double integrand(double *k, size_t dim, void *params);


#endif

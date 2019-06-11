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
functionint init_vegas_dy(std::string order, std::string power="full", std::string process = "qqbar", bool integrated = false);
functionint init_vegas_higgs(std::string order, std::string power="full", std::string process = "gg", bool integrated = false);
results call_vegas(functionint integrand, /*gsl_monte_function G, double *xl, double *xu,*/ lumni_params params, bool verbal = false,  bool high_prec = false);
double integrand(double *k, size_t dim, void *params);


#endif

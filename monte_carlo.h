#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "parameters.h"
#include "k_factors_dy.h"

#ifndef MONTE_H 
#define MONTE_H

struct results{double res; double err;};
void display_results(std::string title, double &result, double &error);
results call_vegas(std::string integrand, double *xl, double *xu, lumni_params params, bool verbal = false,  bool high_prec = false);
double integrand(double *k, size_t dim, void *params);


#endif

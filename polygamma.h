#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include <iostream>

#ifndef POLYGAM_H
#define POLYGAM_H

double Li2(double x);
double Li3(double x);
double S12(double x);
double clenshaw(std::vector<double> coeff, double x);
std::complex<double> Gamma(std::complex<double> x);

#endif

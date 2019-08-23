#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

#ifndef RESUM_H
#define RESUM_H

// the resummation coefficients
std::complex<double> h0(double A1, std::complex<double>lambda);
std::complex<double> h0g(std::complex<double>lambda);
std::complex<double> h0q(std::complex<double>lambda);
std::complex<double> h0qNLP(std::complex<double> N, std::complex<double>lambda);
std::complex<double> h0gNLP(std::complex<double> N, std::complex<double>lambda);
std::complex<double> h1(double A1,double A2,std::complex<double>lambda);
std::complex<double> h1g(std::complex<double>lambda);
std::complex<double> h1q(std::complex<double>lambda);
std::complex<double> h2(double A1,double A2,double A3,std::complex<double>lambda);
std::complex<double> h2g(std::complex<double>lambda);
std::complex<double> h2q(std::complex<double>lambda);
std::complex<double> wideangle(double D2,std::complex<double>lambda);
std::complex<double> NLOmatch_higgs(std::complex<double> N);
std::complex<double> NNLOmatch_higgs(std::complex<double> N);
std::complex<double> NLOmatch_DY(std::complex<double> N);
std::complex<double> NNLOmatch_DY(std::complex<double> N);
double FDY();
#endif

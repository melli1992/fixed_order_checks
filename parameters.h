#include <iostream>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"


#ifndef PARAM_H //need this otherwise compiler error that things were predefined (called a guard)
#define PARAM_H

extern std::complex<double> I;
extern double CMP, phiMP;
extern double tau;
extern double S;
extern double S2;
extern double Q;
extern double Q2;
extern double muF;
extern double muF2;
extern double muR;
extern double muR2;
extern double CF;
extern double CA;
extern double TF;
extern double alphas_muF;
extern double alphas_muR;
extern double alphas_Q;
extern double alphaEM;
extern double GF;
extern double nF;
extern double zeta2;
extern double b0;
extern double zeta3;
extern double pbunits;
extern std::string setname;
extern std::vector<LHAPDF::PDF*> pdfs; //pdf vector
extern double xmin_pdfs, xmax_pdfs; //min x, max x

struct lumni_params {double z; double pT; double xT; double epeta; double emeta; int power; int flavor; int coefficient;};
void print_usage();
void print_defaults();
void read_arguments(int argc, char* argv[]);


#endif

#include <bits/stdc++.h>
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
extern double pbunits;
extern double fbunits;

extern double quarkmasses[2];
extern double mt;
extern double mt2;
extern double mb;
extern double mb2;
extern double mZ2;
extern double mW2;
extern double mH;
extern double mH2;
extern double mHeavy2;
extern double mA2;
extern double GammaH;
extern double Lt;

extern double M_gammaE;
extern double zeta2;
extern double zeta3;
extern double zeta5;

extern double CF;
extern double CA;
extern double TF;
extern double alphas_muF;
extern double alphas_muR;
extern double alphas_Q;
extern double alphaEM;
extern double GF;
extern double nF;

extern double b0;
extern double b1;
extern double b2;
extern double beta0;
//for W+W-, ZZ
extern double sinw;
extern double cosw;
extern double e;
extern double ez;
extern double guL;
extern double guR;
extern double gdL;
extern double gdR;
extern double gVu;
extern double gAu;
extern double gVd;
extern double gAd;
extern double ctt;
//for SUSY dihiggs;

extern double ght;
extern double gHt;
extern double gAt;
extern double ghb;
extern double gHb;
extern double gAb;
extern double lambdaZAh;
extern double lambdaZAH;
extern double lambdahhh;
extern double lambdaHhh;
extern double lambdaHHh;
extern double lambdaHHH;
extern double lambdahAA;
extern double lambdaHAA;


extern double ISLL;
extern double ISNLL;
extern double ISNNLL;
extern double ISNLP;

extern double A1q;
extern double A1g;
extern double A2q;
extern double A2g;
extern double A3q;
extern double A3g;
extern double B1q;
extern double B1g;
extern double C2qgamma;
extern double Dbar2q;
extern double D1DY;
extern double D2DY;
extern double D3DY;
extern double D1higgs;
extern double D2higgs;
extern double D3higgs;

extern std::string setname;
extern std::vector<LHAPDF::PDF*> pdfs; //pdf vector
extern double xmin_pdfs, xmax_pdfs; //min x, max x
extern int use_member; //the member that one needs to use
struct lumni_params {double z; double pT; double xT; double epeta; double emeta; int power; int flavor; int coefficient;};

extern bool DY, higgs, hh, WW, ZZ, diff, full, PF, LO, NLO, NNLO, RES, realPDF, fitPDF;

void print_usage();
void update_defaults(bool printout = true , bool pdfset = true);
void set_SUSY(int set, bool printout = true);
void read_arguments(int argc, char* argv[]);

extern std::unordered_map<double, std::vector<std::vector<double>>> fitcoeff;
#endif

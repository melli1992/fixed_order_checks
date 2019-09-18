#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORDIHIGGS_H
#define KFACTORDIHIGGS_H

//LO part
std::complex<double> Cab(double s, double mQ2);
std::complex<double> Cac(double t, double mH2, double mQ2);
std::complex<double> Cad(double u, double mH2, double mQ2);
std::complex<double> Cbc(double u, double mH2, double mQ2);
std::complex<double> Cbd(double t, double mH2, double mQ2);
std::complex<double> Ccd(double s, double mH2, double mQ2);
std::complex<double> Dabc(double s, double t, double mH2, double mQ2);
std::complex<double> Dacb(double t, double u, double mH2, double mQ2);
std::complex<double> Dbac(double s, double t, double mH2, double mQ2);
std::complex<double> Fdelta(double s, std::complex<double> cab);
std::complex<double> Fbox(double s, double t, double u, std::complex<double> cab, std::complex<double> cac, std::complex<double> cad, std::complex<double> cbc, std::complex<double> cbd, std::complex<double> ccd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> Gbox(double s, double t, double u, std::complex<double> cab, std::complex<double> cac, std::complex<double> cad, std::complex<double> cbc, std::complex<double> cbd, std::complex<double> ccd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
double dihiggs_LO_factor(double scale2, double ctheta1);
double dihiggs_LO_factor_approx(double scale2, double ctheta1);

#endif

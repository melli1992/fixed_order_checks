#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "mellin_functions.h"
#include "k_factors_dy.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"

using namespace std;

///////////////////////////
/// test functions
///////////////////////////
/// full transformation
double vegas_fofx2_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
  //cout << k[1] << " " << pow(k[1],Nint+1.) << endl;
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint-1.)*pow((1.-k[1]),2.));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// derivative approach
double vegas_fofx2_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
  double result = -2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint)/Nint*-2.*(1.-k[1]));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// deformation of the contour
double vegas_fofx2_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*I;
  complex<double> Njac = 1./pow(1.-k[0],2)*I;
  double wr = k[1]/(1.+k[1]);
  double wjac = 1./(pow(1.+k[1],2));
  complex<double> x = exp(wr/Nint);
  //cout << "x " << x << endl;
	double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*wjac*exp(wr)/Nint*pow(1.-x,2));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// direct N-space result
double vegas_fofx2_Nspace(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	return 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*2./(Nint*(Nint*Nint+3.*Nint+2.)));
}

/////////////////////////
/// LO mellin transform
/////////////////////////
/// checked and is correct, returns dsigma/dQ2 at LO
double vegas_sigma0_nomel(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
  double result = LO_factor()*pdf_sum_qqbar_charge_weighted(k[0], tau);
	if (isnan(result)){return 0;}
	else{return result;}
}

double vegas_sigma0_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*LO_factor()*mellin_pdf_sum_qqbar_charge_weighted(k[1],k[2],Nint));
	if (isnan(result)){return 0;}
	else{return result;}
}


double vegas_sigma0_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  double wr1 = k[1]/(1.+k[1]); //note deformation of the contour!
  double wjac1 =  1./(pow(1.+k[1],2));
  complex<double> x1 = exp(wr1/Nint);
  double wr2 =  k[2]/(1.+k[2]);
  double wjac2 =  1./(pow(1.+k[2],2));
  complex<double> x2 = exp(wr2/Nint);
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac1*wjac2*pow(Nint,-2)*LO_factor()*exp(wr1)*exp(wr2)*fit_sum_qqbar_charge_weighted(x1,x2));
  if (isnan(result)){return 0;}
	else{return result;}
}


//these test cases seem to work, although there is a slight deviation
/// find out where the deviation is comming from, is it numerical or true?
double vegas_test_nomel(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
  double result = LO_factor()*1./k[0]*(1.-pow(k[0],2)*(1.-pow(tau/k[0],2)));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_test_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*LO_factor()*mellin_test_pdf(k[1],k[2],Nint));
	if (isnan(result)){return 0;}
	else{return result;}
}

double vegas_test_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  double wr1 = k[1]/(1.+k[1]); //note deformation of the contour!
  double wjac1 =  1./(pow(1.+k[1],2));
  complex<double> x1 = exp(wr1/Nint);
  double wr2 =  k[2]/(1.+k[2]);
  double wjac2 =  1./(pow(1.+k[2],2));
  complex<double> x2 = exp(wr2/Nint);
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac1*wjac2*pow(Nint,-2)*LO_factor()*exp(wr1)*exp(wr2)*(1.-pow(x1,2))*(1.-pow(x2,2)));
  if (isnan(result)){return 0;}
	else{return result;}
}

////////////////////////
/// resummation code
////////////////////////

double vegas_resum_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  double wr1 = k[1]/(1.+k[1]); //note deformation of the contour!
  double wjac1 =  1./(pow(1.+k[1],2));
  complex<double> x1 = exp(wr1/Nint);
  double wr2 =  k[2]/(1.+k[2]);
  double wjac2 =  1./(pow(1.+k[2],2));
  complex<double> x2 = exp(wr2/Nint);
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac1*wjac2*pow(Nint,-2)*exp(2.*LP_LL_q(Nint))*LO_factor()*exp(wr1)*exp(wr2)*fit_sum_qqbar_charge_weighted(x1,x2));
  if (isnan(result)){return 0;}
	else{return result;}
}

////////////////////////
/// resummation factors
////////////////////////
complex<double> LP_LL_q(complex<double>N){
	return alphas_Q/M_PI*CF*(pow(log(N),2)+log(N)/N);
}

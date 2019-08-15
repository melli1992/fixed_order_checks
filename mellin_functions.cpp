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
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint)*(1.-k[1]));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// derivative approach
double vegas_fofx2_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
  double result = -2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint)/Nint*(1.-2.*k[1]));
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
	double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*wjac*exp(wr)/Nint*x*(1.-x));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// direct N-space result
double vegas_fofx2_Nspace(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	return 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*1./((Nint*Nint+3.*Nint+2.)));
}

/////////////////////////
/// LO mellin transform
/////////////////////////
/// checked and is correct, returns dsigma/dQ2 at LO
/// also checked the implementation of the fit_sum_qqbar function, returns the same value within 4.20677e-07 (fake) vs 4.21225e-07 (real) 
double vegas_sigma0_true(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
    double result = LO_factor()*real(pdf_sum_qqbar_charge_weighted(k[0], tau)); // this is the real PDF
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_sigma0_nomel(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
    double result = LO_factor()*real(fit_sum_qqbar_charge_weighted(k[0], tau)); // this is the fitted PDF
	if (isnan(result)){return 0;}
	else{return result;}
}
// this function gives the same result as the one above
double vegas_sigma0_exact(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*LO_factor()*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}
}
// the derivative method, seems very unstable
double vegas_sigma0_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*LO_factor()*mellin_pdf_sum_qqbar_charge_weighted(k[1],k[2],Nint));
	if (isnan(result)){return 0;}
	else{return result;}
}
// deformation, implementation has an error (the PDF convolution is wrong)
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

////////////////////////
/// resummation code
////////////////////////

/// this implementation fully works
double vegas_resum_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
    //double wr1 = k[1]/(1.+k[1]); //note deformation of the contour!
    //double wjac1 =  1./(pow(1.+k[1],2));
    //complex<double> x1 = exp(wr1/Nint);
    //double wr2 =  k[2]/(1.+k[2]);
    //double wjac2 =  1./(pow(1.+k[2],2));
    //complex<double> x2 = exp(wr2/Nint);
    //double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac1*wjac2*pow(Nint,-2)*exp(2.*LP_LL_q(Nint))*LO_factor()*exp(wr1)*exp(wr2)*fit_sum_qqbar_charge_weighted(x1,x2));
   
	complex<double> lambda = alphas_Q*b0*log(Nint*exp(M_gammaE));
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*LO_factor()*exp(2.*(1./alphas_Q*h0q(lambda)+h1q(lambda)+0.*h0qNLP(Nint,lambda)))*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}
}

////////////////////////
/// resummation factors
////////////////////////
complex<double> LP_LL_q(complex<double>N){
	return alphas_Q/M_PI*A1q*(pow(log(N),2)+log(N)/N);
}

// LP LL (quark)
complex<double> h0q(complex<double>lambda){
	return A1q/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}

// LP LL (gluon)
complex<double> h0g(complex<double>lambda){
	return A1g/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}

// NLP LL (quark)
complex<double> h0qNLP(complex<double> N, complex<double>lambda){
	return -A1q/(2.*M_PI*b0)*(log(1.-2.*lambda))/N;
}

// NLP LL (gluon)
complex<double> h0gNLP(complex<double> N, complex<double>lambda){
	return -A1g/(2.*M_PI*b0)*(log(1.-2.*lambda))/N;
}

// LP NLL (quark)
complex<double> h1q(complex<double>lambda){
	return 1./(2.*M_PI*b0)*(-A2q/(M_PI*b0)+A1q*log(Q2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1q*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1q/(M_PI*b0)*log(Q2/muF2);
}

// LP NLL (gluon)
complex<double> h1g(complex<double>lambda){
	return 1./(2.*M_PI*b0)*(-A2g/(M_PI*b0)+A1g*log(Q2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1g*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1g/(M_PI*b0)*log(Q2/muF2);
}

//constants (DY case in MSbar scheme - hep-ph/0508284 Eq. 3.13 and 4.6)
double FDY(){
	return -alphas_Q/(M_PI)*CF*(3./2.*zeta2)+pow(alphas_Q/M_PI,2)*(CA*CF*(607./324.-469./144.*zeta2+1./4.*pow(zeta2,2)-187./72*zeta3)+(-41./162.+35./72.*zeta2+17./36.*zeta3)*nf*CF);
}

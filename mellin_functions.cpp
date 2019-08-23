#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "mellin_functions.h"
#include "resum_functions.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_higgs.h"
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

// first the DY functions
/// checked the implementation of the fit_sum_qqbar function, returns the same value within 4.20677e-07 (fake) vs 4.21225e-07 (real) 
double vegas_DY_true(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
    double result = DY_LO_factor()*real(pdf_sum_qqbar_charge_weighted(k[0], tau)); // this is the real PDF
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
    double result = DY_LO_factor()*real(fit_sum_qqbar_charge_weighted(k[0], tau)); // this is the fitted PDF
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_DY_mellin(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}
}
// the derivative method, seems very unstable
double vegas_DY_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*mellin_pdf_sum_qqbar_charge_weighted(k[1],k[2],Nint));
	if (isnan(result)){return 0;}
	else{return result;}
}
// deformation, implementation has an error (the PDF convolution is wrong)
double vegas_DY_defor(double *k, size_t dim, void *params){
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
    double result = 2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*wjac1*wjac2*pow(Nint,-2)*DY_LO_factor()*exp(wr1)*exp(wr2)*fit_sum_qqbar_charge_weighted(x1,x2));
    if (isnan(result)){return 0;}
	else{return result;}
}

/// //////////////////////////////////////////////
// and the higgs functions
double vegas_higgs_mellin(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()/tau*fit_mellin_pdf_sum_gg(Nint));
	if (isnan(result)){return 0;}
	else{return result;}
}

////////////////////////
/// resummation code
////////////////////////
/// this implementation fully works
double vegas_resum_DY(double *k, size_t dim, void *params){
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
   
	complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*exp(ISNNLL*alphas_muR*wideangle(D2DY,lambda)+2.*(1./alphas_muR*ISLL*h0q(lambda)+ISNLL*h1q(lambda)+ISNNLL*alphas_muR*h2q(lambda)+ISNLP*h0qNLP(Nint,lambda)))*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_resum_DY_expanded_NLO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*NLOmatch_DY(Nint)*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}	
}
double vegas_resum_DY_expanded_NNLO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*DY_LO_factor()*NNLOmatch_DY(Nint)*fit_mellin_pdf_sum_qqbar_charge_weighted(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}	
}

/// higgs function
double vegas_resum_higgs(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  
	complex<double> lambda = alphas_muR*b0*log(Nint*exp(M_gammaE));
	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*exp(ISNNLL*alphas_muR*wideangle(D2higgs,lambda)+2.*(1./alphas_muR*ISLL*h0g(lambda)+ISNLL*h1g(lambda)+ISNNLL*alphas_muR*h2g(lambda)+ISNLP*h0gNLP(Nint,lambda)))*fit_mellin_pdf_sum_gg(Nint-1.));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_resum_higgs_expanded_NLO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*NLOmatch_higgs(Nint)*fit_mellin_pdf_sum_gg(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}	
}
double vegas_resum_higgs_expanded_NNLO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(phiMP*I);
    complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
  	double result =  2.*imag(1./(2*M_PI)*Njac*pow(tau,-Nint)*higgs_LO_factor()*NNLOmatch_higgs(Nint)*fit_mellin_pdf_sum_gg(Nint-1.));
	if (isnan(result)){return 0;}
	else{return result;}	
}

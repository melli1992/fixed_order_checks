#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "k_factors_higgs.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for higgs 
/// split up in LO, NLO and LP, NLP, NNLP, NNNLP, full
///
//////////////////////////////////////////////////////////



//////// 
/// LO
////////

////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in https://arxiv.org/pdf/hep-ph/0207004.pdf, beneath eqn. 43
/// note that we replace v^2 with GF and Q2, S2 like in https://arxiv.org/pdf/0809.4283.pdf
////////////////////////////////////////////////////////////////////////////////////////////////////
double higgs_LO_factor(){
	return pbunits*alphas_Q*alphas_Q*Q2*sqrt(2)*GF/576./M_PI/S2;
}
////////////////////////////////////////////////////////////////
/// constant delta contribution
////////////////////////////////////////////////////////////////
double vegas_higgs_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;
	double x = k[0]*jacobian;
	return higgs_LO_factor()*(jacobian*pdf_sum_gg(k[0]*jacobian,tau));
}


//////// 
/// NLO
////////

///////////////////////////
/// gg channel
///////////////////////////

/////////// integration routines
double vegas_higgs_NLO_gg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_full(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_gg_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_LP(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_gg(tau+k[1]*(1.-tau),tau));
}
double vegas_higgs_NLO_gg_LP_corr(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; //needed to transform the boundary dependent terms
	double x = tau+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_LP(z))*jacobian*(-pdf_sum_gg(x,tau));
}
double vegas_higgs_NLO_gg_NLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_NLP(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_gg_NNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_NNLP(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_gg_NNNLP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_NNNLP(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_gg_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;
	double x = k[0]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_delta())*(jacobian*pdf_sum_gg(k[0]*jacobian,tau));
}

/////// and the NLO functions
double higgs_NLO_gg_LP(double x){
	return (12*alphas_Q*log(1 - x))/(M_PI - M_PI*x);
}
double higgs_NLO_gg_delta(){
	return (alphas_Q*(5.5 + 6*zeta2))/M_PI;
}
double higgs_NLO_gg_NLP(double x){
	return (6.*alphas_Q*(1. - 4.*log(1. - x)))/M_PI;
}
double higgs_NLO_gg_NNLP(double x){
	return (-9*alphas_Q*(-1 + x)*(-1 + 4*log(1 - x)))/M_PI;
}
double higgs_NLO_gg_NNNLP(double x){
	return (-2*alphas_Q*pow(-1 + x,2)*(-7 + 12*log(1 - x)))/M_PI;
}
double higgs_NLO_gg_full(double x){
	return -(alphas_Q*(-11*pow(-1 + x,4) + 24*x*(-2 + 3*x - 2*pow(x,2) + pow(x,3))*log(1 - x) - 12*pow(1 - x + pow(x,2),2)*log(x)))/(2.*M_PI*(-1 + x));
}

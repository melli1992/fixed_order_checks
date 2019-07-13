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
/// this is sigma0*tau!
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
double vegas_higgs_NLO_gg_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double jacobian = 1.;
	double x = k[0]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_delta())*(jacobian*pdf_sum_gg(k[0]*jacobian,tau));
}
double vegas_higgs_NLO_gg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_expansion(z, fp->power))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}

/////// and the NLO functions
double higgs_NLO_gg_LP(double x){
	return (-6*alphas_Q*(log(Q2/muF2) + 2*log(1 - x)))/(M_PI*(-1 + x));
	//return (12*alphas_Q*log(1 - x))/(M_PI - M_PI*x);
}
double higgs_NLO_gg_delta(){
	return (alphas_Q*(5.5 + 6*zeta2))/M_PI;
}
double higgs_NLO_gg_full(double x){
	
	return -(alphas_Q*(-11*pow(-1 + x,4) + 12*(1 - 3*x + 3*pow(x,2) - 2*pow(x,3) + pow(x,4))*log((Q2*pow(-1 + x,2))/(muF2*x)) - 12*x*log(x)))/(2.*M_PI*(-1 + x)*x);

	//return (alphas_Q*(12*(log(Q2/muF2) + 2*log(1 - x)) - (-11*pow(-1 + x,4) + 24*pow(1 - x + pow(x,2),2)*log(1 - x) - 12*pow(1 - x + pow(x,2),2)*log(x))/x))/(2.*M_PI*(-1 + x));
	//return -(alphas_Q*(-11*pow(-1 + x,4) + 24*x*(-2 + 3*x - 2*pow(x,2) + pow(x,3))*log(1 - x) - 12*pow(1 - x + pow(x,2),2)*log(x)))/(2.*M_PI*(-1 + x));
}
double higgs_NLO_gg_expansion(double x, int power){
	if(power==1){
		return (-6*alphas_Q*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI;
	}
	else if(power==2){
		return (-3*alphas_Q*(-1 + x)*(-1 + 4*log(Q2/muF2) + 8*log(1 - x)))/M_PI;
	}
	else if(power==3){
		return (11*alphas_Q*pow(-1 + x,2))/M_PI;
	}
	else if(power==4){
		return (-6*alphas_Q*pow(-1 + x,3)*(log(Q2/muF2) + 2*log(1 - x)))/M_PI;
	}
	else if(power==5){
		return (3*alphas_Q*pow(-1 + x,4)*(7 + 10*log(Q2/muF2) + 20*log(1 - x)))/(5.*M_PI);
	}
	else if(power==6){
		return (alphas_Q*pow(1 - x,5)*(6.3 + 6*log(Q2/muF2) + 12*log(1 - x)))/M_PI;
	}
	else if(power==7){
		return (alphas_Q*pow(-1 + x,6)*(7.757142857142857 + 6*log(Q2/muF2) + 12*log(1 - x)))/M_PI;
	}
	else if(power==8){
		return (alphas_Q*pow(1 - x,7)*(8.892857142857142 + 6*log(Q2/muF2) + 12*log(1 - x)))/M_PI;
	}
	else if(power==9){
		return (alphas_Q*pow(-1 + x,8)*(9.83095238095238 + 6*log(Q2/muF2) + 12*log(1 - x)))/M_PI;
	}
	else if(power==10){
		return (alphas_Q*pow(1 - x,9)*(10.633333333333333 + 6*log(Q2/muF2) + 12*log(1 - x)))/M_PI;
	}
	else{
		return 0;
	}
}

///////////////////////////
/// qg channel
///////////////////////////

/////////// integration routines
double vegas_higgs_NLO_qg_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qg_full(z))*(jacobian*pdf_sum_qg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_qg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qg_expansion(z, fp->power))*(jacobian*pdf_sum_qg(tau/z+k[1]*jacobian,tau/z)/z);
}

///// NLO functions
double higgs_NLO_qg_full(double x){
	return (alphas_Q*(-3 - (-6 + x)*x + 2*(2 + (-2 + x)*x)*log(Q2/muF2) + 4*(2 + (-2 + x)*x)*log(1 - x) - 2*(2 + (-2 + x)*x)*log(x)))/(3.*M_PI*x);
}
double higgs_NLO_qg_expansion(double x, int power){		
	if(power==1){
		return (2*alphas_Q*(1 + log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==2){
		return (-2*alphas_Q*(-1 + x)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==3){
		return (4*alphas_Q*pow(-1 + x,2)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==4){
		return (-4*alphas_Q*pow(-1 + x,3)*(2 + 3*log(Q2/muF2) + 6*log(1 - x)))/(9.*M_PI);
	}
	if(power==5){
		return (alphas_Q*pow(-1 + x,4)*(25 + 24*log(Q2/muF2) + 48*log(1 - x)))/(18.*M_PI);
	}
	if(power==6){
		return (alphas_Q*pow(1 - x,5)*(5.233333333333333 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==7){
		return (alphas_Q*pow(-1 + x,6)*(91 + 60*log(Q2/muF2) + 120*log(1 - x)))/(45.*M_PI);
	}
	if(power==8){
		return (alphas_Q*pow(1 - x,7)*(6.752380952380952 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==9){
		return (alphas_Q*pow(-1 + x,8)*(7.335714285714285 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==10){
		return (alphas_Q*pow(1 - x,9)*(7.843650793650793 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	else{
		return 0;
	}
}


///////////////////////////
/// qg channel
///////////////////////////

/////////// integration routines
double vegas_higgs_NLO_qqbar_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qqbar_full(z))*(jacobian*pdf_sum_qqbar(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_qqbar_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qqbar_expansion(z, fp->power))*(jacobian*pdf_sum_qqbar(tau/z+k[1]*jacobian,tau/z)/z);
}

///// NLO functions
double higgs_NLO_qqbar_full(double x){
	return (-32*alphas_Q*pow(-1 + x,3))/(27.*M_PI*x);
}
double higgs_NLO_qqbar_expansion(double x, int power){		
	if(power==1){
		return 0;
	}
	if(power==2){
		return 0;
	}
	if(power==3){
		return 0;
	}
	if(power==4){
		return (-32*alphas_Q*pow(-1 + x,3))/(27.*M_PI);
	}
	if(power==5){
		return (32*alphas_Q*pow(-1 + x,4))/(27.*M_PI);
	}
	if(power==6){
		return (-32*alphas_Q*pow(-1 + x,5))/(27.*M_PI);
	}
	if(power==7){
		return (32*alphas_Q*pow(-1 + x,6))/(27.*M_PI);
	}
	if(power==8){
		return (-32*alphas_Q*pow(-1 + x,7))/(27.*M_PI);
	}
	if(power==9){
		return (32*alphas_Q*pow(-1 + x,8))/(27.*M_PI);
	}
	if(power==10){
		return (-32*alphas_Q*pow(-1 + x,9))/(27.*M_PI);
	}
	else{
		return 0;
	}
}

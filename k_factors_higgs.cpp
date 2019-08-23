#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_higgs.h"
using namespace std;

////////////////////////////////////////////////////////////
///
/// contains all K factors for higgs 
/// split up in LO, NLO and power corrections
/// also have two routines: either fitted pdfs or real ones
///
////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////
/// LO
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in https://arxiv.org/pdf/hep-ph/0207004.pdf, beneath eqn. 43
/// note that we replace v^2 with GF ( v^2 = 1/(Sqrt[2]*GF)) like in https://arxiv.org/pdf/0809.4283.pdf
/// this is sigma0*tau!
////////////////////////////////////////////////////////////////////////////////////////////////////
double higgs_LO_factor(){
	complex<double> AQtot = 0;
	for(int i = 0; i<2; i++){AQtot = AQtot + AQ(4.*pow(quarkmasses[i],2)/pow(mH,2));}
	AQtot = norm(AQtot);
	AQtot = 1;
	return pbunits*alphas_muR*alphas_muR*Q2*sqrt(2)*GF/576./M_PI/S2*real(AQtot);
}
// this factor is checked, get the same pb as the results in 0809.4283 fig 1 (also changed higgs mass to check it)
complex<double> AQ(double x){
	if(x>=1){ return 3./2.*x*(1.+(1.-x)*pow(asin(1./sqrt(x)),2));}
	if(x<1){ return 3./2.*x*(1.+(1.-x)*-1./4.*pow(log((1.+sqrt(1.-x))/(1.-sqrt(1.-x))-I*M_PI),2));}
}
////////////////////////////////////////////////////////////////
/// constant delta contribution
////////////////////////////////////////////////////////////////
double vegas_higgs_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = higgs_LO_factor()*(pdf_sum_gg(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_LO_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = higgs_LO_factor()*real(fit_sum_gg(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}

/////////////////////////////////////////////////////////////////////////////////////////
/// NLO
/////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////
/// gg channel
///////////////////////////

/////////// integration routines

double vegas_higgs_NLO_gg_zdep(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	
	
	double result = higgs_LO_factor()*(higgs_NLO_gg_reg(z)*jacobian*real(pdf_sum_gg(x,tau/z))/z+higgs_NLO_gg_plus(z)*(jacobian*real(pdf_sum_gg(x,tau/z))/z-(1.-tau)*real(pdf_sum_gg(tau+k[1]*(1.-tau),tau))));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_zdep_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0];
	double jacobian = 1.-tau/z; 
	double x = tau/z+k[1]*jacobian;
	
	double result = higgs_LO_factor()*(higgs_NLO_gg_reg(z)*jacobian*real(fit_sum_gg(x,tau/z))/z+higgs_NLO_gg_plus(z)*(jacobian*real(fit_sum_gg(x,tau/z))/z-(1.-tau)*real(fit_sum_gg(tau+k[1]*(1.-tau),tau))));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_zdepcorr(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	double result =  higgs_LO_factor()*(higgs_NLO_gg_plus(z))*jacobian*(-pdf_sum_gg(x,tau));
	
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_zdepcorr_fit(double *k, size_t dim, void *params){ /*correction from not integrating from 0 to 1 but tau to 1*/
	(void)(dim);
	(void)(params);

	double z = k[0]; // from 0 to tau
	double jacobian = 1.-tau; 
	double x = tau+k[1]*jacobian;
	double result =  higgs_LO_factor()*(higgs_NLO_gg_plus(z))*jacobian*(-real(fit_sum_gg(x,tau)));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_delta(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = higgs_LO_factor()*higgs_NLO_gg_delta()*real(pdf_sum_gg(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_delta_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double result = higgs_LO_factor()*higgs_NLO_gg_delta()*real(fit_sum_gg(k[0],tau));
	if (isnan(result)){return 0;}
	else{return result;}
}
double vegas_higgs_NLO_gg_LP(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_plus(z))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z - (1.-tau)*pdf_sum_gg(tau+k[1]*(1.-tau),tau));
}
double vegas_higgs_NLO_gg_LP_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_plus(z))*(jacobian*real(fit_sum_gg(tau/z+k[1]*jacobian,tau/z))/z - (1.-tau)*(real(fit_sum_gg(tau+k[1]*(1.-tau),tau))));
}
double vegas_higgs_NLO_gg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_expansion(z, fp->power))*(jacobian*pdf_sum_gg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_gg_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_gg_expansion(z, fp->power))*(jacobian*real(fit_sum_gg(tau/z+k[1]*jacobian,tau/z))/z);
}

/// and the NLO functions
//https://arxiv.org/pdf/0809.4283.pdf eqn. 9 with extra 1/x (see eqn. 1)
double higgs_NLO_gg_reg(double x){
	return alphas_muR/M_PI*(6.*(1./x-2.+x-pow(x,2))*log(Q2/muF2)+(11.*pow(-1.+x,4)-24.*(-1.+x)*(-1.+x*(2.+(-1.+x)*x))*log(1.-x) + 12.*pow(1.+(-1.+x)*x,2)*log(x))/(2.*(-1.+x)*x));
}
//https://arxiv.org/pdf/0809.4283.pdf eqn. 7 and 8
double higgs_NLO_gg_plus(double x){
	return alphas_muR/M_PI*((12.*log(1.-x)+6.*log(Q2/muF2))/(1.-x));
}
//https://arxiv.org/pdf/0809.4283.pdf eqn. 7 (see eqn. 1) 
double higgs_NLO_gg_delta(){
	return alphas_muR/M_PI*(11./2.+pow(M_PI,2));
}
double higgs_NLO_gg_expansion(double x, int power){
	if(power==1){
		return (-6*alphas_muR*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI;
	}
	if(power==2){
		return (-3*alphas_muR*(-1 + x)*(-1 + 4*log(Q2/muF2) + 8*log(1 - x)))/M_PI;
	}
	if(power==3){
		return (11*alphas_muR*pow(-1 + x,2))/M_PI;
	}
	if(power==4){
		return (-6*alphas_muR*pow(-1 + x,3)*(log(Q2/muF2) + 2*log(1 - x)))/M_PI;
	}
	if(power==5){
		return (3*alphas_muR*pow(-1 + x,4)*(7 + 10*log(Q2/muF2) + 20*log(1 - x)))/(5.*M_PI);
	}
	if(power==6){
		return (-3*alphas_muR*pow(-1 + x,5)*(21 + 20*log(Q2/muF2) + 40*log(1 - x)))/(10.*M_PI);
	}
	if(power==7){
		return (3*alphas_muR*pow(-1 + x,6)*(181 + 140*log(Q2/muF2) + 280*log(1 - x)))/(70.*M_PI);
	}
	if(power==8){
		return (-3*alphas_muR*pow(-1 + x,7)*(83 + 56*log(Q2/muF2) + 112*log(1 - x)))/(28.*M_PI);
	}
	if(power==9){
		return (alphas_muR*pow(-1 + x,8)*(4129 + 2520*log(Q2/muF2) + 5040*log(1 - x)))/(420.*M_PI);
	}
	if(power==10){
		return -(alphas_muR*pow(-1 + x,9)*(319 + 180*log(Q2/muF2) + 360*log(1 - x)))/(30.*M_PI);
	}
	if(power==-1){
		return higgs_NLO_gg_reg(x)-(-(alphas_muR*(-25878 + 147185*x - 445357*pow(x,2) + 835769*pow(x,3) - 1046335*pow(x,4) + 894569*pow(x,5) - 520159*pow(x,6) + 197543*pow(x,7) - 44323*pow(x,8) + 4466*pow(x,9) + 2520*(-8 + 44*x - 119*pow(x,2) + 210*pow(x,3) - 252*pow(x,4) + 210*pow(x,5) - 120*pow(x,6) + 45*pow(x,7) - 10*pow(x,8) + pow(x,9))*log(Q2/muF2) + 5040*(-8 + 44*x - 119*pow(x,2) + 210*pow(x,3) - 252*pow(x,4) + 210*pow(x,5) - 120*pow(x,6) + 45*pow(x,7) - 10*pow(x,8) + pow(x,9))*log(1 - x)))/(420.*M_PI));
	}
	if(power==-2){
		return -(alphas_muR*(-25878 + 147185*x - 445357*pow(x,2) + 835769*pow(x,3) - 1046335*pow(x,4) + 894569*pow(x,5) - 520159*pow(x,6) + 197543*pow(x,7) - 44323*pow(x,8) + 4466*pow(x,9) + 2520*(-8 + 44*x - 119*pow(x,2) + 210*pow(x,3) - 252*pow(x,4) + 210*pow(x,5) - 120*pow(x,6) + 45*pow(x,7) - 10*pow(x,8) + pow(x,9))*log(Q2/muF2) + 5040*(-8 + 44*x - 119*pow(x,2) + 210*pow(x,3) - 252*pow(x,4) + 210*pow(x,5) - 120*pow(x,6) + 45*pow(x,7) - 10*pow(x,8) + pow(x,9))*log(1 - x)))/(420.*M_PI);
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
double vegas_higgs_NLO_qg_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qg_full(z))*(jacobian*real(fit_sum_qg(tau/z+k[1]*jacobian,tau/z))/z);
}
double vegas_higgs_NLO_qg_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qg_expansion(z, fp->power))*(jacobian*pdf_sum_qg(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_qg_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qg_expansion(z, fp->power))*(jacobian*real(fit_sum_qg(tau/z+k[1]*jacobian,tau/z))/z);
}

///// NLO functions
double higgs_NLO_qg_full(double x){
	return (alphas_muR*(-3 - (-6 + x)*x + 2*(2 + (-2 + x)*x)*log(Q2/muF2) + 4*(2 + (-2 + x)*x)*log(1 - x) - 2*(2 + (-2 + x)*x)*log(x)))/(3.*M_PI*x);
}
double higgs_NLO_qg_expansion(double x, int power){		
	if(power==1){
		return (2*alphas_muR*(1 + log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==2){
		return (-2*alphas_muR*(-1 + x)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==3){
		return (4*alphas_muR*pow(-1 + x,2)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
	}
	if(power==4){
		return (-4*alphas_muR*pow(-1 + x,3)*(2 + 3*log(Q2/muF2) + 6*log(1 - x)))/(9.*M_PI);
	}
	if(power==5){
		return (alphas_muR*pow(-1 + x,4)*(25 + 24*log(Q2/muF2) + 48*log(1 - x)))/(18.*M_PI);
	}
	if(power==6){
		return (alphas_muR*pow(1 - x,5)*(5.233333333333333 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==7){
		return (alphas_muR*pow(-1 + x,6)*(91 + 60*log(Q2/muF2) + 120*log(1 - x)))/(45.*M_PI);
	}
	if(power==8){
		return (alphas_muR*pow(1 - x,7)*(6.752380952380952 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==9){
		return (alphas_muR*pow(-1 + x,8)*(7.335714285714285 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==10){
		return (alphas_muR*pow(1 - x,9)*(7.843650793650793 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
	}
	if(power==-1){
		return higgs_NLO_qg_full(x)-(-(alphas_muR*(-53002 + 332361*x - 1015440*pow(x,2) + 1888740*pow(x,3) - 2342928*pow(x,4) + 1993992*pow(x,5) - 1156176*pow(x,6) + 438240*pow(x,7) - 98190*pow(x,8) + 9883*pow(x,9) + 2520*(-18 + 89*x - 240*pow(x,2) + 420*pow(x,3) - 504*pow(x,4) + 420*pow(x,5) - 240*pow(x,6) + 90*pow(x,7) - 20*pow(x,8) + 2*pow(x,9))*log(Q2/muF2) + 5040*(-18 + 89*x - 240*pow(x,2) + 420*pow(x,3) - 504*pow(x,4) + 420*pow(x,5) - 240*pow(x,6) + 90*pow(x,7) - 20*pow(x,8) + 2*pow(x,9))*log(1 - x)))/(3780.*M_PI));
	}
	if(power==-2){
		return -(alphas_muR*(-53002 + 332361*x - 1015440*pow(x,2) + 1888740*pow(x,3) - 2342928*pow(x,4) + 1993992*pow(x,5) - 1156176*pow(x,6) + 438240*pow(x,7) - 98190*pow(x,8) + 9883*pow(x,9) + 2520*(-18 + 89*x - 240*pow(x,2) + 420*pow(x,3) - 504*pow(x,4) + 420*pow(x,5) - 240*pow(x,6) + 90*pow(x,7) - 20*pow(x,8) + 2*pow(x,9))*log(Q2/muF2) + 5040*(-18 + 89*x - 240*pow(x,2) + 420*pow(x,3) - 504*pow(x,4) + 420*pow(x,5) - 240*pow(x,6) + 90*pow(x,7) - 20*pow(x,8) + 2*pow(x,9))*log(1 - x)))/(3780.*M_PI);
	}
	else{return 0;}

}


///////////////////////////
/// qqbar channel
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
double vegas_higgs_NLO_qqbar_full_fit(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qqbar_full(z))*(jacobian*real(fit_sum_qqbar(tau/z+k[1]*jacobian,tau/z))/z);
}
double vegas_higgs_NLO_qqbar_power(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qqbar_expansion(z, fp->power))*(jacobian*pdf_sum_qqbar(tau/z+k[1]*jacobian,tau/z)/z);
}
double vegas_higgs_NLO_qqbar_power_fit(double *k, size_t dim, void *params){
	(void)(dim);
	struct lumni_params * fp = (struct lumni_params *)params;
	double z = k[0];
	double jacobian = 1.-tau/z; //needed to transform the boundary dependent terms
	double x = tau/z+k[1]*jacobian;
	return higgs_LO_factor()*(higgs_NLO_qqbar_expansion(z, fp->power))*(jacobian*real(fit_sum_qqbar(tau/z+k[1]*jacobian,tau/z))/z);
}

///// NLO functions
double higgs_NLO_qqbar_full(double x){
	return (-32*alphas_muR*pow(-1 + x,3))/(27.*M_PI*x);
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
		return (-32*alphas_muR*pow(-1 + x,3))/(27.*M_PI);
	}
	if(power==5){
		return (32*alphas_muR*pow(-1 + x,4))/(27.*M_PI);
	}
	if(power==6){
		return (-32*alphas_muR*pow(-1 + x,5))/(27.*M_PI);
	}
	if(power==7){
		return (32*alphas_muR*pow(-1 + x,6))/(27.*M_PI);
	}
	if(power==8){
		return (-32*alphas_muR*pow(-1 + x,7))/(27.*M_PI);
	}
	if(power==9){
		return (32*alphas_muR*pow(-1 + x,8))/(27.*M_PI);
	}
	if(power==10){
		return (-32*alphas_muR*pow(-1 + x,9))/(27.*M_PI);
	}
	if(power==-1){
		return higgs_NLO_qqbar_full(x)- (-32*alphas_muR*pow(-1 + x,3)*(7 - 21*x + 35*pow(x,2) - 35*pow(x,3) + 21*pow(x,4) - 7*pow(x,5) + pow(x,6)))/(27.*M_PI);
	}
	if(power==-2){
		return (-32*alphas_muR*pow(-1 + x,3)*(7 - 21*x + 35*pow(x,2) - 35*pow(x,3) + 21*pow(x,4) - 7*pow(x,5) + pow(x,6)))/(27.*M_PI);
	}
	else{
		return 0;
	}
}
